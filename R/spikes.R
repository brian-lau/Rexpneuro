#' @export
read_spike <- function(fname) {
  library(magrittr)
  library(dplyr)
  library(tibble)

  dat <- R.matlab::readMat(fname)

  ## Neuron info
  # http://lukaspuettmann.com/2017/03/13/matlab-struct-to-r-dataframe/
  varNames    <- names(dat$neuron.info[,,1])
  datList     <- dat$neuron.info
  temp = (t(matrix(unlist(datList), nrow = length(varNames))))
  colnames(temp) <- varNames
  neuron_info = as_tibble(temp)
  neuron_info %<>% mutate(tStart = as.numeric(tStart),
                          tEnd = as.numeric(tEnd),
                          channel = as.integer(channel),
                          unit = as.integer(unit),
                          gap = as.integer(gap)
                          ) %>%
    relocate(tStart, .before = tEnd)

  session = temp = unlist(dat$session[, , 1])

  spk = matrix(dat$times, ncol = length(neuron_info$name))
  colnames(spk) <- neuron_info$name
  f <- function(x) unlist(x, use.names = F)
  spk <- apply(spk, c(1, 2), f)

  spktimes <- as_tibble(spk) %>%
    add_column(counter_total_trials = 1:nrow(spk), .before = 1)

  quality <- dat$quality
  colnames(quality) <- neuron_info$name
  quality <- as_tibble(quality) %>%
    add_column(counter_total_trials = 1:nrow(spk), .before = 1)

  session_info <- tibble::tibble(fname_ephys = session["eFname"],
                                 fname_eventide = session["bFname"],
                                 trigger = session["trigger"],
                                 artifact = as.logical(session["artifact"]),
                                 target = session["target"],
                                 grid_x = session["grid.x"],
                                 grid_y = session["grid.y"],
                                 tip_depth = as.numeric(session["depth"])
                                 )

  out <- list(
    session_info = session_info,
    trigger_timestamps = tibble(t = as.vector(dat$event.timestamps)),
    neuron_info = neuron_info,
    spktimes = spktimes,
    quality = quality
  )

  class(out) <- "GNGspikedata"

  return(out)
}

#' @export
count_in_window <- function(x, t1, t2) {
  out = lapply(x, function(xx) {length(which((xx>=t1) & (xx<=t2)))})
  unlist(out)
}

#' @export
bin_times <- function(x, breaks) {
  out = lapply(x, function(xx) { xx[xx>=breaks[1] & xx<=breaks[length(breaks)]] })
  out = lapply(out, function(xx) { tabulate(.bincode(xx, breaks, include.lowest = T), nbins = length(breaks) - 1) })
}

#' @export
align_to <- function(x, group_info) {
  library(purrr)
  x %>% unnest(cols = c(neurons)) %>%
    mutate(across(starts_with("AD"), ~map2(.x, shift, ~.x - .y))) %>%
    select(-shift) %>%
    pivot_longer(!counter_total_trials, values_to = "times")
    #nest(neurons = starts_with("AD")) # re-nest?
}

#' @export
add_covariates <- function(x, group_info, y) {
  x %>% left_join(y %>% filter(session == group_info$session),
                  by = c("session", "counter_total_trials")) %>%
    select(-session)
}

#' @export
drop_trials <- function(x, group_info, y) {
  x %>% semi_join(y %>% filter(session == group_info$session, !is_abort),
                  by = c("session", "counter_total_trials")) %>%
    select(-session)
}

#' @export
drop_neurons <- function(x, group_info, min_trial = 50) {
  temp <- x %>% group_by(uname) %>%
    filter(t == unique(t)[1]) %>%
    mutate(n = n()) %>% slice(1) %>%
    filter(n>min_trial)
  x %>% filter(uname %in% temp$uname)
}

#' @export
prep_for_model <- function(obj,
                           align = "cue_onset_time",
                           t_start = 0,
                           t_end = 0.5,
                           binwidth = t_end,
                           bl_align = "cue_onset_time",
                           bl_t_start = -0.5,
                           bl_t_end = 0,
                           mask_by_quality = TRUE,
                           min_trial = 50,
                           normalize = TRUE
                           ) {
  # Pivot spike quality
  m <- obj$spike_mask %>%
    group_by(session) %>%
    group_modify(function(x,y) x %>%
                   unnest(cols = c(neurons)) %>%
                   pivot_longer(!counter_total_trials, values_to = "mask")) %>%
    ungroup()

  # Align spike times to event
  t <- obj$spike_times %>%
    bind_cols(shift = obj[["trial_data"]][[align]]) %>%
    group_by(session) %>%
    group_modify(align_to) %>%
    ungroup() %>%
    filter(m$mask > 0)

  # Count spikes in baseline window & summarise mean and std
  bl <- obj$spike_times %>%
    bind_cols(shift = obj[["trial_data"]][[bl_align]]) %>%
    group_by(session) %>%
    group_modify(align_to) %>%
    ungroup()  %>%
    filter(m$mask > 0) %>%
    mutate(fr_bl = count_in_window(times, bl_t_start, bl_t_end)/(bl_t_end - bl_t_start)) %>%
    select(-times) %>%
    group_by(session, name) %>%
    summarise(fr_bl_mean = mean(fr_bl), fr_bl_sd = sd(fr_bl), .groups = "drop")

  # Count spikes in windows
  #bw = .25
  breaks = seq(t_start, t_end, by = binwidth)
  mids = breaks[-length(breaks)] + diff(breaks)/2
  c <- t %>%
    mutate(uname = fct_cross(as_factor(session), name), .after = name) %>%
    arrange(session, uname) %>%
    mutate(t = list(mids)) %>%
    mutate(binned =  bin_times(times, breaks)) %>%
    select(-times) %>%
    unnest(cols = c(t, binned)) %>%
    mutate(fr = binned/binwidth)

  # Create normalized firing rates
  c %<>% left_join(bl, by = c("session", "name")) %>%
    mutate(fr_norm1 = fr - fr_bl_mean,
           fr_norm2 = fr/fr_bl_mean,
           fr_norm3 = (fr - fr_bl_mean) / fr_bl_sd)

  # Bind covariates
  c %<>% group_by(session) %>%
    group_modify(add_covariates, obj$trial_data %>% select(session, counter_total_trials, block, gng, direction), .keep = TRUE) %>%
    ungroup()

  # Drop trials
  c %<>% group_by(session) %>%
    group_modify(drop_trials, obj$trial_data %>% select(session, counter_total_trials, is_correct, is_incorrect, is_abort), .keep = TRUE)

  c %<>% group_by(session) %>%
    group_modify(drop_neurons, min_trial = min_trial)

  return(c)
}
