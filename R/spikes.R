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
                                 target = trimws(session["target"]),
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
  # https://stackoverflow.com/questions/13970333/fast-vectorized-function-to-check-if-a-value-is-in-an-interval
  # https://stackoverflow.com/questions/2190756/how-to-count-true-values-in-a-logical-vector
  # f1 <- function(xx) {length(which((xx>=t1) & (xx<=t2)))}
  # f2 <- function(xx) {tabulate(.bincode(xx, c(t1,t2)))}
  # f3 <- function(xx) {sum(!is.na(.bincode(xx,c(t1,t2))))}
  # f4 <- function(xx) {sum(.bincode(xx,c(t1,t2)), na.rm = TRUE)}
  # f5 <- function(xx) {sum(complete.cases(.bincode(xx,c(t1,t2))))}
  # microbenchmark(lapply(x, f1), lapply(x, f2), lapply(x, f3), lapply(x, f4), lapply(x, f5), times = 10)

  out = lapply(x, function(xx) {sum(.bincode(xx,c(t1,t2)), na.rm = TRUE)})
  unlist(out)
}

#' @export
in_window <- function(x, t1, t2) {
  out = lapply(x, function(xx) {xx[!is.na(.bincode(xx,c(t1,t2)))]})
}

#' @export
bin_times <- function(x, breaks) {
  # f1 <- function(xx) { xx[xx>=breaks[1] & xx<=breaks[length(breaks)]] }
  # f2 <- function(xx) { tabulate(.bincode(xx, breaks, include.lowest = T), nbins = length(breaks) - 1) }
  # microbenchmark(lapply(lapply(x, f1), f2), lapply(x, f2), times = 10)
  # all.equal(lapply(lapply(x, f1), f2), lapply(x, f2))
  lapply(x, function(xx) { tabulate(.bincode(xx, breaks, include.lowest = T), nbins = length(breaks) - 1) })
}

#' #' @export
#' align_to <- function(x, group_info) {
#'   #library(furrr)
#'   #albrowser()
#'   x %>% unnest(cols = c(neurons)) %>%
#'     mutate(across(starts_with("AD"), ~map2(.x, shift, ~.x - .y))) %>%
#'     select(-shift) %>%
#'     pivot_longer(!counter_total_trials, values_to = "times")
#'     #nest(neurons = starts_with("AD")) # re-nest?
#' }

#' #' @export
#' add_covariates <- function(x, group_info, y) {
#'   x %>% left_join(
#'     y %>% filter(id == group_info$id, session == group_info$session),
#'                   by = c("id", "session", "counter_total_trials")) %>%
#'     select(-id, -session)
#' }
#'
#' #' @export
#' drop_aborts <- function(x, group_info, y) {
#'   x %>% semi_join(
#'     y %>% filter(id == group_info$id, session == group_info$session, !is_abort),
#'                   by = c("id", "session", "counter_total_trials")) %>%
#'     select(-id, -session)
#' }

#' @export
drop_neurons <- function(x, group_info, min_trial = 50) {
  if ("t" %in% names(x)) {
    temp <- x %>% group_by(uname) %>%
      filter(t == unique(t)[1]) %>%
      mutate(n = n()) %>% slice(1) %>%
      filter(n>min_trial)
  } else {
    temp <- x %>% group_by(uname) %>%
      mutate(n = n()) %>% slice(1) %>%
      filter(n>min_trial)
  }
  x %>% filter(uname %in% temp$uname)
}

#' @export
spks_in_window <- function(obj,
                           align = "cue_onset_time",
                           t_start = -0.5,
                           t_end = 0,
                           binwidth = t_end,
                           mask_by_quality = TRUE
) {
  # Pivot spike quality
  m <- obj$spike_mask %>% unnest_spike_mask()

  # Long format spike times
  t_long <- obj$spike_times %>%
    bind_cols(shift = obj[["trial_data"]][[align]]) %>%
    unnest_spike_times()

  # Align spike times to event
  t <- t_long %>%
    mutate(times = map2(times, shift, ~.x - .y)) %>%
    ungroup() %>%
    select(-shift) %>%
    filter(m$mask > 0) %>%
    mutate(times2 = in_window(times, t_start, t_end)) # Drop spikes outside of window

  # Create unique label for neurons
  t %<>%
    mutate(uname = fct_cross(as_factor(id), as_factor(session), name), .after = name) %>%
    arrange(id, session, uname)

  return(t %>% ungroup())
}

#' @export
unnest_spike_mask <- function(df) {
  df %>%
    group_by(id, session) %>%
    group_modify(function(x, y) x %>%
                   unnest(cols = c(neurons)) %>%
                   pivot_longer(!counter_total_trials, values_to = "mask")) %>%
    ungroup()
}

#' @export
unnest_spike_times <- function(df) {
  df %>%
    group_by(id, session) %>%
    group_modify(function(x, y) x %>%
                   unnest(cols = c(neurons)) %>%
                   pivot_longer(starts_with("AD"), values_to = "times")) %>%
    ungroup()
}

#' @export
get_psth <- function(obj,
                     align = "cue_onset_time",
                     method = "kde", # should add bin
                     t_start = 0,
                     t_end = 0.5,
                     pre_trunc_event = "none", # for censored averaging
                     post_trunc_event = "none",
                     dt = 0.005,
                     h = 0.025,
                     pad = 5,
                     mask_by_quality = TRUE,
                     drop_abort_trials = TRUE,
                     drop_incorrect_trials = TRUE,
                     add_trial_vars = c("is_correct",
                                        "is_incorrect",
                                        "is_abort",
                                        "block",
                                        "gng",
                                        "direction",
                                        "rt"),
                     add_trial_time_vars = c("fix_onset_time",
                                             "fix_acquired_onset_time",
                                             "cue_onset_time",
                                             "target_onset_time",
                                             "liftoff_onset_time"), # Will be shifted by align eventtime
                     min_trial = 50
) {
  # Pivot spike quality
  m <- obj$spike_mask %>% unnest_spike_mask()

  # Long format spike times
  t_long <- obj$spike_times %>%
    bind_cols(shift = obj[["trial_data"]][[align]]) %>%
    bind_cols(pre_trunc = obj[["trial_data"]][[pre_trunc_event]]) %>%
    bind_cols(post_trunc = obj[["trial_data"]][[post_trunc_event]]) %>%
    unnest_spike_times()

  # Align spike times to event
  t <- t_long %>% mutate(times = map2(times, shift, ~.x - .y))
  if (pre_trunc_event == "none") {
    t %<>% mutate(pre_trunc = t_start)
  } else {
    t %<>% mutate(pre_trunc = pre_trunc - shift)
  }
  if (post_trunc_event == "none") {
    t %<>% mutate(post_trunc = t_end)
  } else {
    t %<>% mutate(post_trunc = post_trunc - shift)
  }

  if (mask_by_quality) t %<>% filter(m$mask > 0)

  # Smooth
  x_eval <- seq(t_start, t_end, by = dt)
  t %<>% mutate(psth = map(times, ~smpsth(t = .x, h = h,
                                           from = t_start,
                                           to = t_end,
                                           ngrid = length(x_eval),
                                           pad = pad),
                           )) %>%
    select(-times)

  # Create unique label for neurons
  t %<>%
    mutate(uname = fct_cross(as_factor(id), as_factor(session), name), .after = name) %>%
    arrange(id, session, uname)

  # Add trial, time-independent covariates
  t %<>% ungroup %>%
    left_join(obj$trial_data %>% select(id, session, counter_total_trials, all_of(add_trial_vars)),
                              by = c("id", "session", "counter_total_trials"))

  # Add trial, time-dependent covariates (shift to alignment event)
  t %<>%
    left_join(obj$trial_data %>% select(id, session, counter_total_trials, all_of(add_trial_time_vars)),
              by = c("id", "session", "counter_total_trials")) %>%
    mutate(across(all_of(add_trial_time_vars), ~.x - shift))

  if (drop_abort_trials) t %<>% filter(!is_abort)

  if (drop_incorrect_trials) t %<>% filter(!is_incorrect)

  if (min_trial > 0) t %<>% group_modify(drop_neurons, min_trial = min_trial)

  # Censor psth outside of pre- and post-events
  if ((pre_trunc_event != "none") | (post_trunc_event != "none")) {
    mask_t <- function(x, x_eval, pre, post) {
      ind = is.na(.bincode(x_eval, c(pre,post)))
      x[ind] = NA
      return(x)
    }
    t %<>% ungroup() %>% rowwise() %>% mutate(psth = list(mask_t(psth, x_eval, pre_trunc, post_trunc)))
  }

  return(t %>% ungroup())
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
                           drop_abort_trials = TRUE,
                           drop_incorrect_trials = TRUE,
                           add_trial_vars = c("is_correct",
                                              "is_incorrect",
                                              "is_abort",
                                              "block",
                                              "gng",
                                              "direction",
                                              "rt"),
                           min_trial = 50
                           ) {
  # Pivot spike quality
  m <- obj$spike_mask %>% unnest_spike_mask()

  # Long format spike times
  t_long <- obj$spike_times %>%
    bind_cols(shift = obj[["trial_data"]][[align]],
              bl_shift = obj[["trial_data"]][[bl_align]]) %>%
    unnest_spike_times()

  # Align spike times to event
  t <- t_long %>%
    mutate(times = map2(times, shift, ~.x - .y)) %>%
    select(-shift, -bl_shift)

  if (mask_by_quality) t %<>% filter(m$mask > 0)

  # Count spikes in baseline window & summarise mean and std
  bl <- t_long %>%
    mutate(times = map2(times, bl_shift, ~.x - .y)) %>%
    select(-shift, -bl_shift)

  if (mask_by_quality) bl %<>% filter(m$mask > 0)

  bl %<>% mutate(fr_bl = count_in_window(times, bl_t_start, bl_t_end)/(bl_t_end - bl_t_start)) %>%
    select(-times) %>%
    group_by(id, session, name) %>%
    summarise(fr_bl_mean = mean(fr_bl), fr_bl_sd = sd(fr_bl), .groups = "drop")

  # Create unique label for neurons
  t %<>%
    mutate(uname = fct_cross(as_factor(id), as_factor(session), name), .after = name) %>%
    arrange(id, session, uname)

  # Count spikes in windows
  breaks = seq(t_start, t_end, by = binwidth)
  mids = breaks[-length(breaks)] + diff(breaks)/2
  c <- t %>%
    mutate(t = list(mids)) %>%
    mutate(binned =  bin_times(times, breaks)) %>%
    select(-times) %>%
    unnest(cols = c(t, binned)) %>%
    mutate(fr = binned/binwidth)

  # Create normalized firing rates
  c %<>% left_join(bl, by = c("id", "session", "name")) %>%
    mutate(fr_norm1 = fr - fr_bl_mean,
           fr_norm2 = fr/fr_bl_mean,
           fr_norm3 = (fr - fr_bl_mean) / fr_bl_sd)

  # Add trial, time-independent covariates
  c %<>% ungroup %>%
    left_join(obj$trial_data %>% select(id, session, counter_total_trials, all_of(add_trial_vars)),
              by = c("id", "session", "counter_total_trials"))

  if (drop_abort_trials) c %<>% filter(!is_abort)

  if (drop_incorrect_trials) c %<>% filter(!is_incorrect)

  # # Bind covariates
  # c %<>% group_by(id, session) %>%
  #   group_modify(add_covariates,
  #                obj$trial_data %>%
  #                  select(id, session, counter_total_trials, block, gng, direction),
  #                .keep = TRUE)
  #
  # if (drop_abort_trials) {
  #   c %<>%
  #     group_modify(drop_aborts,
  #                  obj$trial_data %>%
  #                    select(id, session, counter_total_trials, is_abort),
  #                  .keep = TRUE)
  # }

  if (min_trial > 0) c %<>% group_modify(drop_neurons, min_trial = min_trial)

  return(c %>% ungroup())
}

#' @export
regularity <- function(t, R = 0.005) {
  # t <- c(0.0097, 0.0272, 0.0615, 0.0779, 0.1918, 0.2574, 0.4438, 0.4561, 0.7816, 0.9658)
  isi <- diff(sort(t))
  n = length(isi)

  if (n==0) {
    return(tibble(cv = NA, cv2 = NA, lv = NA, lvr = NA))
  }

  cv <- sd(isi)/mean(isi)
  isi1 <- isi[1:(n-1)]
  isi2 <- isi[2:n]
  cv2 <- mean(2*abs(isi2 - isi1)/(isi2 + isi1))
  lv <- 3*mean(((isi1 - isi2)/(isi1 + isi2))^2)
  lvr <- 3*mean((1 - (4*isi1*isi2)/(isi1 + isi2)^2) * (1 + (4*R)/(isi1 + isi2)))

  #return(tibble(cv = cv, cv2 = cv2, lv = lv, lvr = lvr))
  return(list(cv = cv, cv2 = cv2, lv = lv, lvr = lvr))
}

#' @export
smpsth <- function(t, from, to, ngrid = 1000, h = 0.025, pad = 0, ...) {
  if (length(t) != 0) {
    if (pad) {
      dt = (to - from)/ngrid
      from = from - dt*pad
      to = to + dt*pad
      ngrid = ngrid + 2*pad
    }
    y = KernSmooth::bkde(t, bandwidth = h, range.x = c(from, to),
                     gridsize = ngrid, truncate = TRUE)[["y"]]*length(t)
    if (pad) {
      y = y[(pad+1):(length(y)-pad)]
    }
    return(y)
  } else {
    rep(0, ngrid)
  }
}

#' @export
filter_runlengths <- function(x, n, filtered_value = 0) {
  temp = rle(x) %>%
    unclass() %>%
    as.data.frame() %>%
    rowwise() %>%
    mutate(values = ifelse(lengths < n, filtered_value, values))

  unlist(temp %>%
           rowwise() %>%
           mutate(y = list(rep(values,lengths))) %>%
           pull(y))
}
