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

  f <- function(x) as.vector(unlist(x))
  spk <- apply(spk, c(1,2), f)
  spktimes <- as_tibble(spk) %>%
    add_column(counter_total_trials = 1:nrow(spk), .before = 1)

  quality <- dat$quality
  #colnames(quality) <- paste0("q", neuron_info$name)
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

# merge.GNGspikedata <- function() {
# # use eventide filename to find and load spikes
#
#   class(out) <- "GNGmerged"
# }

#' @export
count_in_window <- function(x, t1, t2) {
  out = lapply(x, function(xx) {length(which((xx>=t1) & (xx<=t2)))})
  unlist(out)
}

#' @export
align_to <- function(x, t) {
  unlist(x) - t
}

#normalize <- function(x, )
#' @export
count_spks_in_window <- function(df, align = "cue_onset_time") {
    df2 <- df %>% filter(!is_abort_trial) %>%
      mutate(across(starts_with("AD"), ~map2(.x, cue_onset_time, ~align_to(.x,.y)))) %>%  # Align to event
      select(counter_total_trials, starts_with("AD")) %>%
      pivot_longer(!counter_total_trials) %>%
      mutate(count = count_in_window(value, t1, t2)) %>%
      select(-value) %>%
      left_join(df %>% select(counter_total_trials, condition, rt))

    return(df2)
  }
