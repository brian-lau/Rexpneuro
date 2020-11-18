#' @export
read_matching_spkdata <- function(obj,
                          basedir = getwd()
) {
  library(tidyr)

  # Incoming eventide files to match to corresponding spike data
  fnames_eventide <- obj$info$fname_eventide
  base_name <- purrr::map_chr(stringr::str_split(fnames_eventide, ".txt"), 1)

  # Available spike data
  fnames <- list.files(path = basedir,
                       pattern = glob2rx(paste0("*GNG", "*.mat")))
  available_name <- purrr::map_chr(stringr::str_split(fnames, ".mat"), 1)

  # Which eventide sessions having corresponding spike data?
  ind <- base_name %in% available_name

  # Drop eventide and tracker sessions where there is no spike data
  if (any(ind)) {
    keep_sessions = obj$info$session[ind]
    obj$info %<>% filter(session %in% keep_sessions)
    obj$info$session <- as.integer(droplevels(as.factor(obj$info$session)))
    obj$trial_data %<>% filter(session %in% keep_sessions)
    obj$trial_data$session <- as.integer(droplevels(as.factor(obj$trial_data$session)))
    if (!is.null(obj$tracker_data)) {
      obj$tracker_data %<>% filter(session %in% keep_sessions)
      obj$tracker_data$session <- as.integer(droplevels(as.factor(obj$tracker_data$session)))
    }
  } else {
    # return NULL
  }

  fnames_spk <- paste0(base_name[ind], ".mat")
  spike_list <- purrr::map(fnames_spk, read_spike)

  spike_times <- purrr::map_dfr(spike_list , "spktimes", .id = "session") %>%
    mutate(session = as.integer(session))
  spike_mask <- purrr::map_dfr(spike_list , "quality", .id = "session") %>%
    mutate(session = as.integer(session))

  # How to check matching trials? Drop in eventide and tracker?
  #spk_data %<>% semi_join(obj$trial_data, by = c("session", "counter_total_trials"))
  #diff(obj$trial_data$define_trial_onset_time_absolute)
  #diff(spk_list[[1]]$trigger_timestamps)

  dropped_eventide <- obj$trial_data %>% anti_join(spike_times, by = c("session", "counter_total_trials"))
  obj$trial_data %<>% anti_join(dropped_eventide, by = c("session", "counter_total_trials"))
  if (!is.null(obj$tracker_data)) {
    obj$tracker_data %<>% anti_join(dropped_eventide, by = c("session", "counter_total_trials"))
  }

  dropped_spike <- spike_times %>% anti_join(obj$trial_data, by = c("session", "counter_total_trials"))
  spike_times %<>% anti_join(dropped_spike, by = c("session", "counter_total_trials"))
  spike_mask %<>% anti_join(dropped_spike, by = c("session", "counter_total_trials"))

  info <- purrr::map_dfr(spk_list , "session_info", .id = "session") %>%
    inner_join(purrr::map_dfr(spk_list , "neuron_info", .id = "session") %>% group_by(session) %>% nest(neuron_info = !session)) %>%
    mutate(session = as.integer(session))
  obj$info %<>% inner_join(info)

  # out <- list(
  #   evi = obj,
  #   spk = spkdat,
  # )

  return(out)
}
