#' @export
read_matched_spike_data <- function(obj,
                          basedir = getwd()
) {
  library(tidyr)
  library(forcats)
  library(furrr)

  unames <- unique(obj$info$id)
  # Handle case where we have multiple ids
  if (length(unames) > 1) {
    temp = list()
    # Process each id
    for (i in 1:length(unames)) {
      obj_id <- list(call = obj$call,
                   info = obj$info %>% filter(id == unames[i]),
                   trial_data = obj$trial_data %>% filter(id == unames[i])
      )
      if (!is.null(obj$tracker_data)) {
        obj_id$tracker_data <- obj$tracker_data %>% filter(id == unames[i])
      } else {
        obj_id$tracker_data <- NULL
      }
      temp[[i]] <- read_matched_spike_data(obj_id, basedir = basedir)
    }

    # Bind together
    out <- list(call = obj$call,
                info = map_dfr(temp, ~.x$info),
                trial_data = map_dfr(temp, ~.x$trial_data),
                tracker_data = map_dfr(temp, ~.x$tracker_data),
                spike_times = map_dfr(temp, ~.x$spike_times),
                spike_mask = map_dfr(temp, ~.x$spike_mask),
                dropped = map(temp, ~.x$dropped),
                trial_duration = map_dfr(temp, ~.x$trial_duration)
    )
    class(out) <- "GNGeventide"

    return(out)
  }

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
    keep_session = obj$info$session[ind]
    dropped_session <- obj$info %>% filter(!(session %in% keep_session))
    obj$info %<>% filter(session %in% keep_session)
    obj$info$session <- as.integer(fct_drop(as_factor(obj$info$session)))
    obj$trial_data %<>% filter(session %in% keep_session)
    obj$trial_data$session <- as.integer(fct_drop(as_factor(obj$trial_data$session)))
    if (!is.null(obj$tracker_data)) {
      obj$tracker_data %<>% filter(session %in% keep_session)
      obj$tracker_data$session <- as.integer(fct_drop(as_factor(obj$tracker_data$session)))
    }
  } else {
    # return NULL
  }

  fnames_spk <- paste(basedir, paste0(base_name[ind], ".mat"), sep = .Platform$file.sep)
  spike_list <- furrr::future_map(fnames_spk, read_spike)
  #spike_list <- purrr::map(fnames_spk, read_spike)

  # Add another level of nesting since we can have different #neurons per session
  spike_times <- purrr::map(spike_list, ~.x$spktimes %>% nest(neurons = starts_with("AD"))) %>%
    bind_rows(.id = "session") %>%
    mutate(session = as.integer(session))

  spike_mask <- purrr::map(spike_list, ~.x$quality %>% nest(neurons = starts_with("AD"))) %>%
    bind_rows(.id = "session") %>%
    mutate(session = as.integer(session))

  # Restrict to matching trials (within sessions)
  dropped_eventide <- obj$trial_data %>% anti_join(spike_times, by = c("session", "counter_total_trials"))
  obj$trial_data %<>% anti_join(dropped_eventide, by = c("session", "counter_total_trials"))
  if (!is.null(obj$tracker_data)) {
    obj$tracker_data %<>% anti_join(dropped_eventide, by = c("session", "counter_total_trials"))
  }

  dropped_spike <- spike_times %>% anti_join(obj$trial_data, by = c("session", "counter_total_trials"))
  spike_times %<>% anti_join(dropped_spike, by = c("session", "counter_total_trials"))
  spike_mask %<>% anti_join(dropped_spike, by = c("session", "counter_total_trials"))

  # Update info
  info <- purrr::map_dfr(spike_list , "session_info", .id = "session") %>%
    inner_join(purrr::map_dfr(spike_list , "neuron_info", .id = "session") %>%
                 group_by(session) %>%
                 nest(neuron_info = !session),
               by = "session") %>%
    mutate(session = as.integer(session), target = stringr::str_to_lower(target))
  obj$info %<>% inner_join(info, by = c("session", "fname_eventide"))

  # Calculate trial durations from eventide (time between define trial state)
  trial_dur_eventide <- obj$trial_data %>%
    group_by(session) %>%
    mutate(duration_eventide = dplyr::lead(define_trial_onset_time_absolute) - define_trial_onset_time_absolute) %>%
    select(session, counter_total_trials, define_trial_onset_time_absolute, duration_eventide)

  # Calculate trials durations from plexon
  trial_dur <- purrr::map_dfr(spike_list , "trigger_timestamps", .id = "session") %>%
    mutate(session = as.integer(session)) %>%
    group_by(session) %>%
    mutate(counter_total_trials = row_number(), .after = session) %>%
    mutate(duration_spike = dplyr::lead(t) - t) %>%
    inner_join(trial_dur_eventide, by = c("session", "counter_total_trials")) %>%
    select(-define_trial_onset_time_absolute) %>%
    mutate(dt = duration_eventide - duration_spike) %>%
    mutate(flag_ms = ifelse(dt > 0.002, dt*1000, NA))

  # Pack into output
  obj$spike_times <- spike_times %>%
    mutate(id = obj$trial_data$id, .before = 1)
  obj$spike_mask <- spike_mask %>%
    mutate(id = obj$trial_data$id, .before = 1)
  obj$dropped <- list(session = dropped_session,
                      eventide = dropped_eventide,
                      spike = dropped_spike)
  obj$trial_duration <- trial_dur %>%
    ungroup() %>%
    mutate(id = obj$trial_data$id, .before = 1)

  return(obj)
}
