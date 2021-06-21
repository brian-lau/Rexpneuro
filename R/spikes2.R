#' @export
read_matched_spike_data <- function(obj,
                          basedir = getwd(),
                          resort = TRUE
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
      temp[[i]] <- read_matched_spike_data(obj_id, basedir = basedir, resort = resort)
    }

    # Bind together
    out <- list(call = obj$call,
                info = purrr::map_dfr(temp, ~.x$info),
                trial_data = purrr::map_dfr(temp, ~.x$trial_data),
                tracker_data = purrr::map_dfr(temp, ~.x$tracker_data),
                spike_times = purrr::map_dfr(temp, ~.x$spike_times),
                spike_mask = purrr::map_dfr(temp, ~.x$spike_mask),
                dropped = purrr::map(temp, ~.x$dropped),
                trial_duration = purrr::map_dfr(temp, ~.x$trial_duration)
    )
    class(out) <- "GNGeventide"

    return(out)
  }

  # Incoming eventide files to match to corresponding spike data
  fnames_eventide <- obj$info$fname_eventide
  base_name <- purrr::map_chr(stringr::str_split(fnames_eventide, ".txt"), 1)

  # Available spike data
  if (resort) {
    base_name <- paste0(base_name,'_resort')
    fnames <- list.files(path = basedir,
                         pattern = glob2rx(paste0("*GNG", "*_resort.mat")))
  } else {
    fnames <- list.files(path = basedir,
                         pattern = glob2rx(paste0("*GNG", "*.mat")))
  }
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
  if (resort) {
    spike_list <- furrr::future_map(fnames_spk, read_spike_resort)
  } else {
    spike_list <- furrr::future_map(fnames_spk, read_spike)
  }
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

#' @export
restrict_neurons <- function(obj,
                             good # data frame must include id, fname_ephys, fname_eventide, name, good
) {
  # restrict_fname = paste0(basedir, "STN_good_neuron_list.csv")
  # good = as_tibble(read.csv(restrict_fname, header = T))
  # good %<>% filter(good == "all")
  neuron_info <- obj$info %>%
    unnest(cols = neuron_info, names_repair = 'universal')

  good_info <- neuron_info %>% semi_join(good, by = c("id", "fname_ephys", "fname_eventide", "name"))
  bad_info <- neuron_info %>% anti_join(good, by = c("id", "fname_ephys", "fname_eventide", "name"))

  # un-nest
  spike_times <- obj$spike_times %>% unnest_spike_times()
  # filtering join
  spike_times %<>% semi_join(good_info, by = c("id", "session", "name"))

  # re-nest, this is clearly not efficient...
  pivot_nest <- function(x, y) {
    x %>% pivot_wider(names_from = name, values_from = times) %>% nest(neurons = everything())
  }
  tic()
  spike_times %<>% group_by(id, session, counter_total_trials) %>%
    group_modify(.f = ~pivot_nest(.x))
  toc()

  #   # re-nest
  # spike_times %<>% group_by(id, session, counter_total_trials) %>%
  #   nest() %>%
  #   rename(neurons = data)

  spike_mask <- obj$spike_mask %>% unnest_spike_mask()
  spike_mask %<>% semi_join(good_info, by = c("id", "session", "name"))
  # spike_mask %<>% group_by(id, session, counter_total_trials) %>%
  #   nest() %>%
  #   rename(neurons = data)
  pivot_nest2 <- function(x, y) {
    x %>% pivot_wider(names_from = name, values_from = mask) %>% nest(neurons = everything())
  }
  tic()
  spike_mask %<>% group_by(id, session, counter_total_trials) %>%
    group_modify(.f = ~pivot_nest2(.x))
  toc()


  # spike_times %>% group_by(id, session, counter_total_trials) %>%
  #   pivot_wider(names_from = name, values_from = times) %>%
  #   nest(neurons = starts_with("AD"))

  good_info %<>% nest( data = c(name, channel, channelName, unit, tStart,
                                tEnd, gap)) %>%
    rename(neuron_info = data)

  obj$info <- good_info
  obj$trial_data %<>% semi_join(good_info, by = c("id", "session"))
  if (nrow(obj$tracker_data) != 0) {
    obj$tracker_data %>% semi_join(good_info, by = c("id", "session"))
  }
  obj$spike_times <- spike_times
  obj$spike_mask <- spike_mask
  obj$trial_duration %<>% semi_join(good_info, by = c("id", "session"))

  return(obj)
}
