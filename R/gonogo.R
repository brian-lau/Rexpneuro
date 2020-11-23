#' @export
read_eventide <- function(fname = NULL,
                          name = NULL,
                          basedir = getwd(),
                          start_date = "30012017", # daymonthyear
                          end_date = "30012021",   # daymonthyear
                          min_trials = 1,
                          include_tracker = FALSE,
                          include_spike = FALSE,
                          ...
) {
  library(magrittr)
  library(dplyr)

  if (is.null(names) & is.null(fname)) stop("Specify subject name or filename.")

  if (!is.null(fname)) {
    fnames <- fname
    name <- stringr::str_split(fnames,"_")[[1]][1]
  } else {
    fnames <- list.files(path = basedir,
                         pattern = glob2rx(paste0(name, "_", "GNG", "*.txt")),
                         ignore.case = TRUE)
  }

  d <- purrr::map_chr(stringr::str_split(fnames,"_"), 3)
  t <- purrr::map_chr(stringr::str_split(fnames,"_"), 4)
  t <- purrr::map_chr(stringr::str_split(t,".txt"), 1)
  d <- as.POSIXct(paste(d, t),
                  "%d%m%Y %H-%M", tz = "Europe/Paris")

  ## Eventide trial data
  # Sort by ascending experiment date
  ind <- order(d)
  fnames <- fnames[ind]
  d <- d[ind]

  start_date <- as.POSIXlt(start_date, format = "%d%m%Y", tz = "Europe/Paris")
  end_date <- as.POSIXlt(end_date, format = "%d%m%Y", tz = "Europe/Paris")
  ind <- (start_date <= d) & (d <= end_date)
  fnames <- fnames[ind]
  d <- d[ind]

  ## Eventide tracker data
  if (include_tracker) {
    td <- parse_tracker_filenames(basedir = basedir)
    #dt <- difftime(d, td$date, units = "secs")

    # Find index for matching tracker file by time difference
    dt <- t(matrix(unlist(purrr::map(d, ~(abs(difftime(.x, td$date, units = "secs"))))), ncol = length(d)))
    mind <- which(dt==matrixStats::rowMins(dt),arr.ind=T)
    #map_int(d, ~which.min(abs(difftime(.x, td$date, units = "secs"))))
    dt <- dt[mind]

    td$fnames <- td$fnames[mind[,2]]
    td$date <- td$date[mind[,2]]
    td$dt <- dt
  }

  if (!any(ind)) {
    out <- list(
      call = match.call(),
      name = name,
      info = NULL,
      trial_data = NULL
    )
  } else {
    # Read session data
    dat <- purrr::map(paste(basedir, fnames, sep = .Platform$file.sep), read_eventide_single, ...)
    trial_data = purrr::map_dfr(dat, "trial_data", .id = "session")
    info <- tibble::tibble(session = unique(trial_data$session),
                           version = purrr::map_chr(dat, "version"),
                           date = unlist(purrr::map(dat, "date") %>% purrr::reduce(c)),
                           fname_eventide = fnames
    )

    if (include_tracker) {
      dat_tracker <- purrr::map(paste(basedir, td$fnames, sep = .Platform$file.sep), read_eventide_tracker, ...)
      tracker_data <- purrr::map_dfr(dat_tracker, "tracker_data", .id = "session")
      info %<>% tibble::add_column(date_tracker = td$date,
                                   fname_tracker = td$fnames,
                                   dt = dt
      )

      # Remove extra tracker data (from session terminating before trial_data written)
      tracker_data %<>% group_by(session, counter_total_trials) %>%
        nest() %>% semi_join(trial_data, by = c("session", "counter_total_trials"))

      # Convert to trial time
      if ("define_trial_onset_time_absolute" %in% names(trial_data)) {
        tracker_data$define_trial_onset_time_absolute = trial_data$define_trial_onset_time_absolute
        f <- function(df, t0) {
          df$t <- df$t - t0
          return(df)
        }
        tracker_data %<>% mutate(data = map(data, ~f(.x, define_trial_onset_time_absolute)))
      }
    } else {
      tracker_data <- NULL
    }

    out <- list(
      call = match.call(),
      name = name,
      info = info,
      trial_data = trial_data,
      tracker_data = tracker_data
    )
    #temp = trial_data %>% nest_join(tracker_data)
  }

  class(out) <- "GNGeventide"

  return(out)
}

#' @export
summary.GNGeventide <- function(obj,
                                skim_func = NULL, # see skimr::skim_with
                                summarise_durations = FALSE
) {
  library(dplyr)
  library(skimr)

  cat("Call: ")
  print(obj$call)
  cat("Name:\t\t", obj$name, "\n")
  cat("# sessions:\t", nrow(obj$info), "\n")
  cat("# trials:\t", nrow(obj$trial_data), "\n")

  if (is.null(skim_func)) {
    skim_func <- skim_with(
      base = sfl(missing = n_missing),
      numeric = sfl(mean = mean,
                    std = sd,
                    med = median,
                    mad = mad,
                    hist = function(x) inline_hist(x, 20)),
      append = F
    )
  }

  x <- obj$trial_data %>%
    group_by(condition) %>%
    skim_func(is_correct, is_abort)
  print(x, include_summary = FALSE, width = NULL)

  x <- obj$trial_data %>%
    group_by(block) %>%
    skim_func(is_correct, is_abort)
  print(x, include_summary = FALSE, width = NULL)

  cat("\nFilter by correct trials")
  x <- obj$trial_data %>%
    group_by(block) %>%
    filter(is_correct & (condition != "nogo")) %>%
    skim_func(counter_total_trials, counter_trials_in_block, rt, mt)
  print(x, include_summary = FALSE, width = NULL)

  if (summarise_durations) {
    ind = stringr::str_detect(names(obj$trial_data), "^measured")
    x <- obj$trial_data %>%
      skim_func(names(obj$trial_data)[stringr::str_detect(names(obj$trial_data), "^measured")])
    print(x, include_summary = FALSE, width = NULL)
  }

}

#' @export
read_eventide_single <- function(fname,
                          remove_measured = FALSE,
                          zero_trial_start_time = TRUE
) {
  library(magrittr)
  library(dplyr)
  library(readr)

  # Parse filename

  ## Parse header
  x <- stringr::str_replace_all(readLines(fname, n = 8), ";", "")
  hdr_txt <- stringr::str_split(
    stringr::str_replace_all(x, stringr::fixed(" "), ""),
    ":", simplify = TRUE)
  ind <- hdr_txt[,2] != ""
  hdr_txt <- hdr_txt[ind, ]

  out <- list(
    name = hdr_txt[5,2],
    date = as.POSIXct(paste0(hdr_txt[3,2], " ", hdr_txt[4,2], ":", hdr_txt[4,3]),
                      "%d/%m/%Y %H:%M", tz = "Europe/Paris"),
    version = hdr_txt[1,2]
  )

  ## Read data
  ct <- readr::cols(
    .default = col_double(),
    "Counter Total Trials" = col_integer(),
    "Counter Trials In Block" = col_integer(),
    "Blocked Mode" = col_logical(),
    "Condition Name" = col_character(),
    "Trial Result Str" = col_character(),
    "Is Correct Trial" = col_logical(),
    "Is Incorrect Trial" = col_logical(),
    "Is Abort Trial" = col_logical(),
    "Is Repeat Trial" = col_logical()
  )

  df <- readr::read_csv2(fname, col_names = TRUE, col_types = ct, skip = 8, locale(decimal_mark = ","))

  df %<>% janitor::remove_empty(which = "cols") %>%
    janitor::clean_names()

  # Movement time
  df$mt <- df$tt - df$rt

  # Fill in cue set index for sessions with only one
  if (!("cue_set_index" %in% colnames(df))) {
    df %<>% tibble::add_column(cue_set_index = 0, .after = "block_index")
  }

  cnames <- colnames(df)

  df %<>% rename(cueset = cue_set_index) %>%
    mutate(cueset = factor(cueset, levels = c(0,1), labels = c("old", "new")))

  # convert time to seconds
  msec_to_sec <- function(x, na.rm = FALSE) (x/1000)
  tvars <- c("rt", "rt2", "tt", "mt", "reward_delay")
  tvars <- c(tvars, cnames[stringr::str_ends(cnames, "_duration")],
             cnames[stringr::str_ends(cnames, "_time")])

  df %<>% mutate_at(which(cnames %in% tvars), msec_to_sec) %>%
    tibble::add_column(direction = contra_ipsi_tar(df$tar_x, tolower(out$name)),
                       .after = "cueset") %>%
    tibble::add_column(block = as.factor(df$block_index),
                       .after = "block_index") %>%
    relocate(mt, .after = rt2)

  # set factor levels of trial_result_str, condition_name, direction
  levels(df$block) <- c("con", "mix")

  df$condition_name[df$condition_name=="Go"] = "go"
  df$condition_name[df$condition_name=="Go control"] = "go_con"
  df$condition_name[df$condition_name=="Nogo"] = "nogo"
  df %<>% tibble::add_column(condition = factor(df$condition_name,
                                                levels = c("go_con", "go", "nogo")),
                             .after = "condition_name")
  df$condition_name[df$condition_name=="go_con"] = "go"
  df %<>% rename(gng = condition_name) %>%
    mutate(gng = factor(gng, levels = c("nogo", "go")))

  df %<>% rename(is_correct = is_correct_trial, is_incorrect = is_incorrect_trial,
                 is_abort = is_abort_trial, is_repeat = is_repeat_trial)

  df$trial_result_str[df$trial_result_str=="target touch"] = "target_touch"
  df$trial_result_str[df$trial_result_str=="fixation holding"] = "fixation_holding"
  df$trial_result_str[df$trial_result_str=="Late retouch"] = "err_late_retouch"
  df$trial_result_str[df$trial_result_str=="cue touch aborted"] = "err_cue_touch_aborted"
  df$trial_result_str[df$trial_result_str=="Else-where touch"] = "err_elsewhere_touch"
  df$trial_result_str[df$trial_result_str=="Overall too late"] = "err_overall_too_late"
  df$trial_result_str[df$trial_result_str=="Anticipation"] = "err_anticipation"
  df$trial_result_str[df$trial_result_str=="Early target release"] = "err_early_target_release"

  df %<>% tibble::add_column(event = factor(df$trial_result_str,
                                            levels = c("fixation_holding",
                                                       "target_touch",
                                                       "err_anticipation",
                                                       "err_cue_touch_aborted",
                                                       "err_late_retouch",
                                                       "err_elsewhere_touch",
                                                       "err_early_target_release",
                                                       "err_overall_too_late"
                                            )),
                             .after = "trial_result_str")

  df$cue_duration_programmed <- df$cue_duration
  df$cue_duration <- df$measured_cue_duration

  if (remove_measured) df %<>% select(-starts_with("measured"))

  if (zero_trial_start_time) {
    df %<>% mutate(define_trial_onset_time_absolute = define_trial_onset_time, .before = define_trial_onset_time)
    df %<>% mutate(across(ends_with("_onset_time"), ~purrr::map2_dbl(.x, define_trial_onset_time, ~.x - .y)))
  }

  out$trial_data <- df

  class(out) <- "GNGeventide_single"

  return(out)
}

#' @export
read_eventide_tracker <- function(fname, Fs = 100) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)

  # Parse filename

  ## Parse header
  x <- stringr::str_replace_all(readLines(fname, n = 3), ";", "")
  hdr_txt <- stringr::str_split(
    stringr::str_replace_all(x, stringr::fixed(" "), ""),
    ":", simplify = TRUE, n = 2)
  ind <- hdr_txt[,2] != ""
  hdr_txt <- hdr_txt[ind, ]

  out <- list(
    date = as.POSIXct(hdr_txt[1,2], "%Y.%d.%m%H:%M", tz = "Europe/Paris"),
    version = hdr_txt[2,2]
  )

  ## Read data
  ct <- readr::cols_only(
    `User Field` = col_integer(),
    `Current Event` = col_character(),
    `EventIDE TimeStamp` = col_double(),
    `Gaze CVX` = col_double(), # degrees
    `Gaze CVY` = col_double(),
    #`Gaze X` = col_double(), # pixels
    #`Gaze Y` = col_double(),
    Pressure = col_double(),
    #`Is Touch` = col_logical(),
    X10 = col_logical()
  )

  df <- readr::read_csv2(fname, col_names = TRUE, col_types = ct, skip = 3, locale(decimal_mark = ","))

  df %<>% janitor::remove_empty(which = "cols") %>%
    rename("counter_total_trials" = "User Field",
           "state" = "Current Event",
           "t" = "EventIDE TimeStamp",
           "x" = "Gaze CVX",
           "y" = "Gaze CVY",
           "pressure" = "Pressure")


  df %<>% filter((state == "Fixation") |
                   (state == "Cue") |
                   (state == "Target>Holding Fixation ROI") |
                   (state == "Target>Waiting") |
                   (state == "Target>Target touch") |
                   (state == "Eval") |
                   (state == "Correct>Delay") |
                   (state == "Abort")) %>%
    mutate(x = ifelse(pressure==0, NA, x), y = ifelse(pressure==0, NA, y))

  lev <- c("Fixation", "Cue", "Target>Holding Fixation ROI", "Target>Waiting",
          "Target>Target touch", "Eval", "Correct>Delay", "Abort")
  lev <- lev[lev %in% unique(df$state)]
  df$state <- factor(df$state, levels = lev)

  ## Linearly interpolate to regular grid
  # Create a regular grid
  myseq <- function(from,to,by) tibble(t = seq(from, to, by))
  df2 <- df %>% group_by(counter_total_trials) %>%
    summarise(start = min(t), end = max(t), .groups = "drop") %>%
    group_by(counter_total_trials) %>%
    mutate(t_r = map2(start, end, ~myseq(start, end, 1000/Fs))) %>% # times in msec
    select(-start,-end)

  # Join with original data
  df2 %<>% full_join(df %>% group_by(counter_total_trials) %>% nest(), by = "counter_total_trials")

  # Interpolate
  myapprox <- function(x, y, xout, method = "linear") {
    tibble(r = approx(x, y, xout, ties = min, na.rm = FALSE, method = method)$y)
  }
  df2 %<>% mutate(state = map2(data, t_r, ~myapprox(.x$t, as.integer(.x$state), .y$t, method = "constant")),
                  x = map2(data, t_r, ~myapprox(.x$t, .x$x, .y$t)),
                  y = map2(data, t_r, ~myapprox(.x$t, .x$y, .y$t)),
                  pressure = map2(data, t_r, ~myapprox(.x$t, .x$pressure, .y$t))) %>%
    select(-data) %>%
    unnest(cols = c(t_r, state, x, y, pressure), names_sep = "_") %>%
    rename(t = t_r_t, state = state_r, x = x_r, y = y_r, pressure = pressure_r) %>%
    mutate(pressure = ifelse(is.na(x), 0, pressure), t = t/1000)

  df2$state <- factor(df2$state, labels = lev)

  out$Fs <- Fs
  out$tracker_data <- df2 %>% ungroup()

  class(out) <- "GNGeventide_tracker"

  return(out)
}

contra_ipsi_tar <- function(x, subject) {
  # Contra/Ipsi relative to arm used
  dir = x
  if (subject == "tess") {
    dir[x < 0] = "ipsi"
    dir[x > 0] = "contra"
  } else if (subject == "chanel") {
    dir[x < 0] = "ipsi"
    dir[x > 0] = "contra"
  } else if (subject == "flocky") {
    dir[x > 0] = "ipsi"
    dir[x < 0] = "contra"
  }
  dir = factor(dir, levels = c("ipsi", "contra"))
  return(dir)
}

#' @export
parse_tracker_filenames <- function(basedir = getwd()) {
  fnames <- list.files(path = basedir, pattern =
                         glob2rx(paste0("TrackerLog--ELOTouchTracker--", "*.txt")))

  d <- purrr::map_chr(stringr::str_split(fnames,"--"), 3)
  t <- purrr::map_chr(stringr::str_split(fnames,"--"), 4)
  t <- purrr::map_chr(stringr::str_split(t,".txt"), 1)

  d <- as.POSIXct(paste(d, t),
                    "%Y-%d-%m %H-%M", tz = "Europe/Paris")

  # Sort by ascending experiment date
  ind <- order(d)
  fnames <- fnames[ind]
  d <- d[ind]

  return(list(fnames = fnames, date = d))
}
