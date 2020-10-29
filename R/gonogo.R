# Allow wild-card searches, or date restrictions
read_eventide_multi <- function(name,
                                task = "GNG",
                                start_date,
                                end_date,
                                min_trials
                                ) {
  library(dplyr)
  library(readr)
}

read_eventide <- function(fname) {
  library(dplyr)
  library(readr)

  # Parse filename

  ## Parse header
  hdr_txt <- stringr::str_split(
    stringr::str_replace_all(readLines(fname, n = 8), stringr::fixed(" "), ""),
    ":", simplify = TRUE)
  ind <- hdr_txt[,2] != ""
  hdr_txt <- hdr_txt[ind, ]

  out = list(
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

  df <- df %>%
    janitor::remove_empty(which = "cols") %>%
    janitor::clean_names()

  cnames <- colnames(df)

  df$mt <- df$tt - df$rt
  if ("cue_set_index" %in% cnames) df$cue_set_index <- 0

  # convert time to seconds
  msec_to_sec <- function(x, na.rm = FALSE) (x/1000)
  tvars <- c("rt", "rt2", "tt", "mt", "reward_delay")
  tvars <- c(tvars, cnames[stringr::str_ends(cnames, "_duration")],
             cnames[stringr::str_ends(cnames, "_time")])

  df <- df %>%
    mutate_at(which(cnames %in% tvars), msec_to_sec) %>%
    tibble::add_column(direction = contra_ipsi_tar(df$tar_x, out$name),
                       .after = "cue_set_index") %>%
    relocate(mt, .after = rt2)

  out$trial_data <- df

  class(out) <- "GNG_session_data"

  return(out)
}


# read_eventide <- function(fname) {
#   library(dplyr)
#   library(readr)
#
#   # Parse filename
#
#   # Parse header
#   hdr_txt <- stringr::str_split(
#     stringr::str_replace_all(readLines(fname, n = 8), stringr::fixed(" "), ""),
#                                ":", simplify = TRUE)
#   ind <- hdr_txt[,2] != ""
#   hdr_txt <- hdr_txt[ind, ]
#
#   hdr <- tibble::tibble(Name = hdr_txt[5,2],
#                Date = as.POSIXct(paste0(hdr_txt[3,2], " ", hdr_txt[4,2], ":", hdr_txt[4,3]),
#                                  "%d/%m/%Y %H:%M", tz = "Europe/Paris"),
#                Version = hdr_txt[1,2])
#
#   # Read data
#   ct <- readr::cols(
#     .default = col_double(),
#     "Counter Total Trials" = col_integer(),
#     "Counter Trials In Block" = col_integer(),
#     "Blocked Mode" = col_logical(),
#     "Condition Name" = col_character(),
#     "Trial Result Str" = col_character(),
#     "Is Correct Trial" = col_logical(),
#     "Is Incorrect Trial" = col_logical(),
#     "Is Abort Trial" = col_logical(),
#     "Is Repeat Trial" = col_logical()
#   )
#
#   df <- readr::read_csv2(fname, col_names = TRUE, col_types = ct, skip = 8, locale(decimal_mark = ","))
#
#   df <- hdr %>% full_join(df, by = character()) %>%
#     janitor::remove_empty(which = "cols") %>%
#     janitor::clean_names()
#
#   cnames <- colnames(df)
#
#   df$mt <- df$tt - df$rt
#   if ("cue_set_index" %in% cnames) df$cue_set_index <- 0
#
#   # convert time to seconds
#   msec_to_sec <- function(x, na.rm = FALSE) (x/1000)
#   tvars <- c("rt", "rt2", "tt", "mt", "reward_delay")
#   tvars <- c(tvars, cnames[stringr::str_ends(cnames, "_duration")],
#              cnames[stringr::str_ends(cnames, "_time")])
#
#   df <- df %>%
#     mutate_at(which(cnames %in% tvars), msec_to_sec)
#
#   df$direction <- contra_ipsi_tar(df$tar_x, df$name[1])
#
#   # reorder columns
#
#   return(df)
# }


contra_ipsi_tar <- function(x, subject) {
  # Contra/Ipsi relative to arm used
  dir = x
  if (subject == "Tess") {
    dir[x < 0] = "ipsi"
    dir[x > 0] = "contra"
  } else if (subject == "Chanel") {
    dir[x < 0] = "ipsi"
    dir[x > 0] = "contra"
  } else if (subject == "Flocky") {
    dir[x > 0] = "ipsi"
    dair[x < 0] = "contra"
  }
  dir = as.factor(dir)
  return(dir)
}
