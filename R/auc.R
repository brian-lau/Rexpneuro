#' @export
batch_auc <- function(datafile = "~/ownCloud/ForFarah/pallidum_GNG.Rdata",
                      comparison = "gng",
                      split_factor = "condition", # condition, direction
                      lev = c("go", "nogo"), # c("go", "nogo")
                      keep_conditions = c("nogo", "go"),
                      align = "target_onset_time",
                      drop_after_event = "liftoff_onset_time",
                      t_start = -0.2,
                      t_end = 0.4,
                      psth_par = list(h = 0.02, dt = 0.005, pad = 30),
                      min_trials_per_cond = 5,
                      bootstrap_auc = TRUE,
                      filename = NULL
) {
  library(dplyr)
  library(magrittr)
  library(tictoc)

  load(datafile)

  obj %<>% drop_abort_trials()
  obj %<>% drop_incorrect_trials()

  # Extract neuron information
  neuron_info <- obj$info %>%
    unnest(cols = neuron_info, names_repair = 'universal')

  # Align spikes and behavioral markers to requested event
  events <- shift_events(obj, align = align)

  # In case we want to bin data instead of smoothing
  #spks <- spks_in_window(obj, align = align, t_start = t_start, t_end = t_end, mask_by_quality = TRUE)
  # df <- spks %>%
  #   left_join(events, by = c("id", "session", "counter_total_trials"))

  # PSTHs
  psth_obj <- get_psth(obj, align = align,
                       t_start = t_start, t_end = t_end, h = psth_par$h, dt = psth_par$dt, pad = psth_par$pad,
                       min_trial = 0, add_trial_vars = NULL, add_trial_time_vars = NULL)

  df_psth <- psth_obj$psth_df %>%
    left_join(events, by = c("id", "session", "counter_total_trials"))

  df_psth %<>% filter(condition %in% keep_conditions)

  if (comparison == "block") {
    n_trials <- df_psth %>%
      group_by(uname, condition) %>%
      count() %>%
      pivot_wider(names_from = condition, values_from = n)
    exclude <- n_trials %>%
      mutate(good = (go>min_trials_per_cond) & (go_con>min_trials_per_cond)) %>%
      filter(good %in% c(FALSE,NA)) %>%
      pull(uname)
  } else if (comparison == "gng") {
    n_trials <- df_psth %>%
      group_by(uname, condition) %>%
      count() %>%
      pivot_wider(names_from = condition, values_from = n)
    exclude <- n_trials %>%
      mutate(good = (go>min_trials_per_cond) & (nogo>min_trials_per_cond)) %>%
      filter(good %in% c(FALSE,NA)) %>%
      pull(uname)
  } else if (comparison == "direction") {
    n_trials <- df_psth %>%
      group_by(uname, direction) %>%
      count() %>%
      pivot_wider(names_from = direction, values_from = n)
    exclude <- n_trials %>%
      mutate(good = (ipsi>min_trials_per_cond) & (contra>min_trials_per_cond)) %>%
      filter(good %in% c(FALSE,NA)) %>%
      pull(uname)
  }

  if (exists("exclude")) {
    df_psth %<>% filter(!uname %in% exclude)
  }

  tic()
  df_psth %<>%
    # drop some columns to save memory
    mutate(t = list(psth_obj$t)) %>%
    unnest(cols = c(psth, t))
  toc()

  ## Drop time samples following subsequent event
  if (!is.null(drop_after_event)) {
    if (comparison == "block") {
      df_psth %<>%
        mutate(psth = ifelse(t >= .data[[drop_after_event]], NA, psth))
    } else if (comparison == "gng") {
      df_psth %<>%
        mutate(psth = ifelse(condition == "go", ifelse(t >= .data[[drop_after_event]], NA, psth), psth))
    } else if (comparison == "direction") {
      df_psth %<>%
        mutate(psth = ifelse(t >= .data[[drop_after_event]], NA, psth))
    }
  }

  if (bootstrap_auc) {
    tic()
    df_auc <- run_auc_boot(df_psth, split_factor = split_factor, metric = "psth", lev = lev, dir = ">")
    toc()

    df_auc %<>%
      left_join(neuron_info %>% select(id, session, uname, rel_depth, area, type), by = "uname") %>%
      ungroup()
  } else {
    tic()
    df_auc <- run_auc(df_psth, split_factor = split_factor, metric = "psth", lev = lev, dir = ">")
    toc()
    df_auc %<>%
      left_join(neuron_info %>% select(id, session, uname, rel_depth, area, type), by = "uname") %>%
      select(-roc)
  }

  df_auc %<>% relocate(id, session)

  df_auc %<>% mutate(area_type = interaction(area,type), .after = "type")

  if (!is.null(filename)) {
    df_psth_mean <- df_psth %>%
      select(id, session, counter_total_trials, uname, split_factor, t, psth) %>%
      chop(cols = c(t,psth)) %>%
      group_by(uname, get(split_factor), t) %>%
      summarise(psth_sd = list(matrixStats::colSds(do.call(rbind, psth))),
                psth_mean = list(colMeans(do.call(rbind, psth))) ) %>%
      ungroup()

    df_psth_mean %<>% left_join(df_auc %>% select(uname, area, type, area_type), by = "uname")

    save(df_psth_mean, file = paste0(str_split(filename, ".Rdata")[[1]][[1]], "_psth_mean.Rdata"))
    save(df_auc, file = filename)
  }
}

#' @export
run_auc <- function(df,
                    split_factor = "condition",
                    metric = "psth",
                    lev = c("go", "nogo"),
                    dir = ">"
) {
  library(pROC)

  df_auc <- df %>%
    group_by(uname, t) %>%
    summarise(roc = list(roc(.data[[split_factor]], .data[[metric]], levels = lev, direction = dir))) %>%
    mutate(auc = purrr::map_dbl(roc, ~.x$auc)) %>%
    mutate(ci = purrr::map(roc, ~ci.auc(.x))) %>%
    unnest_wider(col = ci, names_sep = "_") %>%
    rename(ci_low = ci_1, auc_median = ci_2, ci_hi = ci_3) %>%
    ungroup()
}

#' @export
percentile_p <- function(x, h0) {
  half.pval <- mean(x > h0) + 0.5*mean(x == 0.5)
  2*min(c(half.pval, 1 - half.pval))
}

#' @export
run_auc_boot <- function(df,
                          split_factor = "condition",
                          metric = "psth",
                          lev = c("go", "nogo"),
                          dir = ">",
                          nboot = 1000,
                          min_cases = 1
) {
  library(fbroc)

  if (dir == ">") {
    pos_class <- lev[1]
  } else if (dir == "<") {
    pos_class <- lev[2]
  }

  get_auc <- function(x, y, metric, split_factor, pos_class, min_cases, nboot) {
    n_pos <- sum(!is.na(x[[metric]]) & (x[[split_factor]] == pos_class))
    n_neg <- sum(!is.na(x[[metric]]) & (x[[split_factor]] != pos_class))

    if ((n_pos>min_cases) & (n_neg>min_cases)) {
      roc <- boot.roc(x[[metric]], x[[split_factor]]==pos_class, n.boot = nboot)
      temp <- perf(roc, metric = "auc")
      auc <- temp$Observed.Performance
      auc_boot_median <- median(temp$boot.results)
      auc_boot_mean <- mean(temp$boot.results)
      ci_low <- temp$CI.Performance[1]
      ci_hi <- temp$CI.Performance[2]
      p_value <- percentile_p(temp$boot.results, 0.5)
    } else {
      auc <- NA
      auc_boot_median <- NA
      auc_boot_mean <- NA
      ci_low <- NA
      ci_hi <- NA
      p_value <- NA
    }

    data.frame(n_pos=n_pos, n_neg=n_neg, auc=auc,
               auc_boot_median=auc_boot_median, auc_boot_mean=auc_boot_mean,
               ci_low=ci_low, ci_hi=ci_hi, p_value=p_value)
  }

  df_auc <- df %>%
    group_by(uname, t) %>%
    group_modify(~get_auc(.x, metric = metric, split_factor = split_factor,
                          pos_class = pos_class, min_cases = min_cases, nboot = nboot),
                 .keep = TRUE)

}


#' @export
find_auc_change <- function(x, y, min_runlength = 10) {
  ind <- (x$t > 0.01) & (x$t <= 0.5)

  min_val <- min(x$auc[ind])
  max_val <- max(x$auc[ind])

  ind_neg <- as.logical(filter_runlengths(x$ci_hi[ind]<0.5, min_runlength))
  ind_pos <- as.logical(filter_runlengths(x$ci_low[ind]>0.5, min_runlength))

  t_neg = x$t[ind][ind_neg][1]
  t_pos = x$t[ind][ind_pos][1]

  if (is.na(t_neg) & is.na(t_pos)) {
    resp = "none"
    t_change = NA
  } else if (is.na(t_neg)) {
    resp = ">"
    t_change = t_pos
  } else if (is.na(t_pos)) {
    resp = "<"
    t_change = t_neg
  } else if (t_neg > t_pos) {
    resp = ">"
    t_change = t_pos
  } else if (t_neg < t_pos) {
    resp = "<"
    t_change = t_neg
  } else {
    resp = "oops"
    t_change = NA
  }

  data.frame(min_val = min_val, max_val = max_val,
             t_neg = t_neg, t_pos = t_pos, t_change = t_change, resp = resp)
}

#' @export
find_auc_change_by_p <- function(x, y, p_thresh = 0.05, min_runlength = 10) {
  ind <- (x$t > 0.01) & (x$t <= 0.5)

  min_val <- min(x$auc[ind])
  max_val <- max(x$auc[ind])

  ind_sig <- as.logical(filter_runlengths(x$p_value[ind] < p_thresh, min_runlength))

  s_change <- x$auc[ind][ind_sig][1]
  t_change <- x$t[ind][ind_sig][1]

  if (is.na(t_change)) {
    resp = "none"
  } else if (s_change > 0.5) {
    resp = ">"
  } else if (s_change < 0.5) {
    resp = "<"
  } else {
    resp = "oops"
  }

  data.frame(min_val = min_val, max_val = max_val,
             t_change = t_change, resp = resp)
}
