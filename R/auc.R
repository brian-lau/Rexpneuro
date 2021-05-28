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
                         nboot = 1000
) {
  library(fbroc)

  if (dir == ">") {
    pos_class <- lev[1]
  } else if (dir == "<") {
    pos_class <- lev[2]
  }

  df_auc <- df %>%
    group_by(uname, t) %>%
    summarise(roc = list(boot.roc(.data[[metric]], .data[[split_factor]]==pos_class, n.boot = nboot))) %>%
    mutate(auc = purrr::map_dbl(roc, ~.x$auc)) %>%
    mutate(ci = purrr::map(roc, ~perf(.x, metric = "auc"))) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(auc_boot_median = median(ci$boot.results),
           auc_boot_mean = mean(ci$boot.results),
           p_value = percentile_p(ci$boot.results, 0.5),
           ci_low = ci$CI.Performance[1],
           ci_hi = ci$CI.Performance[2])
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
