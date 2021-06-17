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
                      psth_par = list(h = 0.02, dt = 0.005, pad = 30, shift = 0),
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

  # if (! ("area" %in% names(neuron_info))) {
  #   neuron_info %<>% mutate(area = "target", .after = "target")
  # }

  # Align spikes and behavioral markers to requested event
  events <- shift_events(obj, align = align)

  # In case we want to bin data instead of smoothing
  #spks <- spks_in_window(obj, align = align, t_start = t_start, t_end = t_end, mask_by_quality = TRUE)
  # df <- spks %>%
  #   left_join(events, by = c("id", "session", "counter_total_trials"))

  # PSTHs
  psth_obj <- get_psth(obj, align = align,
                       t_start = t_start, t_end = t_end, h = psth_par$h,
                       dt = psth_par$dt, pad = psth_par$pad, shift = psth_par$shift,
                       min_trial = 0, add_trial_vars = NULL, add_trial_time_vars = NULL)

  df_psth <- psth_obj$psth_df %>%
    left_join(events, by = c("id", "session", "counter_total_trials"))

  df_psth %<>% filter(condition %in% keep_conditions)

  n_trials <- df_psth %>%
    group_by(across(all_of(c("uname",split_factor)))) %>%
    count() %>%
    pivot_wider(names_from = split_factor, values_from = n, names_prefix = "n_")

  tic()
  df_psth %<>%
    # drop some columns to save memory?
    left_join(n_trials, by = "uname") %>%
    mutate(t = list(psth_obj$t)) %>%
    unnest(cols = c(psth, t))
  toc()

  # Get an estimate of the pre-cue firing rate as an estimate of baseline for z-scoring
  pre_cue <- spks_in_window(obj, align = "cue_onset_time", t_start = -0.5, t_end = 0)
  # Make sure we have matching neurons
  pre_cue %<>% semi_join(df_psth, by = c("id", "session", "counter_total_trials"))

  firing_rate_mean <- function(x) {
    n <- purrr::map_dbl(x,length)
    mean(n)
  }
  firing_rate_sd <- function(x) {
    n <- purrr::map_dbl(x,length)
    sd(n)
  }
  pre_cue %<>% group_by(uname) %>% summarise(fr_bl_mean = firing_rate_mean(times2),
                                             fr_bl_sd = firing_rate_sd(times2))

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

  tic()
  if (bootstrap_auc) {
    df_auc <- run_auc_boot(df_psth, split_factor = split_factor, metric = "psth", lev = lev, dir = ">")
    df_auc %<>% tidyr::chop(cols = c(t:p_value))
  } else {
    df_auc <- run_auc(df_psth, split_factor = split_factor, metric = "psth", lev = lev, dir = ">")
    df_auc %<>% tidyr::chop(cols = c(t:ci_hi))
  }
  toc()

  df_auc %<>%
    ungroup() %>%
    left_join(neuron_info %>% select(id, session, uname, rel_depth, area, type), by = "uname")

  df_auc %<>% left_join(n_trials, by = "uname")

  df_auc %<>% relocate(id, session, uname, rel_depth, area, type)

  df_auc %<>% mutate(area_type = interaction(area,type), .after = "type")

  if (!is.null(filename)) {
    df_psth_mean <- df_psth %>%
      select(id, session, counter_total_trials, uname, split_factor, t, psth) %>%
      tidyr::chop(cols = c(t,psth)) %>%
      group_by(uname, get(split_factor), t) %>%
      summarise(psth_sd = list(matrixStats::colSds(do.call(rbind, psth))),
                psth_mean = list(colMeans(do.call(rbind, psth))) ) %>%
      ungroup() %>%
      left_join(pre_cue, by = "uname")

    df_psth_mean %<>% left_join(df_auc %>% select(uname, area, type, area_type), by = "uname")

    saveRDS(df_psth_mean, paste0(str_split(filename, ".rds")[[1]][[1]], "_psth_mean.rds"))
    saveRDS(df_auc, filename)
  }
}

#' @export
batch_find_auc_change <- function(df,
                                  t_min = -0.2,
                                  t_max = 0.4,
                                  min_auc_runlength = 15,  # 15 #(for dt - 0.005)
                                  use_bootstrap_p = TRUE,
                                  use_com = F,
                                  p_thresh = 0.05,
                                  p_adjust_method = "none")
{
  df %<>% unnest(cols = c(t, n_pos, n_neg, auc, auc_boot_median, auc_boot_mean, ci_low,
                          ci_hi, p_value))

  if (use_bootstrap_p) {
    df_t <- df %>%
      group_by(uname) %>%
      group_modify(~find_auc_change_by_p(.x, t_min = t_min, t_max = t_max,
                                         min_runlength = min_auc_runlength,
                                         use_com = use_com,
                                         p_thresh = p_thresh,
                                         p_adjust_method = p_adjust_method))
  } else {
    df_t <- df %>%
      group_by(uname) %>%
      group_modify(~find_auc_change(.x, t_min = t_min, t_max = t_max,
                                    min_runlength = min_auc_runlength,
                                    p_thresh = p_thresh,
                                    p_adjust_method = p_adjust_method))
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

  if (dir == ">") {
    pos_class <- lev[1]
  } else if (dir == "<") {
    pos_class <- lev[2]
  }

  get_auc <- function(x, y, metric, split_factor, pos_class, min_cases) {
    n_pos <- sum(!is.na(x[[metric]]) & (x[[split_factor]] == pos_class))
    n_neg <- sum(!is.na(x[[metric]]) & (x[[split_factor]] != pos_class))

    if ((n_pos>min_cases) & (n_neg>min_cases)) {
      roc <- roc(x[[split_factor]], x[[metric]], levels = lev, direction = dir)
      ci <- ci.auc(roc)
      auc <- roc$auc
      ci_low <- temp$CI.Performance[1]
      ci_hi <- temp$CI.Performance[2]
    } else {
      auc <- NA
      ci_low <- NA
      ci_hi <- NA
    }

    data.frame(n_pos=n_pos, n_neg=n_neg, auc=auc,
               ci_low=ci_low, ci_hi=ci_hi)
  }

  df_auc <- df %>%
    group_by(uname, t) %>%
    group_modify(~get_auc(.x, metric = metric, split_factor = split_factor,
                          pos_class = pos_class, min_cases = min_cases, nboot = nboot),
                 .keep = TRUE)
  # df_auc <- df %>%
  #   group_by(uname, t) %>%
  #   summarise(roc = list(roc(.data[[split_factor]], .data[[metric]], levels = lev, direction = dir))) %>%
  #   mutate(auc = purrr::map_dbl(roc, ~.x$auc)) %>%
  #   mutate(ci = purrr::map(roc, ~ci.auc(.x))) %>%
  #   unnest_wider(col = ci, names_sep = "_") %>%
  #   rename(ci_low = ci_1, auc_median = ci_2, ci_hi = ci_3) %>%
  #   ungroup()
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
find_auc_change_by_p <- function(x, y, t_min, t_max,
                                 p_thresh = 0.05, min_runlength = 10,
                                 use_com = F,
                                 p_adjust_method = "none") {
  ind <- (x$t >= t_min) & (x$t <= t_max)

  min_val <- min(x$auc[ind])
  max_val <- max(x$auc[ind])

  p_vec <- p.adjust(x$p_value[ind], method = p_adjust_method)
  ind_sig <- as.logical(filter_runlengths(p_vec <= p_thresh, min_runlength))

  if (!use_com) {
    # # AUC value at change
    s_change <- x$auc[ind][ind_sig][1]
    # Time at change
    t_change <- x$t[ind][ind_sig][1]
  } else {
    # AUC value at change
    s <- x$auc[ind][ind_sig][1:min_runlength]
    # Time at change
    t <- x$t[ind][ind_sig][1:min_runlength]
    t_change <- sum(s*t)/sum(s)
    s_change <- s[1]
  }

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

#' @export
merge_t_change <- function(df_t1, # auc detections epoch 1
                           df_t2, # auc detections epoch 2
                           t_max1 = 0.3,
                           t_min2 = -0.2
) {

  merge_t_change_for_group <- function(df,
                                       z, # placeholder for group_modify
                                       t_max1 = 0.3,
                                       t_min2 = -0.2
  ) {
    t_change.x = df$t_change.x
    resp.x = df$resp.x
    t_change.y = df$t_change.y
    resp.y = df$resp.y

    if (!is.na(t_change.x) & (t_change.x <= t_max1)) {
      # Change in epoch 1,
      t_change1 <- t_change.x
      resp1 <- resp.x
      # Epoch 2 forced to follow, even if sign changes
      t_change2 <- t_min2 - t_change.x/1e6
      #t_change2 <- t_min2 + t_change.x/1e6
      resp2 <- resp.x
      epoch <- 1
    } else if (is.na(t_change.x) & !is.na(t_change.y) & (t_change.y >= t_min2)) {
      # Pinned to left of epoch 2, but not detected in 1
      t_change2 <- t_change.y
      resp2 <- resp.y
      # Epoch 1 forced to follow
      t_change1 <- t_max1 + t_change.y/1e6
      #t_change1 <- t_max1 - t_change.y/1e6
      resp1 <- resp.y
      epoch <- 2
    } else {
      t_change1 = NA
      resp1 = "none"
      t_change2 = NA
      resp2 = "none"
      epoch <- NA
    }

    data.frame(t_change.x = t_change.x, resp.x = resp.x,
               t_change.y = t_change.y, resp.y = resp.y,
               t_change1 = t_change1, resp1 = resp1,
               t_change2 = t_change2, resp2 = resp2,
               epoch = epoch
    )
  }

  df <- df_t1 %>%
    select(uname, t_change, resp) %>%
    left_join(df_t2 %>%
                select(uname, t_change, resp), by = "uname")

  df2 <- df %>%
    group_by(uname) %>%
    group_modify(~ merge_t_change_for_group(.x, t_max1 = t_max1, t_min2 = t_min2),
                 .keep = T)
}

#' @export
plot_auc_heatmap_pallidum <- function(df_auc,
                                      t_start,
                                      t_end,
                                      pos_first = F,
                                      body_width = NULL,
                                      fix_height = NULL,
                                      left_annotate = TRUE,
                                      id_annotate = TRUE,
                                      label_rows = TRUE,
                                      plot_latency = TRUE,
                                      ignore_epoch = NULL
) {
  library(ComplexHeatmap) # BiocManager::install("ComplexHeatmap")
  library(circlize)

  cm = colormaps()
  if (left_annotate) {
    hm_colors <- c(cm["AreaType"], cm["Id"])
  } else {
    hm_colors <- NULL
  }

  # breaks = c(-.5, -.4, -.3, -.15, 0, .15, .3, .4, .5) + 0.5
  # col_fun = circlize::colorRamp2(breaks,
  #                                khroma::colour("BuRd")(9))
  col_fun = circlize::colorRamp2(seq(from=0, to=1, length.out = 27),
                                 khroma::colour("BuRd")(27))
  # col_fun = circlize::colorRamp2(seq(from=0, to=1, length.out = 13),
  #                                khroma::colour("sunset")(13))

  df_auc %<>% tidyr::unnest(cols = c(t, n_pos, n_neg,
                                     auc, auc_boot_median, auc_boot_mean, ci_low,
                                     ci_hi, p_value))
  df_auc %<>% filter((t>=t_start)&(t<=t_end)) %>%
    tidyr::chop(cols = c(t, n_pos, n_neg,
                         auc, auc_boot_median, auc_boot_mean, ci_low,
                         ci_hi, p_value))

  t_vec <- df_auc[1,]$t[[1]]
  t_vec <- t_vec[(t_vec>=t_start) & (t_vec<=t_end)]
  t_names <- round(t_vec, digits = 3)
  t_names[!t_names %in% round(seq(-.2,.7, by = 0.1), digits = 3)] = NA
  t_names <- as.character(t_names)
  t_names[is.na(t_names)] <- ""

  gen_heatmap <- function(df,
                          colormap,
                          tag,
                          flip = F,
                          column_names = NULL,
                          label_rows = F,
                          ann_col = NULL,
                          id_annotate = F,
                          fix_height = NULL,
                          fix_width = NULL,
                          plot_latency = TRUE,
                          ignore_epoch = NULL
  ) {
    if (nrow(df)==0) {
      return(NULL)
    }
    M <- do.call(rbind, df  %>% pull(auc))
    t <- df %>% pull(t_change)
    if ("t_max" %in% names(df)) {
      if (!all(is.na(df %>% pull(t_max)))) {
        t_max <- df %>% pull(t_max)
      }
    }
    if (label_rows) {
      rownames(M) <- df %>% pull(uname)
    }
    if (!is.null(column_names)) {
      colnames(M) <- column_names
    }

    t_vec <- df[1,]$t[[1]]
    if (!all(is.na(t))) {
      t_ind <- t %>% map_int(~which(abs(t_vec-.x)==min(abs(t_vec-.x))))
      if (exists("t_max")) {
        t_max_ind <- t_max %>% map_int(~which(abs(t_vec-.x)==min(abs(t_vec-.x))))
      }
    } else {
      t_ind <- t
    }

    epoch <- df$epoch
    if (flip) {
      M <- apply(M, 2, rev)
      t <- rev(t)
      t_ind <- rev(t_ind)
      if (exists("t_max")) {
        t_max <- rev(t_max)
        t_max_ind <- rev(t_max_ind)
      }
      epoch <- rev(epoch)
    }

    if (!is.null(ann_col)) {
      # row_ha = rowAnnotation(AreaType = df %>% pull(area_type), show_annotation_name = F,
      #                        col = ann_col, show_legend = F,
      #                        simple_anno_size = unit(4, "mm"))
      if (id_annotate) {
        row_ha = rowAnnotation(AreaType = df %>% pull(area_type), Id = df %>% pull(id),
                               show_annotation_name = F, col = ann_col, show_legend = F,
                               simple_anno_size = unit(3, "mm"))
      } else {
        row_ha = rowAnnotation(AreaType = df %>% pull(area_type),
                               show_annotation_name = F, col = ann_col, show_legend = F,
                               simple_anno_size = unit(3, "mm"))
      }
    } else {
      row_ha = NULL
    }

    hm = Heatmap(M, name = tag, cluster_rows = F, cluster_columns = F,
                 col = colormap, left_annotation = row_ha,
                 height = fix_height,
                 width = fix_width,
                 show_heatmap_legend = F,
                 row_names_gp = gpar(fontsize = 3),
                 column_names_gp = gpar(fontsize = 6),
                 column_names_rot = 0,
                 column_names_centered = T)

    if (plot_latency) {
      if (!is.null(ignore_epoch)) {
        t[epoch == ignore_epoch] = NA
        t_ind[epoch == ignore_epoch] = NA
      }
    } else {
      t <- t*NA
      t_ind <- t_ind*NA
    }
    if (!exists("t_max")) {
      t_max = rep(NA,length(t))
      t_max_ind = t_max
    }
    list(hm = hm, t = t, t_ind = t_ind, t_max = t_max, t_max_ind = t_max_ind)
  }

  #fix_height =  NULL #unit(10, "mm") #
  hm_gpe_hfd_neg <- gen_heatmap(df_auc %>% filter(area_type=="gpe.hfd", resp=="<"),
                                tag = "gpe_hfd_neg",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip = ifelse(pos_first, F, T),
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)
  hm_gpe_hfd_pos <- gen_heatmap(df_auc %>% filter(area_type=="gpe.hfd", resp==">"),
                                tag = "gpe_hfd_pos",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip = ifelse(pos_first, T, F),
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)
  hm_gpe_hfd_ns <- gen_heatmap(df_auc %>% filter(area_type=="gpe.hfd", resp=="none"),
                               tag = "gpe_hfd_ns",
                               column_names = t_names,
                               label_rows = label_rows,
                               colormap = col_fun,
                               ann_col = hm_colors,
                               id_annotate = id_annotate,
                               fix_height = fix_height,
                               fix_width = body_width,
                               flip = F,
                               plot_latency = plot_latency,
                               ignore_epoch = ignore_epoch)

  hm_gpe_hfdp_neg <- gen_heatmap(df_auc %>% filter(area_type=="gpe.hfd-p", resp=="<"),
                                 tag = "gpe_hfdp_neg",
                                 column_names = t_names,
                                 label_rows = label_rows,
                                 colormap = col_fun,
                                 ann_col = hm_colors,
                                 id_annotate = id_annotate,
                                 fix_height = fix_height,
                                 fix_width = body_width,
                                 flip = ifelse(pos_first, F, T),
                                 plot_latency = plot_latency,
                                 ignore_epoch = ignore_epoch)
  hm_gpe_hfdp_pos <- gen_heatmap(df_auc %>% filter(area_type=="gpe.hfd-p", resp==">"),
                                 tag = "gpe_hfdp_pos",
                                 column_names = t_names,
                                 label_rows = label_rows,
                                 colormap = col_fun,
                                 ann_col = hm_colors,
                                 id_annotate = id_annotate,
                                 fix_height = fix_height,
                                 fix_width = body_width,
                                 flip = ifelse(pos_first, T, F),
                                 plot_latency = plot_latency,
                                 ignore_epoch = ignore_epoch)
  hm_gpe_hfdp_ns <- gen_heatmap(df_auc %>% filter(area_type=="gpe.hfd-p", resp=="none"),
                                tag = "gpe_hfdp_ns",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip =  F,
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)

  hm_gpe_lfd_neg <- gen_heatmap(df_auc %>% filter(area_type=="gpe.lfd", resp=="<"),
                                tag = "gpe_lfd_neg",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip = ifelse(pos_first, F, T),
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)
  hm_gpe_lfd_pos <- gen_heatmap(df_auc %>% filter(area_type=="gpe.lfd", resp==">"),
                                tag = "gpe_lfd_pos",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip = ifelse(pos_first, T, F),
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)
  hm_gpe_lfd_ns <- gen_heatmap(df_auc %>% filter(area_type=="gpe.lfd", resp=="none"),
                               tag = "gpe_lfd_ns",
                               column_names = t_names,
                               label_rows = label_rows,
                               colormap = col_fun,
                               ann_col = hm_colors,
                               id_annotate = id_annotate,
                               fix_height = fix_height,
                               fix_width = body_width,
                               flip =  F,
                               plot_latency = plot_latency,
                               ignore_epoch = ignore_epoch)

  hm_gpe_lfdb_neg <- gen_heatmap(df_auc %>% filter(area_type=="gpe.lfd-b", resp=="<"),
                                 tag = "gpe_lfdb_neg",
                                 column_names = t_names,
                                 label_rows = label_rows,
                                 colormap = col_fun,
                                 ann_col = hm_colors,
                                 id_annotate = id_annotate,
                                 fix_height = fix_height,
                                 fix_width = body_width,
                                 flip = ifelse(pos_first, F, T),
                                 plot_latency = plot_latency,
                                 ignore_epoch = ignore_epoch)
  hm_gpe_lfdb_pos <- gen_heatmap(df_auc %>% filter(area_type=="gpe.lfd-b", resp==">"),
                                 tag = "gpe_lfdb_pos",
                                 column_names = t_names,
                                 label_rows = label_rows,
                                 colormap = col_fun,
                                 ann_col = hm_colors,
                                 id_annotate = id_annotate,
                                 fix_height = fix_height,
                                 fix_width = body_width,
                                 flip = ifelse(pos_first, T, F),
                                 plot_latency = plot_latency,
                                 ignore_epoch = ignore_epoch)
  hm_gpe_lfdb_ns <- gen_heatmap(df_auc %>% filter(area_type=="gpe.lfd-b", resp=="none"),
                                tag = "gpe_lfdb_ns",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip =  F,
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)

  hm_gpi_hfd_neg <- gen_heatmap(df_auc %>% filter(area_type=="gpi.hfd", resp=="<"),
                                tag = "gpi_hfd_neg",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip = ifelse(pos_first, F, T),
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)
  hm_gpi_hfd_pos <- gen_heatmap(df_auc %>% filter(area_type=="gpi.hfd", resp==">"),
                                tag = "gpi_hfd_pos",
                                column_names = t_names,
                                label_rows = label_rows,
                                colormap = col_fun,
                                ann_col = hm_colors,
                                id_annotate = id_annotate,
                                fix_height = fix_height,
                                fix_width = body_width,
                                flip = ifelse(pos_first, T, F),
                                plot_latency = plot_latency,
                                ignore_epoch = ignore_epoch)
  hm_gpi_hfd_ns <- gen_heatmap(df_auc %>% filter(area_type=="gpi.hfd", resp=="none"),
                               tag = "gpi_hfd_ns",
                               column_names = t_names,
                               label_rows = label_rows,
                               colormap = col_fun,
                               ann_col = hm_colors,
                               id_annotate = id_annotate,
                               fix_height = fix_height,
                               fix_width = body_width,
                               flip =  F,
                               plot_latency = plot_latency,
                               ignore_epoch = ignore_epoch)

  if (pos_first) {
    hm_list = hm_gpe_hfd_pos[[1]] %v% hm_gpe_hfd_neg[[1]] %v% hm_gpe_hfd_ns[[1]] %v%
      hm_gpe_hfdp_pos[[1]] %v% hm_gpe_hfdp_neg[[1]] %v% hm_gpe_hfdp_ns[[1]] %v%
      hm_gpe_lfd_pos[[1]] %v% hm_gpe_lfd_neg[[1]] %v% hm_gpe_lfd_ns[[1]] %v%
      hm_gpe_lfdb_pos[[1]] %v% hm_gpe_lfdb_neg[[1]] %v% hm_gpe_lfdb_ns[[1]] %v%
      hm_gpi_hfd_pos[[1]] %v% hm_gpi_hfd_neg[[1]] %v% hm_gpi_hfd_ns[[1]]
  } else {
    hm_list = hm_gpe_hfd_neg[[1]] %v% hm_gpe_hfd_pos[[1]] %v% hm_gpe_hfd_ns[[1]] %v%
      hm_gpe_hfdp_neg[[1]] %v% hm_gpe_hfdp_pos[[1]] %v% hm_gpe_hfdp_ns[[1]] %v%
      hm_gpe_lfd_neg[[1]] %v% hm_gpe_lfd_pos[[1]] %v% hm_gpe_lfd_ns[[1]] %v%
      hm_gpe_lfdb_neg[[1]] %v% hm_gpe_lfdb_pos[[1]] %v% hm_gpe_lfdb_ns[[1]] %v%
      hm_gpi_hfd_neg[[1]] %v% hm_gpi_hfd_pos[[1]] %v% hm_gpi_hfd_ns[[1]]
  }

  intra_gap = 0.25
  inter_gap = 2
  gaps = c(rep(intra_gap,2), inter_gap,
           rep(intra_gap,2), inter_gap,
           rep(intra_gap,2), inter_gap,
           rep(intra_gap,2), inter_gap,
           rep(intra_gap,2))
  gaps = unit(gaps, "mm")

  if (left_annotate) {
    type_legend <- Legend(labels = names(hm_colors$AreaType),
                          title = "Type", legend_gp = gpar(fill = hm_colors$AreaType))

    if (id_annotate) {
      id_legend <- Legend(labels = names(hm_colors$Id),
                          title = "Id", legend_gp = gpar(fill = hm_colors$Id))
      auc_legend <- list(type_legend, id_legend, Legend(col_fun = col_fun, title = "AUC"))
    } else {
      auc_legend <- list(type_legend, Legend(col_fun = col_fun, title = "AUC"))
    }

  } else {
    auc_legend <- list()
  }

  hm_all <- ComplexHeatmap::draw(hm_list, gap = gaps, merge_legends = F, heatmap_legend_list = auc_legend, newpage = F)

  types = c("gpe_hfd", "gpe_hfdp", "gpe_lfd", "gpe_lfdb", "gpi_hfd")
  resps = c("neg", "pos", "ns")

  for (i in 1:length(types)) {
    for (j in 1:length(resps)) {
      decorate_heatmap_body(paste0(types[i], '_', resps[j]), {
        ind = which(t_names=="0")
        x = ind/length(t_names)
        grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = .5, lty = 1))
      }, slice = 1)

      decorate_heatmap_body(paste0(types[i], '_', resps[j]), {
        varname <- paste0("hm_", types[i], '_', resps[j])
        x = get(varname)[[3]]
        x = x/length(t_names)
        y = ((length(x)-1):0)/(length(x)) + (1/length(x))/2

        grid.points(x, y, pch = 19, size = unit(0.2, "mm"), default.units = "npc")
      }, slice = 1)

      # For plotting maximum time
      # decorate_heatmap_body(paste0(types[i], '_', resps[j]), {
      #   varname <- paste0("hm_", types[i], '_', resps[j])
      #   x = get(varname)[[5]]
      #   x = x/length(t_names)
      #   y = ((length(x)-1):0)/(length(x)) + (1/length(x))/2
      #
      #   grid.points(x, y, pch = 15, size = unit(0.2, "mm"), default.units = "npc")
      # }, slice = 1)
    }
  }
}

