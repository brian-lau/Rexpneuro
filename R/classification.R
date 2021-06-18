#' @export
classify_dirXgng <- function(datafile,
                             align = "cue_onset_time", # must exist in obj$trial_data
                             t_start = -0.5,
                             t_end = 0,
                             min_n = 10,
                             binwidth = t_end - t_start)
{
  load(datafile)
  obj %<>% drop_abort_trials()
  obj %<>% drop_incorrect_trials()

  # Extract neuron information
  neuron_info <- obj$info %>%
    unnest(cols = neuron_info, names_repair = 'universal')

  # Align spikes and behavioral markers to requested event
  events <- shift_events(obj, align = align)

  # Bin data
  # Should allow multiple aligns with different t_starts and ends
  # add variable epoch
  spks <- count_spks_in_window(obj, align = align, t_start = t_start, t_end = t_end,
                               binwidth = binwidth, mask_by_quality = TRUE)
  # Bind in some useful information
  spks2 <- spks %>%
    left_join(neuron_info %>% select(id, session, name, uname, area, type)) %>%
    left_join(events %>% select(id, session, counter_total_trials, direction, gng, condition)) %>%
    mutate(area_type = interaction(area, type)) %>%
    rename(spk_count = n) %>%
    select(-area, -type, -id, -session, -name) %>%
    mutate(cond = interaction(direction, gng)) %>%
    select(-direction, -gng) %>%
    filter(condition != "go_con") %>%
    select(-condition, -counter_total_trials)

  spks3 <- spks2 %>%
    unnest(cols = c(t,spk_count))

  # Move this to loop below
  # Restrict to type
  spks4 <- spks3 %>%
    filter(area_type == "gpe.lfd") %>%
    select(-area_type)

  df_count <- spks4 %>%
    group_by(uname, cond, t) %>%
    count %>%
    ungroup() %>%
    pivot_wider(names_from = cond, values_from = n) %>%
    group_by(uname) %>%
    slice_head(n=1) %>%
    mutate(n = ipsi.nogo + contra.nogo + ipsi.go + contra.go) %>%
    mutate(good = (ipsi.nogo >= min_n) + (contra.nogo >= min_n) +
             (ipsi.go >= min_n) + (contra.go >= min_n)) %>%
    arrange(n)

  keep_names <- df_count %>% filter(good == 4) %>% pull(uname)

  spks4 %<>%
    filter(uname %in% keep_names)

  tvec <- unique(spks4$t)

  res = list()
  for (i in 1:length(tvec)) {
    spks5 <- spks4 %>%
      filter(t == tvec[i]) %>%
      select(-t)

    # Resample trials to be the same
    spks_pseudo <- spks5 %>% group_by(uname, cond) %>% slice_sample(n = min_n, replace = F) %>% ungroup

    dat <- spks_pseudo %>%
      group_by(uname, cond) %>%
      mutate(trial = 1:n()) %>%
      ungroup %>%
      pivot_wider(names_from = uname, values_from = spk_count) %>%
      select(-trial)

    #set.seed(1234)
    #repeats = 1
    ncv_splits = nested_cv(dat, outside = vfold_cv(v = 5, strata = "cond"),
                           inside = vfold_cv(v = 5, strata = "cond"))

    tic()
    res[[i]] = ncv_fit(ncv_splits, model = "rf")
    toc()

    # Predict at all other times
  }
  return(res)

}

#' @export
ncv_fit <- function(ncv_splits, model = "rf", ...) {
  res = list()
  for (i in 1:nrow(ncv_splits)) {
    res[[i]] <- tune_mod_multiclass(ncv_splits[i,], model = model, ...)
  }
  return(res)
}

#' @export
tune_mod_multiclass <- function(ncv_split, model = "rf", grid_size = 10, verbose = F) {
  metrics <- metric_set(roc_auc, accuracy)
  best_metric <- "roc_auc"

  ## Data
  # Extract outer fold
  split_obj <- ncv_split$splits[[1]]
  df_train <- training(split_obj)

  # Extract inner fold(s)
  cv_folds <- ncv_split$inner_resamples[[1]]

  ## Setup and tune model
  recette <- recipe(cond ~ ., data = df_train) #%>%
  #step_nzv(all_predictors()) %>%
  #step_normalize(all_predictors())

  if (model == "glm") {
    mod <- multinom_reg(
      penalty = tune::tune(),
      mixture = tune::tune()) %>%
      set_mode("classification") %>%
      set_engine("glmnet")

    grid <- grid_max_entropy(
      penalty(),
      mixture(),
      size = grid_size)
  } else if (model == "rf") {
    mod <- rand_forest(mtry = tune::tune(),
                       trees = 500,
                       min_n = tune::tune()) %>%
      set_engine("ranger", importance = "permutation",
                 num.threads = 10) %>%
      set_mode("classification")

    grid <- grid_max_entropy(
      mtry(range = c(1, 4)),
      min_n(range = c(2, 15)),
      size = grid_size)
  } else if (model == "svm") {
    mod <- svm_poly(cost = tune(),
                    degree = tune()) %>%
      set_mode("classification") %>%
      set_engine("kernlab")

    grid <- grid_max_entropy(
      cost(),
      degree(),
      size = grid_size)
  } else if (model == "svm-rbf") {
    mod <- svm_rbf(cost = tune(),
                   rbf_sigma = tune()) %>%
      set_mode("classification") %>%
      set_engine("kernlab")

    grid <- grid_max_entropy(
      cost(),
      rbf_sigma(),
      size = grid_size)
  } else if (model == "xgb") {
    mod <- boost_tree(
      trees = 1000,
      tree_depth = tune(), min_n = 10,
      loss_reduction = 0,                     ## first three: model complexity
      sample_size = 0.5, mtry = tune(),         ## randomness
      learn_rate = tune(),                         ## step size
    ) %>%
      set_engine("xgboost") %>%
      set_mode("classification")

    grid <- grid_latin_hypercube(
      tree_depth(),
      #min_n(),
      #loss_reduction(),
      #sample_size = sample_prop(),
      finalize(mtry(c(2,16)), df_train),
      learn_rate(),
      #over_ratio(range = c(1, 1.75)),
      size = grid_size
    )
  }

  wkflw <- workflow() %>%
    add_recipe(recette) %>%
    add_model(mod)

  tune_fits <- tune_grid(
    wkflw,
    resamples = cv_folds,
    grid = grid,
    metrics = metrics,
    control = control_grid(verbose = verbose)
  )

  # Finalize workflow
  final_wkflw <- wkflw %>%
    finalize_workflow(select_best(tune_fits, metric = best_metric))
  # Fit to training data
  fit_train <- final_wkflw %>%
    fit(data = juice(prep(recette)))
  # Fit to outer fold, evaluating on test
  fit_test <- final_wkflw %>%
    last_fit(split = split_obj, metrics = metrics)

  return(list(tune_res = tune_fits, fit_train = fit_train, fit_test = fit_test))

}
