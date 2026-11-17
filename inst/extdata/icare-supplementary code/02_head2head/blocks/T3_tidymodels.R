m <- floor(sqrt(ncol(d_tr) - 1))
rec <- recipe(y ~ ., data = d_tr) |> step_normalize(all_predictors())
specs <- list(
  glmnet = logistic_reg(penalty = 0.01, mixture = 0.5) |> set_engine("glmnet"),
  ranger = rand_forest(mtry = m, trees = 500, min_n = 1) |>
    set_engine("ranger") |> set_mode("classification"),
  rpart  = decision_tree(cost_complexity = 0.01) |>
    set_engine("rpart") |> set_mode("classification"))
wfs <- lapply(specs, function(sp) workflow() |> add_recipe(rec) |> add_model(sp))
set.seed(s)
folds <- vfold_cv(d_tr, v = 5, strata = y)
cv_auc <- sapply(wfs, function(w)
  collect_metrics(fit_resamples(w, folds, metrics = metric_set(roc_auc)))$mean)
fits <- lapply(wfs, fit, data = d_tr)
probs <- lapply(fits, function(f) predict(f, d_va, type = "prob")$.pred_Case)
