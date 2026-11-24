rec <- recipe(y ~ ., data = d_tr) |>
  step_normalize(all_predictors()) |>
  step_select_roc(all_predictors(), outcome = "y", top_p = K)
wf <- workflow() |>
  add_model(logistic_reg() |> set_engine("glm") |> set_mode("classification")) |>
  add_recipe(rec)
fit <- fit(wf, data = d_tr)
prob <- predict(fit, d_va, type = "prob")$.pred_Case
sel <- setdiff(names(bake(extract_recipe(fit), new_data = d_va)), "y")
