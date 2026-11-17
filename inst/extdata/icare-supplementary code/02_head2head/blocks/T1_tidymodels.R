rec <- recipe(y ~ ., data = d_tr) |>
  step_center(all_predictors()) |>
  step_scale(all_predictors())
wf <- workflow() |>
  add_model(logistic_reg() |> set_engine("glm") |> set_mode("classification")) |>
  add_recipe(rec)
fit <- fit(wf, data = d_tr)
prob <- predict(fit, d_va, type = "prob")$.pred_Case
