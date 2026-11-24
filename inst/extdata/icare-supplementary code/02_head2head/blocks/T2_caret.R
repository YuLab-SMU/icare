fn <- caretSBF
fn$filter <- function(score, x, y) rank(score, ties.method = "first") <= K
fit <- sbf(d_tr[, setdiff(names(d_tr), "y")], d_tr$y, method = "glm",
           preProcess = c("center", "scale"),
           trControl = trainControl(method = "none"),
           sbfControl = sbfControl(functions = fn, method = "cv", number = 5))
prob <- predict(fit, d_va[, setdiff(names(d_va), "y")])[, "Case"]
sel <- fit$optVariables
