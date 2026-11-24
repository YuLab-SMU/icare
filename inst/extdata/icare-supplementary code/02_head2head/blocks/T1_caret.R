pp <- preProcess(d_tr[, setdiff(names(d_tr), "y")], method = c("center", "scale"))
tr_s <- predict(pp, d_tr)
va_s <- predict(pp, d_va)
fit <- train(y ~ ., data = tr_s, method = "glm", family = binomial,
             trControl = trainControl(method = "none"))
prob <- predict(fit, va_s, type = "prob")[, "Case"]
