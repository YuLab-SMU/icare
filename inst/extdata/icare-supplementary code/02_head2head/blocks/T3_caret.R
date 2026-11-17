m <- floor(sqrt(ncol(d_tr) - 1))
ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                     summaryFunction = twoClassSummary)
grids <- list(glmnet = data.frame(alpha = 0.5, lambda = 0.01),
              ranger = data.frame(mtry = m, splitrule = "gini", min.node.size = 1),
              rpart  = data.frame(cp = 0.01))
fits <- lapply(names(grids), function(a) {
  set.seed(s)
  train(y ~ ., data = d_tr, method = a, metric = "ROC", trControl = ctrl,
        tuneGrid = grids[[a]], preProcess = c("center", "scale"))
})
names(fits) <- names(grids)
cv_auc <- sapply(fits, function(f) f$results$ROC)
probs <- lapply(fits, function(f) predict(f, d_va, type = "prob")[, "Case"])
