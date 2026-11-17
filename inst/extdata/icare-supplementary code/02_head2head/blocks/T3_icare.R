m <- floor(sqrt(ncol(d_tr) - 1))
mm <- CreateModelObject(data = rbind(d_tr, d_va), group_col = "y")
mm@split.data <- list(training = d_tr, testing = d_va)
mm@split.scale.data <- mm@split.data
mm@filtered.set <- mm@split.data
mm <- ModelTrainAnalysis(mm, methods = c("glmnet", "ranger", "rpart"),
        control = list(method = "cv", number = 5, repeats = 1),
        tune_grids = list(glmnet = data.frame(alpha = 0.5, lambda = 0.01),
                          ranger = data.frame(mtry = m, splitrule = "gini", min.node.size = 1),
                          rpart  = data.frame(cp = 0.01)),
        preProcess = c("center", "scale"), loocv_threshold = 10,
        metric_selection = "ROC", imbalance_handling = "none", seed = s)
cv_auc <- setNames(mm@process.info$cv_metrics$CV_Mean, mm@process.info$cv_metrics$Model)
probs <- lapply(mm@train.models, function(f)
  predict(f, mm@split.scale.data$testing, type = "prob")[, "Case"])
