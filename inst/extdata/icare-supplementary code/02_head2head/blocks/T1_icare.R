ppi <- FitApplyPreprocess(d_tr, d_va, group_col = "y", method = c("center", "scale"))
mm <- CreateModelObject(data = rbind(ppi$training, ppi$testing), group_col = "y")
mm@split.data <- list(training = ppi$training, testing = ppi$testing)
mm@split.scale.data <- mm@split.data
mm@filtered.set <- mm@split.data
mm <- ModelTrainAnalysis(mm, methods = "glm",
                         control = list(method = "cv", number = 5, repeats = 1),
                         tuneLength = 1, preProcess = NULL, loocv_threshold = 10,
                         metric_selection = "ROC", imbalance_handling = "none", seed = s)
prob <- predict(mm@train.models$glm, mm@filtered.set$testing, type = "prob")[, "Case"]
