fn <- ldaSBF
fn$filter <- function(score, x, y) rank(score, ties.method = "first") <= K
ppi <- FitApplyPreprocess(d_tr, d_va, group_col = "y", method = c("center", "scale"))
fs <- FeatureSelectSBF(CreateModelObject(data = ppi$training, group_col = "y"),
                       sbf_func = fn, method = "cv", number = 5, repeats = 1, seed = s)
mm <- CreateModelObject(data = rbind(ppi$training, ppi$testing), group_col = "y")
mm@split.data <- list(training = ppi$training, testing = ppi$testing)
mm <- ApplyFeatureSelection(mm, fs$opt_vars)
mm@split.scale.data <- mm@split.data
mm@filtered.set <- mm@split.data
mm <- ModelTrainAnalysis(mm, methods = "glm",
                         control = list(method = "cv", number = 5, repeats = 1),
                         tuneLength = 1, preProcess = NULL, loocv_threshold = 10,
                         metric_selection = "ROC", imbalance_handling = "none", seed = s)
prob <- predict(mm@train.models$glm, mm@filtered.set$testing, type = "prob")[, "Case"]
sel <- fs$opt_vars
