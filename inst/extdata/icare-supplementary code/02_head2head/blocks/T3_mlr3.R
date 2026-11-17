m <- floor(sqrt(ncol(d_tr) - 1))
task <- as_task_classif(rbind(d_tr, d_va), target = "y", positive = "Case")
tr_ids <- seq_len(nrow(d_tr))
va_ids <- nrow(d_tr) + seq_len(nrow(d_va))
lrns <- list(
  glmnet = lrn("classif.glmnet", alpha = 0.5, s = 0.01, predict_type = "prob"),
  ranger = lrn("classif.ranger", num.trees = 500, mtry = m, min.node.size = 1,
               predict_type = "prob"),
  rpart  = lrn("classif.rpart", predict_type = "prob"))
gls <- lapply(lrns, function(l) as_learner(po("scale") %>>% l))
set.seed(s)
cv_auc <- sapply(gls, function(g)
  unname(resample(task$clone()$filter(tr_ids), g, rsmp("cv", folds = 5))$aggregate(msr("classif.auc"))))
probs <- lapply(gls, function(g) { g$train(task, tr_ids); g$predict(task, va_ids)$prob[, "Case"] })
