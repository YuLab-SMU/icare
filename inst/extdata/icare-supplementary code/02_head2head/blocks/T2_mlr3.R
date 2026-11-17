task <- as_task_classif(rbind(d_tr, d_va), target = "y", positive = "Case")
gl <- as_learner(po("scale") %>>% po("filter", filter = flt("anova"), filter.nfeat = K) %>>%
                   lrn("classif.log_reg", predict_type = "prob"))
gl$train(task, row_ids = seq_len(nrow(d_tr)))
prob <- gl$predict(task, row_ids = nrow(d_tr) + seq_len(nrow(d_va)))$prob[, "Case"]
sel <- gl$model$classif.log_reg$train_task$feature_names
