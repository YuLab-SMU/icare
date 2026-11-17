# ==============================================================================
# Undiagnosed dysglycemia - Step 3 of 3
# Train with iCARE and validate independently on the other NHANES cycle
# ------------------------------------------------------------------------------
# Scenario      : HbA1c >= 6.5% among adults self-reporting no diabetes (DIQ010 = 2)
# Clinical action: Who needs an HbA1c test / referral for confirmatory diagnosis
#
# Workflow follows the official iCARE Module 2 Quick Start (yulab-smu.top/icare):
#   model object -> stratified split -> training-only preprocessing ->
#   ModelTrainAnalysis -> SelectBestModel
# Then, following the official Advanced chapter section 13, the fitted model and
# preprocessing are packaged with ModelDeployment and replayed on the external cycle.
#
# Two deliberate departures from the official default (consistent with this
# project's leakage-control focus):
#   1) Algorithm selection uses TRAINING-CV ROC ranking, not the default
#      test-set AUC maximum (the latter optimistically biases the reported value).
#   2) The official quick start runs feature selection on all data before the
#      split; this script performs no cross-split feature selection. If screening
#      is needed, it must be done inside the training partition only.
#
# Input  : output/cohort.rds (produced by 02_cohort.R)
# Output : output/results.csv, cv_rankings.csv, clinical_utility.csv,
#          roc_external.png, model.rds
# ==============================================================================

suppressMessages({
  library(icare)
  library(caret)
  library(pROC)
})

# ---- 0. Locate this folder ---------------------------------------------------
script_dir <- function() {
  args <- commandArgs(FALSE)
  hit  <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (nzchar(p) && file.exists(p)) return(dirname(normalizePath(p)))
  }
  stop("Source this script from RStudio, or run it with Rscript and its full path.")
}
SCRIPT_DIR <- script_dir()

OUTPUT_DIR <- file.path(SCRIPT_DIR, "output")

# ---- 1. Constants -------------------------------------------------------------
SEED        <- 2026
TRAIN_FRAC  <- 0.70
CV_FOLDS    <- 5
METHODS     <- c("glm", "glmnet", "rf", "nb")
NEG         <- "Negative"
POS         <- "Positive"
GROUP       <- "group"
TARGET_SENS <- c(0.90, 0.95)   # target sensitivities for the rule-out rule

# ---- 2. Helpers ---------------------------------------------------------------
auc_with_ci <- function(y, p) {
  stopifnot(mean(p[y == POS]) > mean(p[y == NEG]))
  r  <- pROC::roc(y, p, levels = c(NEG, POS), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  c(auc = as.numeric(pROC::auc(r)), low = ci[1], high = ci[3])
}

performance_at <- function(y, p, threshold) {
  pred <- ifelse(p >= threshold, POS, NEG)
  tp <- sum(pred == POS & y == POS); tn <- sum(pred == NEG & y == NEG)
  fp <- sum(pred == POS & y == NEG); fn <- sum(pred == NEG & y == POS)
  c(threshold   = threshold,
    sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp), npv         = tn / (tn + fn))
}

# Highest threshold on the internal held-out split still reaching target sensitivity
rule_out_threshold <- function(y, p, target) {
  candidates <- sort(unique(p))
  sens <- vapply(candidates, function(t) mean(p[y == POS] >= t), numeric(1))
  ok <- which(sens >= target)
  if (!length(ok)) return(NA_real_)
  candidates[max(ok)]
}

# ---- 3. Load cohort -----------------------------------------------------------
cohort   <- readRDS(file.path(OUTPUT_DIR, "cohort.rds"))
dev      <- cohort$dev
ext      <- cohort$ext
FEATURES <- cohort$features

# ---- 4. Stratified split (official Module 2 Quick Start style) ---------------
# Only the development cycle is split; the external cycle is untouched until final evaluation
set.seed(SEED)
idx      <- caret::createDataPartition(dev[[GROUP]], p = TRAIN_FRAC, list = FALSE)
train    <- dev[idx, ]
internal <- dev[-idx, ]
cat(sprintf("train %d (positives %d) | internal held-out %d | external %d (positives %d)\n",
            nrow(train), sum(train[[GROUP]] == POS), nrow(internal),
            nrow(ext), sum(ext[[GROUP]] == POS)))

# ---- 5. Training-only collinearity screen (|r| > 0.95) -----------------------
cm   <- cor(train[, FEATURES], use = "pairwise.complete.obs")
drop <- caret::findCorrelation(cm, cutoff = 0.95, exact = TRUE)
FEATURES <- FEATURES[-drop]
cat(sprintf("Collinearity screen (training only): dropped %d, remaining %d features\n",
            length(drop), length(FEATURES)))

# ---- 6. Create the iCARE modelling object ------------------------------------
# The official workflow starts from a Module 1 Stat object via ConvertObject();
# with no Stat object, CreateModelObject(data, group_col) is the equivalent entry point.
model_frame <- function(d) d[, c(FEATURES, GROUP), drop = FALSE]

model_obj <- icare::CreateModelObject(
  data      = rbind(model_frame(train), model_frame(internal)),
  group_col = GROUP
)
model_obj@split.data <- list(
  training = model_frame(train),
  testing  = model_frame(internal)
)
model_obj@filtered.set <- model_obj@split.data

# ---- 7. Train: center/scale fitted on the training partition by iCARE --------
model_obj <- suppressWarnings(icare::ModelTrainAnalysis(
  model_obj,
  methods            = METHODS,
  control            = list(method = "cv", number = CV_FOLDS),
  tuneLength         = 2,
  preProcess         = c("center", "scale"),
  metric_selection   = "ROC",
  task_type          = "binary",
  imbalance_handling = "none",
  loocv_threshold    = 10,
  seed               = SEED
))

cv <- model_obj@process.info$cv_metrics
cv <- cv[order(cv$CV_Mean, decreasing = TRUE), ]
cv_winner <- as.character(cv$Model[1])

# Record what the official default rule (max test-set AUC) would pick, for comparison
test_winner <- as.character(model_obj@all.results$Model[
  which.max(model_obj@all.results$auc)])
cat(sprintf("CV winner = %s (%.4f) | default test-AUC-max = %s -> %s\n",
            cv_winner, cv$CV_Mean[1], test_winner,
            if (cv_winner == test_winner) "same" else "DIFFERENT"))

model_obj <- icare::SelectBestModel(model_obj, custom_selection = cv_winner)

# ---- 8. Official deployment object: replay model + training preProc ----------
preproc <- model_obj@process.info$preprocess_model
stopifnot(!is.null(preproc))

deploy <- icare::ModelDeployment(
  object            = model_obj,
  preproc           = preproc,
  class_labels      = c(NEG, POS),
  model_description = "Undiagnosed dysglycemia triage model"
)

predict_prob <- function(d) {
  deploy$predict_fn(
    newdata        = d[, FEATURES, drop = FALSE],
    selected_model = cv_winner
  )[[POS]]
}

# Self-check: official deployment prediction must match manual replay exactly
manual_prob <- function(d) {
  x <- as.data.frame(predict(preproc, d[, FEATURES, drop = FALSE]))
  as.numeric(predict(model_obj@best.model.result$model,
                     newdata = x, type = "prob")[[POS]])
}
stopifnot(max(abs(predict_prob(internal) - manual_prob(internal))) < 1e-8)

# ---- 9. Thresholds fixed on the internal held-out split only -----------------
p_int   <- predict_prob(internal)
roc_int <- pROC::roc(internal[[GROUP]], p_int,
                     levels = c(NEG, POS), direction = "<", quiet = TRUE)
thr_youden <- as.numeric(pROC::coords(
  roc_int, x = "best", best.method = "youden",
  ret = "threshold", transpose = FALSE)[1])

# ---- 10. External-cycle evaluation (prediction only; nothing re-estimated) ---
p_ext   <- predict_prob(ext)
roc_ext <- pROC::roc(ext[[GROUP]], p_ext,
                     levels = c(NEG, POS), direction = "<", quiet = TRUE)
a_int <- auc_with_ci(internal[[GROUP]], p_int)
a_ext <- auc_with_ci(ext[[GROUP]], p_ext)
perf  <- performance_at(ext[[GROUP]], p_ext, thr_youden)

cat(sprintf("internal AUC %.4f (%.4f-%.4f) | external AUC %.4f (%.4f-%.4f)\n",
            a_int[1], a_int[2], a_int[3], a_ext[1], a_ext[2], a_ext[3]))

results <- data.frame(
  scenario            = "Undiagnosed dysglycemia",
  dev_n               = nrow(dev),  dev_positive = sum(dev[[GROUP]] == POS),
  ext_n               = nrow(ext),  ext_positive = sum(ext[[GROUP]] == POS),
  n_features          = length(FEATURES),
  cv_winner           = cv_winner, cv_roc = cv$CV_Mean[1],
  default_test_winner = test_winner,
  selection_agreement = cv_winner == test_winner,
  internal_auc        = a_int[1], internal_low = a_int[2], internal_high = a_int[3],
  external_auc        = a_ext[1], external_low = a_ext[2], external_high = a_ext[3],
  auc_decay           = a_ext[1] - a_int[1],
  youden_threshold    = thr_youden,
  external_sensitivity = perf["sensitivity"],
  external_specificity = perf["specificity"],
  external_ppv         = perf["ppv"],
  external_npv         = perf["npv"]
)

# Transparent per-algorithm comparison: CV / internal / external.
# All trained models are re-scored on the internal and external partitions so the
# audit table can show whether the CV winner would also win on the test set.
prob_of_model <- function(model, d) {
  x <- as.data.frame(predict(preproc, d[, FEATURES, drop = FALSE]))
  as.numeric(predict(model, newdata = x, type = "prob")[[POS]])
}
cv_table <- do.call(rbind, lapply(names(model_obj@train.models), function(nm) {
  ai <- tryCatch(auc_with_ci(internal[[GROUP]],
                             prob_of_model(model_obj@train.models[[nm]], internal))[1],
                 error = function(e) NA_real_)
  ae <- tryCatch(auc_with_ci(ext[[GROUP]],
                             prob_of_model(model_obj@train.models[[nm]], ext))[1],
                 error = function(e) NA_real_)
  data.frame(model = nm,
             cv_roc = cv$CV_Mean[match(nm, cv$Model)],
             internal_auc = ai, external_auc = ae,
             selected = nm == cv_winner,
             selected_by_default_test_max = nm == test_winner)
}))
rownames(cv_table) <- NULL

# Clinical utility: how many gold-standard tests could be skipped on the external cycle
utility <- do.call(rbind, lapply(TARGET_SENS, function(target) {
  cut <- rule_out_threshold(internal[[GROUP]], p_int, target)
  skip <- p_ext < cut
  data.frame(
    target_sensitivity   = target,
    rule_out_threshold   = cut,
    external_n           = nrow(ext),
    spared_n             = sum(skip),
    spared_pct           = 100 * mean(skip),
    tests_still_needed_pct = 100 * mean(!skip),
    missed_cases         = sum(skip & ext[[GROUP]] == POS),
    false_omission_pct   = 100 * sum(skip & ext[[GROUP]] == POS) / sum(skip),
    external_sensitivity = sum(!skip & ext[[GROUP]] == POS) /
                           sum(ext[[GROUP]] == POS)
  )
}))

# ---- 11. Write results --------------------------------------------------------
write.csv(results,   file.path(OUTPUT_DIR, "results.csv"),           row.names = FALSE)
write.csv(cv_table,  file.path(OUTPUT_DIR, "cv_rankings.csv"),       row.names = FALSE)
write.csv(utility,   file.path(OUTPUT_DIR, "clinical_utility.csv"),  row.names = FALSE)

png(file.path(OUTPUT_DIR, "roc_external.png"),
    width = 1400, height = 1350, res = 200)
plot(roc_ext, col = "#2C7FB8", lwd = 2.6, legacy.axes = TRUE,
     main = "External validation: Undiagnosed dysglycemia",
     xlab = "False positive rate (1 - specificity)",
     ylab = "True positive rate",
     grid = FALSE, auc.polygon = FALSE, print.auc = TRUE,
     print.auc.x = 0.58, print.auc.y = 0.28)
grid(col = "grey88", lty = 3)
legend("bottomright", bty = "n", cex = 0.95,
       sprintf("NHANES 2009-2010 (n = %d)", nrow(ext)))
dev.off()

saveRDS(
  list(deploy = deploy, features = FEATURES,
       cv_winner = cv_winner, youden_threshold = thr_youden,
       rule_out_thresholds = utility$rule_out_threshold,
       results = results, roc_external = roc_ext),
  file = file.path(OUTPUT_DIR, "model.rds")
)

cat("\nClinical utility (gold-standard tests spared):\n")
print(utility[, c("target_sensitivity", "spared_pct",
                  "missed_cases", "false_omission_pct")], row.names = FALSE)
cat("\nDone. Results are in output/\n")
