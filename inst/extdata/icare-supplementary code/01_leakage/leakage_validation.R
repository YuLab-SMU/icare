# ==============================================================================
# 01 - Data-leakage validation of the iCARE pipeline
# ------------------------------------------------------------------------------
# Reviewer question addressed:
#   "Leakage prevention is a core selling point of iCARE, but the manuscript
#    contains no empirical validation."
#
# Idea
#   The whole pipeline is built with iCARE's own functions:
#     CreateModelObject -> FeatureSelectSBF -> ApplyFeatureSelection ->
#     FitApplyPreprocess -> ModelTrainAnalysis
#   The two arms differ in exactly ONE thing: which samples feature screening and
#   preprocessing are allowed to learn from.
#     correct arm : training samples only
#     leaky arm   : all samples, including the future test set
#   Every other step and every other parameter is identical, so the AUC
#   difference is the inflation caused by leakage.
#
# Experiments
#   A  pure-noise data (true AUC = 0.5)    leakage is visible at a glance
#   B  real cohort PMID37633276 (n = 601)  2 x 2 design: preprocessing x screening
#   C  (optional) Module 1 Stat object     fitted on the full cohort vs on training
#   D  optimism bias of selecting the algorithm on the TEST set
#   E  summary numbers and figures
#
# How to run
#   Batch     Rscript 01_leakage/run_leakage.R        (or run_all.R at the root)
#   Stepwise  open this file in RStudio and run block by block; each block is
#             self-contained and leaves its objects in the global environment
#   Fast      QUICK=1 Rscript 01_leakage/run_leakage.R   (5 repetitions per part)
#
# Outputs (git-ignored, written to outputs/01_leakage/)
#   tables/*.csv   figures/*.pdf   sessionInfo.txt
# ==============================================================================

## --- 0. Setup -----------------------------------------------------------------
# `here` locates the repository root from any working directory inside it.
if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

check_packages(required = c("icare", "caret", "pROC", "ggplot2", "readxl"))

# This script follows the icare >= 1.1.0 contract: the group label lives in
# info.data ONLY (see NEWS.md, "Breaking change"). Older versions (<= 1.0.1)
# expect the label inside the feature matrix and fail in Part C with
# "Specified group_col 'group' not found in data columns".
stopifnot("icare >= 1.1.0 is required (restart R after updating)" =
            packageVersion("icare") >= "1.1.0")

suppressPackageStartupMessages({
  library(icare); library(caret); library(pROC)
  library(ggplot2); library(readxl)
})

OUT <- out_dir("01_leakage")     # every file this module writes lands here
dir.create(file.path(OUT, "tables"),  showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE)
set_output_root(OUT)             # icare's own plot dumps follow the same tree

# Repetitions. QUICK=1 gives a fast smoke test with the same code path; the
# published numbers use the defaults.
QUICK   <- env_flag("LEAK_QUICK", env_flag("QUICK"))
R_NOISE <- if (QUICK) 5 else 20   # A: repetitions on noise data
R_REAL  <- if (QUICK) 5 else 20   # B/C: random 70/30 splits of the real cohort
R_SEL   <- if (QUICK) 5 else 30   # D: repetitions of the selection experiment

# ==============================================================================
# A. Pure-noise data: correct vs leaky usage (iCARE functions throughout)
#    120 samples, 500 pure-noise features, random labels -> true AUC = 0.5.
#    Any clearly-above-0.5 test AUC can only come from leakage.
# ==============================================================================

## A1. Generate noise data -------------------------------------------------------
set.seed(1)
n <- 120; p <- 500
X <- as.data.frame(matrix(rnorm(n * p), n, p))            # 500 features unrelated to the label
names(X) <- paste0("V", 1:p)
X$y <- factor(rep(c("Control", "Case"), length.out = n),  # labels: 60 per class, unrelated to features
              levels = c("Control", "Case"))
X <- X[sample(n), ]                                       # shuffle rows
rownames(X) <- paste0("s", 1:n)                           # give each sample an ID (needed for the audit below)
table(X$y)

## A2. Train / test split (70/30) ---------------------------------------------------
idx <- as.vector(createDataPartition(X$y, p = 0.7, list = FALSE))
train_ids <- rownames(X)[idx]                             # 84 training samples
test_ids  <- rownames(X)[-idx]                            # 36 test samples (to be "locked away")

## A3. Create two iCARE modeling objects ----------------------------------------------
obj_train <- CreateModelObject(data = X[idx, ], group_col = "y")  # contains training samples only
obj_all   <- CreateModelObject(data = X,        group_col = "y")  # contains all samples (train + test)
obj_all@split.data <- list(training = X[idx, ], testing = X[-idx, ])  # record the split inside the object

## A4. Sample-level audit: which samples will feature screening "see"? ------------------
# The screening functions read the object's @clean.df, so simply count how
# many test samples it contains:
length(intersect(rownames(obj_train@clean.df), test_ids))  # expected 0  -> no test samples
length(intersect(rownames(obj_all@clean.df),   test_ids))  # expected 36 -> all test samples present

## A5. [CORRECT] Feature screening on the object that holds training samples only ---------
# FeatureSelectSBF = univariate filter (caret::sbf), the same family as the very
# common "pick the significant variables first" step in clinical studies.
sbf_c <- FeatureSelectSBF(obj_train, method = "cv", number = 3, repeats = 1, seed = 1)
feats_correct <- sbf_c$opt_vars
length(feats_correct)

## A6. [LEAKY] Feature screening on the object that holds all samples ------------------------
# Same function, same parameters; only the object differs.
sbf_l <- FeatureSelectSBF(obj_all, method = "cv", number = 3, repeats = 1, seed = 1)
feats_leaky <- sbf_l$opt_vars
length(feats_leaky)

## A7. Identical downstream steps in both arms:
##     apply features -> preprocess -> train glm -> test-set AUC ----------------------------
## --- Correct arm ---
m_c <- ApplyFeatureSelection(obj_all, feats_correct)      # keep the selected features; the (row-level) split is preserved
pp_c <- FitApplyPreprocess(m_c@split.data$training,       # preprocessing parameters are learned from the training set only
                           m_c@split.data$testing,        # and then applied as-is to the test set
                           group_col = "y", method = c("center", "scale"))
m_c@split.scale.data <- list(training = pp_c$training, testing = pp_c$testing)
m_c@filtered.set     <- m_c@split.scale.data              # the slot read by ModelTrainAnalysis
m_c <- ModelTrainAnalysis(m_c, methods = "glm",
                          control = list(method = "cv", number = 5, repeats = 1),
                          tuneLength = 1, preProcess = NULL, loocv_threshold = 10,
                          metric_selection = "ROC", imbalance_handling = "none", seed = 1)
prob_c <- predict(m_c@train.models$glm, m_c@filtered.set$testing, type = "prob")[, "Case"]
auc_correct <- auc_of(m_c@filtered.set$testing$y, prob_c)

## --- Leaky arm (line-by-line identical to the above; only feats_leaky differs) ---
m_l <- ApplyFeatureSelection(obj_all, feats_leaky)
pp_l <- FitApplyPreprocess(m_l@split.data$training, m_l@split.data$testing,
                           group_col = "y", method = c("center", "scale"))
m_l@split.scale.data <- list(training = pp_l$training, testing = pp_l$testing)
m_l@filtered.set     <- m_l@split.scale.data
m_l <- ModelTrainAnalysis(m_l, methods = "glm",
                          control = list(method = "cv", number = 5, repeats = 1),
                          tuneLength = 1, preProcess = NULL, loocv_threshold = 10,
                          metric_selection = "ROC", imbalance_handling = "none", seed = 1)
prob_l <- predict(m_l@train.models$glm, m_l@filtered.set$testing, type = "prob")[, "Case"]
auc_leaky <- auc_of(m_l@filtered.set$testing$y, prob_l)

## A8. Single-run result ------------------------------------------------------------
c(correct = auc_correct, leaky = auc_leaky, inflation = auc_leaky - auc_correct)
# The data contain no signal, so the true AUC is 0.5. The leaky arm should be
# clearly above 0.5; the correct arm should fluctuate around 0.5.

## A9. Repeat R_NOISE times (new data and split each time) to get a distribution ---------
# This loop simply repeats A1-A8 with no new logic; the inner "arm" loop
# replaces the two copied blocks above.
exp_A <- NULL
for (s in 1:R_NOISE) {
  set.seed(s)
  n <- 120; p <- 500
  Xs <- as.data.frame(matrix(rnorm(n * p), n, p)); names(Xs) <- paste0("V", 1:p)
  Xs$y <- factor(rep(c("Control", "Case"), length.out = n), levels = c("Control", "Case"))
  Xs <- Xs[sample(n), ]; rownames(Xs) <- paste0("s", 1:n)
  ix <- as.vector(createDataPartition(Xs$y, p = 0.7, list = FALSE))
  o_tr  <- CreateModelObject(data = Xs[ix, ], group_col = "y")
  o_all <- CreateModelObject(data = Xs,       group_col = "y")
  o_all@split.data <- list(training = Xs[ix, ], testing = Xs[-ix, ])
  res <- c()
  for (arm in c("correct", "leaky")) {
    src   <- if (arm == "correct") o_tr else o_all                    # the only difference between arms
    feats <- FeatureSelectSBF(src, method = "cv", number = 3, repeats = 1, seed = s)$opt_vars
    mm <- ApplyFeatureSelection(o_all, feats)
    pp <- FitApplyPreprocess(mm@split.data$training, mm@split.data$testing, "y", c("center", "scale"))
    mm@split.scale.data <- list(training = pp$training, testing = pp$testing)
    mm@filtered.set <- mm@split.scale.data
    mm <- suppressWarnings(ModelTrainAnalysis(mm, methods = "glm",
                                              control = list(method = "cv", number = 5, repeats = 1), tuneLength = 1,
                                              preProcess = NULL, loocv_threshold = 10, metric_selection = "ROC",
                                              imbalance_handling = "none", seed = s))
    pr <- predict(mm@train.models$glm, mm@filtered.set$testing, type = "prob")[, "Case"]
    res[arm] <- auc_of(mm@filtered.set$testing$y, pr)
    # Sample-level audit: number of test samples seen by screening (must be 0 for the correct arm)
    res[paste0("seen_test_", arm)] <- length(intersect(rownames(src@clean.df), rownames(Xs)[-ix]))
  }
  exp_A <- rbind(exp_A, data.frame(sim = s, correct = res[["correct"]], leaky = res[["leaky"]],
                                   seen_test_correct = res[["seen_test_correct"]],
                                   seen_test_leaky = res[["seen_test_leaky"]]))
}
exp_A$inflation <- exp_A$leaky - exp_A$correct
write.csv(exp_A, file.path(OUT, "tables", "expA_noise_icare.csv"), row.names = FALSE)

## A10. Summary ---------------------------------------------------------------------
colMeans(exp_A[, c("correct", "leaky", "inflation")])
ci95(exp_A$inflation)                                          # 95% CI of the inflation
wilcox.test(exp_A$leaky, exp_A$correct, paired = TRUE)         # paired test
table(exp_A$seen_test_correct)                                 # must be all 0
table(exp_A$seen_test_leaky)                                   # should be all 36


# ==============================================================================
# B. Real cohort (PMID37633276, 601 samples): 2 x 2 leakage design
#    Factor 1  Where preprocessing (center/scale/YeoJohnson) is fitted:
#              training set / all samples
#    Factor 2  Which samples feature screening uses:
#              training set / all samples (including the validation set)
#    "train|train" is the correct workflow; each of the other three arms
#    introduces one kind of leakage.
# ==============================================================================

## B1. Read the data -------------------------------------------------------------------
xlsx <- file.path(ROOT, "data", "PMID37633276_DIA_plasmaproomic.xlsx")
if (!file.exists(xlsx))   # fall back to the copy shipped inside the icare package
  xlsx <- system.file("extdata", "PMID37633276_DIA_plasmaproomic.xlsx", package = "icare")
stopifnot("Workbook PMID37633276_DIA_plasmaproomic.xlsx not found" = nzchar(xlsx) && file.exists(xlsx))
raw_dis <- as.data.frame(read_excel(xlsx, sheet = 3))
raw_dis <- raw_dis[!grepl("^KC", raw_dis$`Case ID`), ]    # drop KC samples (as in the manuscript)
feat_cols <- c("White blood cell count", "Lymphocytes%", "Monocytes%", "Neutrophilic granulocyte%",
               "Eosinophil%", "Basophil%", "Lymphocytes count", "Monocytes count",
               "Neutrophilic granulocyte count", "Eosinophil count", "Basophil count",
               "Red blood cell count", "Hemoglobin", "Hematocrit", "Mean red blood cell volume",
               "Mean cell hemoglobin", "Mean corpuscular hemoglobin concentration",
               "Red blood cell distribution width CV", "Red blood cell distribution width SD",
               "Platelet count", "Mean platelet volume", "PCT", "Platelet volume distribution width")
m <- as.matrix(raw_dis[, feat_cols]); storage.mode(m) <- "numeric"
rownames(m) <- raw_dis$`Case ID`
mat_raw <- as.data.frame(m, check.names = FALSE)           # version with original column names (used in Part C)
dim(mat_raw); sum(is.na(mat_raw))                          # 601 x 23; no missing values in this data

## B2. Build a "features + y" data frame ----------------------------------------------------
dat <- modify_column_names(mat_raw)                        # icare's column-name normalization (removes spaces and %)
dat$y <- factor(ifelse(grepl("^HC", raw_dis$`Case ID`), "Control", "Case"),
                levels = c("Control", "Case"))
table(dat$y)

## B3. One split: walk through the [CORRECT] arm line by line -----------------------------------
set.seed(1)
tr <- as.vector(createDataPartition(dat$y, p = 0.7, list = FALSE))
va <- setdiff(seq_len(nrow(dat)), tr)
dat_tr <- dat[tr, ]; dat_va <- dat[va, ]
num_cols <- setdiff(names(dat), "y")
PP <- c("center", "scale", "YeoJohnson")                   # preprocessing methods (no missing data, so no imputation)

pp_tr <- FitApplyPreprocess(dat_tr, dat_va, group_col = "y", method = PP)   # parameters learned from the training set only
obj_b <- CreateModelObject(data = pp_tr$training, group_col = "y")          # object holds training samples only
sbf_b <- FeatureSelectSBF(obj_b, method = "cv", number = 5, repeats = 1, seed = 1)
feats_b <- sbf_b$opt_vars                                  # feature screening used the training set only
feats_b

mb <- CreateModelObject(data = rbind(pp_tr$training, pp_tr$testing), group_col = "y")
mb@split.data <- list(training = pp_tr$training, testing = pp_tr$testing)
mb <- ApplyFeatureSelection(mb, feats_b)
mb@split.scale.data <- mb@split.data; mb@filtered.set <- mb@split.data    # already preprocessed; place directly
mb <- ModelTrainAnalysis(mb, methods = "glm",
                         control = list(method = "cv", number = 5, repeats = 1), tuneLength = 1,
                         preProcess = NULL, loocv_threshold = 10, metric_selection = "ROC",
                         imbalance_handling = "none", seed = 1)
pb <- predict(mb@train.models$glm, mb@filtered.set$testing, type = "prob")[, "Case"]
auc_of(mb@filtered.set$testing$y, pb)                      # validation AUC of the correct arm

## B4. The other three arms only change "which samples are used as input" (code is in the
##     B5 loop). The logic in brief:
#   Preprocessing leak : FitApplyPreprocess(training = dat (all samples), testing = dat_va)
#                        -> parameters come from the full cohort
#   Screening leak     : CreateModelObject(data = rbind(train, validation))
#                        -> screening sees the validation set
#   Both leak          : both of the above at once

## B5. Repeat over R_REAL random splits and run the full 2 x 2 ---------------------------------
N_NOISE <- 0        # set to 500 to add 500 noise features as a stress test (slow; start with 0)
exp_B <- NULL
for (s in 1:R_REAL) {
  set.seed(s)
  d <- dat
  if (N_NOISE > 0) {                                       # optional: append noise features to amplify leakage
    nz <- as.data.frame(matrix(rnorm(nrow(d) * N_NOISE), nrow(d)))
    names(nz) <- paste0("noise", seq_len(N_NOISE)); rownames(nz) <- rownames(d)
    d <- cbind(nz, d)
  }
  tr_ <- as.vector(createDataPartition(d$y, p = 0.7, list = FALSE)); va_ <- setdiff(seq_len(nrow(d)), tr_)
  d_tr <- d[tr_, ]; d_va <- d[va_, ]; cols <- setdiff(names(d), "y")
  pp_train <- FitApplyPreprocess(d_tr, d_va, "y", PP)      # preprocessing fitted on the training set
  pp_full  <- FitApplyPreprocess(d,    d_va, "y", PP)      # preprocessing fitted on all samples (leak)
  tr_leak <- d_tr; tr_leak[, cols] <- predict(pp_full$preProc, d_tr[, cols])
  prep_sets <- list(train = list(tr = pp_train$training, va = pp_train$testing),
                    full  = list(tr = tr_leak,           va = pp_full$testing))
  out <- c()
  for (prep in c("train", "full")) for (scr in c("train", "all")) {
    D <- prep_sets[[prep]]
    scr_data <- if (scr == "train") D$tr else rbind(D$tr, D$va)       # samples visible to screening
    feats <- FeatureSelectSBF(CreateModelObject(data = scr_data, group_col = "y"),
                              method = "cv", number = 5, repeats = 1, seed = s)$opt_vars
    mm <- CreateModelObject(data = rbind(D$tr, D$va), group_col = "y")
    mm@split.data <- list(training = D$tr, testing = D$va)
    mm <- ApplyFeatureSelection(mm, feats)
    mm@split.scale.data <- mm@split.data; mm@filtered.set <- mm@split.data
    mm <- suppressWarnings(ModelTrainAnalysis(mm, methods = "glm",
                                              control = list(method = "cv", number = 5, repeats = 1), tuneLength = 1,
                                              preProcess = NULL, loocv_threshold = 10, metric_selection = "ROC",
                                              imbalance_handling = "none", seed = s))
    pr <- predict(mm@train.models$glm, mm@filtered.set$testing, type = "prob")[, "Case"]
    out[paste(prep, scr, sep = "|")] <- auc_of(mm@filtered.set$testing$y, pr)
  }
  exp_B <- rbind(exp_B, data.frame(sim = s, t(out), check.names = FALSE))
}
write.csv(exp_B, file.path(OUT, "tables", "expB_real_cohort_2x2.csv"), row.names = FALSE)

## B6. Summary: bias of each arm relative to the correct arm ----------------------------------
arm_names <- c("train|train", "full|train", "train|all", "full|all")
summ_B <- do.call(rbind, lapply(arm_names, function(a) data.frame(
  arm = a, mean_AUC = mean(exp_B[[a]]), lo = ci95(exp_B[[a]])[1], hi = ci95(exp_B[[a]])[2],
  bias_vs_correct = mean(exp_B[[a]] - exp_B[["train|train"]]))))
summ_B
# Interpretation note: the real cohort has only 23 features and a strong signal,
# so the inflation caused by leakage is usually small (possibly invisible). That
# is an honest result in itself -- the harm of leakage grows with many features,
# weak signal and small samples, which is exactly the situation Experiment A shows.


# ==============================================================================
# C. (Optional) Module 1 Stat object: fitted on the full cohort vs training only
#    This is the most realistic leakage path in the tutorials: Module 1 cleaning
#    (imputation / outliers / scaling) is run on the whole cohort first, and only
#    then is the object converted with ConvertObject() and split -- so the
#    preprocessing parameters already contain the test samples.
#    Requires: a full icare installation + icare_stat_replay_helpers.R in the
#    same directory.
# ==============================================================================
source(file.path(ROOT, "R", "stat_replay.R"))   # build_stat / replay_stat_data / transform_rows

## C1. Build the info table needed by build_stat (row names = sample IDs, with a group column) ----
inf <- data.frame(group = factor(ifelse(grepl("^HC", raw_dis$`Case ID`), "Healthy Control", "Case"),
                                 levels = c("Healthy Control", "Case")),
                  row.names = raw_dis$`Case ID`)

## C2. One split: fit the Stat object on the training set vs on the full cohort ------------------
set.seed(1)
tr <- as.vector(createDataPartition(inf$group, p = 0.7, list = FALSE)); va <- setdiff(seq_len(nrow(inf)), tr)
so_train <- build_stat(mat_raw[tr, ], inf[tr, , drop = FALSE])    # correct: parameters from the training set only
so_full  <- build_stat(mat_raw,       inf)                        # leaky: parameters from the full cohort

## C3. Apply both parameter sets as-is to training / validation samples (no re-estimation) ------
Xtr_c <- transform_rows(so_train, mat_raw[tr, ]); Xva_c <- transform_rows(so_train, mat_raw[va, ])
Xtr_l <- transform_rows(so_full,  mat_raw[tr, ]); Xva_l <- transform_rows(so_full,  mat_raw[va, ])

## C4. Remaining steps identical: SBF screening on training samples only -> glm -> validation AUC
y <- factor(ifelse(inf$group == "Case", "Case", "Control"), levels = c("Control", "Case")); names(y) <- rownames(inf)
auc_stat <- c()
for (arm in c("correct", "leaky")) {
  Xtr <- if (arm == "correct") Xtr_c else Xtr_l
  Xva <- if (arm == "correct") Xva_c else Xva_l
  Xtr <- Xtr[, vapply(Xtr, is.numeric, logical(1)), drop = FALSE]; Xva <- Xva[, names(Xtr), drop = FALSE]
  trd <- cbind(Xtr, y = y[rownames(Xtr)]); vad <- cbind(Xva, y = y[rownames(Xva)])
  feats <- FeatureSelectSBF(CreateModelObject(data = trd, group_col = "y"),
                            method = "cv", number = 5, repeats = 1, seed = 1)$opt_vars
  mm <- CreateModelObject(data = rbind(trd, vad), group_col = "y")
  mm@split.data <- list(training = trd, testing = vad)
  mm <- ApplyFeatureSelection(mm, feats)
  mm@split.scale.data <- mm@split.data; mm@filtered.set <- mm@split.data
  mm <- suppressWarnings(ModelTrainAnalysis(mm, methods = "glm",
                                            control = list(method = "cv", number = 5, repeats = 1), tuneLength = 1,
                                            preProcess = NULL, loocv_threshold = 10, metric_selection = "ROC",
                                            imbalance_handling = "none", seed = 1))
  pr <- predict(mm@train.models$glm, mm@filtered.set$testing, type = "prob")[, "Case"]
  auc_stat[arm] <- auc_of(mm@filtered.set$testing$y, pr)
}
auc_stat
# For multiple splits, wrap C2-C4 in a for (s in 1:R_REAL) loop
# (call set.seed(s) at the start of each iteration to redraw the split).


# ==============================================================================
# D. Optimism bias of selecting the algorithm on the TEST set
#    iCARE's ModelTrainAnalysis() stores @all.results as performance on the
#    TEST set, and SelectBestModel(metric = "auc") picks the best model from
#    that table by default.
#    Below, on pure-noise data, we compare two ways of selecting:
#      (1) pick the model with the highest test-set AUC
#      (2) pick the model by training-set CV, then look at its test-set AUC
#    The true AUC is 0.5, so how far the selected model's test AUC exceeds 0.5
#    is the selection bias.
#    IMPORTANT: iCARE computes the AUC stored in @all.results with pROC
#    direction = "auto" (model_05_train.R), which flips any AUC < 0.5 upward.
#    On null data this inflates every number, including the "CV-picked" one
#    (mean 0.556 instead of ~0.5). We therefore recompute each model's test AUC
#    with the direction fixed ("<") via auc_of(), and keep the iCARE values
#    alongside for transparency.
# ==============================================================================

## D1. One run: noise data (10 features only), train 4 algorithms ---------------------------------
set.seed(11)
n <- 120
Xd <- as.data.frame(matrix(rnorm(n * 10), n, 10)); names(Xd) <- paste0("V", 1:10)
Xd$y <- factor(rep(c("Control", "Case"), length.out = n), levels = c("Control", "Case"))
Xd <- Xd[sample(n), ]; rownames(Xd) <- paste0("s", 1:n)
ix <- as.vector(createDataPartition(Xd$y, p = 0.7, list = FALSE))
md <- CreateModelObject(data = Xd, group_col = "y")
md@split.data <- list(training = Xd[ix, ], testing = Xd[-ix, ])
md@split.scale.data <- md@split.data; md@filtered.set <- md@split.data
md <- suppressWarnings(ModelTrainAnalysis(md, methods = c("glm", "rf", "knn", "rpart"),
                                          control = list(method = "cv", number = 5, repeats = 1), tuneLength = 1,
                                          preProcess = NULL, loocv_threshold = 10, metric_selection = "ROC",
                                          imbalance_handling = "none", seed = 11))

## D2. Two tables: one on the test set, one from CV inside the training set -------------------------
te_d <- md@filtered.set$testing
# Test AUC of every model with the ROC direction fixed (iCARE's own value uses direction = "auto")
auc_fixed <- sapply(md@all.results$Model, function(k)
  auc_of(te_d$y, predict(md@train.models[[k]], te_d, type = "prob")[, "Case"]))
cbind(md@all.results[, c("Model", "auc")], auc_fixed)   # 'auc' = iCARE (TEST set, direction auto); 'auc_fixed' = TEST set, direction "<"
md@process.info$cv_metrics               # training-set CV performance; model selection should use this
max(auc_fixed)                           # method (1): highest TEST AUC among models (inflated by selection bias)
cv_pick <- md@process.info$cv_metrics$Model[1]          # method (2): the model chosen by CV
auc_fixed[[cv_pick]]                     # its test-set AUC (unbiased)

## D3. Repeat R_SEL times --------------------------------------------------------------------------
exp_D <- NULL
for (s in 1:R_SEL) {
  set.seed(100 + s); n <- 120
  Xs <- as.data.frame(matrix(rnorm(n * 10), n, 10)); names(Xs) <- paste0("V", 1:10)
  Xs$y <- factor(rep(c("Control", "Case"), length.out = n), levels = c("Control", "Case"))
  Xs <- Xs[sample(n), ]; rownames(Xs) <- paste0("s", 1:n)
  ix <- as.vector(createDataPartition(Xs$y, p = 0.7, list = FALSE))
  mm <- CreateModelObject(data = Xs, group_col = "y")
  mm@split.data <- list(training = Xs[ix, ], testing = Xs[-ix, ])
  mm@split.scale.data <- mm@split.data; mm@filtered.set <- mm@split.data
  mm <- suppressWarnings(ModelTrainAnalysis(mm, methods = c("glm", "rf", "knn", "rpart"),
                                            control = list(method = "cv", number = 5, repeats = 1), tuneLength = 1,
                                            preProcess = NULL, loocv_threshold = 10, metric_selection = "ROC",
                                            imbalance_handling = "none", seed = s))
  cvp <- mm@process.info$cv_metrics$Model[1]
  te  <- mm@filtered.set$testing
  fx  <- sapply(mm@all.results$Model, function(k)
    auc_of(te$y, predict(mm@train.models[[k]], te, type = "prob")[, "Case"]))
  exp_D <- rbind(exp_D, data.frame(sim = s,
                                   pick_by_test      = max(fx),                       # direction fixed (used for the figure / table)
                                   pick_by_cv        = fx[[cvp]],
                                   pick_by_test_auto = max(mm@all.results$auc),       # iCARE value (direction = "auto"), for transparency
                                   pick_by_cv_auto   = mm@all.results$auc[mm@all.results$Model == cvp]))
}
write.csv(exp_D, file.path(OUT, "tables", "expD_model_selection_bias.csv"), row.names = FALSE)
colMeans(exp_D[, c("pick_by_test", "pick_by_cv", "pick_by_test_auto", "pick_by_cv_auto")])
ci95(exp_D$pick_by_test - exp_D$pick_by_cv)
wilcox.test(exp_D$pick_by_test, exp_D$pick_by_cv, paired = TRUE)


# ==============================================================================
# E. Summary numbers and figures (PDF output only)
# ==============================================================================

## E1. Figure A: noise data, correct vs leaky ----------------------------------------------------
long_A <- rbind(data.frame(sim = exp_A$sim, Pipeline = "Correct (screen on train only)", AUC = exp_A$correct),
                data.frame(sim = exp_A$sim, Pipeline = "Leaky (screen on all samples)",   AUC = exp_A$leaky))
fig_A <- ggplot(long_A, aes(Pipeline, AUC, fill = Pipeline)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_line(aes(group = sim), colour = "grey70", linewidth = .3) +
  geom_boxplot(width = .45, alpha = .7, outlier.shape = NA) +
  geom_point(position = position_jitter(width = .03), size = 1.5) +
  scale_fill_manual(values = c("#2c7fb8", "#d7301f")) +
  theme_bw(base_size = 11) + theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(title = "Held-out AUC on pure-noise data (true AUC = 0.5), iCARE pipeline",
       subtitle = sprintf("n = 120, p = 500, %d repeats; mean inflation = %.3f (95%% CI %.3f-%.3f)",
                          R_NOISE, mean(exp_A$inflation), ci95(exp_A$inflation)[1], ci95(exp_A$inflation)[2]),
       x = NULL, y = "Held-out AUC")
fig_A
ggsave(file.path(OUT, "figures", "FigA_noise_icare.pdf"), fig_A, width = 6, height = 4.8)

## E2. Figure B: real cohort 2x2 -------------------------------------------------------------------
summ_B$arm <- factor(summ_B$arm, levels = arm_names,
                     labels = c("Correct\n(prep: train, screen: train)", "Prep leak\n(prep: all)",
                                "Screen leak\n(screen: train+valid)", "Both leak"))
fig_B <- ggplot(summ_B, aes(arm, mean_AUC, colour = arm == levels(arm)[1])) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6, linewidth = 0.8) +
  scale_colour_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "#d7301f"), guide = "none") +
  coord_cartesian(ylim = c(0.90, 1.00)) +                  # fixed range: avoids exaggerating differences of <0.01
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8)) +
  labs(title = "Leakage factorial on the PMID37633276 cohort (iCARE functions)",
       subtitle = sprintf("validation AUC, mean +/- 95%% CI over %d random 70/30 splits", R_REAL),
       x = NULL, y = "Validation AUC")
fig_B
ggsave(file.path(OUT, "figures", "FigB_real_cohort_2x2.pdf"), fig_B, width = 7, height = 4.5)

## E3. Figure D: selecting on the test set vs selecting by CV ----------------------------------------
long_D <- rbind(data.frame(Selection = "By test-set AUC", AUC = exp_D$pick_by_test),
                data.frame(Selection = "By training-set CV", AUC = exp_D$pick_by_cv))
fig_D <- ggplot(long_D, aes(Selection, AUC, fill = Selection)) +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_boxplot(width = .45, alpha = .7, outlier.shape = NA) +
  geom_jitter(width = .06, size = 1.4) + scale_fill_manual(values = c("#d7301f", "#1a9850")) +
  theme_bw(base_size = 11) + theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(title = "Test AUC of the selected model on noise data",
       subtitle = "True AUC = 0.5; 4 algorithms; selecting on the test set inflates the reported AUC",
       x = NULL, y = "Test AUC (ROC direction fixed)")
fig_D
ggsave(file.path(OUT, "figures", "FigD_selection_bias.pdf"), fig_D, width = 5.5, height = 4.6)

## E4. Numbers ready to paste into the manuscript ---------------------------------------------------
cat(sprintf("[A] noise: correct AUC %.3f, leaky %.3f, inflation %.3f (95%% CI %.3f-%.3f); held-out samples seen by screening: correct=%d, leaky>0 in %d/%d runs\n",
            mean(exp_A$correct), mean(exp_A$leaky), mean(exp_A$inflation),
            ci95(exp_A$inflation)[1], ci95(exp_A$inflation)[2],
            sum(exp_A$seen_test_correct), sum(exp_A$seen_test_leaky > 0), nrow(exp_A)))
print(summ_B)
cat(sprintf("[D] noise (direction fixed): AUC of model picked on test set %.3f vs picked by CV %.3f; iCARE auto-direction values %.3f vs %.3f\n",
            mean(exp_D$pick_by_test), mean(exp_D$pick_by_cv),
            mean(exp_D$pick_by_test_auto), mean(exp_D$pick_by_cv_auto)))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))