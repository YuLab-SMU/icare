# ==============================================================================
# 02 - Head-to-head: iCARE vs caret / tidymodels / mlr3
# ------------------------------------------------------------------------------
# Reviewer question addressed:
#   "Compare iCARE with caret / mlr3 / tidymodels pipelines in performance and in
#    the amount of code."
#
# The first version of this experiment compared a single task (z-score + glm) on a
# single dataset with a single metric, and all four frameworks returned
# bit-identical AUCs. That shows equivalence on the simplest possible task but says
# nothing about the tasks where iCARE is meant to help. This version compares the
# four frameworks on FOUR tasks, FOUR datasets and FIVE axes.
#
# TASKS  (each framework's code lives in blocks/<task>_<framework>.R and this
#         script EXECUTES those files, so the lines counted are the lines run)
#   T1  z-score (fitted on training) -> glm                             [basic]
#   T2  z-score -> univariate top-K screening -> glm, all fitted on the training
#       partition only                              [leakage-prone, multi-step]
#   T3  multi-algorithm benchmark: glmnet / ranger / rpart with identical fixed
#       hyper-parameters, 5-fold CV on the training partition, best algorithm by
#       CV AUC, evaluated on the held-out partition
#   T4  scalability: runtime of the T2 pipeline over a grid of n and p
#
# DATASETS
#   UTUC_lab        real, n = 601, p = 23   (clinical laboratory; strong signal)
#   UTUC_prot       real, n = 601, p = 730  (plasma proteome; p > n_train)
#   UTUC_prot_n120  real, n = 120, p = 730  (small-n / high-dim subsample)
#   Noise           synthetic, n = 120, p = 500, random labels (true AUC 0.5; any
#                   AUC clearly above 0.5 would indicate leakage)
#
# AXES
#   1  discrimination (AUC, direction fixed) and calibration (Brier score)
#   2  paired equivalence against iCARE (TOST, margin DELTA, same splits)
#   3  agreement of the selected features (T2) / selected algorithm (T3)
#   4  runtime, including scaling with n and p (T4)
#   5  code volume: executable lines, function calls, distinct functions, packages
#
# FAIRNESS NOTES (also reported in the Methods)
#   * All frameworks see the SAME training / validation split for each seed.
#   * iCARE's ModelTrainAnalysis() always runs CV on the training partition; in
#     T1/T2 the other frameworks fit a single model with no CV, so iCARE does more
#     work for its lines. T3 uses 5-fold CV everywhere.
#   * iCARE is built on caret, so identical numbers to caret are expected by
#     construction; the informative comparisons are tidymodels and mlr3.
#   * T3 hyper-parameters are fixed and identical (glmnet alpha 0.5 / lambda 0.01;
#     ranger 500 trees, mtry = floor(sqrt(p)), minimum node size 1; rpart cp 0.01),
#     so AUC differences reflect the framework, not tuning defaults.
#   * T2 screening: caret / mlr3 / iCARE rank features by the same ANOVA F-test and
#     keep the top K (the selected sets are checked for identity). tidymodels has
#     no built-in supervised filter; the closest standard option is
#     colino::step_select_roc (ranks by univariate ROC AUC), so its selected set may
#     differ -> reported as set overlap, not hidden.
#   * Runtime is wall-clock of the block only (AUC / Brier excluded), measured after
#     one warm-up pass; frameworks are timed sequentially, not in parallel.
#
# How to run
#   Rscript 02_head2head/framework_head2head.R             ~30-45 min, full settings
#   QUICK=1 Rscript 02_head2head/framework_head2head.R     3 repetitions, small grid
#   H2H_DATASETS=Noise QUICK=1 Rscript ...                 smoke test, one dataset
#
# Outputs (git-ignored, written to outputs/02_head2head/)
#   tables/*.csv   figures/*.pdf,*.png   expE_sessionInfo.txt
# ==============================================================================

## --- 0. Setup -----------------------------------------------------------------
# `here` locates the repository root from any working directory inside it.
if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

QUICK <- env_flag("H2H_QUICK", env_flag("QUICK"))

need <- c("icare", "caret", "tidymodels", "mlr3", "mlr3pipelines", "mlr3learners",
          "mlr3filters", "pROC", "readxl", "ggplot2", "dplyr", "glmnet", "ranger", "rpart")
check_packages(required = need, optional = "colino")

# colino adds supervised filter steps to recipes (tidymodels has none built in)
HAS_COLINO <- requireNamespace("colino", quietly = TRUE)
if (!HAS_COLINO)
  message("Package 'colino' not found -> tidymodels is skipped in T2/T4 ",
          "(install.packages('colino')).")

suppressMessages({
  library(icare); library(caret); library(tidymodels)
  if (HAS_COLINO) library(colino)
  library(mlr3); library(mlr3pipelines); library(mlr3learners); library(mlr3filters)
  library(pROC); library(readxl); library(ggplot2); library(dplyr)
})
lgr::get_logger("mlr3")$set_threshold("warn")

# Every artefact lands in outputs/02_head2head; icare's plot dumps follow the same tree.
OUT    <- out_dir("02_head2head")
BLOCKS <- file.path(ROOT, "02_head2head", "blocks")
for (d in c("tables", "figures"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
invisible(utils::capture.output(set_output_root(OUT)))

FW      <- c("caret", "tidymodels", "mlr3", "icare")     # icare is the reference
ALGOS   <- c("glmnet", "ranger", "rpart")                # T3 algorithms
R_H2H   <- if (QUICK) 3 else 20                          # random 70/30 splits per dataset
R_SCALE <- if (QUICK) 1 else 3                           # repetitions per cell in T4
SCALE_N <- if (QUICK) c(200, 1000) else c(200, 1000, 2000)
SCALE_P <- if (QUICK) c(100, 1000) else c(100, 1000, 2000)
DELTA   <- 0.01                                          # equivalence margin on AUC

if (nzchar(Sys.getenv("H2H_REPS")))                       # override for longer/shorter runs
  R_H2H <- as.integer(Sys.getenv("H2H_REPS"))

FW_COL <- c(caret = "#1b9e77", tidymodels = "#d95f02", mlr3 = "#7570b3", icare = "#e7298a")
## E2. Data ---------------------------------------------------------------------
xlsx <- file.path(ROOT, "data", "PMID37633276_DIA_plasmaproomic.xlsx")
if (!file.exists(xlsx))   # fall back to the copy shipped inside the icare package
  xlsx <- system.file("extdata", "PMID37633276_DIA_plasmaproomic.xlsx", package = "icare")
stopifnot("Workbook PMID37633276_DIA_plasmaproomic.xlsx not found" =
            nzchar(xlsx) && file.exists(xlsx))

## --- UTUC_lab: 23 routine laboratory features (identical to Experiment B) -------
raw_dis <- as.data.frame(read_excel(xlsx, sheet = 3))
raw_dis <- raw_dis[!grepl("^KC", raw_dis$`Case ID`), ]
feat_cols <- c("White blood cell count", "Lymphocytes%", "Monocytes%", "Neutrophilic granulocyte%",
               "Eosinophil%", "Basophil%", "Lymphocytes count", "Monocytes count",
               "Neutrophilic granulocyte count", "Eosinophil count", "Basophil count",
               "Red blood cell count", "Hemoglobin", "Hematocrit", "Mean red blood cell volume",
               "Mean cell hemoglobin", "Mean corpuscular hemoglobin concentration",
               "Red blood cell distribution width CV", "Red blood cell distribution width SD",
               "Platelet count", "Mean platelet volume", "PCT", "Platelet volume distribution width")
m <- as.matrix(raw_dis[, feat_cols]); storage.mode(m) <- "numeric"
rownames(m) <- raw_dis$`Case ID`
D_lab <- modify_column_names(as.data.frame(m, check.names = FALSE))
D_lab$y <- factor(ifelse(grepl("^HC", raw_dis$`Case ID`), "Control", "Case"),
                  levels = c("Control", "Case"))

## --- UTUC_prot: plasma proteome of the same 601 patients (Table S1E) ------------
# Label-free preprocessing only (uses no outcome information, so it cannot leak):
#   keep proteins quantified in >= 80% of all samples; log2(x + 1); a protein that
#   was not detected (NA) is set to 0, i.e. "below detection" -- a constant, not a
#   fitted parameter. Everything that is fitted (scaling, screening, models) is
#   done inside each framework on the training partition only.
S1B <- as.data.frame(read_excel(xlsx, sheet = "Table S1B"))
S1E <- as.data.frame(read_excel(xlsx, sheet = "Table S1E"))
S1B <- S1B[S1B$`Plasma proteome ID` %in% colnames(S1E), ]        # drops tissue-only rows
pm  <- as.matrix(S1E[, -1]); storage.mode(pm) <- "numeric"
rownames(pm) <- make.unique(S1E[[1]])
pm  <- t(pm)[S1B$`Plasma proteome ID`, ]; rownames(pm) <- S1B$`Case ID`
pm  <- log2(pm[, colMeans(!is.na(pm)) >= 0.8] + 1); pm[is.na(pm)] <- 0
D_prot <- modify_column_names(as.data.frame(pm, check.names = FALSE))
D_prot$y <- factor(ifelse(grepl("^HC", rownames(pm)), "Control", "Case"),
                   levels = c("Control", "Case"))

## --- UTUC_prot_n120: balanced small-n subsample of the proteome (n << p) --------
set.seed(2024)
D_small <- D_prot[c(sample(which(D_prot$y == "Control"), 60),
                    sample(which(D_prot$y == "Case"),    60)), ]

## --- Noise: regenerated for every repetition (true AUC = 0.5) -------------------
make_noise <- function(s, n = 120, p = 500) {
  set.seed(s)
  X <- as.data.frame(matrix(rnorm(n * p), n, p)); names(X) <- paste0("V", 1:p)
  X$y <- factor(rep(c("Control", "Case"), length.out = n), levels = c("Control", "Case"))
  X <- X[sample(n), ]; rownames(X) <- paste0("s", 1:n); X
}

## --- Dataset registry: K = number of features kept by the T2 screening step ----
DS <- list(
  UTUC_lab       = list(data = D_lab,   K = 10, tasks = c("T1", "T2", "T3")),
  UTUC_prot      = list(data = D_prot,  K = 30, tasks = c("T2", "T3")),
  UTUC_prot_n120 = list(data = D_small, K = 20, tasks = c("T2", "T3")),
  Noise          = list(data = NULL,    K = 10, tasks = c("T2", "T3")))
## Optional subset for a smoke test: H2H_DATASETS=Noise (or Noise,UTUC_lab)
.sel <- trimws(strsplit(Sys.getenv("H2H_DATASETS", ""), ",")[[1]])
.sel <- .sel[nzchar(.sel)]
if (length(.sel)) {
  unknown <- setdiff(.sel, names(DS))
  if (length(unknown))
    stop("Unknown dataset(s) in H2H_DATASETS: ", paste(unknown, collapse = ", "),
         "\n  available: ", paste(names(DS), collapse = ", "), call. = FALSE)
  DS <- DS[.sel]
}

for (nm in names(DS)) if (!is.null(DS[[nm]]$data))
  cat(sprintf("%-15s n = %3d, p = %3d, cases = %d\n", nm, nrow(DS[[nm]]$data),
              ncol(DS[[nm]]$data) - 1, sum(DS[[nm]]$data$y == "Case")))
cat(sprintf("%-15s n = 120, p = 500, cases = 60 (regenerated per repetition)\n", "Noise"))

## E3. Harness: run one framework block on one split ----------------------------
# A block (blocks/<task>_<framework>.R) reads d_tr, d_va, s (seed) and K from its
# environment and must leave `prob` (T1/T2) or `cv_auc` + `probs` (T3) behind.
run_block <- function(task, fw, d_tr, d_va, s, K) {
  env <- new.env(parent = globalenv())
  env$d_tr <- d_tr; env$d_va <- d_va; env$s <- s; env$K <- K
  file <- file.path(BLOCKS, sprintf("%s_%s.R", task, fw))
  set.seed(s)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)   # swallow iCARE's plots
  t0  <- proc.time()[["elapsed"]]
  err <- tryCatch({
    utils::capture.output(suppressWarnings(suppressMessages(eval(parse(file), envir = env))))
    NULL
  }, error = function(e) conditionMessage(e))
  list(env = env, sec = proc.time()[["elapsed"]] - t0, err = err)
}

mk_row <- function(task, ds, algo, fw, s, auc = NA, brier = NA, cv = NA, sec = NA,
                   n_sel = NA, note = NA_character_, err = NA_character_)
  data.frame(task = task, dataset = ds, algo = algo, framework = fw, seed = s,
             auc = auc, brier = brier, cv_auc = cv, sec = sec, n_sel = n_sel,
             note = note, err = err, stringsAsFactors = FALSE)

# Turn the objects left by a T1/T2 block into one result row
score_simple <- function(task, ds, fw, s, d_va, r) {
  if (!is.null(r$err)) return(mk_row(task, ds, "-", fw, s, sec = r$sec, err = r$err))
  p <- as.numeric(r$env$prob)
  if (length(p) != nrow(d_va) || anyNA(p))
    return(mk_row(task, ds, "-", fw, s, sec = r$sec, err = "prob has wrong length or NA"))
  mk_row(task, ds, "-", fw, s, auc = auc_of(d_va$y, p), brier = brier_of(d_va$y, p),
         sec = r$sec, n_sel = if (is.null(r$env$sel)) NA else length(r$env$sel))
}

# Turn the objects left by a T3 block into rows: one per algorithm + "CV-picked"
score_bench <- function(ds, fw, s, d_va, r) {
  if (!is.null(r$err)) return(mk_row("T3", ds, "CV-picked", fw, s, sec = r$sec, err = r$err))
  cv <- unlist(r$env$cv_auc); names(cv) <- sub("\\..*$", "", names(cv))   # strip mlr3 suffix
  ok <- setequal(names(cv), ALGOS) && setequal(names(r$env$probs), ALGOS)
  if (!ok) return(mk_row("T3", ds, "CV-picked", fw, s, sec = r$sec,
                         err = "block did not return all three algorithms"))
  cv  <- cv[ALGOS]
  auc <- vapply(ALGOS, function(a) auc_of(d_va$y, as.numeric(r$env$probs[[a]])), numeric(1))
  bri <- vapply(ALGOS, function(a) brier_of(d_va$y, as.numeric(r$env$probs[[a]])), numeric(1))
  pick <- ALGOS[which.max(cv)]
  rbind(do.call(rbind, lapply(seq_along(ALGOS), function(i)
          mk_row("T3", ds, ALGOS[i], fw, s, auc = auc[i], brier = bri[i], cv = cv[i]))),
        mk_row("T3", ds, "CV-picked", fw, s, auc = auc[pick], brier = bri[pick],
               cv = max(cv), sec = r$sec, note = pick))
}

skip_fw <- function(task, fw) fw == "tidymodels" && task %in% c("T2", "T4") && !HAS_COLINO

## E4. Warm-up (lazy loading, JIT) so that the first timed run is not penalised ----
cat("\nWarm-up ...\n")
local({
  set.seed(1)
  ix <- as.vector(createDataPartition(D_lab$y, p = 0.7, list = FALSE))
  for (tk in c("T1", "T2", "T3")) for (fw in FW) if (!skip_fw(tk, fw))
    invisible(run_block(tk, fw, D_lab[ix, ], D_lab[-ix, ], 1, 10))
})

## E5. Main loop: T1-T3 over datasets and random splits --------------------------
res <- list(); SEL <- list()
for (ds in names(DS)) {
  cfg <- DS[[ds]]
  cat(sprintf("\n=== %s ===\n", ds))
  for (s in seq_len(R_H2H)) {
    dat <- if (ds == "Noise") make_noise(s) else cfg$data
    set.seed(s)
    tr   <- as.vector(createDataPartition(dat$y, p = 0.7, list = FALSE))
    d_tr <- dat[tr, ]; d_va <- dat[-tr, ]
    line <- sprintf("seed %2d:", s)
    for (tk in cfg$tasks) for (fw in FW) {
      if (skip_fw(tk, fw)) next
      r <- run_block(tk, fw, d_tr, d_va, s, cfg$K)
      row <- if (tk == "T3") score_bench(ds, fw, s, d_va, r) else score_simple(tk, ds, fw, s, d_va, r)
      res[[length(res) + 1]] <- row
      if (tk == "T2" && is.null(r$err) && !is.null(r$env$sel))
        SEL[[paste(ds, s, fw, sep = "|")]] <- r$env$sel
      pick <- row[row$algo %in% c("-", "CV-picked"), ]
      line <- paste0(line, sprintf(" %s/%s %.3f", tk, substr(fw, 1, 2), pick$auc[1]))
      if (!is.null(r$err)) warning(sprintf("[%s %s %s seed %d] %s", ds, tk, fw, s, r$err), call. = FALSE)
    }
    cat(line, "\n")
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(OUT, "tables", "expE_runs_long.csv"), row.names = FALSE)

# Failures are recorded, not hidden
n_err <- sum(!is.na(res$err))
cat(sprintf("\nRuns with errors: %d of %d\n", n_err, nrow(res)))
if (n_err) print(unique(res[!is.na(res$err), c("task", "dataset", "framework", "err")]))
ok <- res[is.na(res$err) & !is.na(res$auc), ]

# Backward-compatible wide table of the original comparison (T1 on UTUC_lab)
t1 <- ok[ok$task == "T1" & ok$dataset == "UTUC_lab", ]
wide <- reshape(t1[, c("seed", "framework", "auc")], idvar = "seed",
                timevar = "framework", direction = "wide")
names(wide) <- sub("^auc\\.", "", names(wide))
write.csv(wide[, c("seed", intersect(FW, names(wide)))],
          file.path(OUT, "tables", "expE_framework_head2head.csv"), row.names = FALSE)

## E6. T4 scalability: runtime of the T2 pipeline over a grid of n and p ----------
sim_data <- function(n, p, s) {
  set.seed(s)
  X <- as.data.frame(matrix(rnorm(n * p), n, p)); names(X) <- paste0("V", 1:p)
  y <- factor(rep(c("Control", "Case"), length.out = n), levels = c("Control", "Case"))
  X[, 1:10] <- X[, 1:10] + 0.6 * (y == "Case")           # 10 informative features
  X$y <- y; X[sample(n), ]
}
cat("\n=== T4 scalability (T2 pipeline, K = 10) ===\n")
scal <- list()
for (n in SCALE_N) for (p in SCALE_P) for (s in seq_len(R_SCALE)) {
  dat  <- sim_data(n, p, s)
  set.seed(s); tr <- as.vector(createDataPartition(dat$y, p = 0.7, list = FALSE))
  line <- sprintf("n %4d p %4d rep %d:", n, p, s)
  for (fw in FW) {
    if (skip_fw("T4", fw)) next
    r   <- run_block("T2", fw, dat[tr, ], dat[-tr, ], s, 10)
    auc <- if (is.null(r$err)) auc_of(dat[-tr, ]$y, as.numeric(r$env$prob)) else NA
    scal[[length(scal) + 1]] <- data.frame(n = n, p = p, rep = s, framework = fw,
                                           sec = r$sec, auc = auc,
                                           err = if (is.null(r$err)) NA_character_ else r$err)
    line <- paste0(line, sprintf(" %s %.1fs", substr(fw, 1, 2), r$sec))
  }
  cat(line, "\n")
}
scal <- do.call(rbind, scal)
write.csv(scal, file.path(OUT, "tables", "expE_scaling.csv"), row.names = FALSE)

## E7. Code volume of every block ------------------------------------------------
# executable lines = non-blank, non-comment lines; calls = function-call tokens;
# packages = attached packages that provide the distinct functions called
# (approximate: base/stats/utils/methods are not counted).
code_stats <- function(file) {
  ln <- readLines(file); ln <- ln[!grepl("^\\s*#", ln) & nzchar(trimws(ln))]
  pd <- utils::getParseData(parse(file, keep.source = TRUE))
  fn <- pd$text[pd$token == "SYMBOL_FUNCTION_CALL"]
  pk <- unique(unlist(lapply(unique(fn), function(f) {
    w <- utils::find(f); w <- sub("^package:", "", w[grepl("^package:", w)])
    setdiff(w[1], c("base", "stats", "utils", "methods", "graphics", "grDevices", NA))
  })))
  if (any(pd$token == "SPECIAL" & pd$text == "%>>%")) pk <- union(pk, "mlr3pipelines")
  data.frame(executable_lines = length(ln), function_calls = length(fn),
             distinct_functions = length(unique(fn)), n_packages = length(pk),
             packages = paste(sort(pk), collapse = "+"), stringsAsFactors = FALSE)
}
blocks <- expand.grid(task = c("T1", "T2", "T3"), framework = FW, stringsAsFactors = FALSE)
blocks$file <- file.path(BLOCKS, sprintf("%s_%s.R", blocks$task, blocks$framework))
blocks <- blocks[file.exists(blocks$file), ]
code <- cbind(blocks[, c("task", "framework")],
              do.call(rbind, lapply(blocks$file, code_stats)))
code$framework <- factor(code$framework, levels = FW)
code <- code[order(code$task, code$framework), ]
write.csv(code, file.path(OUT, "tables", "expE_code_complexity.csv"), row.names = FALSE)
write.csv(data.frame(framework = code$framework[code$task == "T1"],
                     executable_lines = code$executable_lines[code$task == "T1"]),
          file.path(OUT, "tables", "expE_loc_counts.csv"), row.names = FALSE)
print(code[, c("task", "framework", "executable_lines", "function_calls",
               "distinct_functions", "n_packages")], row.names = FALSE)

## E8. Summaries -----------------------------------------------------------------
# 8a. performance, per task x dataset x algorithm x framework
perf <- ok |>
  group_by(task, dataset, algo, framework) |>
  summarise(n = n(), AUC_mean = mean(auc), AUC_sd = sd(auc), Brier_mean = mean(brier),
            sec_median = median(sec, na.rm = TRUE), .groups = "drop")
write.csv(perf, file.path(OUT, "tables", "expE_summary_performance.csv"), row.names = FALSE)

# 8b. paired equivalence of each framework against iCARE (same splits, TOST)
eq <- do.call(rbind, lapply(split(ok, list(ok$task, ok$dataset, ok$algo), drop = TRUE), function(g) {
  ref <- g[g$framework == "icare", c("seed", "auc")]
  do.call(rbind, lapply(setdiff(FW, "icare"), function(f) {
    x  <- g[g$framework == f, c("seed", "auc")]
    mg <- merge(ref, x, by = "seed", suffixes = c(".icare", ".fw"))
    if (nrow(mg) < 3) return(NULL)
    d  <- mg$auc.fw - mg$auc.icare
    se <- sd(d) / sqrt(length(d)); q <- qt(0.95, length(d) - 1)
    data.frame(task = g$task[1], dataset = g$dataset[1], algo = g$algo[1], framework = f,
               n_pairs = length(d), mean_diff = mean(d), lo90 = mean(d) - q * se,
               hi90 = mean(d) + q * se, max_abs_diff = max(abs(d)),
               equivalent = (mean(d) - q * se > -DELTA) & (mean(d) + q * se < DELTA))
  }))
}))
rownames(eq) <- NULL
write.csv(eq, file.path(OUT, "tables", "expE_equivalence_vs_icare.csv"), row.names = FALSE)
cat(sprintf("\nEquivalence vs iCARE (TOST, 90%% CI of paired AUC difference within +/-%.3f): %d of %d comparisons equivalent\n",
            DELTA, sum(eq$equivalent), nrow(eq)))
print(eq[!eq$equivalent, ], row.names = FALSE)

# 8c. T2: do the frameworks keep the same features?
sel_tab <- do.call(rbind, lapply(names(DS), function(ds) do.call(rbind, lapply(seq_len(R_H2H), function(s) {
  ref <- SEL[[paste(ds, s, "caret", sep = "|")]]
  if (is.null(ref)) return(NULL)
  do.call(rbind, lapply(setdiff(FW, "caret"), function(f) {
    o <- SEL[[paste(ds, s, f, sep = "|")]]
    if (is.null(o)) return(NULL)
    data.frame(dataset = ds, seed = s, framework = f, n_sel = length(o),
               jaccard_vs_caret = jaccard(ref, o), identical = setequal(ref, o))
  }))
}))))
if (!is.null(sel_tab)) {
  sel_sum <- sel_tab |> group_by(dataset, framework) |>
    summarise(n_sel = mean(n_sel), mean_jaccard = mean(jaccard_vs_caret),
              identical_sets = mean(identical), .groups = "drop")
  write.csv(sel_sum, file.path(OUT, "tables", "expE_T2_selection_overlap.csv"), row.names = FALSE)
  cat("\nT2 feature-set agreement with caret:\n"); print(as.data.frame(sel_sum), row.names = FALSE)
}

# 8d. T3: which algorithm does each framework pick, and does it agree with iCARE?
pk <- ok[ok$task == "T3" & ok$algo == "CV-picked", c("dataset", "seed", "framework", "note", "cv_auc")]
if (nrow(pk)) {
  pk_sum <- do.call(rbind, lapply(split(pk, pk$dataset), function(g) {
    ref <- g[g$framework == "icare", c("seed", "note")]
    do.call(rbind, lapply(setdiff(FW, "icare"), function(f) {
      mg <- merge(ref, g[g$framework == f, c("seed", "note")], by = "seed")
      if (!nrow(mg)) return(NULL)
      data.frame(dataset = g$dataset[1], framework = f, same_pick_as_icare = mean(mg$note.x == mg$note.y))
    }))
  }))
  pick_freq <- as.data.frame(table(dataset = pk$dataset, framework = pk$framework, picked = pk$note))
  write.csv(pk_sum, file.path(OUT, "tables", "expE_T3_pick_agreement.csv"), row.names = FALSE)
  write.csv(pick_freq, file.path(OUT, "tables", "expE_T3_pick_frequency.csv"), row.names = FALSE)
  cat("\nT3 agreement of the CV-selected algorithm with iCARE:\n"); print(pk_sum, row.names = FALSE)
}

# 8e. Runtime summary (T1-T3, median seconds) and T4 scaling
rt <- ok[ok$algo %in% c("-", "CV-picked"), ] |>
  group_by(task, dataset, framework) |>
  summarise(median_sec = median(sec), iqr_sec = IQR(sec), .groups = "drop")
write.csv(rt, file.path(OUT, "tables", "expE_runtime_summary.csv"), row.names = FALSE)
scal_sum <- scal[is.na(scal$err), ] |> group_by(n, p, framework) |>
  summarise(median_sec = median(sec), mean_auc = mean(auc), .groups = "drop")
write.csv(scal_sum, file.path(OUT, "tables", "expE_scaling_summary.csv"), row.names = FALSE)

## E9. Figures -------------------------------------------------------------------
ok$framework   <- factor(ok$framework, levels = FW)
ok$ds_label    <- factor(ok$dataset, levels = names(DS))
th <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                       strip.background = element_rect(fill = "grey95"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, "figures", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, "figures", paste0(name, ".pdf")), p, width = w, height = h)
}

# FigE1: held-out AUC of every framework (T1, T2, T3 = CV-picked algorithm)
d1 <- ok[ok$algo %in% c("-", "CV-picked"), ]
d1$panel <- factor(d1$task, levels = c("T1", "T2", "T3"),
                   labels = c("T1: scale + glm", "T2: scale + screen + glm", "T3: CV-picked of 3 algorithms"))
p1 <- ggplot(d1, aes(framework, auc, colour = framework)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = NA) +
  geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.6) +
  facet_grid(panel ~ ds_label, scales = "free_y") +
  scale_colour_manual(values = FW_COL, guide = "none") +
  labs(x = NULL, y = "Held-out AUC (direction fixed)",
       title = "Held-out AUC of the four frameworks on the same splits",
       subtitle = sprintf("%d random 70/30 splits per dataset; Noise panel: true AUC = 0.5", R_H2H)) + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p1, "FigE1_auc_by_framework", 11, 7.5)

# FigE2: paired difference to iCARE with the equivalence region
eq$label <- paste(eq$task, eq$dataset, ifelse(eq$algo == "-", "", eq$algo), sep = " | ")
eq$label <- factor(sub(" \\| $", "", eq$label), levels = rev(unique(sub(" \\| $", "", eq$label))))
eq$framework <- factor(eq$framework, levels = FW)
p2 <- ggplot(eq, aes(mean_diff, label, colour = framework)) +
  annotate("rect", xmin = -DELTA, xmax = DELTA, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = lo90, xmax = hi90), position = position_dodge(width = 0.7), size = 0.3) +
  scale_colour_manual(values = FW_COL, name = "vs iCARE") +
  labs(x = "Mean paired AUC difference to iCARE (90% CI)", y = NULL,
       title = "Equivalence to iCARE",
       subtitle = sprintf("grey band = equivalence margin +/-%.3f AUC; CI inside the band = equivalent (TOST)", DELTA)) + th
save_fig(p2, "FigE2_equivalence_vs_icare", 8.5, max(5, 0.32 * nlevels(eq$label) + 1.5))

# FigE3: runtime
d3 <- d1; d3$panel <- d1$panel
p3 <- ggplot(d3, aes(framework, sec, colour = framework)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = NA) +
  geom_jitter(width = 0.12, height = 0, size = 0.8, alpha = 0.5) +
  scale_y_log10() + facet_grid(panel ~ ds_label, scales = "free_y") +
  scale_colour_manual(values = FW_COL, guide = "none") +
  labs(x = NULL, y = "Wall-clock time per run (s, log scale)",
       title = "Runtime per run (one split, training + prediction)",
       subtitle = "iCARE always runs 5-fold CV inside ModelTrainAnalysis(); caret/tidymodels/mlr3 in T1/T2 fit once") + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p3, "FigE3_runtime", 11, 7.5)

# FigE4: code volume
cd <- code; cd$task <- factor(cd$task, levels = c("T1", "T2", "T3"),
                              labels = c("T1: scale + glm", "T2: scale + screen + glm", "T3: 3-algorithm CV benchmark"))
cd_long <- rbind(data.frame(cd[, c("task", "framework")], metric = "Executable lines", value = cd$executable_lines),
                 data.frame(cd[, c("task", "framework")], metric = "Function calls", value = cd$function_calls))
p4 <- ggplot(cd_long, aes(framework, value, fill = framework)) +
  geom_col(width = 0.7) + geom_text(aes(label = value), vjust = -0.3, size = 3) +
  facet_grid(metric ~ task, scales = "free_y") +
  scale_fill_manual(values = FW_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL, title = "Code volume of each framework's implementation",
       subtitle = "Counted from the block files that were executed (blank lines and comments excluded)") + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p4, "FigE4_code_volume", 10, 5.5)

# FigE5: scaling
ss <- scal_sum; ss$framework <- factor(ss$framework, levels = FW)
p5 <- ggplot(ss, aes(n, median_sec, colour = framework, group = framework)) +
  geom_line() + geom_point(size = 1.8) +
  scale_x_log10(breaks = SCALE_N) + scale_y_log10() +
  facet_wrap(~ p, labeller = label_both, nrow = 1) +
  scale_colour_manual(values = FW_COL, name = NULL) +
  labs(x = "Number of samples n (log scale)", y = "Median run time (s, log scale)",
       title = "Scalability of the T2 pipeline (scale + top-10 screen + glm)",
       subtitle = sprintf("simulated data, 10 informative features, %d repetitions per cell", R_SCALE)) + th
save_fig(p5, "FigE5_scaling", 10, 3.8)

## E10. Headline numbers for the manuscript ---------------------------------------
cat("\n================ HEADLINE ================\n")
t1max <- max(abs(t1$auc[t1$framework == "icare"] -
                   t1$auc[t1$framework == "caret"]), na.rm = TRUE)
cat(sprintf("T1 (UTUC_lab): max |AUC(iCARE) - AUC(caret)| over splits = %.2e\n", t1max))
nz <- perf[perf$dataset == "Noise" & perf$algo %in% c("-", "CV-picked"), ]
cat("Noise data (true AUC 0.5), mean held-out AUC (all frameworks must stay near 0.5):\n")
print(as.data.frame(nz[, c("task", "framework", "AUC_mean", "AUC_sd")]), row.names = FALSE)
cat("\nMedian runtime (s):\n"); print(as.data.frame(rt), row.names = FALSE)

writeLines(c(capture.output(sessionInfo()), "",
             sprintf("icare %s | caret %s | tidymodels %s | mlr3 %s | mlr3pipelines %s | colino %s",
                     packageVersion("icare"), packageVersion("caret"), packageVersion("tidymodels"),
                     packageVersion("mlr3"), packageVersion("mlr3pipelines"),
                     if (HAS_COLINO) packageVersion("colino") else "not installed")),
           file.path(OUT, "expE_sessionInfo.txt"))
cat("\nDone. Tables: ", file.path(OUT, "tables"), "  Figures: ", file.path(OUT, "figures"), "\n", sep = "")
