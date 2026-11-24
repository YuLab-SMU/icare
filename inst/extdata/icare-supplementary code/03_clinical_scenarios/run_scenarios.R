#!/usr/bin/env Rscript
# ==============================================================================
# 03_clinical_scenarios/run_scenarios.R - run the three NHANES scenarios
# ------------------------------------------------------------------------------
# Each scenario is a folder with the same three steps:
#   01_download.R  download the NHANES tables needed (skipped when data/ is cached)
#   02_cohort.R    merge, apply the eligibility rules, build the development and
#                  the independent external cohort
#   03_model.R     train with iCARE, validate on the other cycle
#
#   Rscript 03_clinical_scenarios/run_scenarios.R                 all three
#   Rscript 03_clinical_scenarios/run_scenarios.R 02_albuminuria  just one
#
# Steps run in separate processes on purpose: 03_model.R is memory hungry and a
# failure in one scenario must not take the others down.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

args <- commandArgs(trailingOnly = TRUE)
root <- file.path(ROOT, "03_clinical_scenarios")

dirs <- list.dirs(root, recursive = FALSE, full.names = FALSE)
dirs <- sort(dirs[grepl("^[0-9]{2}_", dirs)])
if (length(args)) {
  unknown <- setdiff(args, dirs)
  if (length(unknown))
    stop("Unknown scenario(s): ", paste(unknown, collapse = ", "),
         "\n  available: ", paste(dirs, collapse = ", "), call. = FALSE)
  dirs <- args
}

rscript <- file.path(R.home("bin"), "Rscript")
steps   <- c("01_download.R", "02_cohort.R", "03_model.R")

for (d in dirs) {
  cat(sprintf("\n================ %s ================\n", d))
  for (step in steps) {
    f <- file.path(root, d, step)
    if (!file.exists(f)) { warning("missing step: ", f, call. = FALSE); next }
    cat(sprintf("[%s] %s\n", d, step))
    status <- system2(rscript, shQuote(f))
    if (status != 0)
      stop("Step failed: ", d, "/", step, " (exit ", status, ")", call. = FALSE)
  }
}

res <- file.path(root, dirs, "output", "results.csv")
res <- res[file.exists(res)]
if (length(res)) {
  cat("\n=== headline results ===\n")
  for (f in res) {
    x <- utils::read.csv(f, check.names = FALSE)
    cat(sprintf("%-55s internal AUC %.3f (%.3f-%.3f)  external AUC %.3f (%.3f-%.3f)\n",
                basename(dirname(dirname(f))),
                x$internal_auc, x$internal_low, x$internal_high,
                x$external_auc, x$external_low, x$external_high))
  }
}
cat("\nDone.\n")
