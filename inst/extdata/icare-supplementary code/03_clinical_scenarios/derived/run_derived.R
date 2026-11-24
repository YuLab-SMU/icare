#!/usr/bin/env Rscript
# ==============================================================================
# 03_clinical_scenarios/derived/run_derived.R
# ------------------------------------------------------------------------------
# The three scenario folders produce one fitted model each (03_model.R). The
# scripts in this folder read those artefacts and add the analyses reported in
# the secondary tables of the paper. Nothing is refitted, so every number here
# is consistent with the main results by construction.
#
#   01_calibration.R          external calibration: CITL and slope
#   02_nri_idi_external.R     NRI / IDI and the external confusion matrices
#   03_threshold_positions.R  where the Youden threshold sits in the risk distribution
#   04_auc_gain_delong.R      paired DeLong test for the AUC gain (Fig. S13)
#   05_composite_figures.R    Fig. 2 / Fig. 3 composite panels
#
#   Rscript 03_clinical_scenarios/derived/run_derived.R
#
# Requires the scenario outputs, i.e. run 03_clinical_scenarios/run_scenarios.R first.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

need <- file.path(ROOT, "03_clinical_scenarios",
                  c("01_undiagnosed_dysglycemia", "02_albuminuria", "03_undiagnosed_hypertension"),
                  "output", "model.rds")
if (!all(file.exists(need)))
  stop("Scenario outputs missing:\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "),
       "\n  run 03_clinical_scenarios/run_scenarios.R first.", call. = FALSE)

steps <- c("01_calibration.R", "02_nri_idi_external.R", "03_threshold_positions.R",
           "04_auc_gain_delong.R", "05_composite_figures.R")
here_dir <- file.path(ROOT, "03_clinical_scenarios", "derived")
rscript  <- file.path(R.home("bin"), "Rscript")

t0 <- Sys.time()
for (s in steps) {
  cat(sprintf("\n=== %s ===\n", s))
  status <- system2(rscript, shQuote(file.path(here_dir, s)))
  if (status != 0) stop("Step failed: ", s, " (exit ", status, ")", call. = FALSE)
}
cat(sprintf("\nDerived analyses done in %.1f min - outputs in %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            out_dir("03_clinical_scenarios", create = FALSE)))
