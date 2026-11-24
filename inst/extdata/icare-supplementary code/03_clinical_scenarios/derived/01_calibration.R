# ==============================================================================
# External calibration for the three scenarios (derived from saved fit objects)
# ------------------------------------------------------------------------------
# Uses the already-fitted ModelDeployment in each folder's output/model.rds and
# the held-out external cohort in output/cohort.rds. No model is re-trained, so
# predictions are identical to those in the main results.
#
# For each scenario, writes <scenario>/output/calibration.csv with:
#   observed_event_rate, mean_predicted_risk, calibration_in_the_large (CITL),
#   calibration_slope (slope of logit(p) ~ event; 1 = perfectly calibrated).
# Run: Rscript 03_clinical_scenarios/derived/01_calibration.R   (after the three
#      scenarios have been run; reads <scenario>/output/{cohort,model}.rds)
#      Output: outputs/03_clinical_scenarios/tables/scenario_calibration.csv
# ==============================================================================

suppressMessages({
  library(icare)
  library(caret)
  library(pROC)
})

## --- Locate the repository ----------------------------------------------------
if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

SCEN <- file.path(ROOT, "03_clinical_scenarios")   # the three scenario folders
TAB  <- out_dir("03_clinical_scenarios", "tables")


FOLDERS <- c("01_undiagnosed_dysglycemia", "02_albuminuria",
             "03_undiagnosed_hypertension")
POS <- "Positive"; NEG <- "Negative"; G <- "group"

out <- list()
for (folder in FOLDERS) {
  out_dir <- file.path(SCEN, folder, "output")
  co <- readRDS(file.path(out_dir, "cohort.rds"))
  mo <- readRDS(file.path(out_dir, "model.rds"))
  ext <- co$ext
  features <- mo$features

  p <- mo$deploy$predict_fn(
    newdata = ext[, features, drop = FALSE],
    selected_model = mo$cv_winner)[[POS]]
  y <- as.integer(ext[[G]] == POS)

  # Calibration-in-the-large and logistic calibration slope on the external cycle
  slope <- coef(glm(
    y ~ qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6)),
    family = binomial))[2]

  row <- data.frame(
    scenario = folder,
    observed_event_rate = mean(y),
    mean_predicted_risk = mean(p),
    calibration_in_the_large = mean(p) - mean(y),
    calibration_slope = as.numeric(slope)
  )
  write.csv(row, file.path(out_dir, "calibration.csv"), row.names = FALSE)
  out[[folder]] <- row
  cat(sprintf("%-34s CITL %+.4f | slope %.3f\n", folder,
              row$calibration_in_the_large, row$calibration_slope))
}

all_cal <- do.call(rbind, out)
write.csv(all_cal, file.path(TAB, "scenario_calibration.csv"),
          row.names = FALSE)
cat("\nSaved output/calibration.csv for each scenario and tables/scenario_calibration.csv\n")
