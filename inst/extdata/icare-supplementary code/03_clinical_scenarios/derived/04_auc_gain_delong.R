# ==============================================================================
# Paired DeLong test for the AUC gain shown in Figure S13
# ------------------------------------------------------------------------------
# Figure S13 compares, on the independent 2009-2010 cycle, a demographics-only
# reference model (age + sex, fitted on the 2007-2008 training partition) with
# the full routine laboratory model. The figure reports AUC with DeLong 95% CIs
# but not the test for the DIFFERENCE. This script adds it, plus a bootstrap
# check, so the legend/text can quote a p value that matches the figure.
#
# Nothing is refitted here beyond the reference glm and the internal 70/30 split,
# both reproduced with the same seed (2026) as 03_model.R / 05_nri_idi_external.R.
#
# Run    : Rscript 03_clinical_scenarios/derived/04_auc_gain_delong.R
# Output : outputs/03_clinical_scenarios/tables/auc_gain_delong.csv
# ==============================================================================

suppressMessages({
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
LABELS  <- c("Undiagnosed dysglycemia", "Albuminuria", "Undiagnosed hypertension")
POS <- "Positive"; NEG <- "Negative"; G <- "group"

rows <- list()
for (i in seq_along(FOLDERS)) {
  co <- readRDS(file.path(SCEN, FOLDERS[i], "output", "cohort.rds"))
  mo <- readRDS(file.path(SCEN, FOLDERS[i], "output", "model.rds"))
  dev <- co$dev; ext <- co$ext; FEATURES <- mo$features

  p_new <- mo$deploy$predict_fn(newdata = ext[, FEATURES, drop = FALSE],
                                selected_model = mo$cv_winner)[[POS]]

  set.seed(2026)
  idx <- caret::createDataPartition(dev[[G]], p = 0.70, list = FALSE)
  tr  <- dev[idx, ]
  ref_fit <- glm(as.formula(sprintf("%s ~ RIDAGEYR + RIAGENDR", G)),
                 data = tr, family = binomial())
  p_ref <- as.numeric(predict(ref_fit, newdata = ext, type = "response"))

  y <- ext[[G]]
  r_new <- roc(y, p_new, levels = c(NEG, POS), direction = "<", quiet = TRUE)
  r_ref <- roc(y, p_ref, levels = c(NEG, POS), direction = "<", quiet = TRUE)

  ci_new <- as.numeric(ci.auc(r_new, method = "delong"))
  ci_ref <- as.numeric(ci.auc(r_ref, method = "delong"))
  dt <- roc.test(r_new, r_ref, method = "delong", paired = TRUE)

  # bootstrap the gain itself (participant-level resampling, the same pairs)
  set.seed(2026)
  yy <- as.integer(y == POS); B <- 1000
  gain <- replicate(B, {
    s <- sample.int(length(yy), replace = TRUE)
    if (length(unique(yy[s])) < 2) return(NA_real_)
    a1 <- as.numeric(auc(roc(yy[s], p_new[s], levels = c(0, 1),
                            direction = "<", quiet = TRUE)))
    a0 <- as.numeric(auc(roc(yy[s], p_ref[s], levels = c(0, 1),
                            direction = "<", quiet = TRUE)))
    a1 - a0
  })
  gain_ci <- quantile(gain, c(.025, .975), na.rm = TRUE)

  cat(sprintf("\n== %s  n=%d (events %d)\n", LABELS[i], length(y), sum(y == POS)))
  cat(sprintf("   demographics only  AUC %.3f (%.3f-%.3f)\n",
              ci_ref[2], ci_ref[1], ci_ref[3]))
  cat(sprintf("   full lab panel     AUC %.3f (%.3f-%.3f)\n",
              ci_new[2], ci_new[1], ci_new[3]))
  cat(sprintf("   gain %+.3f  DeLong paired p = %.3g  |  bootstrap 95%% CI %+.3f to %+.3f\n",
              ci_new[2] - ci_ref[2], dt$p.value, gain_ci[1], gain_ci[2]))

  rows[[i]] <- data.frame(
    scenario = LABELS[i], external_n = length(y), events = sum(y == POS),
    auc_demographics = ci_ref[2], demo_low = ci_ref[1], demo_high = ci_ref[3],
    auc_full_panel = ci_new[2], full_low = ci_new[1], full_high = ci_new[3],
    auc_gain = ci_new[2] - ci_ref[2],
    gain_boot_low = as.numeric(gain_ci[1]), gain_boot_high = as.numeric(gain_ci[2]),
    delong_p = dt$p.value)
}

out <- do.call(rbind, rows)
write.csv(out, file.path(TAB, "auc_gain_delong.csv"), row.names = FALSE)
cat("\nWritten: ", file.path(TAB, "auc_gain_delong.csv"), "\n", sep = "")
print(out[, c("scenario", "auc_demographics", "auc_full_panel", "auc_gain",
              "delong_p")], digits = 4, row.names = FALSE)
