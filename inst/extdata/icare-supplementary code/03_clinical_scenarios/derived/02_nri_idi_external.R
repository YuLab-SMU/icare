# ==============================================================================
# External-validation NRI / IDI and confusion matrices for the three scenarios
# ------------------------------------------------------------------------------
# Question: on the INDEPENDENT 2009-2010 cycle, how much does the full routine
# laboratory panel add over a demographics-only reference model (age + sex)?
#
#   Reference model : glm(Positive ~ age + sex), fitted on the development
#                     cohort only (training partition), no leakage.
#   New model       : the already-fitted iCARE model saved in output/model.rds
#                     (full laboratory / CBC panel + demographics).
#   Evaluation      : the whole 2009-2010 cycle, never used for fitting.
#
# Uses the official icare NRI/IDI toolchain:
#   NRI_IDI_Analysis() / CalculateCategoryNRI() / PlotIDICurve()
# Continuous NRI, categorical NRI (risk thresholds), and IDI are reported with
# bootstrap 95% CIs. Confusion matrices use the Youden threshold fixed on the
# internal held-out set in 03_model.R (nothing re-estimated on the external set).
#
# Run: Rscript 05_nri_idi_external.R  (from any working directory)
# Output: <scenario>/output/nri_idi_results.csv
#         <scenario>/output/confusion_matrix.csv
#         <scenario>/output/nri_idi/  (official icare NRI/IDI plots)
#         figures/ (confusion-matrix and IDI plot objects for the composite)
# ==============================================================================

suppressMessages({
  library(icare)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(dplyr)
})

## --- Locate the repository ----------------------------------------------------
if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

SCEN <- file.path(ROOT, "03_clinical_scenarios")   # the three scenario folders
TAB  <- out_dir("03_clinical_scenarios", "tables")
CACHE <- out_dir("03_clinical_scenarios", "cache")

set.seed(2026)

FOLDERS <- c("01_undiagnosed_dysglycemia", "02_albuminuria",
             "03_undiagnosed_hypertension")
SHORT <- c("Undiagnosed dysglycemia", "Albuminuria", "Undiagnosed hypertension")
POS <- "Positive"; NEG <- "Negative"; G <- "group"
RISK_THRESHOLDS <- c(0.02, 0.10, 0.50, 0.95)   # clinical risk categories for categorical NRI
BOOT <- 400

# ---- helpers -----------------------------------------------------------------
auc_ci <- function(y, p) {
  r <- pROC::roc(y, p, levels = c(NEG, POS), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  c(as.numeric(pROC::auc(r)), ci[1], ci[3])
}

confusion <- function(y, p, threshold) {
  pred <- factor(ifelse(p >= threshold, POS, NEG), levels = c(NEG, POS))
  cm <- table(predicted = pred, actual = y)
  tp <- cm[POS, POS]; tn <- cm[NEG, NEG]
  fp <- cm[POS, NEG]; fn <- cm[NEG, POS]
  list(tp = as.integer(tp), tn = as.integer(tn),
       fp = as.integer(fp), fn = as.integer(fn),
       sensitivity = as.numeric(tp / (tp + fn)),
       specificity = as.numeric(tn / (tn + fp)),
       ppv = as.numeric(tp / (tp + fp)),
       npv = as.numeric(tn / (tn + fn)))
}

results <- list()
cm_results <- list()
plot_objs <- list()

for (i in seq_along(FOLDERS)) {
  folder <- FOLDERS[i]
  label  <- SHORT[i]
  out_dir <- file.path(SCEN, folder, "output")
  nri_dir <- file.path(out_dir, "nri_idi")
  dir.create(nri_dir, recursive = TRUE, showWarnings = FALSE)

  co <- readRDS(file.path(out_dir, "cohort.rds"))
  mo <- readRDS(file.path(out_dir, "model.rds"))
  dev <- co$dev
  ext <- co$ext
  FEATURES <- mo$features

  # ---- full iCARE model probabilities on the external cycle -----------------
  p_new <- mo$deploy$predict_fn(
    newdata = ext[, FEATURES, drop = FALSE],
    selected_model = mo$cv_winner)[[POS]]

  # ---- demographics-only reference model, fitted on development data --------
  # Reproduce the same internal 70/30 split used in 03_model.R for consistency
  set.seed(2026)
  idx <- caret::createDataPartition(dev[[G]], p = 0.70, list = FALSE)
  tr <- dev[idx, ]
  ref_formula <- stats::as.formula(sprintf("%s ~ RIDAGEYR + RIAGENDR", G))
  ref_fit <- stats::glm(ref_formula, data = tr, family = stats::binomial())
  p_ref <- as.numeric(stats::predict(ref_fit, newdata = ext, type = "response"))

  y <- ext[[G]]

  # ---- official icare NRI/IDI analysis (writes its 6 plots) ------------------
  set.seed(2026)
  nri_analysis <- tryCatch(
    icare::NRI_IDI_Analysis(
      ref_prob = p_ref, new_prob = p_new, truth = y, positive = POS,
      labels = c("Demographics only", "Full lab panel"),
      outcome_labels = c("Negative", "Positive"),
      risk_thresholds = RISK_THRESHOLDS,
      show_ci = TRUE, n_boot = BOOT, conf_int = TRUE, conf_level = 0.95,
      save_plots = TRUE, save_dir = nri_dir, seed = 2026),
    error = function(e) { message("NRI_IDI_Analysis failed: ", conditionMessage(e)); NULL })

  # ---- direct category NRI (to capture the detailed result object) -----------
  set.seed(2026)
  cat_nri <- icare::CalculateCategoryNRI(
    truth = y, ref_prob = p_ref, new_prob = p_new,
    risk_thresholds = RISK_THRESHOLDS,
    conf_int = TRUE, n_boot = BOOT, conf_level = 0.95, seed = 2026)

  # ---- continuous NRI and IDI (self-contained, transparent formulas) ---------
  # Continuous NRI: using the same threshold scale for both models.
  events <- which(y == POS); nonevents <- which(y == NEG)
  nri_ev <- mean(p_new[events] > p_ref[events]) - mean(p_new[events] < p_ref[events])
  nri_ne <- mean(p_new[nonevents] < p_ref[nonevents]) -
            mean(p_new[nonevents] > p_ref[nonevents])
  nri_cont <- nri_ev + nri_ne
  idi <- (mean(p_new[events]) - mean(p_ref[events])) -
         (mean(p_new[nonevents]) - mean(p_ref[nonevents]))

  # ---- discrimination --------------------------------------------------------
  auc_ref <- auc_ci(y, p_ref)
  auc_new <- auc_ci(y, p_new)

  # ---- confusion matrix at the internal Youden threshold --------------------
  thr <- mo$youden_threshold
  cm <- confusion(y, p_new, thr)
  cm_df <- data.frame(
    scenario = label, threshold = thr,
    true_negative = cm$tn, false_positive = cm$fp,
    false_negative = cm$fn, true_positive = cm$tp,
    sensitivity = cm$sensitivity, specificity = cm$specificity,
    ppv = cm$ppv, npv = cm$npv, external_n = nrow(ext),
    external_positive = sum(y == POS))
  write.csv(cm_df, file.path(out_dir, "confusion_matrix.csv"), row.names = FALSE)
  cm_results[[folder]] <- cm_df

  # ---- assemble result row ---------------------------------------------------
  # .boot_ci returns a data.frame: metric, estimate, lower, upper
  get_ci <- function(ci_df, key) {
    if (is.null(ci_df) || !is.data.frame(ci_df)) return(c(NA_real_, NA_real_))
    row <- ci_df[ci_df$metric == key, ]
    if (!nrow(row)) return(c(NA_real_, NA_real_))
    c(row$lower[1], row$upper[1])
  }
  nri_tot_lohi <- get_ci(cat_nri$ci, "nri_total")
  idi_ci <- if (!is.null(nri_analysis)) nri_analysis$idi$ci else NULL
  # PlotIDICurve labels the metric "IDI" (upper case) in its CI data.frame
  idi_lohi <- if (!is.null(idi_ci)) {
    out <- get_ci(idi_ci, "IDI")
    if (anyNA(out)) out <- get_ci(idi_ci, "idi")
    out
  } else c(NA_real_, NA_real_)
  # Continuous-NRI CI (computed the same way as icare's category NRI bootstrap)
  nri_cont_ci <- get_ci(cat_nri$ci, "nri_total")
  if (!is.null(nri_analysis) && !is.null(nri_analysis$nri)) {
    cont_ci <- get_ci(nri_analysis$nri$ci, "nri_total")
    if (!anyNA(cont_ci)) nri_cont_ci <- cont_ci
  }

  results[[folder]] <- data.frame(
    scenario = label,
    external_n = nrow(ext),
    external_positive = sum(y == POS),
    auc_demographics = auc_ref[1], auc_demo_low = auc_ref[2], auc_demo_high = auc_ref[3],
    auc_full_panel = auc_new[1], auc_full_low = auc_new[2], auc_full_high = auc_new[3],
    auc_gain = auc_new[1] - auc_ref[1],
    nri_continuous = nri_cont,
    nri_events = nri_ev, nri_nonevents = nri_ne,
    nri_category_total = cat_nri$nri_total,
    nri_category_low = nri_tot_lohi[1], nri_category_high = nri_tot_lohi[2],
    idi = idi,
    idi_low = idi_lohi[1], idi_high = idi_lohi[2],
    youden_threshold = thr,
    sensitivity = cm$sensitivity, specificity = cm$specificity,
    ppv = cm$ppv, npv = cm$npv,
    stringsAsFactors = FALSE)

  # ---- build ggplot objects for the composite figure --------------------------
  # Confusion matrix as a labelled 2x2 tile plot
  cm_long <- expand.grid(
    actual = c(NEG, POS), predicted = c(NEG, POS),
    stringsAsFactors = FALSE)
  cm_long$n <- c(cm$tn, cm$fn, cm$fp, cm$tp)
  cm_plot <- ggplot(cm_long, aes(x = actual, y = predicted, fill = n)) +
    geom_tile(color = "grey30", linewidth = 0.6) +
    geom_text(aes(label = n), size = 6, colour = ifelse(cm_long$n > 400, "white", "grey15")) +
    scale_fill_gradient(low = "#EAF3DE", high = "#2E8B57", name = "n") +
    labs(title = label,
         subtitle = sprintf("Youden thr = %.3f | Sens %.3f | Spec %.3f",
                            thr, cm$sensitivity, cm$specificity),
         x = "Actual", y = "Predicted") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey30"))
  plot_objs[[paste0("cm_", folder)]] <- cm_plot

  # IDI curve from the official icare plot object (if available)
  if (!is.null(nri_analysis)) {
    idi_plot <- nri_analysis$idi$plot +
      labs(title = sprintf("IDI = %+.4f", nri_analysis$idi$idi)) +
      theme(plot.title = element_text(face = "bold", size = 11))
    plot_objs[[paste0("idi_", folder)]] <- idi_plot
  }

  # Category-NRI transition bar plot (up/down movement by cases and non-cases)
  nri_bar <- data.frame(
    group = factor(c("Cases\nup", "Cases\ndown", "Non-cases\nup", "Non-cases\ndown"),
                   levels = c("Cases\nup", "Cases\ndown", "Non-cases\nup", "Non-cases\ndown")),
    frac = c(cat_nri$events_up / max(1, length(events)),
             -cat_nri$events_down / max(1, length(events)),
             cat_nri$nonevents_up / max(1, length(nonevents)),
             -cat_nri$nonevents_down / max(1, length(nonevents))))
  nri_bar$grp <- ifelse(grepl("Cases", nri_bar$group), "Cases", "Non-cases")
  bar_plot <- ggplot(nri_bar, aes(x = group, y = frac, fill = grp)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    scale_fill_manual(values = c(Cases = "#D94841", `Non-cases` = "#2C7FB8")) +
    scale_y_continuous(labels = scales::percent_format(),
                       limits = c(-max(abs(nri_bar$frac)) * 1.2,
                                   max(abs(nri_bar$frac)) * 1.2)) +
    labs(title = sprintf("Category NRI = %+.3f", cat_nri$nri_total),
         subtitle = "Reclassification across clinical risk categories",
         x = NULL, y = "Fraction of group") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(size = 8),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
  plot_objs[[paste0("nri_", folder)]] <- bar_plot

  # AUC comparison bar (demographics vs full panel)
  auc_bar <- data.frame(
    model = factor(c("Demographics\nonly", "Full lab\npanel"),
                   levels = c("Demographics\nonly", "Full lab\npanel")),
    auc = c(auc_ref[1], auc_new[1]),
    lo = c(auc_ref[2], auc_new[2]), hi = c(auc_ref[3], auc_new[3]))
  auc_plot <- ggplot(auc_bar, aes(x = model, y = auc, colour = model)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.7) +
    geom_point(size = 3.2, show.legend = FALSE) +
    annotate("segment", x = 1, xend = 2, y = auc_ref[1], yend = auc_new[1],
             linetype = 3, linewidth = 0.5, colour = "grey60") +
    scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
    scale_colour_manual(values = c("#9AA5B1", "#2E8B57"), guide = "none") +
    labs(title = sprintf("AUC gain = %+.3f", auc_new[1] - auc_ref[1]),
         subtitle = "External cycle (2009-2010), DeLong 95% CI",
         x = NULL, y = "AUC") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(size = 8),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
  plot_objs[[paste0("auc_", folder)]] <- auc_plot

  cat(sprintf("%-34s AUC demo %.3f -> full %.3f (+%.3f) | NRI_cat %+.4f | IDI %+.4f | CM TP=%d FN=%d TN=%d FP=%d\n",
              folder, auc_ref[1], auc_new[1], auc_new[1] - auc_ref[1],
              cat_nri$nri_total, idi, cm$tp, cm$fn, cm$tn, cm$fp))
}

res_df <- do.call(rbind, results)
write.csv(res_df, file.path(TAB, "nri_idi_external.csv"),
          row.names = FALSE)
# Dedicated confusion-matrix table (one row per scenario, 2x2 + metrics)
cm_all <- do.call(rbind, cm_results)
write.csv(cm_all, file.path(TAB, "external_confusion_matrices.csv"),
          row.names = FALSE)
saveRDS(plot_objs, file.path(CACHE, "nri_idi_plot_objects.rds"))
cat("\nSaved NRI/IDI results, confusion matrices and plot objects.\n")
print(res_df[, c("scenario", "auc_demographics", "auc_full_panel", "auc_gain",
                 "nri_category_total", "idi")], row.names = FALSE)
