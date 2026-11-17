# ==============================================================================
# Composite figure: external-validation confusion matrices + NRI/IDI + AUC gain
# ------------------------------------------------------------------------------
# Two publication-ready composite panels (fully English):
#
#   Fig. 2  3 scenarios x 4 columns
#       a) confusion matrix at the internal Youden threshold
#       b) external ROC: full lab panel vs demographics-only reference
#       c) AUC gain (DeLong 95% CI)
#       d) categorical NRI reclassification (cases / non-cases, up / down)
#
#   Fig. 3  1 x 3 official icare IDI curves (one per scenario)
#
# All evaluation on the independent 2009-2010 NHANES cycle; thresholds fixed on
# the internal held-out set in 03_model.R; nothing re-estimated here.
#
# Run: Rscript 03_clinical_scenarios/derived/05_composite_figures.R
#      (after 02_nri_idi_external.R, which saves the plot objects it reads)
# ==============================================================================

suppressMessages({
  library(icare)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

## --- Locate the repository ----------------------------------------------------
if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

SCEN <- file.path(ROOT, "03_clinical_scenarios")   # the three scenario folders
TAB  <- out_dir("03_clinical_scenarios", "tables")
FIG  <- out_dir("03_clinical_scenarios", "figures")
CACHE <- out_dir("03_clinical_scenarios", "cache")

FOLDERS <- c("01_undiagnosed_dysglycemia", "02_albuminuria",
             "03_undiagnosed_hypertension")
SHORT <- c("Undiagnosed dysglycemia", "Albuminuria", "Undiagnosed hypertension")
POS <- "Positive"; NEG <- "Negative"; G <- "group"

plots <- readRDS(file.path(CACHE, "nri_idi_plot_objects.rds"))
nri_df <- read.csv(file.path(TAB, "nri_idi_external.csv"),
                    stringsAsFactors = FALSE)

panel_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        legend.position = "none")

rows <- vector("list", length(FOLDERS))
idi_rows <- vector("list", length(FOLDERS))

for (i in seq_along(FOLDERS)) {
  folder <- FOLDERS[i]
  label <- SHORT[i]
  out_dir <- file.path(SCEN, folder, "output")
  co <- readRDS(file.path(out_dir, "cohort.rds"))
  mo <- readRDS(file.path(out_dir, "model.rds"))
  ext <- co$ext
  y <- ext[[G]]

  p_new <- mo$deploy$predict_fn(ext[, mo$features, drop = FALSE],
                                selected_model = mo$cv_winner)[[POS]]

  set.seed(2026)
  idx <- caret::createDataPartition(co$dev[[G]], p = 0.70, list = FALSE)
  tr <- co$dev[idx, ]
  ref_fit <- glm(as.formula(sprintf("%s ~ RIDAGEYR + RIAGENDR", G)),
                 data = tr, family = binomial())
  p_ref <- as.numeric(predict(ref_fit, newdata = ext, type = "response"))

  row <- nri_df[nri_df$scenario == label, ]

  # (a) confusion matrix -------------------------------------------------------
  cm <- plots[[paste0("cm_", folder)]] +
    labs(title = "a) Confusion matrix") + panel_theme

  # (b) external ROC comparison -----------------------------------------------
  roc_new <- pROC::roc(y, p_new, levels = c(NEG, POS), direction = "<", quiet = TRUE)
  roc_ref <- pROC::roc(y, p_ref, levels = c(NEG, POS), direction = "<", quiet = TRUE)
  auc_new <- as.numeric(auc(roc_new)); auc_ref <- as.numeric(auc(roc_ref))
  roc_df <- rbind(
    data.frame(fpr = 1 - roc_new$specificities, tpr = roc_new$sensitivities,
               model = "Full lab panel"),
    data.frame(fpr = 1 - roc_ref$specificities, tpr = roc_ref$sensitivities,
               model = "Demographics only"))
  roc_df$model <- factor(roc_df$model,
                         levels = c("Full lab panel", "Demographics only"))
  roc_plot <- ggplot(roc_df, aes(fpr, tpr, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, colour = "grey55") +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("Full lab panel" = "#2E8B57",
                                   "Demographics only" = "#9AA5B1")) +
    coord_equal() + xlim(0, 1) + ylim(0, 1) +
    labs(title = "b) External ROC", x = "False positive rate",
         y = "True positive rate",
         subtitle = sprintf("AUC %.3f vs %.3f", auc_new, auc_ref)) +
    panel_theme +
    theme(legend.position = "inside", legend.justification = c(1, 0),
          legend.location = "panel", legend.key.size = unit(0.12, "in"),
          legend.text = element_text(size = 7)) +
    guides(colour = guide_legend(title = NULL, nrow = 2))

  # (c) AUC gain ---------------------------------------------------------------
  auc_plot <- plots[[paste0("auc_", folder)]] +
    labs(title = "c) AUC gain") + panel_theme

  # (d) categorical NRI --------------------------------------------------------
  nri_plot <- plots[[paste0("nri_", folder)]] +
    labs(title = "d) Categorical NRI") + panel_theme

  rows[[i]] <- (cm | roc_plot | auc_plot | nri_plot) +
    plot_annotation(
      title = label,
      subtitle = sprintf("n = %d, %d positives | independent 2009-2010 cycle",
                         row$external_n, row$external_positive),
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.01),
                    plot.subtitle = element_text(size = 9, colour = "grey40", hjust = 0.01)))

  # IDI curve (official icare plot). Its default caption contains a stray CJK
  # token, so we replace it with an English caption at the figure level only.
  idi <- plots[[paste0("idi_", folder)]] +
    labs(caption = NULL,
         subtitle = sprintf("IDI %+.4f (95%% CI %.4f to %.4f)",
                            row$idi, row$idi_low, row$idi_high)) + panel_theme
  idi_rows[[i]] <- idi +
    plot_annotation(title = label,
                    theme = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.01)))
}

# ---- assemble & save ---------------------------------------------------------
composite <- wrap_plots(rows, ncol = 1) +
  plot_annotation(
    title = "External validation: routine laboratory panel vs demographics alone",
    subtitle = paste0("Three NHANES test-triage scenarios | reference = age + sex, new = full routine lab panel | ",
                      "independent 2009-2010 cycle | decision thresholds fixed on internal held-out data"),
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.02),
                  plot.subtitle = element_text(size = 10, colour = "grey35", hjust = 0.02))) &
  panel_theme

ggsave(file.path(FIG, "Fig2_external_confusion_NRI_composite.png"), composite,
       width = 17, height = 13.5, dpi = 220, bg = "white")
ggsave(file.path(FIG, "Fig2_external_confusion_NRI_composite.pdf"), composite,
       width = 17, height = 13.5, bg = "white")

idi_composite <- wrap_plots(idi_rows, ncol = 3) +
  plot_annotation(
    title = "Integrated discrimination improvement (IDI) on the external cycle",
    subtitle = paste0("New = full routine laboratory panel; reference = age + sex | bands = bootstrap 95% CI | ",
                      "IS: gain in mean risk for cases; IP: change in mean risk for non-cases (IDI = IS - IP)"),
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.02),
                  plot.subtitle = element_text(size = 9, colour = "grey35", hjust = 0.02))) &
  panel_theme
ggsave(file.path(FIG, "Fig3_IDI_curves_composite.png"), idi_composite,
       width = 17, height = 5.6, dpi = 220, bg = "white")
ggsave(file.path(FIG, "Fig3_IDI_curves_composite.pdf"), idi_composite,
       width = 17, height = 5.6, bg = "white")

cat("Saved composite figures to", FIG, "\n")
