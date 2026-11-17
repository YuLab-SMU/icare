# ==============================================================================
# Composite figure for the leakage-validation experiment (Fig. S12 of the paper)
# ------------------------------------------------------------------------------
# Inputs (written by 01_leakage/leakage_validation.R):
#   outputs/01_leakage/tables/expA_noise_icare.csv
#   outputs/01_leakage/tables/expB_real_cohort_2x2.csv
#   outputs/01_leakage/tables/expD_model_selection_bias.csv
#
# Outputs:
#   outputs/01_leakage/figures/Figure_leakage_validation.{pdf,png}
#   outputs/01_leakage/tables/figure_numbers.csv   (every number shown in the figure)
#
# Run: Rscript 01_leakage/make_composite_figure.R   (after leakage_validation.R)
#
# Layout 2 x 2:
#   A  pure-noise data, correct vs leaky screening (held-out AUC)
#   B  real cohort, 2 x 2 factorial (preprocessing x screening)
#   C  model selection on the test set vs by training-set CV
#   D  magnitude of the bias contributed by each leakage route (95% CI)
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))
check_packages(required = c("ggplot2", "patchwork", "scales"))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(scales)
})

OUT  <- out_dir("01_leakage")
TAB  <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

missing <- c("expA_noise_icare.csv", "expB_real_cohort_2x2.csv",
             "expD_model_selection_bias.csv")
missing <- missing[!file.exists(file.path(TAB, missing))]
if (length(missing))
  stop("Missing input table(s): ", paste(missing, collapse = ", "),
       "\n  run 01_leakage/leakage_validation.R first.", call. = FALSE)

exp_A <- read.csv(file.path(TAB, "expA_noise_icare.csv"))
exp_B <- read.csv(file.path(TAB, "expB_real_cohort_2x2.csv"), check.names = FALSE)
exp_D <- read.csv(file.path(TAB, "expD_model_selection_bias.csv"))

R_NOISE <- nrow(exp_A); R_REAL <- nrow(exp_B); R_SEL <- nrow(exp_D)
ARM <- c("train|train", "full|train", "train|all", "full|all")
stopifnot(all(ARM %in% names(exp_B)))

# --- shared theme ------------------------------------------------------------------------
blue <- "#2c7fb8"; red <- "#d7301f"; green <- "#1a9850"; grey <- "#4d4d4d"
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5, colour = grey),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8.5),
        plot.margin = margin(6, 8, 6, 6))

# ==============================================================================
# Panel A -- pure-noise data: correct vs leaky feature screening
# ==============================================================================
long_A <- rbind(
  data.frame(sim = exp_A$sim, Pipeline = "Correct\n(screening: training only)", AUC = exp_A$correct),
  data.frame(sim = exp_A$sim, Pipeline = "Leaky\n(screening: all samples)",    AUC = exp_A$leaky))
long_A$Pipeline <- factor(long_A$Pipeline, levels = unique(long_A$Pipeline))

mA <- mean(exp_A$correct); mL <- mean(exp_A$leaky); dA <- exp_A$leaky - exp_A$correct
ciA <- ci95(dA)

pA <- ggplot(long_A, aes(Pipeline, AUC)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40", linewidth = .4) +
  geom_line(aes(group = sim), colour = "grey75", linewidth = .25) +
  geom_boxplot(aes(fill = Pipeline), width = .42, alpha = .65, outlier.shape = NA, linewidth = .4) +
  geom_point(aes(colour = Pipeline), position = position_jitter(width = .04, height = 0),
             size = 1.1, alpha = .8) +
  scale_fill_manual(values = c(blue, red), guide = "none") +
  scale_colour_manual(values = c(blue, red), guide = "none") +
  annotate("text", x = 1.5, y = 0.985, hjust = .5, size = 3, colour = red,
           label = sprintf("inflation +%s (95%% CI %s to %s)", fmt(mean(dA)), fmt(ciA[1]), fmt(ciA[2]))) +
  coord_cartesian(ylim = c(0.25, 1.02)) +
  labs(title = "A  Pure-noise data (true AUC = 0.50)",
       subtitle = sprintf("n = 120, p = 500, %d random splits; mean held-out AUC", R_NOISE),
       x = NULL, y = "Held-out AUC") + th

# ==============================================================================
# Panel B -- real cohort: 2 x 2 factorial (preprocessing x screening)
# ==============================================================================
summ_B <- do.call(rbind, lapply(ARM, function(a) {
  ci <- ci95(exp_B[[a]]); b <- exp_B[[a]] - exp_B[["train|train"]]
  data.frame(arm = a, mean_AUC = mean(exp_B[[a]]), lo = ci[1], hi = ci[2],
             bias = mean(b), b_lo = ci95(b)[1], b_hi = ci95(b)[2]) }))
labs_B <- c("Correct\ntrain | train", "Preprocessing leak\nall | train",
            "Screening leak\ntrain | all", "Both leaks\nall | all")
summ_B$arm_lab <- factor(labs_B, levels = labs_B)
summ_B$correct <- summ_B$arm == ARM[1]

pB <- ggplot(summ_B, aes(arm_lab, mean_AUC)) +
  geom_pointrange(aes(ymin = lo, ymax = hi, colour = correct), size = .55, linewidth = .8) +
  geom_text(aes(label = sprintf("%s\n(bias %+.3f)", fmt(mean_AUC), bias)),
            vjust = -.9, size = 2.7, colour = grey, lineheight = .95) +
  scale_colour_manual(values = c(`TRUE` = blue, `FALSE` = red), guide = "none") +
  coord_cartesian(ylim = c(0.935, 0.965)) +
  labs(title = "B  Real cohort (n = 601, 23 routine indices)",
       subtitle = sprintf("mean \u00b1 95%% CI over %d random 70/30 splits; axis range fixed at 0.935-0.965", R_REAL),
       x = NULL, y = "Validation AUC") + th

# ==============================================================================
# Panel C -- model selection on the test set vs by training-set CV
# ==============================================================================
long_C <- rbind(
  data.frame(Selection = "By\ntest-set AUC", AUC = exp_D$pick_by_test),
  data.frame(Selection = "By\ntraining-set CV", AUC = exp_D$pick_by_cv))
long_C$Selection <- factor(long_C$Selection, levels = unique(long_C$Selection))
dC <- exp_D$pick_by_test - exp_D$pick_by_cv; ciC <- ci95(dC)

pC <- ggplot(long_C, aes(Selection, AUC)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40", linewidth = .4) +
  geom_boxplot(aes(fill = Selection), width = .42, alpha = .65, outlier.shape = NA, linewidth = .4) +
  geom_point(aes(colour = Selection), position = position_jitter(width = .05, height = 0),
             size = 1.1, alpha = .8) +
  scale_fill_manual(values = c(red, green), guide = "none") +
  scale_colour_manual(values = c(red, green), guide = "none") +
  annotate("text", x = 1.5, y = 0.815, hjust = .5, size = 3, colour = red,
           label = sprintf("selection bias +%s (95%% CI %s to %s)", fmt(mean(dC)), fmt(ciC[1]), fmt(ciC[2]))) +
  coord_cartesian(ylim = c(0.38, 0.83)) +
  labs(title = "C  Algorithm chosen on the test set (noise data)",
       subtitle = sprintf("true AUC = 0.50; 4 algorithms; %d splits; test AUC recomputed with ROC direction fixed", R_SEL),
       x = NULL, y = "Test-set AUC of the selected model") + th

# ==============================================================================
# Panel D -- magnitude of the bias contributed by each leakage route
# ==============================================================================
bB <- exp_B[["full|train"]] - exp_B[["train|train"]]; ciB <- ci95(bB)
bias_tab <- data.frame(
  route = c("Screening on all samples\n(noise, p = 500)",
            "Selection on the test set\n(noise)",
            "Preprocessing on all samples\n(real cohort)"),
  est = c(mean(dA), mean(dC), mean(bB)),
  lo  = c(ciA[1], ciC[1], ciB[1]),
  hi  = c(ciA[2], ciC[2], ciB[2]))
bias_tab$col <- c(red, red, green)
bias_tab$route <- factor(bias_tab$route, levels = rev(bias_tab$route))
bias_tab$lab_x <- bias_tab$hi + 0.010

pD <- ggplot(bias_tab, aes(est, route)) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = .4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = .18, linewidth = .7, colour = grey, orientation = "y") +
  geom_point(aes(colour = col), size = 3) +
  geom_text(aes(x = lab_x, label = sprintf("%+.3f", est)), hjust = 0,
            size = 2.9, colour = grey) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(-0.02, 0.37), breaks = seq(0, 0.30, .10)) +
  labs(title = "D  Bias contributed by each leakage route",
       subtitle = "mean bias vs the leakage-free pipeline, with 95% CI",
       x = "Bias in reported AUC", y = NULL) +
  th + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 8, lineheight = .95))

# ==============================================================================
# Assemble
# ==============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = NULL,
    caption = "Panels A-C: every comparison keeps the modelling pipeline identical and changes only which samples a step is allowed to see.") &
  theme(plot.caption = element_text(size = 7.5, colour = grey, hjust = 0))

ggsave(file.path(FIG, "Figure_leakage_validation.pdf"), fig, width = 10, height = 8.2)
ggsave(file.path(FIG, "Figure_leakage_validation.png"), fig, width = 10, height = 8.2, dpi = 400, bg = "white")

# ==============================================================================
# Every number used above, in one table (for the legend and the Results text)
# ==============================================================================
nums <- data.frame(
  quantity = c("A_correct_mean", "A_leaky_mean", "A_inflation", "A_inflation_CI_lo", "A_inflation_CI_hi",
               "A_runs_with_leaky_gt0", "A_runs", "A_test_samples_seen_correct", "A_test_samples_seen_leaky",
               paste0("B_", ARM, "_mean"), paste0("B_", ARM, "_bias"),
               "C_selection_on_test_mean", "C_selection_by_cv_mean", "C_selection_bias",
               "C_selection_bias_CI_lo", "C_selection_bias_CI_hi"),
  value = c(mean(exp_A$correct), mean(exp_A$leaky), mean(dA), ciA[1], ciA[2],
            sum(exp_A$seen_test_leaky > 0), nrow(exp_A), sum(exp_A$seen_test_correct), sum(exp_A$seen_test_leaky),
            summ_B$mean_AUC, summ_B$bias,
            mean(exp_D$pick_by_test), mean(exp_D$pick_by_cv), mean(dC), ciC[1], ciC[2]))
write.csv(nums, file.path(TAB, "figure_numbers.csv"), row.names = FALSE)
print(nums, digits = 3)
cat("\n[pA/pB/pC/pD written to", normalizePath(FIG), "]\n")
