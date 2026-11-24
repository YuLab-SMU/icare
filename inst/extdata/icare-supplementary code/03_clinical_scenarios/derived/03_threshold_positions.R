# ==============================================================================
# Diagnostic: where does the Youden threshold actually sit?
# For each scenario: distribution of predicted risk in the external cycle,
# split by true outcome, with the operating threshold marked.
# Output: outputs/03_clinical_scenarios/figures/FigS_threshold_positions.png/.pdf
#         outputs/03_clinical_scenarios/tables/threshold_audit.csv
# ==============================================================================
suppressMessages({ library(icare); library(caret); library(pROC); library(ggplot2); library(patchwork) })
## --- Locate the repository ----------------------------------------------------
if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

SCEN_DIR <- file.path(ROOT, "03_clinical_scenarios")   # the three scenario folders
TAB  <- out_dir("03_clinical_scenarios", "tables")
FIG  <- out_dir("03_clinical_scenarios", "figures")

ordn <- function(x) { x <- round(x); suf <- if (x %% 100 %in% 11:13) "th" else c("th","st","nd","rd","th","th","th","th","th","th")[x %% 10 + 1]; paste0(x, suf) }

SCEN <- c("01_undiagnosed_dysglycemia", "02_albuminuria", "03_undiagnosed_hypertension")
TITLE <- c("Undiagnosed dysglycemia\n(prevalence 3.2%)",
           "Albuminuria\n(prevalence 10.2%)",
           "Undiagnosed hypertension\n(prevalence 9.1%)")

plots <- list(); rows <- list()
for (i in seq_along(SCEN)) {
  s <- SCEN[i]
  coh <- readRDS(file.path(SCEN_DIR, s, "output", "cohort.rds")); mdl <- readRDS(file.path(SCEN_DIR, s, "output", "model.rds"))
  p <- mdl$deploy$predict_fn(newdata = coh$ext[, mdl$features, drop = FALSE],
                             selected_model = mdl$cv_winner)[["Positive"]]
  y <- coh$ext$group; thr <- mdl$youden_threshold
  eps <- 1e-4
  d <- data.frame(risk = pmax(p, eps), outcome = y)
  pct <- 100 * mean(p <= thr); flag <- 100 * mean(p >= thr)
  ys <- function(t) { pr <- ifelse(p >= t, 1, 0)
    mean(pr[y == "Positive"] == 1) + mean(pr[y == "Negative"] == 0) - 1 }
  grid <- exp(seq(log(min(p[p > 0])), log(max(p)), length.out = 200))
  you <- sapply(grid, ys)
  best <- grid[which.max(you)]

  rows[[s]] <- data.frame(scenario = s, threshold = thr,
                          pct_of_patients_below = pct, flag_rate_pct = flag,
                          youden_index = ys(thr), best_on_external = best,
                          youden_optimum_external = max(you),
                          youden_at_flat_lo = grid[max(which(you > max(you) - 0.01)[1])],
                          youden_at_flat_hi = grid[max(which(you > max(you) - 0.01))])

  plots[[i]] <- ggplot(d, aes(risk, fill = outcome)) +
    geom_histogram(bins = 60, position = "stack", alpha = .85, colour = NA) +
    geom_vline(xintercept = thr, linetype = 2, linewidth = .6, colour = "#d7301f") +
    annotate("text", x = min(d$risk), y = Inf, hjust = 0, vjust = 1.4, size = 3,
             colour = "#d7301f",
             label = sprintf("Youden %.3f\n%s pct of risk\nflags %.0f%% of the cohort",
                             thr, ordn(pct), flag)) +
    scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
    scale_fill_manual(values = c(Negative = "#9ecae1", Positive = "#d7301f"),
                      name = "True status") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.30))) +
    labs(title = TITLE[i], x = "Predicted risk (log scale)", y = "Patients") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = i,
          plot.title = element_text(size = 10, face = "bold"))
}

fig <- (plots[[1]] | plots[[2]] | plots[[3]]) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = paste("Dashed line = the Youden threshold fixed on the internal held-out split and applied unchanged to the external cycle.",
                                  "It is not a probability of 0.5: with 3-10% prevalence the whole risk distribution sits below 0.2,",
                                  "and the threshold falls at the 60th-80th percentile of that distribution.", sep = "\n")) &
  theme(plot.caption = element_text(size = 7.5, hjust = 0))

ggsave(file.path(FIG, "FigS_threshold_positions.png"), fig, width = 11, height = 4.6, dpi = 400, bg = "white")
ggsave(file.path(FIG, "FigS_threshold_positions.pdf"), fig, width = 11, height = 4.6)

print(do.call(rbind, rows), digits = 4, row.names = FALSE)
write.csv(do.call(rbind, rows), file.path(TAB, "threshold_audit.csv"), row.names = FALSE)
