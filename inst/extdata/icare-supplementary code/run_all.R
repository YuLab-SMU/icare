#!/usr/bin/env Rscript
# ==============================================================================
# run_all.R - reproduce every analysis in the paper
# ------------------------------------------------------------------------------
#   Rscript run_all.R                       all three modules, full settings
#   Rscript run_all.R leakage               one module
#   Rscript run_all.R leakage head2head     any subset, in the order given
#   QUICK=1 Rscript run_all.R               smoke test (minutes, not hours)
#
# Modules: leakage | head2head | scenarios
#
# Each module runs in its own R process, so a crash in one does not contaminate
# the others and memory is released between them. Progress and timings are
# printed to the console and appended to outputs/run_all.log.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

args  <- commandArgs(trailingOnly = TRUE)
steps <- c(
  leakage   = "01_leakage/run_leakage.R",
  head2head = "02_head2head/framework_head2head.R",
  scenarios = "03_clinical_scenarios/run_scenarios.R")

mods <- if (length(args)) args else names(steps)
unknown <- setdiff(mods, names(steps))
if (length(unknown))
  stop("Unknown module(s): ", paste(unknown, collapse = ", "),
       "\n  available: ", paste(names(steps), collapse = ", "), call. = FALSE)

log_file <- file.path(ROOT, "outputs", "run_all.log")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

say <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
}

say("run_all: %s%s", paste(mods, collapse = ", "),
    if (env_flag("QUICK")) "   [QUICK = 1]" else "")

rscript <- file.path(R.home("bin"), "Rscript")
t_all   <- Sys.time()
failed  <- character(0)

for (m in mods) {
  f <- file.path(ROOT, steps[[m]])
  say("\n=== %s (%s) ===", m, basename(steps[[m]]))
  t0 <- Sys.time()
  status <- system2(rscript, shQuote(f))
  if (status != 0) {
    failed <- c(failed, m)
    say("[%s] FAILED (exit %d) after %.1f min", m, status,
        as.numeric(difftime(Sys.time(), t0, units = "mins")))
    next
  }
  say("[%s] finished in %.1f min", m,
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

say("\ntotal %.1f min - outputs under %s",
    as.numeric(difftime(Sys.time(), t_all, units = "mins")),
    file.path(ROOT, "outputs"))

if (length(failed)) {
  say("FAILED module(s): %s", paste(failed, collapse = ", "))
  quit(status = 1)
}
