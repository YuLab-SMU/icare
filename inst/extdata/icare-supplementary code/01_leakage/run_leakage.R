#!/usr/bin/env Rscript
# ==============================================================================
# 01_leakage/run_leakage.R - batch entry point for the leakage module
# ------------------------------------------------------------------------------
# The experiments themselves live in leakage_validation.R, which is written to be
# read and run block by block. It is fully sequential and every experiment sets
# its own seed, so sourcing it from a clean session reproduces the same numbers.
#
#   Rscript 01_leakage/run_leakage.R              full run (20 / 20 / 30 repeats)
#   QUICK=1 Rscript 01_leakage/run_leakage.R      smoke test (5 repeats, ~2 min)
#
# Then draws the 2 x 2 composite figure of Fig. S12 from the tables just written.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required to locate the repository root: install.packages('here')",
       call. = FALSE)
source(here::here("R", "setup.R"))

t0 <- Sys.time()
cat("[01_leakage] experiments ...\n")
source(here::here("01_leakage", "leakage_validation.R"), echo = FALSE, max.deparse.length = 200)

cat("[01_leakage] composite figure ...\n")
source(here::here("01_leakage", "make_composite_figure.R"), echo = FALSE)

cat(sprintf("[01_leakage] done in %.1f min - outputs in %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            out_dir("01_leakage", create = FALSE)))
