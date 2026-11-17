# ==============================================================================
# Albuminuria - Step 1 of 3
# Download the NHANES source tables needed for this scenario with nhanesA
# ------------------------------------------------------------------------------
# Scenario      : Urine albumin/creatinine ratio UACR >= 30 mg/g (KDIGO A2+)
# Clinical action: Who needs a urine albumin/creatinine test to detect early CKD
# Cycles        : 2007-2008 (_E) development / 2009-2010 (_F) independent external validation
#
# This scenario downloads ONLY the tables listed below; unrelated tables are never touched:
#   DEMO     Demographics: age, sex, pregnancy status
#   ALB_CR   Urine albumin / creatinine: used to compute UACR (label of this scenario)
#   BIOPRO   Biochemistry panel: 23 raw LBXS* values
#   CBC      Complete blood count: 20 items
#   BMX      Anthropometry: BMI, waist circumference
#   BPX      Blood-pressure readings: 4 SBP / DBP measurements
#
# Input  : none (online; CDC NHANES XPT files)
# Output : data/<table>_<cycle>.rds
# Next   : 02_cohort.R
#
# Note: this folder ships with the downloaded cache, so re-running will not access the network.
#       Delete any rds in data/ (or the whole data/ folder) to re-download with nhanesA.
# ==============================================================================

suppressMessages(library(nhanesA))

# ---- 0. Locate this folder ---------------------------------------------------
script_dir <- function() {
  args <- commandArgs(FALSE)
  hit  <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (nzchar(p) && file.exists(p)) return(dirname(normalizePath(p)))
  }
  stop("Source this script from RStudio, or run it with Rscript and its full path.")
}
SCRIPT_DIR <- script_dir()

DATA_DIR <- file.path(SCRIPT_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Table list for this scenario -----------------------------------------
# One file per cycle, named <table>_E / <table>_F
TABLE_NAMES <- c("DEMO", "ALB_CR", "BIOPRO", "CBC", "BMX", "BPX")
CYCLES      <- c("E", "F")

# ---- 2. Download each table (skip if already cached) --------------------------
for (cy in CYCLES) {
  for (name in TABLE_NAMES) {
    tb   <- paste0(name, "_", cy)
    path <- file.path(DATA_DIR, paste0(tb, ".rds"))

    if (file.exists(path)) {
      cat(sprintf("  cached      %-12s\n", tb))
      next
    }

    cat(sprintf("  downloading %-12s ... ", tb))
    tbl <- tryCatch(
      nhanesA::nhanes(tb, translated = FALSE),
      error = function(e) { message(conditionMessage(e)); NULL }
    )
    if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) {
      stop("Download failed: ", tb, " (check network or table name)")
    }
    saveRDS(tbl, path)
    cat(sprintf("OK  %d rows x %d cols\n", nrow(tbl), ncol(tbl)))
  }
}

cat("\nScenario data ready: ", length(TABLE_NAMES) * length(CYCLES),
    " tables -> data/\nNext step: 02_cohort.R\n")
