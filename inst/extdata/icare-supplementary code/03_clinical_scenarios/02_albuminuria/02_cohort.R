# ==============================================================================
# Albuminuria - Step 2 of 3
# Merge the NHANES source tables and build the development / external cohorts
# ------------------------------------------------------------------------------
# Scenario      : Urine albumin/creatinine ratio UACR >= 30 mg/g (KDIGO A2+)
# Label         : uacr >= 30
# Clinical action: Who needs a urine albumin/creatinine test to detect early CKD
#
# Design (test-triage): the label is a gold-standard test result that a clinician
# would otherwise have to order; features are only routine, low-cost measurements.
# The gold-standard variables never enter the feature set.
#
# Input  : data/*.rds (produced by 01_download.R)
# Output : output/cohort.rds, output/cohort_flow.csv
# Next   : 03_model.R
# ==============================================================================

suppressMessages(library(dplyr))

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

DATA_DIR   <- file.path(SCRIPT_DIR, "data")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Constants: cycles, class names, feature panel -------------------------
DEV_CYCLE <- "E"    # 2007-2008, development cohort
EXT_CYCLE <- "F"    # 2009-2010, independent external validation cohort
NEG       <- "Negative"
POS       <- "Positive"
GROUP     <- "group"

DEMO    <- c(
  "RIDAGEYR", "RIAGENDR"
  )
ANTHRO  <- c(
  "BMXBMI", "BMXWAIST"
  )
BP_MEAN <- c("mean_sbp", "mean_dbp")
BIO     <- c(
  "LBXSAL", "LBXSATSI", "LBXSASSI", "LBXSAPSI",
  "LBXSBU", "LBXSCA", "LBXSCH", "LBXSC3SI",
  "LBXSCR", "LBXSGL", "LBXSGTSI", "LBXSIR",
  "LBXSLDSI", "LBXSPH", "LBXSTB", "LBXSTP",
  "LBXSTR", "LBXSUA", "LBXSNASI", "LBXSKSI",
  "LBXSCLSI", "LBXSOSSI", "LBXSGB"
  )
CBC     <- c(
  "LBXWBCSI", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT",
  "LBXEOPCT", "LBXBAPCT", "LBDLYMNO", "LBDMONO",
  "LBDNENO", "LBDEONO", "LBDBANO", "LBXRBCSI",
  "LBXHGB", "LBXHCT", "LBXMCVSI", "LBXMCHSI",
  "LBXMC", "LBXRDW", "LBXPLTSI", "LBXMPSI"
  )

# Variables explicitly excluded for this scenario (label components or derived composites)
# uacr / URXUMA / URXUCR are used only to define the label and never
#   enter the candidate feature pool; serum creatinine LBXSCR and BUN
#   LBXSBU are routine blood tests and are retained as predictors.
#   LBXSOSSI is a derived composite removed from all three scenarios.
EXCLUDE <- c(
  "LBXSOSSI"
  )

# Candidate pool -> features actually used for modelling
POOL     <- c(DEMO, ANTHRO, BP_MEAN, BIO, CBC)
FEATURES <- setdiff(POOL, EXCLUDE)

# ---- 2. Helpers ---------------------------------------------------------------
read_table <- function(name, cycle) {
  path <- file.path(DATA_DIR, paste0(name, "_", cycle, ".rds"))
  if (!file.exists(path)) stop("Missing ", path, "; run 01_download.R first")
  readRDS(path)
}

# Standardized clinic BP: mean of valid readings; 0 or out-of-range = not measured
mean_blood_pressure <- function(bpx) {
  sys <- bpx[, grep("^BPXSY[1-4]$", names(bpx), value = TRUE), drop = FALSE]
  dia <- bpx[, grep("^BPXDI[1-4]$", names(bpx), value = TRUE), drop = FALSE]
  sys[sys <= 0 | sys > 300] <- NA
  dia[dia <= 0 | dia > 200] <- NA

  out <- data.frame(
    SEQN           = bpx$SEQN,
    mean_sbp       = rowMeans(sys, na.rm = TRUE),
    mean_dbp       = rowMeans(dia, na.rm = TRUE),
    n_bp_readings  = rowSums(!is.na(sys))
  )
  out$mean_sbp[!is.finite(out$mean_sbp)] <- NA
  out$mean_dbp[!is.finite(out$mean_dbp)] <- NA
  out
}

# ---- 3. Build the cohort for one cycle ---------------------------------------
build_cycle <- function(cycle) {
  # 3.1 Read tables, keeping only the columns used by this scenario
  demo <- read_table("DEMO", cycle)[, c("SEQN", DEMO, "RIDEXPRG")]
  bio  <- read_table("BIOPRO", cycle)[, c("SEQN", BIO)]
  cbc  <- read_table("CBC", cycle)[, c("SEQN", CBC)]
  bmx  <- read_table("BMX", cycle)[, c("SEQN", ANTHRO)]
  bpx  <- read_table("BPX", cycle)
  bp   <- mean_blood_pressure(bpx)
  alb <- read_table("ALB_CR", cycle)[, c("SEQN", "URXUMA", "URXUCR")]
  # UACR (mg/g): URXUMA is mg/L, URXUCR is mg/dL = 10 mg/L
  uacr <- data.frame(
    SEQN = alb$SEQN,
    uacr = 100 * alb$URXUMA / alb$URXUCR
  )
  uacr$uacr[!is.finite(uacr$uacr)] <- NA


  # 3.2 Left-join on SEQN (DEMO is the master table)
  d <- demo |>
    left_join(bio, by = "SEQN") |>
    left_join(cbc, by = "SEQN") |>
    left_join(bmx, by = "SEQN") |>
    left_join(bp,  by = "SEQN") |>
    left_join(uacr, by = "SEQN") |>
    as.data.frame()

  # 3.3 General eligibility: adults, not pregnant
  d <- d[!is.na(d$RIDAGEYR) & d$RIDAGEYR >= 18 & !(d$RIDEXPRG %in% 1), ]

  # 3.4 Scenario-specific eligibility
  n_eligible_before <- nrow(d)
  d <- d[!is.na(!is.na(d$uacr)) & (!is.na(d$uacr)), ]

  # 3.5 Define the binary label
  is_positive <- d$uacr >= 30
  d <- d[!is.na(is_positive), ]
  d[[GROUP]] <- factor(ifelse(d$uacr >= 30, POS, NEG), levels = c(NEG, POS))

  # 3.6 Complete cases on the feature columns only (label is not checked for missingness)
  n_before_complete <- nrow(d)
  d <- d[complete.cases(d[, FEATURES]), ]

  # 3.7 Cohort flow, for verifying sample sizes in the report
  attr(d, "flow") <- data.frame(
    cycle              = cycle,
    adults_nonpregnant = n_eligible_before,
    scenario_eligible  = n_before_complete,
    complete_cases     = nrow(d),
    positive           = sum(d[[GROUP]] == POS),
    negative           = sum(d[[GROUP]] == NEG)
  )
  d
}

# ---- 4. Build both cycles -----------------------------------------------------
cat("Building cohorts ...\n")
dev <- build_cycle(DEV_CYCLE)
ext <- build_cycle(EXT_CYCLE)

flow <- rbind(attr(dev, "flow"), attr(ext, "flow"))
write.csv(flow, file.path(OUTPUT_DIR, "cohort_flow.csv"), row.names = FALSE)

cat("\nCohort flow:\n")
print(flow, row.names = FALSE)
cat("\nNumber of features: ", length(FEATURES), "\n")
cat("Excluded variables: ", paste(EXCLUDE, collapse = ", "), "\n")

# ---- 5. Save ------------------------------------------------------------------
saveRDS(
  list(
    dev = dev, ext = ext,
    features = FEATURES, exclude = EXCLUDE,
    label_text = "uacr >= 30",
    scenario_title = "Albuminuria",
    clinical_action = "Who needs a urine albumin/creatinine test to detect early CKD",
    dev_cycle = DEV_CYCLE, ext_cycle = EXT_CYCLE
  ),
  file = file.path(OUTPUT_DIR, "cohort.rds")
)
cat("\nSaved output/cohort.rds\nNext step: 03_model.R\n")
