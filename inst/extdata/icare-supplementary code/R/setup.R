# ==============================================================================
# R/setup.R - shared configuration and helpers
# ------------------------------------------------------------------------------
# Sourced by every entry point in this repository. It provides
#
#   icare_code_root()   the repository root (the folder that holds R/setup.R)
#   out_dir()           <root>/outputs/<module>, created on demand
#   check_packages()    one readable error instead of a chain of library() failures
#   auc_of()            AUC with the ROC direction fixed
#   brier_of()          Brier score
#   jaccard()           overlap of two feature sets
#   ci95()              normal-approximation 95% CI, for display only
#   fmt()               fixed-width number formatting
#
# No analysis decision lives here. Seeds, folds, thresholds and equivalence
# margins stay in the script that reports them, so every number in the paper can
# be traced back to one file and one line.
# ==============================================================================

## --- Repository root ----------------------------------------------------------
# Walks up from `start` until it finds R/setup.R, so the same code works when the
# working directory is the repository root, a module folder, or a subfolder of
# one (RStudio, Rscript, or an interactive session).
icare_code_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "R", "setup.R"))) return(d)
    up <- dirname(d)
    if (identical(up, d))
      stop("R/setup.R not found above '", start, "'.\n",
           "  Run this script from inside the repository, i.e. from the folder ",
           "that contains run_all.R.", call. = FALSE)
    d <- up
  }
}

## --- Output folders -----------------------------------------------------------
# The repository root, resolved once when this file is sourced. `here` uses the
# .here marker at the root, so this is correct whether the working directory is
# the root, a module folder, or a subfolder of one.
ROOT <- if (requireNamespace("here", quietly = TRUE)) here::here() else icare_code_root()

# All generated material goes to <root>/outputs/<module>[/...]; that whole tree is
# git-ignored, so a clone never carries stale results around.
out_dir <- function(module, ..., create = TRUE) {
  p <- do.call(file.path, c(list(ROOT, "outputs", module), list(...)))
  if (create) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

## --- Packages -----------------------------------------------------------------
# required : missing -> stop, with the install.packages() line to run
# optional : missing -> message, and the calling script skips the steps that
#            need it (documented in each module's README)
check_packages <- function(required, optional = character()) {
  has  <- function(p) requireNamespace(p, quietly = TRUE)
  miss <- required[!vapply(required, has, logical(1))]
  if (length(miss))
    stop("Missing package(s): ", paste(miss, collapse = ", "),
         "\n  install with: install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE)
  miss_opt <- optional[!vapply(optional, has, logical(1))]
  if (length(miss_opt))
    message("Note: optional package(s) not installed: ",
            paste(miss_opt, collapse = ", "),
            " - the steps that depend on them are skipped.")
  invisible(TRUE)
}

## --- Metrics ------------------------------------------------------------------
# AUC with the ROC direction fixed ("<").
# pROC's direction = "auto" flips every AUC below 0.5 upward, which turns a true
# 0.5 on null data into ~0.55. That single default would hide the very bias this
# repository is built to measure, so direction is always passed explicitly.
auc_of <- function(y, p) {
  as.numeric(pROC::auc(pROC::roc(response = y, predictor = p,
                                 levels = levels(y), direction = "<", quiet = TRUE)))
}

brier_of <- function(y, p) mean((p - as.numeric(y == "Case"))^2)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# Normal-approximation 95% CI of a mean. Descriptive only: every test that
# supports a claim in the paper is a proper paired test (wilcox.test / DeLong).
ci95 <- function(x) mean(x) + c(-1, 1) * 1.96 * sd(x) / sqrt(length(x))

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

## --- Small utilities ----------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# TRUE when the environment variable `var` is set to 1/true/yes
env_flag <- function(var, default = FALSE) {
  v <- tolower(trimws(Sys.getenv(var, "")))
  if (!nzchar(v)) return(default)
  v %in% c("1", "true", "yes", "y")
}

# "source/file.R" -> "file" (used for progress lines)
script_tag <- function(path) sub("\\.R$", "", basename(path))
