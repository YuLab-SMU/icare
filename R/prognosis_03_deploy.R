# ==============================================================================
# prognosis_03_deploy.R
# PrognosiX Deployment Module
#
# Provides:
#   1. Stat_to_PrognosiX()      -- robust Stat -> PrognosiX conversion
#   2. run_prognosis_pipeline() -- end-to-end survival analysis pipeline
#   3. Prog_deploy_dispatcher() -- risk prediction for new samples
#   4. New_Prog_Manager()       -- lightweight deployment manager object
#   5. launch_prog_deploy_app() -- interactive Shiny prediction terminal
# ==============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b


# ==============================================================================
# 1.  Stat -> PrognosiX Conversion
# ==============================================================================

#' Convert a Stat Object to a PrognosiX Object
#'
#' Fixes all common upstream data issues so the resulting object is immediately
#' ready for mlr3-based survival analysis:
#' \itemize{
#'   \item Special characters in numeric columns (\code{">60"}, \code{"<40"})
#'         are stripped and coerced to \code{numeric}.
#'   \item \strong{All \code{character} feature columns are converted to
#'         \code{factor}.}  mlr3 \code{TaskSurv} objects reject \code{character}
#'         features; this conversion is the root fix for the
#'         \emph{"unsupported feature types: character"} error.
#'   \item Rows with non-positive / non-finite time are removed.
#'   \item Rows where status is not \{0, 1\} are removed.
#'   \item Missing feature values are omitted or imputed.
#'   \item Row-name mismatches between \code{clean.data} and \code{info.data}
#'         are resolved by intersection.
#' }
#'
#' @param stat_obj    A \code{Stat} S4 object from \code{CreateStatObject()}.
#' @param time_col    Name of the survival time column (character string).
#' @param status_col  Name of the event status column (character string;
#'                    1 = event occurred, 0 = censored).
#' @param na_action   How to handle missing feature values:
#'   \describe{
#'     \item{\code{"omit"}}{(default) Delete any row that contains an NA in
#'           a feature column.}
#'     \item{\code{"impute_median"}}{Replace numeric NAs with the column
#'           median; replace categorical NAs with the column mode.}
#'   }
#' @param min_events  Minimum event count; a \code{warning()} is issued when
#'                    fewer events are present (default 20).
#' @param verbose     Print a step-by-step conversion log (default \code{TRUE}).
#'
#' @return A \code{PrognosiX} S4 object ready for \code{run_prognosis_pipeline()}.
#'
#' @export
#' @examples
#' \dontrun{
#' library(survival)
#' data("veteran")
#' veteran$celltype <- as.character(veteran$celltype)   # simulate dirty data
#'
#' stat_obj <- CreateStatObject(raw.data   = veteran,
#'                               clean.data = veteran,
#'                               group_col  = "celltype")
#'
#' prog_obj <- Stat_to_PrognosiX(stat_obj,
#'                                time_col   = "time",
#'                                status_col = "status")
#' }
Stat_to_PrognosiX <- function(stat_obj,
                              time_col,
                              status_col,
                              na_action  = c("omit", "impute_median"),
                              min_events = 20,
                              verbose    = TRUE) {
  
  na_action <- match.arg(na_action)
  .log <- function(fmt, ...) if (verbose) message(sprintf(fmt, ...))
  
  # 1. Validate and extract ---------------------------------------------------
  if (!inherits(stat_obj, "Stat"))
    stop("[Stat_to_PrognosiX] 'stat_obj' must be a 'Stat' S4 object.")
  
  .log("[1/6] Extracting data from Stat object...")
  core_data <- if (nrow(stat_obj@clean.data) > 0) {
    stat_obj@clean.data
  } else if (nrow(stat_obj@raw.data) > 0) {
    warning("[Stat_to_PrognosiX] clean.data is empty; using raw.data.")
    stat_obj@raw.data
  } else {
    stop("[Stat_to_PrognosiX] Both clean.data and raw.data are empty.")
  }
  
  info_data <- stat_obj@info.data
  if (nrow(info_data) > 0) {
    common <- intersect(rownames(core_data), rownames(info_data))
    if (length(common) == 0) {
      warning("[Stat_to_PrognosiX] Row names do not match between clean.data and info.data. info.data ignored.")
      info_data <- data.frame(row.names = rownames(core_data))
    } else {
      core_data <- core_data[common, , drop = FALSE]
      info_data <- info_data[common, , drop = FALSE]
    }
  }
  for (col in c(time_col, status_col))
    if (col %in% colnames(info_data) && !(col %in% colnames(core_data)))
      core_data[[col]] <- info_data[[col]]
  
  # 2. Verify required columns ------------------------------------------------
  .log("[2/6] Verifying time / status columns...")
  miss <- setdiff(c(time_col, status_col), colnames(core_data))
  if (length(miss) > 0)
    stop(sprintf("[Stat_to_PrognosiX] Column(s) not found: %s\n  Available: %s",
                 paste(miss, collapse = ", "),
                 paste(colnames(core_data), collapse = ", ")))
  
  # 3. Type coercion ---------------------------------------------------------
  .log("[3/6] Coercing column types...")
  # Strip ">", "<", spaces etc. then force numeric
  core_data[[time_col]]   <- suppressWarnings(
    as.numeric(gsub("[^0-9.-]", "", as.character(core_data[[time_col]]))))
  core_data[[status_col]] <- suppressWarnings(
    as.numeric(as.character(core_data[[status_col]])))
  
  bad_t <- !is.finite(core_data[[time_col]]) | core_data[[time_col]] <= 0
  if (any(bad_t)) {
    warning(sprintf("[Stat_to_PrognosiX] Removing %d row(s) with time <= 0 or NA.", sum(bad_t)))
    core_data <- core_data[!bad_t, , drop = FALSE]
  }
  bad_s <- is.na(core_data[[status_col]]) | !(core_data[[status_col]] %in% c(0, 1))
  if (any(bad_s)) {
    warning(sprintf("[Stat_to_PrognosiX] Removing %d row(s) where status not in {0,1}.", sum(bad_s)))
    core_data <- core_data[!bad_s, , drop = FALSE]
  }
  
  # ROOT FIX: character -> factor
  # mlr3 TaskSurv rejects character feature columns.
  # Converting to factor also lets surv_get_learner() apply its automatic
  # encoding pipeline for learners that cannot handle factors.
  feat_cols <- setdiff(colnames(core_data), c(time_col, status_col))
  for (col in feat_cols) {
    if (is.character(core_data[[col]])) {
      core_data[[col]] <- factor(core_data[[col]])
      .log("    [fix] %-14s  character -> factor  (%d levels)",
           col, nlevels(core_data[[col]]))
    }
  }
  
  # 4. Missing values ---------------------------------------------------------
  .log("[4/6] Missing values (strategy: %s)...", na_action)
  if (na_action == "omit") {
    n_before  <- nrow(core_data)
    core_data <- core_data[complete.cases(core_data[, feat_cols, drop = FALSE]), , drop = FALSE]
    removed   <- n_before - nrow(core_data)
    if (removed > 0)
      warning(sprintf("[Stat_to_PrognosiX] Removed %d/%d row(s) with NA in features.",
                      removed, n_before))
  } else {
    for (col in feat_cols) {
      n_na <- sum(is.na(core_data[[col]]))
      if (n_na > 0) {
        fill <- if (is.numeric(core_data[[col]])) {
          median(core_data[[col]], na.rm = TRUE)
        } else {
          names(sort(table(core_data[[col]]), decreasing = TRUE))[1]
        }
        core_data[[col]][is.na(core_data[[col]])] <- fill
        .log("    [impute] %-14s  %d NA -> %s", col, n_na, fill)
      }
    }
  }
  
  # 5. Event count ------------------------------------------------------------
  .log("[5/6] Checking event count...")
  n_ev <- sum(core_data[[status_col]] == 1)
  n_to <- nrow(core_data)
  .log("    N = %d | Events = %d | Censoring = %.1f%%",
       n_to, n_ev, (1 - n_ev / n_to) * 100)
  if (n_ev < min_events)
    warning(sprintf("[Stat_to_PrognosiX] Only %d events (< min_events=%d). Results may be unstable.",
                    n_ev, min_events))
  
  # 6. Build PrognosiX --------------------------------------------------------
  .log("[6/6] Building PrognosiX object...")
  new_info  <- core_data[, c(time_col, status_col), drop = FALSE]
  new_clean <- core_data[, feat_cols, drop = FALSE]
  
  prog_obj  <- CreatePrognosiXObject(
    clean.data = new_clean,
    info.data  = new_info,
    time_col   = time_col,
    status_col = status_col)
  
  .log("[OK] Done. Features: %d | Samples: %d", ncol(new_clean), nrow(new_clean))
  prog_obj
}


# ==============================================================================
# 2.  Complete Prognosis Analysis Pipeline
# ==============================================================================

#' Run the Complete Prognosis Analysis Pipeline
#'
#' Chains all prognosis analysis steps into a single call:
#' \enumerate{
#'   \item \strong{Feature filtering}  -- univariate Cox (\code{surv_filter_features_clinical})
#'   \item \strong{Algorithm benchmark} -- multi-model CV (\code{surv_run_algorithm_benchmark})
#'   \item \strong{Tuning + training}   -- best algorithm tuned and fitted (\code{surv_train_and_tune})
#'   \item \strong{KM risk curves}      -- stratified survival plots (\code{surv_plot_risk_km})
#'   \item \strong{Time-dependent AUC}  -- dynamic accuracy (\code{surv_plot_time_dependent_auc})
#'   \item \strong{Nomogram}            -- clinical scoring chart (requires \pkg{rms})
#'   \item \strong{SHAP}                -- feature explanation (requires \pkg{survex})
#'   \item \strong{Save}               -- all artefacts written to \code{output_dir}
#' }
#'
#' @param object         \code{PrognosiX} or \code{Stat} S4 object.
#' @param time_col       Survival time column (required for \code{Stat} input).
#' @param status_col     Event status column (required for \code{Stat} input).
#' @param learner_ids    mlr3 learner IDs to benchmark.  Run
#'   \code{surv_list_available_learners()} to see every option.
#'   Default: \code{c("surv.coxph", "surv.cv_glmnet", "surv.ranger")}.
#' @param p_threshold    Univariate Cox p-value cut-off (default \code{0.05}).
#' @param tuning_budget  Number of hyperparameter evaluations per model
#'   (default \code{30}).  Use \code{50}-\code{100} for publication results.
#' @param cutoff_method  Risk stratification cut-point:
#'   \code{"median"}, \code{"tertile"}, \code{"quartile"}, or \code{"p_optimize"}.
#' @param time_points    Nomogram prediction horizons (default \code{c(1,3,5)}).
#' @param output_dir     Root folder for all outputs.
#' @param seed           Random seed (default \code{2025}).
#' @param run_shap       Run SHAP explanation (default \code{FALSE}).
#' @param run_nomogram   Draw nomogram (default \code{TRUE}).
#' @param val_data       Optional external validation data frame.
#' @param subgroup_vars  Optional subgroup variables for forest-plot analysis.
#'
#' @return Invisible named list:
#'   \code{prog_obj}, \code{best_learner}, \code{filter_result},
#'   \code{benchmark_table}, \code{km_plot}, \code{auc_data}, \code{output_dir}.
#' @export
#' @examples
#' \dontrun{
#' # --- Minimal call (PrognosiX object) ---
#' res <- run_prognosis_pipeline(prog_obj)
#'
#' # --- From a Stat object with custom learners ---
#' res <- run_prognosis_pipeline(
#'   object        = stat_obj,
#'   time_col      = "OS_time",
#'   status_col    = "OS_status",
#'   learner_ids   = c("surv.coxph", "surv.ranger", "surv.cv_glmnet"),
#'   tuning_budget = 50,
#'   cutoff_method = "tertile"
#' )
#' }
run_prognosis_pipeline <- function(
    object,
    time_col      = NULL,
    status_col    = NULL,
    learner_ids   = c("surv.coxph", "surv.cv_glmnet", "surv.ranger"),
    p_threshold   = 0.05,
    tuning_budget = 30,
    cutoff_method = c("median", "tertile", "quartile", "p_optimize"),
    time_points   = c(1, 3, 5),
    output_dir    = NULL,
    seed          = 2025,
    run_shap      = FALSE,
    run_nomogram  = TRUE,
    val_data      = NULL,
    subgroup_vars = NULL) {
  
  cutoff_method <- match.arg(cutoff_method)
  set.seed(seed)
  
  if (is.null(output_dir)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path("./icare_output", "m4",
                            paste0("Prognosis_Pipeline_", ts))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Helper: create numbered step sub-directory
  .sdir <- function(n, nm) {
    p <- file.path(output_dir, sprintf("Step%02d_%s", n, nm))
    dir.create(p, recursive = TRUE, showWarnings = FALSE); p
  }
  # Helper: save ggplot safely
  .sp <- function(plt, path, w = 8, h = 6)
    tryCatch(ggplot2::ggsave(path, plt, width = w, height = h),
             error = function(e) NULL)
  
  message("\n", strrep("=", 65))
  message(" PrognosiX Pipeline  |  seed=", seed,
          "  |  output: ", output_dir)
  message(strrep("=", 65))
  
  # Step 0: object conversion -------------------------------------------------
  if (inherits(object, "Stat")) {
    if (is.null(time_col) || is.null(status_col))
      stop("[run_prognosis_pipeline] time_col and status_col required for Stat input.")
    message("\n[Step 0] Stat -> PrognosiX...")
    prog_obj <- Stat_to_PrognosiX(object, time_col, status_col, "omit")
  } else if (inherits(object, "PrognosiX")) {
    prog_obj <- object
  } else {
    stop("[run_prognosis_pipeline] object must be Stat or PrognosiX.")
  }
  saveRDS(prog_obj, file.path(output_dir, "prog_obj_initial.rds"))
  
  # Step 1: univariate Cox feature filtering ----------------------------------
  d1 <- .sdir(1, "Feature_Selection")
  message("\n[Step 1] Univariate Cox feature filtering (p < ", p_threshold, ")...")
  filter_result <- surv_filter_features_clinical(prog_obj, p_threshold = p_threshold)
  write.csv(filter_result$table, file.path(d1, "Univariate_Cox_Results.csv"),
            row.names = FALSE)
  .sp(filter_result$plot, file.path(d1, "Univariate_Feature_Selection.pdf"))
  
  task_filtered  <- filter_result$task
  selected_feats <- task_filtered$feature_names
  message(sprintf("  [OK] %d feature(s): %s",
                  length(selected_feats), paste(selected_feats, collapse = ", ")))
  
  if (length(selected_feats) == 0) {
    message("  [!] 0 features at p<", p_threshold, " -- relaxing to p<0.2 ...")
    filter_result  <- surv_filter_features_clinical(prog_obj, p_threshold = 0.2)
    task_filtered  <- filter_result$task
    selected_feats <- task_filtered$feature_names
    if (length(selected_feats) == 0)
      stop("Still 0 features at p<0.2. Check your data.")
  }
  prog_obj@univariate.analysis <- filter_result
  prog_obj@survival.var        <- list(selected = selected_feats)
  
  # Step 2: multi-algorithm benchmark ----------------------------------------
  d2 <- .sdir(2, "Algorithm_Benchmark")
  message("\n[Step 2] Algorithm benchmark (", paste(learner_ids, collapse = ", "), ")...")
  
  # CRITICAL: wrap via surv_get_learner() so factor/character columns are
  # handled by the automatic encoding pipeline and do NOT cause the
  # "unsupported feature types: character" error in mlr3.
  lrn_list <- Filter(Negate(is.null), lapply(learner_ids, function(lid) {
    tryCatch(surv_get_learner(lid, task_filtered),
             error = function(e) { warning("Skipping ", lid, ": ", e$message); NULL })
  }))
  
  bmr_result <- surv_run_algorithm_benchmark(task_filtered, learners_list = lrn_list)
  perf_tab   <- bmr_result$table
  best_id    <- as.character(
    perf_tab[order(-perf_tab$surv.cindex), ][1, "learner_id"])
  
  write.csv(perf_tab, file.path(d2, "Benchmark_Leaderboard.csv"), row.names = FALSE)
  .sp(bmr_result$plot, file.path(d2, "Algorithm_Benchmark.pdf"), w = 9, h = 5)
  message(sprintf("  [OK] Winner: %s  (C-index = %.4f)",
                  best_id, perf_tab[perf_tab$learner_id == best_id, "surv.cindex"]))
  
  # Step 3: hyperparameter tuning ---------------------------------------------
  d3 <- .sdir(3, "Final_Model")
  message("\n[Step 3] Tuning ", best_id, " (budget=", tuning_budget, " evals)...")
  tune_res     <- surv_train_and_tune(task_filtered, best_id,
                                      tuning_budget = tuning_budget, seed = seed)
  best_learner <- tune_res$learner
  saveRDS(best_learner, file.path(d3, "best_learner.rds"))
  write.csv(
    data.frame(Algorithm  = best_id,
               CV_CIndex  = round(tune_res$cv_performance, 4),
               N_Features = length(selected_feats),
               Params     = paste(names(tune_res$best_params),
                                  unlist(tune_res$best_params),
                                  sep = "=", collapse = "; ")),
    file.path(d3, "Best_Model_Summary.csv"), row.names = FALSE)
  message(sprintf("  [OK] CV C-index = %.4f", tune_res$cv_performance))
  
  # Step 4: KM risk stratification --------------------------------------------
  d4 <- .sdir(4, "Risk_KM")
  message("\n[Step 4] KM risk stratification (", cutoff_method, ")...")
  km_plot <- tryCatch({
    p <- surv_plot_risk_km(best_learner, task_filtered,
                           cutoff_method = cutoff_method, risk_table = TRUE)
    pdf(file.path(d4, paste0("KM_", cutoff_method, ".pdf")), width = 8, height = 7)
    print(p); dev.off(); p
  }, error = function(e) {
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    message("  [!] KM skipped: ", e$message); NULL
  })
  
  # Step 5: time-dependent AUC ------------------------------------------------
  d5 <- .sdir(5, "Time_AUC")
  message("\n[Step 5] Time-dependent AUC...")
  auc_data <- tryCatch({
    adf <- surv_plot_time_dependent_auc(best_learner, task_filtered)
    write.csv(adf, file.path(d5, "Time_AUC_Data.csv"), row.names = FALSE)
    .sp(last_plot(), file.path(d5, "Time_Dependent_AUC.pdf"))
    adf
  }, error = function(e) {
    message("  [!] AUC skipped (need risksetROC): ", e$message); NULL
  })
  
  # Step 6: nomogram (optional) -----------------------------------------------
  if (run_nomogram && requireNamespace("rms", quietly = TRUE)) {
    d6 <- .sdir(6, "Nomogram")
    message("\n[Step 6] Nomogram...")
    tryCatch({
      pdf(file.path(d6, "Clinical_Nomogram.pdf"), width = 13, height = 7)
      surv_generate_nomogram(task_filtered,
                             selected_features = head(selected_feats, 6),
                             time_points = time_points)
      dev.off(); message("  [OK] Nomogram saved.")
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
      message("  [!] Nomogram failed: ", e$message)
    })
  }
  
  # Step 7: SHAP (optional) ---------------------------------------------------
  if (run_shap) {
    d7 <- .sdir(7, "SHAP")
    message("\n[Step 7] SHAP explanation...")
    tryCatch({
      sr <- surv_explain_shap(best_learner, task_filtered)
      .sp(sr$plot, file.path(d7, "SHAP_Importance.pdf"))
      message("  [OK] SHAP saved.")
    }, error = function(e)
      message("  [!] SHAP skipped (need survex): ", e$message))
  }
  
  # Step 8: write back and save -----------------------------------------------
  message("\n[Step 8] Saving results to PrognosiX object...")
  prog_obj@best.model <- list(
    learner_id  = best_id,
    learner     = best_learner,
    best_params = tune_res$best_params,
    cv_cindex   = tune_res$cv_performance,
    features    = selected_feats)
  prog_obj@subgroup.risk <- list(benchmark_table = perf_tab)
  saveRDS(prog_obj, file.path(output_dir, "prog_obj_final.rds"))
  
  message("\n", strrep("=", 65))
  message(" [OK] Pipeline complete -> ", output_dir)
  message(strrep("=", 65))
  
  invisible(list(
    prog_obj        = prog_obj,
    best_learner    = best_learner,
    filter_result   = filter_result,
    benchmark_table = perf_tab,
    km_plot         = km_plot,
    auc_data        = auc_data,
    output_dir      = output_dir))
}


# ==============================================================================
# PrognosiX Deployment Module
# ==============================================================================
# Functions for deploying trained PrognosiX models via Shiny app,
# with support for custom thresholds (binary, tertile, multi-group),
# full UI text and color theming, and a production-ready manager.
# ---- 1. Theme customization ----
#' Set custom theme for the Prognosis Terminal
#' @param primary_color Main accent color (borders, buttons, headers)
#' @param background_color Main background color
#' @param sidebar_color Sidebar and header background
#' @param box_background Box background color
#' @param text_color General text color
#' @param label_color Label text color
#' @param run_button_gradient_start Left side of button gradient
#' @param run_button_gradient_end Right side of button gradient
#' @param risk_high_color Text color for "High Risk"
#' @param risk_medium_color Text color for "Medium Risk"
#' @param risk_low_color Text color for "Low Risk"
#' @param table_header_color Table header text color
#' @param table_row_hover_color Hover background for table rows
#' @param font_family CSS font family (e.g., "Arial, sans-serif")
#' @param font_size_base Base font size in px
#' @return Invisibly stores the theme in options("prog_app_theme_css")
#' @export
#' 
#' @examples
#' \dontrun{
#' set_prog_app_theme(primary_color = "#2c7fb8", background_color = "#f8f9fa")
#' }
set_prog_app_theme <- function(
    primary_color = "#2c7fb8",
    background_color = "#f8f9fa",
    sidebar_color = "#e9ecef",
    box_background = "#ffffff",
    text_color = "#212529",
    label_color = "#2c7fb8",
    run_button_gradient_start = "#2c7fb8",
    run_button_gradient_end = "#1d4e6e",
    risk_high_color = "#d9534f",
    risk_medium_color = "#f0ad4e",
    risk_low_color = "#5cb85c",
    table_header_color = "#2c7fb8",
    table_row_hover_color = "#212529",
    font_family = NULL,
    font_size_base = 14
) {
  css <- sprintf("
    .content-wrapper, .right-side { background: %s !important; }
    .skin-black .main-header .logo, .skin-black .main-header .navbar,
    .skin-black .main-sidebar { background: %s !important; }
    .box { background: %s !important; color: %s !important;
           border-top: 3px solid %s !important; border-radius: 8px; }
    .box-header .box-title { color: %s !important; font-weight: bold; }
    label { color: %s !important; font-weight: bold; font-size: %dpx; }
    .form-control { background: %s !important; color: %s !important;
                    border: 1px solid #ced4da !important; }
    .btn-run { background: linear-gradient(90deg, %s, %s) !important;
               color: #fff !important; font-size: 18px !important; font-weight: bold !important;
               padding: 12px 40px !important; border-radius: 8px !important; border: none !important;
               margin: 16px 0; box-shadow: 0 2px 6px rgba(0,0,0,0.1); }
    .risk-high { font-size: 48px; font-weight: 900; color: %s; }
    .risk-med { font-size: 48px; font-weight: 900; color: %s; }
    .risk-low { font-size: 48px; font-weight: 900; color: %s; }
    .score-lbl { font-size: 14px; color: %s; font-family: monospace; }
    .badge-box { background: %s; border: 1px solid %s; border-radius: 6px;
                 padding: 8px 16px; margin: 4px; display: inline-block; color: %s; }
    table.dataTable { background: %s !important; color: %s !important; }
    table.dataTable thead th { border-bottom: 1px solid %s !important; color: %s !important; }
    table.dataTable tbody tr { background: %s !important; }
    table.dataTable tbody tr:hover { background: %s !important; }
    h4 { color: %s; font-weight: bold; border-left: 3px solid %s; padding-left: 10px; }
    p { color: %s; }
  ",
                 background_color, sidebar_color,
                 box_background, text_color, primary_color,
                 primary_color,
                 label_color, font_size_base,
                 box_background, text_color,
                 run_button_gradient_start, run_button_gradient_end,
                 risk_high_color, risk_medium_color, risk_low_color,
                 text_color, box_background, sidebar_color, primary_color,
                 box_background, text_color, sidebar_color, table_header_color,
                 box_background, table_row_hover_color,
                 primary_color, primary_color, text_color
  )
  if (!is.null(font_family)) {
    css <- paste0("body, .box, label, .form-control, .btn-run, .score-lbl, .badge-box, table, h4, p { font-family: ", font_family, " !important; }\n", css)
  }
  options("prog_app_theme_css" = css)
  invisible(css)
}

#' Predefined theme: advanced grey (light background, dark text)
#' @examples
#' \dontrun{
#' use_app_theme_grey()
#' }
#' @export
use_app_theme_grey <- function() {
  set_prog_app_theme()
}

#' Predefined theme: dark (original dark background)
#' @export
#' @examples
#' \dontrun{
#' use_app_theme_dark()
#' }
use_app_theme_dark <- function() {
  set_prog_app_theme(
    primary_color = "#58a6ff",
    background_color = "#0d1117",
    sidebar_color = "#161b22",
    box_background = "#1c2128",
    text_color = "#e6edf3",
    label_color = "#58a6ff",
    run_button_gradient_start = "#1f6feb",
    run_button_gradient_end = "#388bfd",
    risk_high_color = "#f85149",
    risk_medium_color = "#d29922",
    risk_low_color = "#3fb950",
    table_header_color = "#58a6ff",
    table_row_hover_color = "#21262d",
    font_family = NULL,
    font_size_base = 14
  )
}

#' Predefined theme: light (clean white, blue accents)
#' @examples
#' \dontrun{
#' use_app_theme_light()
#' }
#' @export
use_app_theme_light <- function() {
  set_prog_app_theme(
    primary_color = "#007bff",
    background_color = "#ffffff",
    sidebar_color = "#f8f9fa",
    box_background = "#ffffff",
    text_color = "#212529",
    label_color = "#007bff",
    run_button_gradient_start = "#007bff",
    run_button_gradient_end = "#0056b3",
    risk_high_color = "#dc3545",
    risk_medium_color = "#fd7e14",
    risk_low_color = "#28a745",
    table_header_color = "#007bff",
    table_row_hover_color = "#f1f3f5",
    font_family = NULL,
    font_size_base = 14
  )
}

# ---- 2. Text customization ----
#' Set custom text for the Prognosis Terminal
#' @param ... Named arguments for text keys (see default list)
#' @examples
#' \dontrun{
#' set_prog_app_text(title = "My Prognosis App", prediction_portal = "Risk Calculator")
#' }
#' @export
set_prog_app_text <- function(...) {
  default_text <- list(
    title = "Prognosis Terminal",
    prediction_portal = "Prediction Portal",
    model_info = "Model Info",
    documentation = "Documentation",
    overview_title = "Project Overview",
    abstract = "Abstract",
    reference = "Reference",
    abstract_text = "Prognostic risk stratification from a validated survival model.",
    citation_text = "Icare R package - PrognosiX framework",
    input_box_title = "1. Input Samples",
    risk_strat_label = "Risk Stratification",
    median_choice = "Median -- High/Low (2 groups)",
    tertile_choice = "Tertile -- Low/Med/High (3 groups)",
    custom_choice = "Custom thresholds (e.g., 30, 60)",
    custom_thresholds_label = "Thresholds (comma-separated)",
    show_scores_check = "Show raw risk scores",
    batch_help = "Batch CSV: first column = SampleID; remaining columns = features.",
    batch_tab = "Batch Upload (CSV)",
    single_tab = "Single Sample",
    sample_id_label = "Sample ID:",
    upload_button_label = "Upload .csv",
    download_template_label = "Download Template",
    calculate_button = "CALCULATE RISK",
    results_title = "2. Results",
    risk_group_heading = "Risk Group",
    sample_table_heading = "Sample Table",
    export_csv_button = "Export CSV",
    model_summary_title = "Model Summary",
    algorithm_label = "Algorithm",
    cv_cindex_label = "CV C-index",
    training_n_label = "Training N",
    events_label = "Events",
    selected_features_title = "Selected Features",
    variable_glossary_title = "Variable Glossary",
    feature_col = "Feature",
    description_col = "Description",
    units_col = "Units"
  )
  user_text <- list(...)
  final_text <- utils::modifyList(default_text, user_text)
  options("prog_app_text" = final_text)
  invisible(final_text)
}

#' Get Custom Application Text
#'
#' Retrieves a specific text string from the global application text options.
#' This is an internal helper used by the PrognosiX deployment Shiny app to
#' fetch user-customized labels and messages.
#'
#' @param key Character string specifying the text key to retrieve
#'   (e.g., `"title"`, `"calculate_button"`). If `NULL`, the entire
#'   text list is returned.
#'
#' @return If `key` is provided, returns the corresponding text string,
#'   or an empty string if the key does not exist. If `key` is `NULL`,
#'   returns the full named list of all application texts.
#'
#' @keywords internal
#' @noRd
get_prog_app_text <- function(key = NULL) {
  txt <- getOption("prog_app_text", list())
  if (is.null(key)) return(txt)
  return(txt[[key]] %||% "")
}

# ---- 3. Core prediction functions ----
#' Robust Prediction for PrognosiX Objects with Imputation
#'
#' Predict prognostic risk scores while automatically handling data type mismatches and missing values.
#'
#' @param prog_obj A \code{PrognosiX} object containing a trained survival model.
#' @param newdata A data frame (or data.table) with the same features as used for training.
#' @param impute Logical. Should missing values in \code{newdata} be imputed using the
#'   reference data (median for numeric, mode for factor/character)? Default is \code{TRUE}.
#'
#' @return A numeric vector of risk scores (crank) for each observation in \code{newdata}.
#'   Higher values indicate higher predicted risk.
#'
#' @details
#' The function first extracts the training task from the \code{PrognosiX} object to obtain
#' the exact feature types (integer, numeric, factor) and factor levels. It then coerces
#' the columns of \code{newdata} to these types. If \code{impute = TRUE}, any remaining
#' \code{NA} values are filled using the reference data stored in \code{prog_obj@clean.data}
#' (median for numeric, mode for categorical). Dummy \code{time} and \code{status} columns
#' are added to satisfy \code{TaskSurv} requirements, and the learner's \code{predict}
#' method is called to obtain risk scores.
#'
#' This function resolves the \code{Mlr3ErrorInput} issue that often occurs when
#' the original \code{predict_prognosix} is used with a learner like \code{surv.ranger}.
#'
#' @examples
#' \dontrun{
#' # Assume 'prog' is a trained PrognosiX object and 'new_data' contains the features
#' risk <- predict_prognosix_robust(prog, new_data, impute = TRUE)
#' head(risk)
#' }
#'
#' @importFrom mlr3proba TaskSurv
#' @importFrom data.table as.data.table
#' @export
predict_prognosix <- function(prog_obj, newdata, impute = TRUE) {
  if (!inherits(prog_obj, "PrognosiX")) {
    stop("prog_obj must be a PrognosiX object.")
  }
  if (length(prog_obj@best.model) == 0) {
    stop("No best.model found in the PrognosiX object.")
  }
  
  # Extract training task and features
  train_task <- surv_extract_task(prog_obj)
  features <- prog_obj@best.model$features
  train_task$select(features)
  
  # Check that all required features are present in newdata
  newdata <- as.data.frame(newdata)
  missing_feats <- setdiff(features, colnames(newdata))
  if (length(missing_feats) > 0) {
    stop("Missing features: ", paste(missing_feats, collapse = ", "))
  }
  
  # Subset to training features
  new_clean <- newdata[, features, drop = FALSE]
  
  # Reference data for imputation (original cleaned data)
  ref <- prog_obj@clean.data[, features, drop = FALSE]
  
  # Coerce column types based on training task and optionally impute
  for (feat in features) {
    # Determine target type from training task
    target_type <- train_task$feature_types[train_task$feature_types$id == feat, "type"]
    
    # Convert to the exact storage type expected by the learner
    if (target_type == "integer") {
      new_clean[[feat]] <- as.integer(new_clean[[feat]])
    } else if (target_type == "numeric") {
      new_clean[[feat]] <- as.numeric(new_clean[[feat]])
    } else if (target_type == "factor") {
      # Obtain factor levels from the training task
      lev <- train_task$levels(feat)[[1]]
      new_clean[[feat]] <- factor(as.character(new_clean[[feat]]), levels = lev)
    } else {
      # fallback for other types (character, logical, etc.)
      new_clean[[feat]] <- new_clean[[feat]]
    }
    
    # Impute missing values if requested
    if (impute && anyNA(new_clean[[feat]])) {
      if (target_type %in% c("numeric", "integer")) {
        fill <- median(ref[[feat]], na.rm = TRUE)
        new_clean[[feat]][is.na(new_clean[[feat]])] <- fill
      } else if (target_type == "factor") {
        fill <- names(sort(table(ref[[feat]]), decreasing = TRUE))[1]
        new_clean[[feat]][is.na(new_clean[[feat]])] <- fill
        # Ensure factor levels remain correct after replacement
        new_clean[[feat]] <- factor(new_clean[[feat]], levels = lev)
      } else {
        # For character or other types, use mode
        fill <- names(sort(table(ref[[feat]]), decreasing = TRUE))[1]
        new_clean[[feat]][is.na(new_clean[[feat]])] <- fill
      }
    }
  }
  
  # Add dummy time and status columns (required by TaskSurv)
  new_clean$time <- 1L
  new_clean$status <- 0L
  
  # Create prediction task
  pred_task <- mlr3proba::TaskSurv$new(
    id = "pred_task",
    backend = data.table::as.data.table(new_clean),
    time = "time",
    event = "status"
  )
  
  # Predict risk scores
  learner <- prog_obj@best.model$learner
  if ("crank" %in% learner$predict_types) {
    learner$predict_type <- "crank"
  }
  risk_scores <- learner$predict(pred_task)$crank
  
  return(risk_scores)
}

#' Predict risk groups with arbitrary thresholds
#' @param prog_obj A `PrognosiX` object
#' @param newdata Data frame of new samples
#' @param cutoff_method `"median"`, `"tertile"`, or `"custom"`
#' @param custom_cutoffs Numeric vector (length >=1). For length >2, labels become "Low Risk", "Medium Risk 1", ..., "High Risk"
#' @param return_scores Logical, include raw risk scores
#' @return Data frame with SampleID, risk_group, and optionally risk_score
#' @examples
#' \dontrun{
#' # Requires a trained PrognosiX object
#' # risk_df <- predict_risk_groups(prog_obj, new_data, cutoff_method = "median")
#' }
#' @export
predict_risk_groups <- function(prog_obj, newdata, 
                                cutoff_method = c("median", "tertile", "custom"),
                                custom_cutoffs = NULL,
                                return_scores = TRUE) {
  cutoff_method <- match.arg(cutoff_method)
  risk_scores <- predict_prognosix(prog_obj, newdata, impute = TRUE)
  
  train_data <- prog_obj@survival.data
  train_feats <- prog_obj@best.model$features %||% colnames(prog_obj@clean.data)
  train_scores <- predict_prognosix(prog_obj, train_data[, train_feats, drop = FALSE], impute = FALSE)
  
  cuts <- switch(cutoff_method,
                 median = median(train_scores, na.rm = TRUE),
                 tertile = quantile(train_scores, probs = c(1/3, 2/3), na.rm = TRUE),
                 custom = {
                   if (is.null(custom_cutoffs))
                     stop("custom_cutoffs required")
                   sort(custom_cutoffs)
                 })
  n_groups <- length(cuts) + 1
  if (n_groups == 2) {
    risk_group <- ifelse(risk_scores > cuts, "High Risk", "Low Risk")
  } else {
    labels <- c("Low Risk", paste("Medium Risk", 1:(n_groups-2)), "High Risk")
    breaks <- c(-Inf, cuts, Inf)
    risk_group <- as.character(cut(risk_scores, breaks = breaks, labels = labels))
  }
  result <- data.frame(SampleID = rownames(newdata), risk_group = risk_group, stringsAsFactors = FALSE)
  if (return_scores) result$risk_score <- round(risk_scores, 6)
  return(result)
}

# ---- 4. Deployment dispatcher and manager ----
#' Deployment dispatcher (legacy compatibility)
#' @param prog_train_obj A `PrognosiX` object
#' @param newdata Data frame of new samples
#' @param cutoff_method One of `"median"`, `"tertile"`, `"custom"`
#' @param custom_cutoffs Numeric vector for custom thresholds
#' @param return_scores Logical
#' @examples
#' \dontrun{
#' # Legacy compatibility wrapper
#' # result <- Prog_deploy_dispatcher(prog_obj, new_data, cutoff_method = "tertile")
#' }
#' @export
Prog_deploy_dispatcher <- function(prog_train_obj,
                                   newdata,
                                   cutoff_method = c("median", "tertile", "custom"),
                                   custom_cutoffs = NULL,
                                   return_scores = TRUE) {
  cutoff_method <- match.arg(cutoff_method)
  predict_risk_groups(prog_train_obj, newdata, cutoff_method, custom_cutoffs, return_scores)
}

#' Create a deployment manager
#' @param prog_train_obj A `PrognosiX` object with a trained model
#' @return An S3 object of class `"Prog_Manager"`
#' @examples
#' \dontrun{
#' # Requires a trained PrognosiX object
#' # mgr <- New_Prog_Manager(prog_obj)
#' # result <- mgr$prog_predict(new_data)
#' }
#' @export
New_Prog_Manager <- function(prog_train_obj) {
  if (!inherits(prog_train_obj, "PrognosiX"))
    stop("[New_Prog_Manager] Input must be a PrognosiX object.")
  if (length(prog_train_obj@best.model) == 0)
    stop("[New_Prog_Manager] No best.model. Run a training pipeline first.")
  
  feats <- prog_train_obj@best.model$features %||% colnames(prog_train_obj@clean.data)
  info <- list(
    algorithm  = prog_train_obj@best.model$learner_id,
    cv_cindex  = prog_train_obj@best.model$cv_cindex,
    features   = feats,
    n_features = length(feats),
    n_train    = nrow(prog_train_obj@survival.data),
    n_events   = sum(prog_train_obj@survival.data[[prog_train_obj@status_col]] == 1)
  )
  mgr <- list(
    trained_obj = prog_train_obj,
    model_info  = info,
    prog_predict = function(newdata,
                            cutoff_method = c("median", "tertile", "custom"),
                            custom_cutoffs = NULL,
                            return_scores = TRUE) {
      cutoff_method <- match.arg(cutoff_method)
      predict_risk_groups(prog_train_obj, newdata, cutoff_method, custom_cutoffs, return_scores)
    }
  )
  class(mgr) <- "Prog_Manager"
  cat("\n-- PrognosiX Manager ---------------------------------\n")
  cat(sprintf("  Algorithm  : %s\n",   info$algorithm))
  cat(sprintf("  CV C-index : %.4f\n", info$cv_cindex %||% NA))
  cat(sprintf("  Features   : %d  (%s%s)\n", info$n_features,
              paste(head(info$features, 4), collapse = ", "),
              if (info$n_features > 4) ", ..." else ""))
  cat(sprintf("  Training N : %d (events: %d)\n", info$n_train, info$n_events))
  cat("------------------------------------------------------\n\n")
  return(mgr)
}

# ---- 5. Shiny application ----
#' Launch the interactive Prognosis Terminal
#' @param prog_manager A `Prog_Manager` object
#' @param title Browser title (overrides custom text)
#' @param var_dict Data frame with columns Feature, Description, Units
#' @param project_info List with elements `abstract` and `citation`
#' @return Runs the Shiny app
#' @examples
#' \dontrun{
#' # Requires a Prog_Manager object
#' # mgr <- New_Prog_Manager(prog_obj)
#' # launch_prog_deploy_app(mgr)
#' }
#' @export
launch_prog_deploy_app <- function(prog_manager,
                                   title = NULL,
                                   var_dict = NULL,
                                   project_info = NULL) {
  if (!inherits(prog_manager, "Prog_Manager"))
    stop("Input must be a Prog_Manager object.")
  for (pkg in c("shiny", "shinydashboard", "DT"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("'%s' required: install.packages('%s')", pkg, pkg))
  # Apply default theme if none set
  if (is.null(getOption("prog_app_theme_css"))) use_app_theme_grey()
  
  txt <- get_prog_app_text()
  if (is.null(title)) title <- txt$title %||% "Prognosis Terminal"
  .get_text <- function(key, def) { val <- txt[[key]]; if (is.null(val)) def else val }
  
  info <- prog_manager$model_info
  train_data <- prog_manager$trained_obj@clean.data
  req_vars <- info$features
  default_vals <- sapply(train_data[, req_vars, drop = FALSE],
                         function(x) round(median(as.numeric(x), na.rm = TRUE), 3))
  
  ui <- dashboardPage(
    skin = "black",
    dashboardHeader(title = span(icon("heartbeat"), " ", title)),
    dashboardSidebar(sidebarMenu(
      menuItem(.get_text("prediction_portal", "Prediction Portal"), tabName = "portal", icon = icon("desktop")),
      menuItem(.get_text("model_info", "Model Info"), tabName = "model", icon = icon("chart-bar")),
      menuItem(.get_text("documentation", "Documentation"), tabName = "docs", icon = icon("info-circle"))
    )),
    dashboardBody(
      tags$head(tags$style(HTML(getOption("prog_app_theme_css")))),
      tabItems(
        tabItem(tabName = "portal",
                fluidRow(box(width = 12, title = .get_text("overview_title", "Project Overview"),
                             column(8, h4(.get_text("abstract", "Abstract")),
                                    p(project_info$abstract %||% .get_text("abstract_text", "Prognostic risk stratification from a validated survival model."))),
                             column(4, h4(.get_text("reference", "Reference")),
                                    p(project_info$citation %||% .get_text("citation_text", "Icare R package - PrognosiX framework")))
                )),
                fluidRow(box(title = .get_text("input_box_title", "1. Input Samples"), width = 12,
                             column(3,
                                    selectInput("cutoff_m", .get_text("risk_strat_label", "Risk Stratification"),
                                                choices = {
                                                  lbl_med <- .get_text("median_choice", "Median -- High/Low (2 groups)")
                                                  lbl_ter <- .get_text("tertile_choice", "Tertile -- Low/Med/High (3 groups)")
                                                  lbl_cus <- .get_text("custom_choice", "Custom thresholds (e.g., 0.3, 0.6)")
                                                  setNames(c("median", "tertile", "custom"), c(lbl_med, lbl_ter, lbl_cus))
                                                }),
                                    conditionalPanel(condition = "input.cutoff_m == 'custom'",
                                                     textInput("custom_thresholds", .get_text("custom_thresholds_label", "Thresholds (comma-separated)"),
                                                               value = "0.3, 0.6", placeholder = "e.g., 0.2, 0.5, 0.8")),
                                    checkboxInput("show_sc", .get_text("show_scores_check", "Show raw risk scores"), TRUE),
                                    hr(),
                                    helpText(.get_text("batch_help", "Batch CSV: first column = SampleID; remaining columns = features."))
                             ),
                             column(9,
                                    tabsetPanel(id = "input_mode",
                                                tabPanel(.get_text("batch_tab", "Batch Upload (CSV)"), br(),
                                                         fileInput("up_file", .get_text("upload_button_label", "Upload .csv"), accept = ".csv"),
                                                         downloadButton("dl_tpl", .get_text("download_template_label", "Download Template"), class = "btn-xs btn-info")),
                                                tabPanel(.get_text("single_tab", "Single Sample"), br(),
                                                         fluidRow(column(4, textInput("sid", .get_text("sample_id_label", "Sample ID:"), "SAMPLE_001")),
                                                                  lapply(seq_along(req_vars), function(i)
                                                                    column(4, numericInput(paste0("f_", req_vars[i]), req_vars[i],
                                                                                           value = default_vals[[i]], step = 0.01))))
                                                )
                                    )
                             )
                )),
                fluidRow(column(12, align = "center",
                                actionButton("go", .get_text("calculate_button", "CALCULATE RISK"), icon = icon("play-circle"), class = "btn-run"))),
                fluidRow(box(title = .get_text("results_title", "2. Results"), width = 12,
                             column(4, align = "center", h4(.get_text("risk_group_heading", "Risk Group")),
                                    uiOutput("risk_ui"), hr(), uiOutput("score_ui")),
                             column(8,
                                    div(style = "display:flex;justify-content:space-between;align-items:center",
                                        h4(.get_text("sample_table_heading", "Sample Table")),
                                        downloadButton("dl_res", .get_text("export_csv_button", "Export CSV"), class = "btn-success btn-xs")),
                                    DT::dataTableOutput("res_tbl"))
                ))
        ),
        tabItem(tabName = "model",
                fluidRow(box(title = .get_text("model_summary_title", "Model Summary"), width = 12,
                             fluidRow(
                               column(3, div(class="badge-box", icon("brain"), " ", .get_text("algorithm_label", "Algorithm"), br(), tags$b(info$algorithm))),
                               column(3, div(class="badge-box", icon("chart-line"), " ", .get_text("cv_cindex_label", "CV C-index"), br(), tags$b(round(info$cv_cindex %||% NA, 4)))),
                               column(3, div(class="badge-box", icon("users"), " ", .get_text("training_n_label", "Training N"), br(), tags$b(info$n_train))),
                               column(3, div(class="badge-box", icon("flag"), " ", .get_text("events_label", "Events"), br(), tags$b(info$n_events)))
                             ),
                             br(), h4(.get_text("selected_features_title", "Selected Features")),
                             DT::dataTableOutput("feat_tbl")
                ))
        ),
        tabItem(tabName = "docs",
                fluidRow(box(title = .get_text("variable_glossary_title", "Variable Glossary"), width = 12,
                             DT::dataTableOutput("doc_tbl")))
        )
      )
    )
  )
  
  server <- function(input, output, session) {
    parsed_data <- eventReactive(input$go, {
      if (input$input_mode == "Batch Upload (CSV)") {
        req(input$up_file)
        read.csv(input$up_file$datapath, row.names = 1, check.names = FALSE)
      } else {
        vals <- setNames(sapply(req_vars, function(v) input[[paste0("f_", v)]]), req_vars)
        df <- as.data.frame(t(vals)); rownames(df) <- input$sid; df
      }
    })
    preds <- reactive({
      req(parsed_data())
      custom_cutoffs <- NULL
      if (input$cutoff_m == "custom") {
        thr_str <- gsub(" ", "", input$custom_thresholds)
        custom_cutoffs <- as.numeric(strsplit(thr_str, ",")[[1]])
        if (any(is.na(custom_cutoffs))) showNotification("Invalid custom thresholds", type = "error")
      }
      prog_manager$prog_predict(parsed_data(),
                                cutoff_method = input$cutoff_m,
                                custom_cutoffs = custom_cutoffs,
                                return_scores = input$show_sc)
    })
    output$risk_ui <- renderUI({
      validate(need(preds(), "Awaiting input..."))
      df <- preds()
      if (nrow(df) == 1) {
        g <- df$risk_group[1]
        cl <- if (grepl("High", g)) "risk-high" else if (grepl("Med", g)) "risk-med" else "risk-low"
        tagList(div(class = "score-lbl", "Sample: ", df$SampleID[1]), div(class = cl, g))
      } else {
        tbl <- table(df$risk_group)
        tagList(div(class = "score-lbl", sprintf("Batch: %d samples processed", nrow(df))), br(),
                lapply(names(tbl), function(g) {
                  cl <- if (grepl("High", g)) "risk-high" else if (grepl("Med", g)) "risk-med" else "risk-low"
                  div(span(class = cl, style = "font-size:28px", g),
                      span(style = "color:#8b949e;font-size:13px;margin-left:8px", paste0("n=", tbl[[g]])))
                }))
      }
    })
    output$score_ui <- renderUI({
      validate(need(preds(), ""))
      df <- preds()
      if (nrow(df) == 1 && "risk_score" %in% names(df))
        div(class = "score-lbl", sprintf("Risk Score: %.4f", df$risk_score[1]))
    })
    output$res_tbl <- DT::renderDataTable({
      req(preds())
      DT::datatable(preds(), rownames = FALSE, options = list(dom = "tp", pageLength = 10)) %>%
        DT::formatStyle("risk_group",
                        color = DT::styleEqual(c("High Risk","Medium Risk","Low Risk"),
                                               c("#f85149","#d29922","#3fb950")), fontWeight = "bold")
    })
    output$dl_res <- downloadHandler(filename = function() paste0("Prognosis_Results_", Sys.Date(), ".csv"),
                                     content = function(f) write.csv(preds(), f, row.names = FALSE))
    output$dl_tpl <- downloadHandler(filename = "Input_Template.csv",
                                     content = function(f) write.csv(head(train_data[, req_vars, drop = FALSE], 5), f))
    output$feat_tbl <- DT::renderDataTable(
      DT::datatable(data.frame(`#` = seq_along(req_vars), Feature = req_vars,
                               Median_Training = round(default_vals, 3), check.names = FALSE),
                    rownames = FALSE, options = list(dom = "t", pageLength = 30)))
    output$doc_tbl <- DT::renderDataTable({
      df <- var_dict %||% data.frame(Feature = req_vars,
                                     Description = paste("Predictor:", req_vars), Units = "--", stringsAsFactors = FALSE)
      colnames(df) <- c(.get_text("feature_col", "Feature"),
                        .get_text("description_col", "Description"),
                        .get_text("units_col", "Units"))
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 30))
    })
  }
  shinyApp(ui, server)
}
