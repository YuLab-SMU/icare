# =============================================================================
# model_feature_selection.R
# Comprehensive Feature Selection Module for Train_Model objects
# Integrates RFE, GA, SA, and SBF methods from caret with automated 
# cross-validation and visualization.
# =============================================================================

# 0. Package check ------------------------------------------------------------
.check_fs_packages <- function() {
  required <- c("caret", "ggplot2", "wesanderson", "ggprism", "dplyr", "tidyr",
                "reshape2", "gridExtra", "doParallel", "foreach")
  missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Install them before using feature selection.")
  }
  invisible(TRUE)
}

# -- Internal helpers ---------------------------------------------------------
.get_output_dir <- function(...) {
  get_output_dir(...)
}

.safe_dir <- function(save_dir) {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  invisible(save_dir)
}

.extract_xy <- function(object) {
  # Extract predictors (x) and outcome (y) from a Train_Model object
  if (inherits(object, "Train_Model")) {
    cd <- object@clean.df
    gc <- object@group_col
  } else if (inherits(object, "Stat")) {
    cd <- object@clean.data
    gc <- object@group_col
  } else {
    stop("Object must be Train_Model or Stat.")
  }
  x <- cd[, setdiff(colnames(cd), gc), drop = FALSE]
  y <- cd[[gc]]
  list(x = x, y = y, group_col = gc)
}

# -- 1. RFE: Recursive Feature Elimination ------------------------------------

#' Recursive Feature Elimination (RFE) for Train_Model objects
#'
#' @param object   A Train_Model or Stat S4 object.
#' @param sizes    Numeric vector of subset sizes to evaluate.
#' @param rfe_func A caret RFE function list (e.g., \code{rfFuncs}).
#' @param method   External resampling: "cv", "repeatedcv", "boot".
#' @param number   Number of folds / resampling iterations.
#' @param repeats  For repeatedcv only.
#' @param metric   Evaluation metric: "Accuracy", "Kappa", "ROC".
#' @param allowParallel Logical.
#' @param save_plot Logical. Save RFE profile plot?
#' @param save_dir  Directory to save outputs.
#' @param seed      Random seed.
#' @return A list with elements: \code{result} (rfe object), 
#'   \code{opt_vars} (optimal features), \code{plot} (ggplot).
#' @export
FeatureSelectRFE <- function(object,
                             sizes    = NULL,
                             rfe_func = caret::rfFuncs,
                             method   = "repeatedcv",
                             number   = 5,
                             repeats  = 2,
                             metric   = "Accuracy",
                             allowParallel = FALSE,
                             save_plot = FALSE,
                             save_dir  = NULL,
                             seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  if (is.null(sizes))
    sizes <- unique(round(seq(2, min(30, ncol(xy$x)), length.out = 10)))
  
  ctrl <- caret::rfeControl(
    functions = rfe_func,
    method    = method,
    number    = number,
    repeats   = repeats,
    verbose   = TRUE,
    allowParallel = allowParallel
  )
  
  cat(sprintf("Running RFE: %d sizes, %s (%d-fold, %d repeats)...\n",
              length(sizes), method, number, repeats))
  
  rfe_res <- caret::rfe(xy$x, xy$y, sizes = sizes, rfeControl = ctrl, metric = metric)
  
  opt_vars <- caret::predictors(rfe_res)
  cat(sprintf("RFE selected: %d features\n", length(opt_vars)))
  
  # Build plot
  p <- ggplot2::ggplot(rfe_res, metric = metric) + 
    ggplot2::theme_bw() +
    ggplot2::labs(title = "RFE - Recursive Feature Elimination",
                  subtitle = paste0("Optimal: ", length(opt_vars), " variables"))
  print(p)
  
  if (save_plot) {
    .safe_dir(save_dir)
    ggplot2::ggsave(file.path(save_dir, "rfe_profile.pdf"), p,
                    width = 7, height = 5, device = "pdf")
  }
  
  invisible(list(result = rfe_res, opt_vars = opt_vars, plot = p))
}


# -- 2. GA: Genetic Algorithm Feature Selection -------------------------------

#' Genetic Algorithm Feature Selection for Train_Model objects
#'
#' @param object   A Train_Model or Stat object.
#' @param iters    Number of GA generations.
#' @param popSize  Population size.
#' @param ga_func  A caret GA function list (e.g., \code{caretGA}, \code{rfGA}).
#' @param method   External resampling method.
#' @param number   Number of folds.
#' @param repeats  Repeats.
#' @param metric   Internal fitness metric.
#' @param allowParallel,genParallel Logical.
#' @param save_plot Logical.
#' @param save_dir  Directory.
#' @param seed      Random seed.
#' @return List with \code{result} (gafs object), \code{opt_vars}, \code{plot}.
#' @export
FeatureSelectGA <- function(object,
                            iters    = 10,
                            popSize  = 20,
                            ga_func  = caret::caretGA,
                            method   = "repeatedcv",
                            number   = 5,
                            repeats  = 2,
                            metric   = "Accuracy",
                            allowParallel = FALSE,
                            genParallel   = FALSE,
                            save_plot = FALSE,
                            save_dir  = NULL,
                            seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  ctrl <- caret::gafsControl(
    functions      = ga_func,
    method         = method,
    number         = number,
    repeats        = repeats,
    verbose        = TRUE,
    allowParallel  = allowParallel,
    genParallel    = genParallel
  )
  
  cat(sprintf("Running GA: %d iters, popSize %d, %s...\n",
              iters, popSize, method))
  
  ga_res <- caret::gafs(xy$x, xy$y, iters = iters, popSize = popSize,
                        gafsControl = ctrl, metric = metric)
  
  opt_vars <- ga_res$optVariables
  cat(sprintf("GA selected: %d features\n", length(opt_vars)))
  
  p <- plot(ga_res) + ggplot2::theme_bw() +
    ggplot2::labs(title = "GA - Genetic Algorithm Feature Selection")
  print(p)
  
  if (save_plot) {
    .safe_dir(save_dir)
    ggplot2::ggsave(file.path(save_dir, "ga_trace.pdf"), p,
                    width = 7, height = 5, device = "pdf")
  }
  
  invisible(list(result = ga_res, opt_vars = opt_vars, plot = p))
}


# -- 3. SA: Simulated Annealing Feature Selection -----------------------------

#' Simulated Annealing Feature Selection for Train_Model objects
#'
#' @param object   A Train_Model or Stat object.
#' @param iters    Number of SA iterations.
#' @param sa_func  A caret SA function list (e.g., \code{caretSA}, \code{rfSA}).
#' @param method   External resampling method.
#' @param number   Number of folds.
#' @param repeats  Repeats.
#' @param metric   Internal fitness metric.
#' @param improve  SA restart after `improve` iters without improvement.
#' @param allowParallel Logical.
#' @param save_plot Logical.
#' @param save_dir  Directory.
#' @param seed      Random seed.
#' @return List with \code{result} (safs object), \code{opt_vars}, \code{plot}.
#' @export
FeatureSelectSA <- function(object,
                            iters    = 25,
                            sa_func  = caret::caretSA,
                            method   = "repeatedcv",
                            number   = 5,
                            repeats  = 2,
                            metric   = "Accuracy",
                            improve  = 3L,
                            allowParallel = FALSE,
                            save_plot = FALSE,
                            save_dir  = NULL,
                            seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  ctrl <- caret::safsControl(
    functions      = sa_func,
    method         = method,
    number         = number,
    repeats        = repeats,
    improve        = improve,
    verbose        = TRUE,
    allowParallel  = allowParallel
  )
  
  cat(sprintf("Running SA: %d iters, %s...\n", iters, method))
  
  sa_res <- caret::safs(xy$x, xy$y, iters = iters,
                        safsControl = ctrl, metric = metric)
  
  opt_vars <- sa_res$optVariables
  cat(sprintf("SA selected: %d features\n", length(opt_vars)))
  
  p <- plot(sa_res) + ggplot2::theme_bw() +
    ggplot2::labs(title = "SA - Simulated Annealing Feature Selection")
  print(p)
  
  if (save_plot) {
    .safe_dir(save_dir)
    ggplot2::ggsave(file.path(save_dir, "sa_trace.pdf"), p,
                    width = 7, height = 5, device = "pdf")
  }
  
  invisible(list(result = sa_res, opt_vars = opt_vars, plot = p))
}


# -- 4. SBF: Selection By Univariate Filter -----------------------------------

#' Univariate Filter Feature Selection for Train_Model objects
#'
#' @param object   A Train_Model or Stat object.
#' @param sbf_func A caret SBF function list (e.g., \code{rfSBF}, \code{caretSBF}).
#' @param method   External resampling method.
#' @param number   Number of folds.
#' @param repeats  Repeats.
#' @param metric   Evaluation metric.
#' @param allowParallel Logical.
#' @param save_plot Logical (saves variable count barplot).
#' @param save_dir  Directory.
#' @param seed      Random seed.
#' @return List with \code{result} (sbf object), \code{opt_vars}, \code{plot}.
#' @export
FeatureSelectSBF <- function(object,
                             sbf_func = caret::rfSBF,
                             method   = "repeatedcv",
                             number   = 5,
                             repeats  = 2,
                             metric   = "Accuracy",
                             allowParallel = FALSE,
                             save_plot = FALSE,
                             save_dir  = NULL,
                             seed      = 825) {
  .check_fs_packages()
  set.seed(seed)
  
  if (is.null(save_dir)) 
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  xy <- .extract_xy(object)
  
  if (allowParallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() })
  }
  
  ctrl <- caret::sbfControl(
    functions = sbf_func,
    method    = method,
    number    = number,
    repeats   = repeats,
    verbose   = TRUE
  )
  
  cat(sprintf("Running SBF: %s...\n", method))
  
  sbf_res <- caret::sbf(xy$x, xy$y, sbfControl = ctrl, metric = metric)
  
  opt_vars <- sbf_res$optVariables
  cat(sprintf("SBF selected: %d features\n", length(opt_vars)))
  
  p <- plot(sbf_res) + ggplot2::theme_bw() +
    ggplot2::labs(title = "SBF - Univariate Filter Feature Selection")
  print(p)
  
  if (save_plot) {
    .safe_dir(save_dir)
    ggplot2::ggsave(file.path(save_dir, "sbf_filter.pdf"), p,
                    width = 7, height = 5, device = "pdf")
  }
  
  invisible(list(result = sbf_res, opt_vars = opt_vars, plot = p))
}


# -- 5. Built-In Importance (from any caret train object) ---------------------
#' Explain Model Performance (ROC, Lift, Boxplot)
#'
#' @param explainer DALEX explainer object.
#' @param geom Plot type: "roc", "lift", or "boxplot".
#' @param save_plots Logical, save plot to file.
#' @param save_dir Directory to save the plot.
#' @param plot_width,plot_height Dimensions in inches.
#' @param ... Additional arguments passed to `ggplot2::ggsave`.
#' @return Invisibly, the DALEX model_performance object.
#' @export
ExplainModelPerformance <- function(explainer,
                                    geom        = c("roc", "lift", "boxplot"),
                                    save_plots  = FALSE,
                                    save_dir    = "ModelExplain",
                                    plot_width  = 6,
                                    plot_height = 5,
                                    ...) {
  if (!requireNamespace("DALEX", quietly = TRUE))
    stop("Package 'DALEX' is required.")
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Package 'pROC' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggprism", quietly = TRUE))
    stop("Package 'ggprism' is required.")
  
  geom <- match.arg(geom)
  
  cat("-- Model Performance ------------------------------------------------------\n")
  mp <- DALEX::model_performance(explainer)
  print(mp)
  auc_val <- mp$measures$auc[1]
  
  # ---- ROC curve (custom ggplot) ----
  if (geom == "roc") {
    probs <- explainer$predict_function(explainer$model, explainer$data)
    true  <- explainer$y
    if (is.factor(true)) true <- as.numeric(true) - 1
    if (!is.numeric(true)) true <- as.numeric(as.factor(true)) - 1
    roc_obj <- pROC::roc(true, probs, levels = c(0, 1), direction = "auto", quiet = TRUE)
    roc_df <- data.frame(
      Sensitivity = roc_obj$sensitivities,
      Specificity = roc_obj$specificities
    )
    
    p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - Specificity, y = Sensitivity)) +
      ggplot2::geom_line(color = "#006d2c", linewidth = 1.2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
      ggplot2::labs(
        title    = paste("Model Performance -", explainer$label),
        subtitle = paste("AUC =", round(auc_val, 3)),
        x        = "1 - Specificity",
        y        = "Sensitivity"
      ) +
      ggplot2::coord_equal() +
      ggprism::theme_prism(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    
    print(p)
    
    if (save_plots) {
      if (is.null(save_dir) || length(save_dir) == 0 || !nzchar(save_dir)) {
        save_dir <- "ModelExplain"
      }
      save_dir <- normalizePath(save_dir, mustWork = FALSE)
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      }
      file_path <- file.path(save_dir, "performance_roc.pdf")
      ggplot2::ggsave(filename = file_path, plot = p,
                      width = plot_width, height = plot_height, ...)
      cat("ROC plot saved to:", file_path, "\n")
    }
  } else {
    # Lift or boxplot using DALEX plot
    p <- plot(mp, geom = geom) +
      ggplot2::labs(title = paste("Model Performance -", explainer$label)) +
      ggprism::theme_prism(base_size = 13) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    
    print(p)
    
    if (save_plots) {
      if (is.null(save_dir) || length(save_dir) == 0 || !nzchar(save_dir)) {
        save_dir <- "ModelExplain"
      }
      save_dir <- normalizePath(save_dir, mustWork = FALSE)
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      }
      file_path <- file.path(save_dir, paste0("performance_", geom, ".pdf"))
      ggplot2::ggsave(filename = file_path, plot = p,
                      width = plot_width, height = plot_height, ...)
      cat("Plot saved to:", file_path, "\n")
    }
  }
  
  invisible(mp)
}


# -- 6. Unified Wrapper: run all methods and benchmark ------------------------
#' Comprehensive Feature Selection with Method Benchmarking
#'
#' Runs up to four caret-based selection methods on the same `Train_Model`
#' object, returns the union or intersection of selected features, and
#' generates a comparative Upset plot.
#'
#' @param object   A Train_Model or Stat S4 object.
#' @param methods  Character vector: "rfe","ga","sa","sbf" (default all four).
#'                Can also pass a named list of pre-run selection results.
#' @param combine  How to combine across methods: "union", "intersect", or "none".
#' @param rfe_args List of extra args passed to \code{FeatureSelectRFE}.
#' @param ga_args  List of extra args passed to \code{FeatureSelectGA}.
#' @param sa_args  List of extra args passed to \code{FeatureSelectSA}.
#' @param sbf_args List of extra args passed to \code{FeatureSelectSBF}.
#' @param upset_plot Logical. Draw an Upset plot of selected features.
#' @param save_plot  Logical.
#' @param save_dir   Directory.
#' @param seed       Random seed.
#'
#' @return An invisible list with:
#'   \describe{
#'   \item{results}{Named list of selection results per method.}
#'   \item{selected_features}{Character vector of final features.}
#'   \item{feature_matrix}{Binary matrix used for Upset.}
#'   \item{upset_plot}{ggplot (if \code{upset_plot = TRUE}).}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' fs <- FeatureSelectionPipeline(
#'   object   = model_obj,
#'   methods  = c("rfe", "ga"),
#'   combine  = "union",
#'   rfe_args = list(method = "cv", number = 5),
#'   ga_args  = list(iters = 10, popSize = 20)
#' )
#' model_obj <- SelectFeatures(model_obj, 
#'                              features = c(fs$selected_features, model_obj@group_col))
#' }
FeatureSelectionPipeline <- function(object,
                                     methods     = c("rfe", "ga", "sa", "sbf"),
                                     combine     = c("union", "intersect", "none"),
                                     rfe_args    = list(),
                                     ga_args     = list(),
                                     sa_args     = list(),
                                     sbf_args    = list(),
                                     upset_plot  = TRUE,
                                     save_plot   = FALSE,
                                     save_dir    = NULL,
                                     seed        = 825) {
  .check_fs_packages()
  combine <- match.arg(combine)
  
  if (is.null(save_dir))
    save_dir <- .get_output_dir("m2", "feature_selection")
  
  result_list <- list()
  
  # Run each requested method
  if ("rfe" %in% methods) {
    cat("\n===== RFE =====\n")
    result_list$RFE <- do.call(FeatureSelectRFE,
      c(list(object = object, save_dir = save_dir, seed = seed), rfe_args))
  }
  if ("ga" %in% methods) {
    cat("\n===== GA =====\n")
    result_list$GA <- do.call(FeatureSelectGA,
      c(list(object = object, save_dir = save_dir, seed = seed), ga_args))
  }
  if ("sa" %in% methods) {
    cat("\n===== SA =====\n")
    result_list$SA <- do.call(FeatureSelectSA,
      c(list(object = object, save_dir = save_dir, seed = seed), sa_args))
  }
  if ("sbf" %in% methods) {
    cat("\n===== SBF =====\n")
    result_list$SBF <- do.call(FeatureSelectSBF,
      c(list(object = object, save_dir = save_dir, seed = seed), sbf_args))
  }
  
  # Compile feature lists
  feature_lists <- lapply(result_list, `[[`, "opt_vars")
  
  # Combine
  if (combine == "union") {
    selected <- unique(unlist(feature_lists))
  } else if (combine == "intersect") {
    selected <- Reduce(intersect, feature_lists)
  } else {
    selected <- unlist(feature_lists)
  }
  
  cat(sprintf("\n%s of features across methods: %d\n", combine, length(selected)))
  
  # Upset plot
  p_upset <- NULL
  if (upset_plot && length(feature_lists) >= 2) {
    if (requireNamespace("UpSetR", quietly = TRUE)) {
      mat <- UpSetR::fromList(feature_lists)
      p_upset <- UpSetR::upset(mat, sets = names(feature_lists),
                               order.by = "freq",
                               text.scale = 1.2,
                               mainbar.y.label = "Number Intersected",
                               sets.x.label   = "Number Selected")
      print(p_upset)
      if (save_plot) {
        .safe_dir(save_dir)
        pdf(file.path(save_dir, "feature_upset.pdf"), width = 8, height = 5)
        print(p_upset)
        dev.off()
      }
    } else {
      cat("Install 'UpSetR' for the Upset plot.\n")
      # Fallback: simple venn-style barplot
      all_vars <- unique(unlist(feature_lists))
      bin_mat <- sapply(feature_lists, function(v) as.integer(all_vars %in% v))
      rownames(bin_mat) <- all_vars
      p_upset <- bin_mat
    }
  }
  
  invisible(list(
    results          = result_list,
    selected_features = selected,
    feature_matrix   = if (exists("bin_mat")) bin_mat else NULL,
    upset_plot       = p_upset
  ))
}


# -- 7. Quick train-and-select: built-in importance from multiple models ------
#' Built-in Feature Selection Using Model Importance
#'
#' Trains one or more classification models (e.g., Random Forest, GBM) with
#' cross-validation and extracts built-in variable importance scores. Features
#' are ranked by importance, and a final set is selected as the union or
#' intersection of the top \code{top_n} features from each model.
#'
#' @param object An object containing the training data. Must be compatible with
#'   the internal extractor \code{.extract_xy()} (e.g., a fitted model object or
#'   a data container with predictors \code{x} and response \code{y}).
#' @param models Character vector of caret model names that support built-in
#'   variable importance (e.g., \code{"rf"}, \code{"gbm"}). Only those with
#'   available importance will be used. Default: \code{c("rf", "gbm")}.
#' @param method Resampling method passed to \code{\link[caret]{trainControl}}.
#'   Default: \code{"repeatedcv"}.
#' @param number Number of folds or resampling iterations. Default: \code{5}.
#' @param top_n Integer specifying the number of top important features to
#'   retain from each model. Default: \code{15}.
#' @param combine Character string; either \code{"union"} (take the union of all
#'   selected features) or \code{"intersect"} (take the intersection). Default:
#'   \code{"union"}.
#' @param seed Random seed for reproducibility. Default: \code{825}.
#'
#' @return An invisible list with the following components:
#' \describe{
#'   \item{importance_table}{A data frame with one row per feature. The first
#'     column \code{Feature} gives the feature name (row names are identical).
#'     For each model used, there is a column with the model name, containing
#'     \code{"Yes"} if the feature was among the top \code{top_n} features for
#'     that model, or \code{"-"} otherwise. A final column \code{Selected} marks
#'     the final selected features with \code{"[OK]"}.}
#'   \item{selected_features}{A character vector of the feature names selected
#'     in the final set (union or intersection).}
#'   \item{per_model}{A named list, where each element corresponds to a model
#'     and contains a character vector of the top \code{top_n} feature names
#'     selected by that model.}
#' }
#'
#' @details
#' The function first checks which requested models actually provide built-in
#' importance (via \code{caret::varImp}). Models without support are skipped
#' with a message. For each supported model, it is trained using the specified
#' resampling scheme, and variable importance is extracted. The top
#' \code{top_n} features are retained. Finally, either the union or
#' intersection of these per-model selections is returned.
#'
#' @note
#' This function requires the \pkg{caret} package and the respective modelling
#' packages (e.g., \pkg{randomForest}, \pkg{gbm}) to be installed.
#'
#' @seealso \code{\link[caret]{varImp}}, \code{\link[caret]{train}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_model' is a pre-processed data container with x and y
#' result <- FeatureSelectBuiltin(my_model,
#'                                models = c("rf", "gbm"),
#'                                top_n = 10,
#'                                combine = "intersect")
#' print(result$importance_table)
#' selected <- result$selected_features
#' }
#'
#' @export
FeatureSelectBuiltin <- function(object,
                                 models  = c("rf", "gbm"),
                                 method  = "repeatedcv",
                                 number  = 5,
                                 top_n   = 15,
                                 combine = "union",
                                 seed    = 825) {
  .check_fs_packages()
  set.seed(seed)
 # object<-model_obj
  xy <- .extract_xy(object)
  levels(xy$y) <- make.names(levels(xy$y))
  df <- cbind(xy$x, group = xy$y)
  
  support_df <- check_varImp_availability(models)
  unsupported <- support_df$Model[!support_df$Has_BuiltIn]
  if (length(unsupported) > 0) {
    cat("Note: The following models do NOT have built-in varImp and will be skipped:\n")
    cat(paste("  -", unsupported), sep = "\n")
  }
  
  imp_list <- list()
  for (m in models) {
    if (m %in% unsupported) next
    
    model_info <- caret::getModelInfo(m, regex = FALSE)[[1]]
    needed_pkgs <- unique(c(model_info$library))
    missing_pkg <- needed_pkgs[!sapply(needed_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkg) > 0) {
      warning("Skipping model '", m, "' because package(s) missing: ",
              paste(missing_pkg, collapse = ", "))
      next
    }
    
    cat(sprintf("Training %s for importance...\n", m))
    fit <- tryCatch({
      caret::train(group ~ ., data = df, method = m,
                   trControl = caret::trainControl(
                     method = method, number = number,
                     classProbs = TRUE),
                   metric = "Accuracy")
    }, error = function(e) {
      warning("Training failed for model '", m, "': ", e$message)
      return(NULL)
    })
    
    if (is.null(fit)) next
    
    imp <- tryCatch({
      caret::varImp(fit, scale = TRUE)$importance
    }, error = function(e) {
      warning("Could not extract importance for model '", m, "': ", e$message)
      return(NULL)
    })
    
    if (is.null(imp)) next
    
    feats <- head(rownames(imp)[order(-imp[, 1])], top_n)
    imp_list[[m]] <- feats
  }
  
  if (length(imp_list) == 0) stop("No models could be trained or provided importance.")
  
  selected <- if (combine == "union") unique(unlist(imp_list))
  else Reduce(intersect, imp_list)
  
  all_feats <- unique(unlist(imp_list))
  comb_imp <- data.frame(Feature = all_feats, row.names = all_feats)
  for (m in names(imp_list)) {
    comb_imp[[m]] <- ifelse(all_feats %in% imp_list[[m]], "Yes", "-")
  }
  comb_imp$Selected <- ifelse(all_feats %in% selected, "[OK]", "")
  comb_imp <- comb_imp[order(-rowSums(comb_imp[, names(imp_list)] == "Yes")), ]
  
  cat(sprintf("Built-in (%s): %d features\n", combine, length(selected)))
  invisible(list(importance_table = comb_imp, selected_features = selected,
                 per_model = imp_list))
}


# -- 8. Helper: Update Train_Model with selected features ---------------------

#' Apply Selected Features to a Train_Model Object
#'
#' Subsets the Train_Model to keep only the chosen features, updating all
#' relevant slots. Downstream slots (split, models) are reset.
#'
#' @param object   A Train_Model object.
#' @param features Character vector of feature names to keep.
#' @return An updated Train_Model object.
#' @export
ApplyFeatureSelection <- function(object, features) {
  if (!inherits(object, "Train_Model"))
    stop("Object must be Train_Model.")
  
  gc <- object@group_col
  keep <- unique(c(features, gc))
  keep <- intersect(keep, colnames(object@clean.df))
  
  object@clean.df <- object@clean.df[, keep, drop = FALSE]
  object@data.df  <- object@data.df[, intersect(keep, colnames(object@data.df)), drop = FALSE]
  
  # Reset downstream slots
  object@split.data       <- list()
  object@split.scale.data <- list()
  object@train.models     <- list()
  object@all.results      <- list()
  object@best.model.result <- list()
  
  cat(sprintf("Applied feature selection: %d features kept.\n", length(keep) - 1))
  cat("Downstream slots (split, models, results) have been reset.\n")
  return(object)
}


#' Check Built-in Variable Importance Availability for caret Models
#'
#' @param model_names Character vector of caret model names.
#' @return A data frame with model name and whether varImp is supported.
#' @export
check_varImp_availability <- function(model_names) {
  result <- data.frame(
    Model      = model_names,
    Has_BuiltIn = sapply(model_names, function(m) {
      info <- caret::getModelInfo(m, regex = FALSE)[[1]]
      if (is.null(info)) return(FALSE)
      # Check if model has a built-in varImp method
      !is.null(info$varImp) || m %in% c(
        "ada", "AdaBag", "AdaBoost.M1", "adaboost", "bagEarth", "bagEarthGCV", 
        "bagFDA", "bagFDAGCV", "bartMachine", "blasso", "BstLm", "bstSm", 
        "C5.0", "C5.0Cost", "C5.0Rules", "C5.0Tree", "cforest", "chaid", 
        "ctree", "ctree2", "cubist", "deepboost", "earth", "enet", 
        "evtree", "extraTrees", "fda", "gamboost", "gbm_h2o", "gbm", 
        "gcvEarth", "glmnet_h2o", "glmnet", "glmStepAIC", "J48", "JRip", 
        "lars", "lars2", "lasso", "LMT", "LogitBoost", "M5", "M5Rules", 
        "msaenet", "nodeHarvest", "OneR", "ordinalNet", "ordinalRF", 
        "ORFlog", "ORFpls", "ORFridge", "ORFsvm", "pam", "parRF", "PART", 
        "penalized", "PenalizedLDA", "qrf", "ranger", "Rborist", "relaxo", 
        "rf", "rFerns", "rfRules", "rotationForest", "rotationForestCp", 
        "rpart", "rpart1SE", "rpart2", "rpartCost", "rpartScore", 
        "rqlasso", "rqnc", "RRF", "RRFglobal", "sdwd", "smda", 
        "sparseLDA", "spikeslab", "wsrf", "xgbDART", "xgbLinear", "xgbTree"
      )
    }),
    row.names = NULL
  )
  return(result)
}

