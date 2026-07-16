# =============================================================================
# model_hyperparameter_tuning.R
# Flexible Hyperparameter Tuning Module
# =============================================================================

# 0. Package Check -----------------------------------------------------------
.check_tune_pkgs <- function() {
  required <- c("caret", "rBayesianOptimization", "dplyr", "ggplot2", "ggprism")
  missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Install them before using this module.")
  }
  invisible(TRUE)
}


# =============================================================================
# Step 1: Inspect tunable parameters for any caret model
# =============================================================================

#' Inspect Tunable Hyperparameters for a caret Model
#'
#' Returns a data frame listing all tuneable parameters and their default
#' ranges (when available), so the user can review and modify them before
#' passing to \code{BuildTuningBounds}.
#'
#' @param method Character. A caret model name (e.g., "rf", "xgbTree").
#' @return A data frame with columns: \code{parameter}, \code{label}, 
#'   \code{class}, and \code{default_range} (a character hint).
#' @export
InspectHyperParams <- function(method) {
  .check_tune_pkgs()
  
  info <- caret::getModelInfo(method, regex = FALSE)[[1]]
  if (is.null(info)) stop("Model '", method, "' not found in caret.")
  params <- info$parameters
  # Add a suggestive default range based on parameter type and common practice
  params$default_range <- sapply(params$parameter, function(p) {
    switch(p,
           mtry               = "[1, floor(sqrt(n_features))]",
           n.trees            = "[50, 500]",
           nrounds            = "[50, 300]",
           max_depth          = "[2, 10]",
           interaction.depth  = "[1, 9]",
           eta                = "[0.01, 0.3]",
           shrinkage          = "[0.001, 0.1]",
           gamma              = "[0, 5]",
           colsample_bytree   = "[0.4, 1]",
           min_child_weight   = "[1, 10]",
           subsample          = "[0.5, 1]",
           sigma              = "[0.001, 0.1]",
           C                  = "[0.1, 10]",
           alpha              = "[0, 1]",
           lambda             = "[0.0001, 0.5]",
           cp                 = "[0.0001, 0.01]",
           degree             = "[1, 3]",
           nprune             = "[2, 20]",
           paste0("[0, 1]  # please set manually")
    )
  })
  
  cat("\n========================================\n")
  cat("  Hyperparameters for:", method, "\n")
  cat("========================================\n")
  print(params[, c("parameter", "label", "class", "default_range")], row.names = FALSE)
  cat("\nUse these parameters to build bounds with BuildTuningBounds().\n")
  cat("You can modify the default ranges as needed.\n\n")
  
  invisible(params)
}


# =============================================================================
# Step 2: Build a user-defined bounds list from inspected parameters
# =============================================================================

#' Build Bounds List for Bayesian Optimization
#'
#' Converts a user-provided named list of \code{c(lower, upper)} vectors into
#' the format required by \code{rBayesianOptimization::BayesianOptimization}.
#'
#' @param ... Named arguments, each a numeric vector of length 2 
#'   \code{c(lower, upper)}.
#' @return A named list of bounds.
#' @export
#'
#' @examples
#' \dontrun{
#' bounds <- BuildTuningBounds(
#'   mtry        = c(2, 10),
#'   n.trees     = c(50, 500),
#'   shrinkage   = c(0.001, 0.1) )
#'   }
BuildTuningBounds <- function(...) {
  bounds <- list(...)
  
  # Validate
  for (nm in names(bounds)) {
    if (length(bounds[[nm]]) != 2 || !is.numeric(bounds[[nm]]))
      stop("Each parameter must be a numeric vector of length 2: c(lower, upper). ",
           "Problem with: ", nm)
    if (bounds[[nm]][1] >= bounds[[nm]][2])
      stop("Lower bound must be < upper bound for: ", nm)
  }
  
  cat("Built bounds for", length(bounds), "parameters:", 
      paste(names(bounds), collapse = ", "), "\n")
  
  return(bounds)
}


# =============================================================================
# Step 3: Run Bayesian Optimization with user-defined bounds
# =============================================================================
#' Bayesian Optimization for Model Fine-Tuning
#'
#' @description This function performs hyperparameter tuning using Bayesian Optimization 
#' via the \code{rBayesianOptimization} package and \code{caret}. It is specifically 
#' optimized for classification tasks using the ROC metric.
#'
#' @param model_obj An S4 object of class 'Train_Model' containing the data and configuration.
#' @param method A string specifying the caret model method (e.g., "rf", "xgbTree").
#' @param bounds A named list defining the parameter search space (e.g., \code{list(mtry = c(1L, 10L))}).
#' @param init_points Integer. Number of initial random points for Bayesian exploration.
#' @param n_iter Integer. Number of iterations for Bayesian Optimization.
#' @param cv_folds Integer. Number of cross-validation folds.
#' @param metric A string specifying the optimization metric (default: "ROC").
#' @param summaryFun Function to calculate performance metrics (default: \code{caret::twoClassSummary}).
#' @param use_scaled Logical. Whether to use scaled training data.
#' @param sampling A string for sampling methods (e.g., "smote", "up", "down"), or NULL.
#' @param class_weights Logical. Whether to apply inverse class frequency weights.
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Whether to print detailed training logs.
#'
#' @importFrom caret train trainControl getModelInfo twoClassSummary
#' @importFrom rBayesianOptimization BayesianOptimization
#' @importFrom stats as.formula
#' @return Returns the updated \code{model_obj} with fine-tuned results.
#' @export
FineTuneModel <- function(model_obj,
                          method,
                          bounds,
                          init_points    = 15,
                          n_iter          = 30,
                          cv_folds        = 5,
                          metric          = "ROC",
                          summaryFun      = caret::twoClassSummary,
                          use_scaled      = FALSE,
                          sampling        = NULL,
                          class_weights   = FALSE,
                          seed            = 123,
                          verbose         = TRUE) {
  
  # --- 1. Environment Setup ---
  set.seed(seed)
  if (!inherits(model_obj, "Train_Model")) {
    stop("model_obj must be an object of class 'Train_Model'.")
  }
  
  # --- 2. Data Retrieval ---
  if (use_scaled) {
    train_data <- model_obj@split.scale.data$training
    if (is.null(train_data)) stop("Scaled training data not found in model_obj@split.scale.data$training.")
  } else {
    train_data <- model_obj@filtered.set$training
  }
  
  if (is.null(train_data)) stop("Training dataset is empty.")
  
  # --- 3. Target Variable Processing ---
  # Ensure the target column is a factor (required for ROC/Classification)
  gc <- model_obj@group_col
  if (!is.factor(train_data[[gc]])) {
    train_data[[gc]] <- as.factor(train_data[[gc]])
  }
  
  # Ensure factor levels are valid R variable names (e.g., "X0", "X1" instead of "0", "1")
  # This prevents caret::twoClassSummary from crashing
  levels(train_data[[gc]]) <- make.names(levels(train_data[[gc]]))
  
  n_features <- ncol(train_data) - 1
  if (verbose) {
    cat(">>> Tuning model method:", method, "\n")
    cat(">>> Parameter bounds:", paste(names(bounds), collapse = ", "), "\n")
  }
  
  # --- 4. Parameter Protection (Random Forest mtry) ---
  if (method == "rf" && "mtry" %in% names(bounds)) {
    bounds$mtry[1] <- 1
    bounds$mtry[2] <- min(bounds$mtry[2], n_features)
    if (verbose) cat(">>> mtry range adjusted to: [", bounds$mtry[1], ",", bounds$mtry[2], "]\n")
  }
  
  # --- 5. Class Weights Calculation ---
  wts <- NULL
  if (class_weights) {
    tab <- table(train_data[[gc]])
    wts <- as.numeric(1 / tab[as.character(train_data[[gc]])])
  }
  
  # --- 6. Integer Parameter Identification ---
  # Fetch model metadata from caret to identify integer-class parameters
  model_info <- tryCatch(
    caret::getModelInfo(method, regex = FALSE)[[1]],
    error = function(e) stop("Model method '", method, "' not found. Check if the required package is installed.")
  )
  param_df <- model_info$parameters
  int_params <- param_df$parameter[param_df$class == "integer"]
  
  # Manually include common integer parameters often mislabeled in metadata
  extra_int  <- c("n.trees", "nrounds", "n.minobsinnode", "min_child_weight", "mtry", "max_depth")
  int_params <- unique(c(int_params, intersect(extra_int, names(bounds))))
  
  # --- 7. Objective Function for Bayesian Optimization ---
  obj_func <- function(...) {
    params <- list(...)
    
    # Round parameters that must be integers
    for (p in intersect(names(params), int_params)) {
      params[[p]] <- round(params[[p]])
    }
    
    tune_grid <- do.call(data.frame, params)
    
    # Configure caret trainControl
    ctrl <- caret::trainControl(
      method          = "cv",
      number          = cv_folds,
      classProbs      = TRUE,
      summaryFunction = summaryFun,
      sampling        = sampling,
      verboseIter     = FALSE
    )
    
    # Run training with error handling
    res <- tryCatch({
      mod <- caret::train(
        as.formula(paste(gc, "~ .")),
        data      = train_data,
        method    = method,
        tuneGrid  = tune_grid,
        trControl = ctrl,
        weights   = wts,
        metric    = metric
      )
      
      score <- max(mod$results[[metric]], na.rm = TRUE)
      list(Score = score, Pred = 0)
      
    }, error = function(e) {
      if (verbose) message("\n[Training Error]: ", e$message)
      return(list(Score = -1e6, Pred = 0)) # Return penalty score on failure
    })
    
    return(res)
  }
  
  # --- 8. Bayesian Optimization Execution ---
  if (verbose) cat("\n>>> Starting Bayesian Optimization process...\n")
  opt_res <- rBayesianOptimization::BayesianOptimization(
    FUN         = obj_func,
    bounds      = bounds,
    init_points = init_points,
    n_iter      = n_iter,
    acq         = "ucb", 
    kappa       = 2.576, # Typical value for exploration/exploitation balance
    verbose     = verbose
  )
  
  # --- 9. Final Model Retraining with Best Parameters ---
  if (verbose) cat("\n>>> Retraining final model with optimal parameters...\n")
  best_params <- as.list(opt_res$Best_Par)
  for (p in intersect(names(best_params), int_params)) {
    best_params[[p]] <- round(best_params[[p]])
  }
  final_grid <- do.call(data.frame, best_params)
  
  ctrl_final <- caret::trainControl(
    method          = "cv",
    number          = cv_folds,
    classProbs      = TRUE,
    summaryFunction = summaryFun,
    sampling        = sampling
  )
  
  final_model <- caret::train(
    as.formula(paste(gc, "~ .")),
    data      = train_data,
    method    = method,
    tuneGrid  = final_grid,
    trControl = ctrl_final,
    weights   = wts,
    metric    = metric
  )
  
  # --- 10. Store and Return Results ---
  model_obj@best.model.result$fine_tuned_model <- final_model
  model_obj@best.model.result$tuning_result      <- opt_res
  
  if (verbose) {
    cat(">>> Fine-tuning completed. Best Parameters Found:\n")
    print(final_grid)
  }
  
  return(model_obj)
}

# =============================================================================
# Helper: Plot tuning history
# =============================================================================

#' Plot Bayesian Optimization History
#'
#' Shows the best score found so far across iterations.
#'
#' @param model_obj A \code{Train_Model} object that has been fine-tuned.
#' @param save_plot Logical.
#' @param save_dir  Output directory.
#' @return A ggplot.
#' @export
PlotTuningHistory <- function(model_obj, save_plot = FALSE, save_dir = NULL) {
  opt_res <- model_obj@best.model.result$tuning_result
  if (is.null(opt_res)) stop("No tuning result found. Run FineTuneModel first.")
  
  # $History is a data.table with columns "Round", "mtry", "Value"
  hist_dt <- opt_res$History
  if (!is.data.frame(hist_dt) || nrow(hist_dt) == 0) stop("No tuning history available.")
  
  scores <- hist_dt$Value
  valid_scores <- scores[is.finite(scores)]
  if (length(valid_scores) == 0) stop("All tuning scores are -Inf/NA.")
  
  best_curve <- cummax(valid_scores)   # cumulative maximum
  
  hist_df <- data.frame(
    Iteration  = seq_along(best_curve),
    Best_Score = best_curve
  )
  
  p <- ggplot2::ggplot(hist_df, ggplot2::aes(x = Iteration, y = Best_Score)) +
    ggplot2::geom_line(color = "#b2e2e2", linewidth = 1) +
    ggplot2::geom_point(color = "#006d2c", size = 2) +
    ggplot2::labs(title = "Bayesian Optimization History",
                  x = "Iteration", y = "Best ROC") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuning_history.pdf"), plot = p,
                    width = 6, height = 4, dpi = 300)
    cat("Plot saved to:", file.path(save_dir, "tuning_history.pdf"), "\n")
  }
  return(p)
}


# Internal theme helper (shared with viz_functions.R)
.pub_theme <- function(base_size = 13) {
  ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.title    = ggplot2::element_text(face = "bold"),
      legend.title  = ggplot2::element_text(face = "bold"),
      strip.text    = ggplot2::element_text(face = "bold")
    )
}

.save_plot <- function(p, dir, filename, width, height, format = "pdf") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, paste0(tools::file_path_sans_ext(filename), ".", format))
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
  cat("Plot saved to:", path, "\n")
}

# =============================================================================
# Tuned Model Visualisation Functions
# =============================================================================
#' Plot ROC Curve for Tuned Model
#'
#' @description Generates an ROC curve for the fine-tuned model and optionally 
#' compares it with the original best model. Automatically handles factor level 
#' alignment for numeric-origin groups.
#'
#' @param tuned_model A caret \code{train} object (the fine-tuned model).
#' @param original_best A caret \code{train} object for comparison. Default is NULL.
#' @param test_data A data frame containing the test set.
#' @param group_col Character. Name of the grouping/target column.
#' @param palette_name Character. Wesanderson palette name. Default "Darjeeling1".
#' @param save_plot Logical. Whether to save the plot as a PDF.
#' @param save_dir Character. Output directory for the saved plot.
#' @param width,height Numeric. Plot dimensions in inches.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @importFrom pROC roc auc
#' @importFrom ggplot2 ggplot aes geom_line geom_abline scale_color_manual labs coord_equal theme element_text ggsave
#' @importFrom wesanderson wes_palette
#' @importFrom ggprism theme_prism
PlotTunedROC <- function(tuned_model,
                         original_best = NULL,
                         test_data,
                         group_col,
                         palette_name = "Darjeeling1",
                         save_plot = FALSE,
                         save_dir = NULL,
                         width = 7,
                         height = 6) {
  
  if (!inherits(tuned_model, "train")) stop("tuned_model must be a caret train object.")
  
  # Align test labels with the model's make.names() transformation
  test_labels <- as.factor(test_data[[group_col]])
  levels(test_labels) <- make.names(levels(test_labels))
  
  # Helper to compute ROC data
  .compute_roc <- function(model, model_label) {
    # Get probabilities for the second class
    probs <- stats::predict(model, newdata = test_data, type = "prob")[, 2]
    
    roc_obj <- pROC::roc(test_labels, probs, 
                         levels = levels(test_labels), 
                         direction = "auto", quiet = TRUE)
    
    auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
    data.frame(
      Sensitivity = roc_obj$sensitivities,
      Specificity = roc_obj$specificities,
      Model = paste0(model_label, " (AUC = ", auc_val, ")")
    )
  }
  
  df_tuned <- .compute_roc(tuned_model, "Tuned")
  
  if (!is.null(original_best)) {
    df_best <- .compute_roc(original_best, "Original Best")
    roc_df <- rbind(df_tuned, df_best)
  } else {
    roc_df <- df_tuned
  }
  
  n_models <- length(unique(roc_df$Model))
  cols <- wesanderson::wes_palette(palette_name, max(2, n_models), type = "discrete")
  
  p <- ggplot2::ggplot(roc_df, 
                       ggplot2::aes(x = 1 - .data$Specificity, y = .data$Sensitivity, color = .data$Model)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(title = "ROC Curve: Tuned vs Original",
                  x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::coord_equal() +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = c(0.75, 0.25))
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuned_ROC.pdf"), plot = p,
                    width = width, height = height, dpi = 300)
  }
  
  return(p)
}

#' Plot Confusion Matrix for Tuned Model
#'
#' @description Generates a heatmap-style confusion matrix. Percentages are 
#' calculated per actual class (column-wise).
#'
#' @param tuned_model A caret \code{train} object.
#' @param test_data A data frame containing the test set.
#' @param group_col Character. Name of the grouping column.
#' @param palette Character vector of two colors for the gradient.
#' @param save_plot Logical.
#' @param save_dir Character.
#' @param width,height Numeric.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient labs
PlotTunedConfusion <- function(tuned_model,
                               test_data,
                               group_col,
                               palette = c("#b2e2e2", "#006d2c"),
                               save_plot = FALSE,
                               save_dir = NULL,
                               width = 5,
                               height = 4.5) {
  
  if (!inherits(tuned_model, "train")) stop("tuned_model must be a caret train object.")
  
  # Align labels
  truth <- as.factor(test_data[[group_col]])
  levels(truth) <- make.names(levels(truth))
  
  pred <- stats::predict(tuned_model, newdata = test_data, type = "raw")
  
  cf <- table(Predicted = factor(pred, levels = levels(truth)),
              Actual = truth)
  
  cf_df <- as.data.frame(cf)
  # Calculate percentages per Actual class
  cf_df <- cf_df %>%
    dplyr::group_by(.data$Actual) %>%
    dplyr::mutate(Pct = round(.data$Freq / sum(.data$Freq) * 100, 1)) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(cf_df, ggplot2::aes(x = .data$Actual, y = .data$Predicted, fill = .data$Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(.data$Freq, "\n(", .data$Pct, "%)")),
                       size = 4.5, fontface = "bold") +
    ggplot2::scale_fill_gradient(low = palette[1], high = palette[2]) +
    ggplot2::labs(title = "Confusion Matrix: Tuned Model",
                  x = "Actual Class", y = "Predicted Class", fill = "Count") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuned_confusion.pdf"), plot = p,
                    width = width, height = height, dpi = 300)
  }
  
  return(p)
}

#' Plot Calibration Curve for Tuned Model
#'
#' @description Comprehensive calibration analysis for the fine-tuned model. 
#' Includes binned observations, a LOESS smoother, and key statistics (Brier score, 
#' Slope, Intercept, and Eavg).
#'
#' @param tuned_model A caret \code{train} object (the fine-tuned model).
#' @param test_data A data frame containing the test set.
#' @param group_col Character. Name of the grouping column.
#' @param n_bins Integer. Number of probability bins. Default 10.
#' @param palette Character. Color for the points and smooth line.
#' @param base_size Numeric. Base font size for the plot.
#' @param se Logical. Show standard error on the smoother.
#' @param show_stats Logical. Whether to display calibration metrics on the plot.
#' @param save_plot Logical. Save the plot to a file.
#' @param save_dir Character. Directory to save the plot.
#' @param width,height Numeric. Plot dimensions in inches.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @importFrom stats predict glm binomial coef
#' @importFrom dplyr group_by summarise n ungroup
#' @importFrom ggplot2 ggplot aes geom_abline geom_point geom_smooth scale_size_continuous labs theme annotate ggsave xlim ylim
#' @importFrom ggprism theme_prism
PlotTunedCalibration <- function(tuned_model,
                                 test_data,
                                 group_col,
                                 n_bins     = 10,
                                 palette    = "#006d2c",
                                 base_size  = 13,
                                 se         = FALSE,
                                 show_stats = TRUE,
                                 save_plot  = FALSE,
                                 save_dir   = NULL,
                                 width      = 6,
                                 height     = 5.5) {
  
  if (!inherits(tuned_model, "train")) stop("tuned_model must be a caret train object.")
  
  # ---------- 1. Factor Level Alignment ----------
  truth <- as.factor(test_data[[group_col]])
  # Ensure labels match the 'make.names' transformation used during tuning
  levels(truth) <- make.names(levels(truth))
  levels_true <- levels(truth)
  
  # We assume the second level is the positive class
  pos_level <- levels_true[2]
  truth_numeric <- as.integer(truth) - 1L 
  
  # ---------- 2. Probability Extraction ----------
  prob_mat <- stats::predict(tuned_model, newdata = test_data, type = "prob")
  # Use the column matching the positive level name
  probs <- if (pos_level %in% colnames(prob_mat)) {
    prob_mat[, pos_level]
  } else {
    prob_mat[, 2]
  }
  
  # ---------- 3. Metrics Calculation ----------
  # Brier Score
  brier <- mean((truth_numeric - probs)^2)
  
  # Calibration Slope & Intercept (Logit-based)
  cal_df <- data.frame(truth = truth_numeric, prob = probs)
  cal_df$prob_clip <- pmax(pmin(cal_df$prob, 1 - 1e-6), 1e-6)
  
  cal_glm <- suppressWarnings(
    stats::glm(truth ~ log(prob_clip/(1 - prob_clip)), family = stats::binomial(), data = cal_df)
  )
  intercept <- stats::coef(cal_glm)[1]
  slope     <- stats::coef(cal_glm)[2]
  
  # Binned Data for Eavg
  cal_df$bin <- cut(probs, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  cal_sum <- cal_df %>%
    dplyr::group_by(.data$bin) %>%
    dplyr::summarise(
      mean_pred = mean(.data$prob),
      obs_rate  = mean(.data$truth),
      n         = dplyr::n(),
      .groups   = "drop"
    )
  e_avg <- mean(abs(cal_sum$obs_rate - cal_sum$mean_pred), na.rm = TRUE)
  
  # ---------- 4. Visualization ----------
  p <- ggplot2::ggplot(cal_sum, ggplot2::aes(x = .data$mean_pred, y = .data$obs_rate)) +
    # Ideal calibration line
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    # Binned points
    ggplot2::geom_point(ggplot2::aes(size = .data$n), color = palette, alpha = 0.7) +
    # LOESS smoother
    ggplot2::geom_smooth(method = "loess", se = se, color = palette, linewidth = 1, 
                         fill = paste0(palette, "30"), formula = y ~ x) +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "Samples (n)") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(title = "Calibration Curve: Tuned Model",
                  subtitle = paste0("Model: ", tuned_model$method, " | Bins: ", n_bins),
                  x = "Mean Predicted Probability", y = "Observed Proportion") +
    ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  # Add Stats Annotation
  if (show_stats) {
    stats_text <- paste0(
      "Brier: ", round(brier, 4), "\n",
      "Intercept: ", round(intercept, 3), "\n",
      "Slope: ", round(slope, 3), "\n",
      "Eavg: ", round(e_avg, 4)
    )
    p <- p + ggplot2::annotate("text", x = 0.05, y = 0.95, label = stats_text, 
                               hjust = 0, vjust = 1, size = base_size * 0.28, 
                               family = "mono", fontface = "bold")
  }
  
  # ---------- 5. Save and Return ----------
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- getwd()
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "tuned_calibration.pdf"), plot = p,
                    width = width, height = height, dpi = 300)
    cat(">>> Calibration plot saved to:", file.path(save_dir, "tuned_calibration.pdf"), "\n")
  }
  
  return(p)
}

