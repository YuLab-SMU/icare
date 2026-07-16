#' Benchmark Binary Classification Models with Preprocessing Pipelines
#'
#' Nested cross-validation benchmarking framework for preprocessing
#' methodology studies in binary classification.
#'
#' @param data A data frame, or an S4 object with slots \code{clean.df} and
#'   \code{group_col} (e.g. \code{Train_Model} or \code{Stat}).
#' @param group_col Character string specifying binary outcome column.
#'   Ignored if \code{data} is an S4 object that already stores the group column.
#' @param algorithms Character vector of caret model names.
#' @param impute_methods Character vector of imputation methods.
#'   Supported: \code{"none"}, \code{"median"}, \code{"knn"}, \code{"bag"}.
#' @param norm_methods List of normalization methods.
#' @param outer_folds Number of outer CV folds.
#' @param inner_cv_number Number of folds for inner CV.
#' @param inner_repeats Number of repeats for inner CV.
#' @param tuneLength Number of tuning parameter combinations.
#' @param seed Random seed.
#' @param verbose Logical; if \code{TRUE} prints a concise progress message
#'   for each algorithm.
#'
#' @return A data frame of class \code{"PreprocessingBenchmark"}.
#'   Columns include \code{Algorithm}, \code{Imputation}, \code{Normalization},
#'   \code{Mean_AUC}, \code{SD_AUC}, \code{CI_lower}, \code{CI_upper},
#'   \code{n_success}, and \code{failure_rate}.
#'   An attribute \code{error_log} is attached with any failures encountered.
#'
#' @details
#' The framework avoids near-zero variance filtering and correlation filtering
#' to preserve orthogonality in preprocessing comparisons. All preprocessing
#' operations are learned on the training folds and applied to the corresponding
#' test folds.
#'
#' @references
#' Kuhn M. (2008). Building Predictive Models in R Using the
#' caret Package. \emph{Journal of Statistical Software}, 28(5), 1-26.
#' @export
#' @examples
#' \dontrun{
#' library(mlbench)
#' data(Sonar)
#' result <- PreprocessingBenchmark(
#'   data = Sonar,
#'   group_col = "Class",
#'   algorithms = c("rf", "gbm", "xgbTree", "glmnet"),
#'   impute_methods = c("none", "median"),
#'   norm_methods = list("none", c("center","scale"))
#' )
#' print(result)
#' attr(result, "error_log")
#' }
PreprocessingBenchmark <- function(
    data,
    group_col = NULL,
    algorithms = c("rf", "ranger", "gbm",  "glmnet", "svmRadial",
                   "knn", "nnet"),
    impute_methods = c("none", "median", "knn", "bag"),
    norm_methods = list("none", c("center", "scale"),
                        c("center", "scale", "YeoJohnson")),
    outer_folds = 5,
    inner_cv_number = 5,
    inner_repeats = 2,
    tuneLength = 3,
    seed = 123,
    verbose = TRUE
) {
  
  # =====================================================
  # Handle S4 train/stat objects
  # =====================================================
  if (methods::is(data, "Train_Model")) {
    df <- data@clean.df
    if (is.null(group_col)) group_col <- data@group_col
  } else if (methods::is(data, "Stat")) {
    df <- data@clean.data
    if (is.null(group_col)) group_col <- data@group_col
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("data must be a data.frame, Train_Model, or Stat object.")
  }
  
  if (is.null(group_col)) stop("group_col must be provided.")
  if (!group_col %in% names(df)) stop("group_col not found in data.")
  
  # =====================================================
  # Required packages
  # =====================================================
  required_pkgs <- c("caret", "pROC", "dplyr")
  missing <- required_pkgs[!sapply(required_pkgs, requireNamespace,
                                   quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "))
  }
  
  # =====================================================
  # Convert predictors to numeric
  # =====================================================
  feature_cols <- setdiff(names(df), group_col)
  for (col in feature_cols) {
    if (is.factor(df[[col]])) {
      df[[col]] <- as.numeric(df[[col]])
    } else if (is.character(df[[col]])) {
      df[[col]] <- as.numeric(as.factor(df[[col]]))
    }
  }
  
  # =====================================================
  # Safe binary outcome
  # =====================================================
  safe_levels <- function(y) {
    y <- as.factor(y)
    if (length(levels(y)) != 2) stop("Only binary classification supported")
    levels(y) <- make.names(levels(y))
    y
  }
  df[[group_col]] <- safe_levels(df[[group_col]])
  positive_class <- levels(df[[group_col]])[2]
  
  # =====================================================
  # Inner trainControl
  # =====================================================
  inner_control <- caret::trainControl(
    method = "repeatedcv",
    number = inner_cv_number,
    repeats = inner_repeats,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  # =====================================================
  # Outer folds
  # =====================================================
  set.seed(seed)
  folds <- caret::createFolds(df[[group_col]], k = outer_folds, list = TRUE)
  
  results_list <- list()
  error_log <- list()
  
  # =====================================================
  # Main loop over algorithms, imputations, normalizations
  # =====================================================
  for (algo in algorithms) {
    if (verbose) cat("Processing", algo, "...\n")
    
    for (imp in impute_methods) {
      for (norm in norm_methods) {
        norm_name <- if (identical(norm, "none")) {
          "none"
        } else {
          paste(norm, collapse = "_")
        }
        
        aucs <- c()
        success_count <- 0
        
        for (k in seq_along(folds)) {
          test_idx <- folds[[k]]
          train_df <- df[-test_idx, , drop = FALSE]
          test_df  <- df[ test_idx, , drop = FALSE]
          
          train_x <- train_df[, feature_cols, drop = FALSE]
          test_x  <- test_df[, feature_cols, drop = FALSE]
          train_y <- train_df[[group_col]]
          test_y  <- test_df[[group_col]]
          
          # ---- Imputation ----
          if (imp == "none") {
            keep_idx <- complete.cases(train_x)
            train_x <- train_x[keep_idx, , drop = FALSE]
            train_y <- train_y[keep_idx]
            keep_idx <- complete.cases(test_x)
            test_x <- test_x[keep_idx, , drop = FALSE]
            test_y <- test_y[keep_idx]
          } else if (imp %in% c("median", "knn", "bag")) {
            method_map <- c(median = "medianImpute",
                            knn    = "knnImpute",
                            bag    = "bagImpute")
            pp_imp <- tryCatch(
              caret::preProcess(train_x, method = method_map[imp]),
              error = function(e) NULL
            )
            if (is.null(pp_imp)) next
            train_x <- predict(pp_imp, train_x)
            test_x  <- predict(pp_imp, test_x)
          } else {
            next
          }
          
          # ---- Remove constant / invalid columns ----
          keep_cols <- sapply(train_x, function(x) {
            is.numeric(x) && !all(is.na(x)) && length(unique(x)) > 1
          })
          train_x <- train_x[, keep_cols, drop = FALSE]
          test_x  <- test_x[, keep_cols, drop = FALSE]
          if (ncol(train_x) == 0) next
          
          # ---- Normalization ----
          if (!identical(norm, "none")) {
            pp_norm <- tryCatch(
              caret::preProcess(train_x, method = norm),
              error = function(e) NULL
            )
            if (is.null(pp_norm)) next
            train_x <- predict(pp_norm, train_x)
            test_x  <- predict(pp_norm, test_x)
          }
          
          # ---- Final data frames ----
          train_final <- data.frame(train_x, Class = train_y)
          test_final  <- data.frame(test_x,  Class = test_y)
          
          # ---- Model training (suppress caret messages) ----
          fit <- tryCatch({
            suppressWarnings(suppressMessages(
              caret::train(
                Class ~ .,
                data = train_final,
                method = algo,
                metric = "ROC",
                maximize = TRUE,
                trControl = inner_control,
                tuneLength = tuneLength
              )
            ))
          }, error = function(e) {
            error_log[[length(error_log) + 1]] <<-
              data.frame(
                Algorithm = algo,
                Imputation = imp,
                Normalization = norm_name,
                Fold = k,
                Error = as.character(e$message),
                stringsAsFactors = FALSE
              )
            NULL
          })
          if (is.null(fit)) next
          
          # ---- Predict probabilities ----
          probs <- tryCatch({
            pred <- predict(fit, newdata = test_final, type = "prob")
            pred[, positive_class]
          }, error = function(e) {
            error_log[[length(error_log) + 1]] <<-
              data.frame(
                Algorithm = algo,
                Imputation = imp,
                Normalization = norm_name,
                Fold = k,
                Error = paste("Prediction:", e$message),
                stringsAsFactors = FALSE
              )
            NULL
          })
          if (is.null(probs)) next
          
          # ---- AUC ----
          auc_val <- tryCatch({
            roc_obj <- pROC::roc(
              response  = test_y,
              predictor = probs,
              levels    = levels(test_y),
              direction = "<",
              quiet     = TRUE
            )
            as.numeric(roc_obj$auc)
          }, error = function(e) {
            error_log[[length(error_log) + 1]] <<-
              data.frame(
                Algorithm = algo,
                Imputation = imp,
                Normalization = norm_name,
                Fold = k,
                Error = paste("AUC:", e$message),
                stringsAsFactors = FALSE
              )
            NULL
          })
          if (is.null(auc_val)) next
          
          aucs <- c(aucs, auc_val)
          success_count <- success_count + 1
        } # end folds
        
        # ---- Summarize ----
        if (length(aucs) > 0) {
          mean_auc <- mean(aucs)
          sd_auc   <- sd(aucs)
          n_success <- length(aucs)
          se <- sd_auc / sqrt(n_success)
          ci_lower <- mean_auc - qt(0.975, df = n_success - 1) * se
          ci_upper <- mean_auc + qt(0.975, df = n_success - 1) * se
          
          results_list[[length(results_list) + 1]] <- data.frame(
            Algorithm    = algo,
            Imputation   = imp,
            Normalization = norm_name,
            Mean_AUC     = mean_auc,
            SD_AUC       = sd_auc,
            CI_lower     = ci_lower,
            CI_upper     = ci_upper,
            n_success    = n_success,
            failure_rate = (outer_folds - n_success) / outer_folds,
            stringsAsFactors = FALSE
          )
        }
      } # norm
    } # imp
  } # algo
  
  # =====================================================
  # Final output
  # =====================================================
  results_df <- dplyr::bind_rows(results_list)
  attr(results_df, "error_log") <- dplyr::bind_rows(error_log)
  class(results_df) <- c("PreprocessingBenchmark", class(results_df))
  
  if (verbose) cat("Benchmark completed.\n")
  return(results_df)
}
#' Forest Plot for Preprocessing Benchmark Results
#'
#' Creates a faceted pointrange plot showing mean AUC and 95% confidence
#' intervals for each algorithm, grouped by imputation and normalization methods.
#'
#' @param benchmark_result Output from \code{PreprocessingBenchmark()}.
#' @param palette_name Wes Anderson palette name (default \code{"Zissou1"}).
#' @param remove_na Logical; remove rows with NA Mean_AUC (recommended).
#' @param save_plot Logical; save plot to PDF.
#' @param save_dir Directory to save plot (required if \code{save_plot = TRUE}).
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' result <- PreprocessingBenchmark(...)
#' PlotBenchmarkForest(result)
#'
#' # Save to file
#' PlotBenchmarkForest(result, save_plot = TRUE, save_dir = "~/plots")
#' }
#'
#' @export
PlotBenchmarkForest <- function(benchmark_result,
                                palette_name = "Zissou1",
                                remove_na = TRUE,
                                save_plot = FALSE,
                                save_dir = NULL) {
  
  required_pkgs <- c("ggplot2", "dplyr", "wesanderson", "ggprism")
  missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) stop("Missing packages for plotting: ", paste(missing, collapse = ", "))
  if (!inherits(benchmark_result, "PreprocessingBenchmark") && !is.data.frame(benchmark_result)) {
    stop("Input must be the output from PreprocessingBenchmark()")
  }
  
  df <- as.data.frame(benchmark_result)
  
  if (remove_na) {
    df <- df[!is.na(df$Mean_AUC), ]
    if (nrow(df) == 0) stop("No valid results after removing NAs.")
  } else {
    df <- df[!is.na(df$Algorithm) & !is.na(df$Imputation) & !is.na(df$Normalization), ]
  }
  
  if (nrow(df) == 0) stop("No data available for plotting.")
  
  # Confidence interval label
  if (!"CI_lower" %in% colnames(df) || !"CI_upper" %in% colnames(df)) {
    df$CI_lower <- df$Mean_AUC
    df$CI_upper <- df$Mean_AUC
    ci_label <- "No CI"
  } else {
    ci_label <- "Percentile 95% CI"
  }
  
  # Facet medians
  facet_stats <- df %>%
    dplyr::group_by(Imputation, Normalization) %>%
    dplyr::summarise(Facet_Median = median(Mean_AUC, na.rm = TRUE), .groups = "drop")
  
  # Order imputation levels (aligned with possible methods)
  df$Imputation <- factor(df$Imputation, levels = c("none", "median", "knn", "bag"))
  
  n_algos <- length(unique(df$Algorithm))
  my_cols <- wesanderson::wes_palette(palette_name, n = n_algos, type = "continuous")
  
  global_median <- median(df$Mean_AUC, na.rm = TRUE)
  
  # Order algorithms by median AUC
  algo_order <- df %>%
    group_by(Algorithm) %>%
    summarise(med = median(Mean_AUC, na.rm = TRUE)) %>%
    arrange(desc(med)) %>%
    pull(Algorithm)
  df$Algorithm <- factor(df$Algorithm, levels = algo_order)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Mean_AUC, y = Algorithm, color = Algorithm)) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = CI_lower, xmax = CI_upper),
                             linewidth = 0.85, alpha = 0.9) +
    ggplot2::facet_grid(Imputation ~ Normalization, scales = "free_y") +
    ggplot2::geom_vline(xintercept = global_median, linetype = "dashed",
                        color = "gray40", alpha = 0.7) +
    ggplot2::scale_color_manual(values = my_cols) +
    ggplot2::geom_text(data = facet_stats,
                       ggplot2::aes(x = -Inf, y = Inf,
                                    label = paste0("Median: ", round(Facet_Median, 3))),
                       vjust = 1.6, hjust = -0.05, size = 3.3,
                       color = "black", fontface = "bold", inherit.aes = FALSE) +
    ggplot2::annotate("text", x = global_median, y = -Inf,
                      label = paste("Global Median:", round(global_median, 3)),
                      vjust = -0.8, hjust = 0.5, size = 3.8, color = "gray30") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::labs(x = paste("Mean AUC +/-", ci_label),
                  y = NULL,
                  title = "Preprocessing Benchmark - AUC Performance") +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  print(p)
  
  if (save_plot) {
    if (is.null(save_dir)) stop("save_dir must be specified when save_plot = TRUE")
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(save_dir, "Preprocessing_Benchmark_Forest.pdf")
    ggplot2::ggsave(filename, plot = p, width = 14, height = 10, dpi = 300)
    cat("Plot saved:", filename, "\n")
  }
  
  invisible(p)
}
