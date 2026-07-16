
# ============================================================================
# Internal Utility Functions
# These are helper functions used across modules
# All start with dot (.) to indicate they are internal
# ============================================================================
#' @keywords internal
.check_class <- function(object, allowed) {
  if (!inherits(object, allowed))
    stop("Object must be one of: ", paste(allowed, collapse = ", "),
         ". Got: ", class(object)[1])
}

#' @keywords internal
.subset_df <- function(df, rows = NULL, cols = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!is.null(rows)) df <- df[rownames(df) %in% rows, , drop = FALSE]
  if (!is.null(cols)) df <- df[, colnames(df) %in% cols,  drop = FALSE]
  df
}

#' @keywords internal
.safe_rbind <- function(df1, df2) {
  if (is.null(df1) || nrow(df1) == 0) return(df2)
  if (is.null(df2) || nrow(df2) == 0) return(df1)
  shared <- intersect(colnames(df1), colnames(df2))
  if (length(shared) == 0) {
    warning("No shared columns; returning first data frame only.")
    return(df1)
  }
  if (length(shared) < max(ncol(df1), ncol(df2)))
    message(length(shared), " shared columns kept.")
  out <- rbind(df1[, shared, drop = FALSE], df2[, shared, drop = FALSE])
  rownames(out) <- make.unique(rownames(out))
  out
}

#' @keywords internal
.clean_symbol_values <- function(data) {
  for(col in colnames(data)) {
    if(is.character(data[[col]])) {
      if(any(grepl("^[<>]", data[[col]]), na.rm = TRUE)) {
        data[[col]] <- gsub("[<>]", "", data[[col]])
        data[[col]] <- as.numeric(data[[col]])
        warning(paste("Removed >/< symbols from column", col))
      }
    }
  }
  return(data)
}

#' @keywords internal
.prepare_data <- function(data, data_name) {
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame())
  }
  
  if (anyDuplicated(rownames(data))) {
    warning(paste("Duplicate row names in", data_name, "made unique."))
    rownames(data) <- make.unique(rownames(data))
  }
  
  if (anyDuplicated(colnames(data))) {
    warning(paste("Duplicate column names in", data_name, "made unique."))
    colnames(data) <- make.unique(colnames(data))
  }
  
  if (is.null(colnames(data))){
    stop(paste(data_name, "is missing column names."))
  }
  
  return(data)
}

#' @keywords internal
.ensure_numeric_data <- function(data, data_name, convert_factors = TRUE) {
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame())
  }
  
  row_names <- rownames(data)
  
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x)) {
      if (convert_factors) {
        levels <- levels(x)
        if (length(levels) == 2) {
          return(as.numeric(x) - 1)
        } else {
          warning("Factor with >2 levels in ", data_name)
          return(as.numeric(x))
        }
      } else {
        return(x)
      }
    } else if (is.character(x)) {
      if (convert_factors) {
        x <- as.factor(x)
        levels <- levels(x)
        if (length(levels) == 2) {
          return(as.numeric(x) - 1)
        } else {
          return(as.numeric(x))
        }
      } else {
        return(x)
      }
    } else if (is.numeric(x)) {
      return(x)
    } else {
      warning("Column in ", data_name, " cannot be converted to numeric.")
      return(x)
    }
  }))
  
  rownames(data) <- row_names
  
  # Remove all-NA columns
  na_columns <- colnames(data)[apply(data, 2, function(x) all(is.na(x)))]
  if (length(na_columns) > 0) {
    warning("Removing all-NA columns from ", data_name, ": ", 
            paste(na_columns, collapse = ", "))
    data <- data[, !colnames(data) %in% na_columns]
  }
  
  return(data)
}

#' Remove constant (zero-variance) columns from a data frame or matrix
#'
#' Identifies and removes columns that have only one distinct value (including NA).
#'
#' @param data A data frame or matrix.
#' @param na.rm Logical. If TRUE, NA values are ignored when checking constancy
#'   (i.e., a column with only one non-NA value and the rest NA is considered constant).
#'   If FALSE (default), NA is treated as a distinct value.
#' @param verbose Logical. Print messages about removed columns (default TRUE).
#' @param allow_empty Logical. If TRUE, allow returning an empty data frame when
#'   all columns are constant (default FALSE, which raises an error).
#'
#' @return A data frame or matrix (same class as input) with constant columns removed.
#' @export
#'
#' @examples
#' df <- data.frame(a = 1:5, b = rep(2, 5), c = c(1, NA, 1, 1, 1))
#' remove_constant_columns(df)                 # removes column 'b'
#' remove_constant_columns(df, na.rm = TRUE)   # also removes column 'c'
#' remove_constant_columns(df, verbose = FALSE)
remove_constant_columns <- function(data,
                                    na.rm   = FALSE,
                                    verbose = TRUE,
                                    allow_empty = FALSE) {
  
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix.")
  }
  
  # Preserve original class for return
  original_class <- class(data)
  
  # Convert to data frame for uniform handling
  if (is.matrix(data)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE)
    was_matrix <- TRUE
  } else {
    df <- data
    was_matrix <- FALSE
  }
  
  if (ncol(df) == 0) {
    if (verbose) message("No columns to check.")
    return(data)
  }
  
  # Check each column for constancy
  is_constant <- sapply(df, function(col) {
    if (na.rm) {
      non_na <- col[!is.na(col)]
      if (length(non_na) == 0) {
        # All NA -> considered constant
        return(TRUE)
      }
      length(unique(non_na)) == 1
    } else {
      length(unique(col)) == 1
    }
  })
  
  # If all columns are constant and empty not allowed, error
  if (all(is_constant) && !allow_empty) {
    stop("All columns are constant. Set allow_empty = TRUE to return an empty object.")
  }
  
  # Remove constant columns
  if (any(is_constant)) {
    removed <- names(is_constant)[is_constant]
    if (verbose) {
      cat("Removing constant columns:", paste(removed, collapse = ", "), "\n")
    }
    df <- df[, !is_constant, drop = FALSE]
  } else {
    if (verbose) cat("No constant columns found.\n")
  }
  
  # Restore original format (matrix or data frame)
  if (was_matrix && ncol(df) > 0) {
    result <- as.matrix(df)
    rownames(result) <- rownames(data)
    colnames(result) <- colnames(df)
  } else {
    result <- df
    if (is.data.frame(result) && !is.data.frame(data)) {
      result <- as.data.frame(result, stringsAsFactors = FALSE)
      rownames(result) <- rownames(data)
    }
  }
  
  return(result)
}


#' Example Stat Object for Testing
#'
#' A \code{Stat} S4 object created from the Bacteremia public dataset.
#' Used for demonstration and testing of statistical analysis functions.
#'
#' @details
#' This object was created from the first 1000 rows of the 
#' \href{https://zenodo.org/records/7554815}{Bacteremia public dataset}.
#' The \code{BloodCulture} column is used as the grouping variable.
#'
#' The object contains:
#' \itemize{
#'   \item \code{@raw.data}: Raw data with BloodCulture as grouping column
#'   \item \code{@clean.data}: Processed numeric data
#'   \item \code{@info.data}: Metadata including BloodCulture group
#'   \item \code{@group_col}: "BloodCulture"
#' }
#'
#' @format A \code{Stat} S4 object with slots:
#' \describe{
#'   \item{raw.data}{Original data frame}
#'   \item{clean.data}{Cleaned numeric data matrix}
#'   \item{info.data}{Metadata data frame}
#'   \item{group_col}{Character, grouping column name ("BloodCulture")}
#'   \item{...}{Additional slots for analysis results}
#' }
#'
#' @source \url{https://zenodo.org/records/7554815}
#' @keywords datasets
"stat_obj_test"


#' Example PrognosiX Object for Testing
#'
#' A \code{PrognosiX} S4 object converted from \code{stat_obj_test}.
#' Used for demonstration and testing of survival/prognostic analysis functions.
#'
#' @details
#' This object was created by converting \code{stat_obj_test} using
#' \code{CreatePrognosiXObject()}, with \code{"time"} and \code{"status"}
#' columns extracted from the clinical metadata.
#'
#' The object contains:
#' \itemize{
#'   \item \code{@survival.data}: Combined data with time and status
#'   \item \code{@clean.data}: Numeric feature matrix
#'   \item \code{@info.data}: Metadata with time and status columns
#'   \item \code{@time_col}: "time"
#'   \item \code{@status_col}: "status"
#' }
#'
#' @format A \code{PrognosiX} S4 object with slots:
#' \describe{
#'   \item{clean.data}{Numeric feature matrix}
#'   \item{info.data}{Metadata with time and status}
#'   \item{survival.data}{Combined survival data}
#'   \item{time_col}{Character, time column name}
#'   \item{status_col}{Character, status column name}
#'   \item{...}{Additional slots for survival analysis results}
#' }
#'
#' @source Derived from \code{stat_obj_test}, which originates from
#'   \url{https://zenodo.org/records/7554815}
#' @keywords datasets
"pro_obj_test"


#' Example Subtyping Object for Testing
#'
#' A \code{Subtyping} S4 object converted from \code{stat_obj_test}.
#' Used for demonstration and testing of clustering/subtyping functions.
#'
#' @details
#' This object was created by converting \code{stat_obj_test} using
#' \code{CreateSubtypingObject()} and contains clustering results
#' from methods such as K-means, LPA, and NMF.
#'
#' The object contains:
#' \itemize{
#'   \item \code{@clean.data}: Numeric feature matrix
#'   \item \code{@info.data}: Metadata including clustering labels
#'   \item \code{@clustered.data}: Data with assigned cluster labels
#'   \item \code{@visualization.results}: t-SNE and UMAP embeddings
#' }
#'
#' @format A \code{Subtyping} S4 object with slots:
#' \describe{
#'   \item{clean.data}{Numeric feature matrix}
#'   \item{info.data}{Metadata data frame}
#'   \item{clustered.data}{Data frame with cluster assignments}
#'   \item{visualization.results}{List of dimensionality reduction results}
#'   \item{...}{Additional slots for clustering results}
#' }
#'
#' @source Derived from \code{stat_obj_test}, which originates from
#'   \url{https://zenodo.org/records/7554815}
#' @keywords datasets
"subtype_obj_test"


#' Example Train_Model Object for Testing
#'
#' A \code{Train_Model} S4 object converted from \code{stat_obj_test}.
#' Used for demonstration and testing of machine learning model training
#' and evaluation functions.
#'
#' @details
#' This object was created by converting \code{stat_obj_test} using
#' \code{CreateModelObject()} or \code{ModelTrainAnalysis()}, and contains
#' trained models for binary classification.
#'
#' The object contains:
#' \itemize{
#'   \item \code{@train.models}: List of trained models (glm, rf, gbm, etc.)
#'   \item \code{@split.data}: Train/test split data
#'   \item \code{@group_col}: "BloodCulture"
#'   \item \code{@split.scale.data}: Scaled train/test data
#' }
#'
#' @format A \code{Train_Model} S4 object with slots:
#' \describe{
#'   \item{train.models}{List of trained caret models}
#'   \item{split.data}{List of train/test data splits}
#'   \item{split.scale.data}{Scaled train/test data splits}
#'   \item{group_col}{Character, grouping column name}
#'   \item{...}{Additional slots for model evaluation results}
#' }
#'
#' @source Derived from \code{stat_obj_test}, which originates from
#'   \url{https://zenodo.org/records/7554815}
#' @keywords datasets
"train_obj_test"


#' Global variables used in non-standard evaluation
#' @keywords internal
#' @noRd
utils::globalVariables(c(
  ".",
  "::<-",
  "AUC",
  "AUC_mean",
  "AUC_se",
  "Actual",
  "Algorithm",
  "Best_Score",
  "C_Index",
  "CI_lower",
  "CI_upper",
  "Category",
  "Class",
  "Cluster",
  "Component",
  "Confidence",
  "Count",
  "Dataset",
  "Dim1",
  "Dim2",
  "Dimension 1",
  "Dimension 2",
  "Estimate",
  "FPR",
  "Facet_Median",
  "Feature",
  "Feature_Removed",
  "Fill",
  "Freq",
  "Frequency",
  "G",
  "Generation",
  "Group",
  "HR",
  "HR_95CI",
  "HR_label",
  "ID",
  "Importance",
  "Imputation",
  "Index",
  "Iteration",
  "Label",
  "Lower",
  "Mean",
  "Mean_AUC",
  "Method",
  "Method1",
  "Method2",
  "Metric",
  "Missing_Percentage",
  "Model",
  "NetBenefit",
  "New",
  "Normalization",
  "NumFeatures",
  "OUT_DIR",
  "Overall",
  "Overlap",
  "PC1",
  "PC2",
  "P_Value",
  "P_label",
  "P_value",
  "Parameter",
  "Pct",
  "Performance",
  "Performance_Drop",
  "Predicted",
  "Predicted_Subtype",
  "Probability",
  "Proportion",
  "Ref",
  "Risk",
  "SD",
  "SE",
  "Sample",
  "Score",
  "Selected",
  "Selected_Factor",
  "Selected_Status",
  "Sensitivity",
  "Sig",
  "Specificity",
  "Status",
  "Strategy",
  "Subgroup",
  "TPR",
  "Threshold",
  "Time",
  "TimePoint",
  "Upper",
  "Value",
  "Variable",
  "Variables",
  "bin",
  "change",
  "cluster_name",
  "contribution",
  "correct",
  "dif",
  "dropout_loss",
  "feat_scaled",
  "feature",
  "feature_value",
  "feature_value_num",
  "fill_col",
  "fill_group",
  "fpr",
  "group",
  "groups",
  "id",
  "logFC",
  "log_val",
  "lower",
  "mean_abs",
  "mean_pred",
  "mean_shap",
  "med",
  "mlr_learners",
  "n_features",
  "neg_log10p",
  "new_sens",
  "new_spec",
  "nri_type",
  "nri_value",
  "obs_rate",
  "observation",
  "observed",
  "outcome",
  "perc",
  "predicted",
  "prob",
  "rainbow",
  "ref_sens",
  "ref_spec",
  "sample_id",
  "se_loss",
  "shap",
  "shap_value",
  "significance",
  "sil_width",
  "stratum",
  "tAUC",
  "target_group",
  "threshold",
  "times",
  "tpr",
  "truth",
  "upper",
  "value",
  "variable",
  "x",
  "y",
  "y_position",
  "yhat"
))

#' Match factor levels of two data frames
#' @keywords internal
match_factor_levels <- function(data, ref) {
  for (col in intersect(colnames(data), colnames(ref))) {
    if (is.factor(ref[[col]])) {
      data[[col]] <- factor(data[[col]], levels = levels(ref[[col]]))
    }
  }
  data
}