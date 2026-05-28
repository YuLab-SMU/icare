
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

