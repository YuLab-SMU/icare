#' Convert variable types
#'
#' @param data A data frame or an S4 object (Stat, Train_Model, Subtyping, PrognosiX).
#' @param variable_types A list containing variable type information.
#' @param save_dir Directory to save the cleaned data.
#' @param save_data Logical. Whether to save the data.
#' @param group_col The name of the group column.
#' @param csv_filename The name of the CSV file.
#' @export
#' @examples
#' \dontrun{
#' types <- diagnose_variable_type(stat_obj_test@clean.data, group_col = "SWAB")
#' converted <- convert_variables(stat_obj_test@clean.data, types, 
#' group_col = "SWAB", save_data = FALSE)
#' }
convert_variables <- function(data,
                              variable_types,
                              save_dir = NULL,
                              save_data = F,
                              group_col = NULL,
                              csv_filename = "clean_data.csv") {
  if (is.null(save_dir)) save_dir <- get_output_dir("StatObject", "Data")
  
  is_s4 <- inherits(data, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  if (is_s4) {
    df <- if (inherits(data, "Stat")) data@clean.data else
          if (inherits(data, "Subtyping")) data@clean.data else
          if (inherits(data, "PrognosiX")) data@clean.data else
          if (inherits(data, "Train_Model")) data@clean.df
  } else {
    stopifnot(is.data.frame(data))
    df <- data
  }
  
  all_factor_cols <- names(df)
  
  if (!is.null(group_col) && group_col %in% all_factor_cols) {
    all_factor_cols <- all_factor_cols[all_factor_cols != group_col]
  }
  for (col in all_factor_cols) {
    if (col %in% variable_types$categorical_vars) {
      df[[col]] <- factor(df[[col]])
      cat("Converted", col, "to factor.\n")
    } else {
      df[[col]] <- as.numeric(df[[col]])
      cat("Converted ", col, " to numeric.\n")
    }
  }
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path <- file.path(save_dir, csv_filename)
    write.csv(df, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  
  if (is_s4) {
    if (inherits(data, "Stat")) data@clean.data <- df else
    if (inherits(data, "Subtyping")) data@clean.data <- df else
    if (inherits(data, "PrognosiX")) data@clean.data <- df else
    if (inherits(data, "Train_Model")) data@clean.df <- df
    return(data)
  }
  
  return(df)
}

#' One hot encode
#'
#' @param data A data frame or an S4 object (Stat, Train_Model, Subtyping, PrognosiX).
#' @param group_col The name of the group column.
#' @param max_unique_values Maximum unique values for encoding.
#' @param save_dir Directory to save.
#' @param save_data Logical.
#' @param csv_filename Filename.
#' @export
#' @examples
#' \dontrun{
#' encoded <- one_hot_encode(stat_obj_test@clean.data, group_col = "SWAB", save_data = FALSE)
#' }
one_hot_encode <- function(data,
                           group_col = "group",
                           max_unique_values = 5,
                           save_dir = NULL,
                           save_data = TRUE,
                           csv_filename = "clean_data.csv") {
  if (is.null(save_dir)) save_dir <- get_output_dir("StatObject", "Data")
  
  is_s4 <- inherits(data, c("Stat", "Train_Model", "Subtyping", "PrognosiX"))
  if (is_s4) {
    df <- if (inherits(data, "Stat")) data@clean.data else
          if (inherits(data, "Subtyping")) data@clean.data else
          if (inherits(data, "PrognosiX")) data@clean.data else
          if (inherits(data, "Train_Model")) data@clean.df
  } else {
    if (!is.data.frame(data)) stop("Input must be a data frame or valid S4 object")
    df <- data
  }
  
  if (length(group_col) == 0 || !is.character(group_col) || !(group_col %in% colnames(df))) {
    cat("Group column is not valid, setting to NULL.\n")
    group_col <- NULL
  }
  
  if (!is.numeric(max_unique_values) || max_unique_values <= 0) {
    stop("max_unique_values must be a positive numeric value")
  }
  
  variable_types <- diagnose_variable_type(df, group_col = group_col, max_unique_values = max_unique_values)
  vars_to_encode <- variable_types$vars_to_encode
  
  encoded_data <- df
  row_names <- rownames(df)
  
  for (var in vars_to_encode) { 
    unique_values <- unique(df[!is.na(df[, var]), var])
    cat("Encoding variable:", var, "with unique values:", paste(unique_values, collapse = ", "), "\n")
    for (value in unique_values) {
      col_name <- paste(var, value, sep = "_")
      encoded_data[, col_name] <- as.integer(df[, var] == value)
    }
    encoded_data[, var] <- NULL
  }
  
  if (!is.null(group_col) && group_col %in% names(encoded_data)) {
    group <- encoded_data[[group_col]]
    encoded_data[[group_col]] <- NULL
    encoded_data <- cbind(encoded_data, group)
    colnames(encoded_data)[ncol(encoded_data)] <- group_col
  }
  
  rownames(encoded_data) <- row_names
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    full_path <- file.path(save_dir, csv_filename)
    write.csv(encoded_data, file = full_path, row.names = FALSE)
    cat("Cleaned data saved to:", full_path, "\n")  
  }
  
  if (is_s4) {
    if (inherits(data, "Stat")) data@clean.data <- encoded_data else
    if (inherits(data, "Subtyping")) data@clean.data <- encoded_data else
    if (inherits(data, "PrognosiX")) {
      data@clean.data <- encoded_data
      time_col <- data@time_col
      status_col <- data@status_col
      if (!is.null(time_col) && !is.null(status_col) && 
          time_col %in% colnames(data@info.data) && status_col %in% colnames(data@info.data)) {
        data@survival.data <- cbind(encoded_data, data@info.data[, c(time_col, status_col), drop = FALSE])
        data@survival.data <- data@survival.data[complete.cases(data@survival.data[, c(time_col, status_col)]), ]
      }
    } else
    if (inherits(data, "Train_Model")) data@clean.df <- encoded_data
    return(data)
  }
  
  return(encoded_data)
}
