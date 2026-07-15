#' Plot Missing Data Distribution
#'
#' This function generates density plots to visualize the distribution of missing data in the given dataset.
#' It provides a variable-wise and sample-wise missing data distribution plot, as well as a combined plot.
#' The function also saves the plots as PNG files if specified and returns a summary of missing data statistics.
#
#' @import ggplot2
#' @import here
#' @import wesanderson
#' @import stats
#' @import here 
#' @importFrom ggprism theme_prism
#' @importFrom grDevices pdf
#' @param data A data.frame containing the data to be analyzed. Missing values should be represented as `<NA>` or `NA`.
#' @param palette_name A character string specifying the palette name for the plot colors. Default is `"Royal1"`.
#'                   Supported palettes are from the `wesanderson` package (e.g., "Royal1", "Zissou1", "GrandBudapest1").
#' @param alpha A numeric value between 0 and 1 specifying the transparency level of the density plots. Default is `0.9`.
#' @param save_plots A logical value indicating whether to save the plots as PNG files. Default is `TRUE`.
#' @param save_dir A character string specifying the directory where the plots will be saved. Default is `here("StatObject")`,
#'                which uses the `here` package to generate the path.
#' @param plot_width A numeric value specifying the width of the saved plots. Default is `5`.
#' @param plot_height A numeric value specifying the height of the saved plots. Default is `5`.
#' @param base_size A numeric value specifying the base size for text and elements in the plot. Default is `14`.
#'
#' @returns A list containing the following:
#'   - `var_plot_obj`: The ggplot object for the variable-wise missing data distribution plot.
#'   - `sample_plot_obj`: The ggplot object for the sample-wise missing data distribution plot.
#'   - `combined_plot`: The ggplot object for the combined missing data distribution plot.
#'   - `total_missing_values`: The total number of missing values in the data.
#'   - `total_variables`: The total number of variables (columns) in the data.
#'   - `total_samples`: The total number of samples (rows) in the data.
#'   - `variables_with_missing`: The number of variables with at least one missing value.
#'   - `samples_with_missing`: The number of samples with at least one missing value.
#'   - `max_missing_variable`: The highest percentage of missing values in any variable (column).
#'   - `min_missing_variable`: The lowest percentage of missing values in any variable (column).
#'   - `max_missing_sample`: The highest percentage of missing values in any sample (row).
#'   - `min_missing_sample`: The lowest percentage of missing values in any sample (row).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("stat_obj_test")
#' missing_info <- plot_missing_data(data = stat_obj_test@raw.data, save_plots =FALSE)
#' # Generate missing data plots with customized parameters
#' missing_info <- plot_missing_data(data = stat_obj_test@raw.data, save_plots =FALSE,
#'                               palette_name = "Zissou1", alpha = 0.7, plot_width = 6, plot_height = 6)
#'}
plot_missing_data <- function(data,
                              palette_name = 'Royal1',
                              alpha = 0.9,
                              save_plots = TRUE,
                              save_dir = NULL,
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 14,
                              save_data = TRUE,
                              var_filename = "var_missing_data.csv",
                              sample_filename = "sample_missing_data.csv") {
  if (is.null(save_dir)) save_dir <- get_output_dir("StatObject", "missing_info")

  colors <- wes_palette(n = 3, name = palette_name, type = "discrete")
  colors <- as.list(colors)

  primary_color <- colors[[1]]
  secondary_color <- colors[[2]]
  tertiary_color <- colors[[3]]

  data[data == '<NA>'] <- NA

  var_missing_percentage <- colMeans(is.na(data)) * 100
  sample_missing_percentage <- rowMeans(is.na(data)) * 100

  if (all(var_missing_percentage == 0) && all(sample_missing_percentage == 0)) {
    cat("No missing values in the data.")
    return(list(var_plot_obj = NULL, sample_plot_obj = NULL, combined_plot = NULL))
  }

  if (mean(var_missing_percentage) < 5 && mean(sample_missing_percentage) < 5) {
    cat("The percentage of missing data is low (below 5%).")
  }

  var_missing_df <- data.frame(Variable = names(var_missing_percentage), Missing_Percentage = var_missing_percentage)
  sample_missing_df <- data.frame(Sample = 1:nrow(data), Missing_Percentage = sample_missing_percentage)

  var_plot_obj <- ggplot(data = var_missing_df, aes(x = Missing_Percentage, fill = "Variable-wise")) +
    geom_density(alpha = alpha, color = NA) +
    labs(x = "Missing Percentage", y = "Density", title = "Variable Missing Data") +
    scale_fill_manual(values = primary_color) +
    theme_prism(base_size = base_size) +
    theme(legend.position = "top")

  sample_plot_obj <- ggplot(data = sample_missing_df, aes(x = Missing_Percentage, fill = "Sample-wise")) +
    geom_density(alpha = alpha, color = NA) +
    labs(x = "Missing Percentage", y = "Density", title = "Sample Missing Data") +
    scale_fill_manual(values = secondary_color) +
    theme_prism(base_size = base_size) +
    theme(legend.position = "none")

  combined_plot <- var_plot_obj +
    geom_density(data = sample_missing_df, aes(x = Missing_Percentage, fill = "Sample-wise"), alpha = alpha, color = NA) +
    labs(y = "Density", title = "Overall Missing Data Distribution") +
    scale_fill_manual(values = c(primary_color, secondary_color)) +
    guides(fill = guide_legend(title = NULL)) +
    theme_prism(base_size = base_size)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    cat("\nAll plots saved successfully to: \n",save_dir,"\n")
    
    ggsave(filename = file.path(save_dir, "variable_missing_data_plot.pdf"), plot = var_plot_obj,
           width = plot_width, height = plot_height, device = "pdf")
    ggsave(filename = file.path(save_dir, "sample_missing_data_plot.pdf"), plot = sample_plot_obj,
           width = plot_width, height = plot_height, device = "pdf")
    ggsave(filename = file.path(save_dir, "combined_missing_data_plot.pdf"), plot = combined_plot,
           width = plot_width, height = plot_height, device = "pdf")
  }
  
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)}
    full_path1 <- file.path(save_dir, var_filename)
    write.csv(var_missing_df, file = full_path1, row.names = FALSE)
    full_path2 <- file.path(save_dir, sample_filename)
    write.csv(sample_missing_df, file = full_path2, row.names = FALSE)
    cat("Saved variable missing rate data to:", full_path1, "\n")
    cat("Saved sample missing rate data to:", full_path2, "\n")
  }
  
  missing_info <- list(
    var_plot_obj = var_plot_obj,
    sample_plot_obj = sample_plot_obj,
    combined_plot = combined_plot,
    total_missing_values = sum(is.na(data)),
    total_variables = ncol(data),
    total_samples = nrow(data),
    variables_with_missing = sum(var_missing_percentage > 0),
    samples_with_missing = sum(sample_missing_percentage > 0),
    max_missing_variable = max(var_missing_percentage),
    min_missing_variable = min(var_missing_percentage),
    max_missing_sample = max(sample_missing_percentage),
    min_missing_sample = min(sample_missing_percentage)
  )
  print(combined_plot)
  return(missing_info)
}


#' Plot Missing Data for a Stat Object or Data Frame
#'
#' This function generates and saves missing data plots (variable-wise and sample-wise distributions) for
#' an object of class 'Stat' or a data frame. It also updates the 'Stat' object with missing data information
#' and returns the object or the missing data summary, depending on the input type.
#'
#' @import ggplot2
#' @import here
#' @import wesanderson
#' @import stats
#' @import methods
#' @import here 
#' @importFrom ggprism theme_prism
#' @importFrom grDevices pdf
#' @param object An object of class 'Stat' or a data frame. If it is a 'Stat' object, the function uses the
#'               `raw.data` slot for missing data analysis. If it is a data frame, it directly performs the analysis.
#' @param palette_name A character string specifying the palette name for the plot colors. Default is `"Royal1"`.
#' @param alpha A numeric value between 0 and 1 specifying the transparency level of the density plots. Default is `0.9`.
#' @param save_plots A logical value indicating whether to save the plots as PNG files. Default is `TRUE`.
#' @param save_dir A character string specifying the directory where the plots will be saved. Default is `"here('StatObject', 'missing_info')"` (within the `StatObject` folder).
#' @param plot_width A numeric value specifying the width of the saved plots. Default is `5`.
#' @param plot_height A numeric value specifying the height of the saved plots. Default is `5`.
#'
#' @returns If the input is a 'Stat' object, returns the updated 'Stat' object with the missing data information stored
#'          in the `process.info` slot. If the input is a data frame, returns the missing data summary as a list.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("stat_obj_test")
#' # Generate missing data plots for a Stat object
#'updated_stat <- state_plot_missing_data(stat_obj_test, save_plots = TRUE)
#' # Generate missing data plots for a data frame
#'missing_info <- state_plot_missing_data(stat_obj_test@raw.data, save_plots = FALSE)
#'}
state_plot_missing_data <- function(
    object,
    palette_name = 'Royal1',
    alpha = 0.9,
    save_plots = TRUE,
    save_dir = NULL,
    plot_width = 5,
    plot_height = 5,
    save_data = TRUE,
    var_filename = "var_missing_data.csv",
    sample_filename = "sample_missing_data.csv") {
  if (is.null(save_dir)) save_dir <- get_output_dir("StatObject", "missing_info")

  if (inherits(object, "Stat")) {
    data <- slot(object, "raw.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input.")
  }

  missing_info <- plot_missing_data(data,
                                    palette_name = palette_name,
                                    alpha = alpha,
                                    save_plots = save_plots,
                                    save_dir = save_dir,
                                    plot_width = plot_width,
                                    plot_height = plot_height,
                                    save_data = save_data,
                                    var_filename = var_filename,
                                    sample_filename = sample_filename)

  print(missing_info)

  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@process.info[["missing_count"]] <- missing_info
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'process.info' slot updated.\n")
    return(object)
  }
  return(missing_info)
}

#' Extract Raw Data from Stat Object
#'
#' This function extracts the 'raw.data' slot from an object of class 'Stat'.
#' If the object is not of class 'Stat' or does not contain a 'raw.data' slot,
#' it will return NULL.
#'
#' @param object An object of class 'Stat' which contains a slot named 'raw.data'.
#'               This should be a valid Stat object.
#'
#' @returns Returns the 'raw.data' slot of the Stat object if it exists.
#'          If the object does not have a 'raw.data' slot or the slot is empty,
#'          it returns NULL.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("stat_obj_test")
#' # Assuming 'stat_object' is a valid Stat object
#' raw_data <- ExtractRawData(stat_obj_test)
#' }
ExtractRawData <- function(object) {
  data <- tryCatch(slot(object, "raw.data"), error = function(e) NULL)
  return(data)
}


#' Diagnose Variable Types in Data
#'
#' Analyzes a data frame to classify variables into numeric and categorical
#' types based on both their storage mode and cardinality. Unlike the original
#' version, this function respects numeric columns and only converts them to
#' categorical upon explicit request.
#'
#' The function excludes a specified grouping column (if present) from analysis.
#' Categorical variables (character, factor, logical) are always recognized as
#' such. Numeric columns are kept as numeric by default, unless
#' `treat_low_card_numeric_as_categorical = TRUE` is set, in which case numeric
#' columns with few distinct values become categorical.
#'
#' @param data A data frame.
#' @param group_col A character string naming the grouping column to be excluded
#'   (default is `"group"`). Set to `NULL` to skip exclusion.
#' @param max_unique_values Numeric, maximum number of unique values a column
#'   may have to be considered "low cardinality" when
#'   `treat_low_card_numeric_as_categorical = TRUE`. Default is 5.
#' @param encode_threshold Numeric, minimum number of unique values for a
#'   categorical variable to be flagged for encoding. Default is 10.
#' @param treat_low_card_numeric_as_categorical Logical, whether to treat
#'   low‑cardinality numeric columns (e.g. 0/1, 1–5 ratings) as categorical
#'   variables. Default is `FALSE`.
#'
#' @return A list with three components:
#' \describe{
#' \item{numeric_vars}{Character vector of variable names classified as numeric.}
#' \item{categorical_vars}{Character vector of variable names classified as
#'   categorical (includes characters, factors, logicals, and optionally
#'   low‑cardinality numerics).}
#' \item{vars_to_encode}{Character vector of categorical variables that have
#'   more than `encode_threshold` unique values and may need encoding.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example data
#' df <- data.frame(
#'   group = rep(1:2, each = 5),
#'   age = c(25, 30, 35, 40, 45, 50, 55, 60, 65, 70),
#'   rating = c(1,2,1,2,3,2,1,3,2,1),   # low‑cardinality numeric
#'   city = c("A","B","A","C","D","E","F","G","H","I"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Default: rating stays numeric
#' res <- diagnose_variable_type(df, group_col = "group")
#' print(res)
#'
#' # Force low‑cardinality numerics to categorical
#' res2 <- diagnose_variable_type(df, group_col = "group",
#'                                treat_low_card_numeric_as_categorical = TRUE)
#' print(res2)
#' }
diagnose_variable_type <- function(data,
                                   group_col = "group",
                                   max_unique_values = 5,
                                   encode_threshold = 10,
                                   treat_low_card_numeric_as_categorical = FALSE) {
  # ----- Input validation -------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.null(group_col)) {
    if (length(group_col) != 1L || !is.character(group_col)) {
      stop("`group_col` must be a single character string or NULL.", call. = FALSE)
    }
  }
  if (!is.numeric(max_unique_values) || length(max_unique_values) != 1L || max_unique_values <= 0) {
    stop("`max_unique_values` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(encode_threshold) || length(encode_threshold) != 1L || encode_threshold <= 0) {
    stop("`encode_threshold` must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.logical(treat_low_card_numeric_as_categorical) || length(treat_low_card_numeric_as_categorical) != 1L) {
    stop("`treat_low_card_numeric_as_categorical` must be a logical scalar.", call. = FALSE)
  }
  
  # ----- Initialise result vectors ----------------------------------------
  numeric_vars <- character()
  categorical_vars <- character()
  vars_to_encode <- character()
  
  # ----- Loop over columns ------------------------------------------------
  for (col in names(data)) {
    # Skip grouping column
    if (!is.null(group_col) && col == group_col) next
    
    col_data <- data[[col]]
    n_unique <- length(unique(col_data))
    is_numeric <- is.numeric(col_data)  # covers integer and numeric
    
    if (is_numeric) {
      if (treat_low_card_numeric_as_categorical && n_unique <= max_unique_values) {
        # Numeric column treated as categorical because of low cardinality
        categorical_vars <- c(categorical_vars, col)
        if (n_unique > encode_threshold) {
          vars_to_encode <- c(vars_to_encode, col)
        }
      } else {
        # Normal numeric column
        numeric_vars <- c(numeric_vars, col)
      }
    } else {
      # Non‑numeric columns (character, factor, logical, etc.) are categorical
      categorical_vars <- c(categorical_vars, col)
      if (n_unique > encode_threshold) {
        vars_to_encode <- c(vars_to_encode, col)
      }
    }
  }
  
  # ----- Return -----------------------------------------------------------
  list(
    numeric_vars = numeric_vars,
    categorical_vars = categorical_vars,
    vars_to_encode = vars_to_encode
  )
}

#' Diagnose Variable Types for 'Stat' Objects or Data Frames
#'
#' This function analyzes the variable types (numeric, categorical, and those
#' needing encoding) of a data frame or an object of class "Stat". If the input
#' is a "Stat" object, it extracts the raw data and group column from the object.
#' If the input is a data frame, it directly uses the provided data for diagnosis.
#' It updates the "Stat" object with the diagnosed variable types.
#'
#' @param object An object of class "Stat" or a data frame. If of class "Stat",
#'   the raw data and group column will be extracted. If a data frame, the function
#'   operates directly on it.
#' @param group_col A character string specifying the grouping column (default "group").
#'   This column is excluded from analysis if present. When `object` is a "Stat"
#'   object, the group column stored in the object overrides this argument.
#' @param max_unique_values Numeric, maximum number of unique values a column may have
#'   to be considered "low cardinality" when `treat_low_card_numeric_as_categorical = TRUE`.
#'   Default is 5.
#' @param encode_threshold Numeric, minimum number of unique values for a categorical
#'   variable to be flagged for encoding. Default is 10.
#' @param treat_low_card_numeric_as_categorical Logical, whether to treat low‑cardinality
#'   numeric columns (e.g. 0/1, 1–5 ratings) as categorical variables. Default is `FALSE`.
#'
#' @returns If input is a "Stat" object, the updated object with the diagnosed variable
#'   types stored in the `variable.types` slot. If input is a data frame, a list containing:
#'   \describe{
#'   \item{numeric_vars}{Character vector of numeric variables.}
#'   \item{categorical_vars}{Character vector of categorical variables.}
#'   \item{vars_to_encode}{Character vector of categorical variables that have more than
#'     `encode_threshold` unique values and may need encoding.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Diagnose variables in a "Stat" object
#' stat_obj <- stat_diagnose_variable_type(stat_obj_test, group_col = "group")
#' print(stat_obj)
#'
#' # Example 2: Diagnose variables in a data frame
#' result <- stat_diagnose_variable_type(stat_obj_test@raw.data,
#'                                       treat_low_card_numeric_as_categorical = TRUE)
#' print(result)
#' }
stat_diagnose_variable_type <- function(object,
                                        group_col = "group",
                                        max_unique_values = 5,
                                        encode_threshold = 10,
                                        treat_low_card_numeric_as_categorical = FALSE) {
  
  # ----- Input validation -------------------------------------------------
  if (!inherits(object, "Stat") && !is.data.frame(object)) {
    stop("Input must be an object of class 'Stat' or a data frame", call. = FALSE)
  }
  
  # ----- Extract data and group column ------------------------------------
  if (inherits(object, "Stat")) {
    # For Stat objects, group column is taken from the object's slot
    group_col_slot <- slot(object, "group_col")
    if (length(group_col_slot) > 0 && !is.null(group_col_slot)) {
      group_col <- group_col_slot
    } else {
      group_col <- NULL
    }
    data <- slot(object, "raw.data")
  } else {
    # For data frames, use provided group_col
    data <- object
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input", call. = FALSE)
  }
  
  # ----- Call the core diagnosis function --------------------------------
  variable_types <- diagnose_variable_type(
    data = data,
    group_col = group_col,
    max_unique_values = max_unique_values,
    encode_threshold = encode_threshold,
    treat_low_card_numeric_as_categorical = treat_low_card_numeric_as_categorical
  )
  
  # ----- Report summary ---------------------------------------------------
  cat("Diagnosed variable types:\n")
  cat("Numeric variables:", length(variable_types$numeric_vars), "\n")
  cat("Categorical variables:", length(variable_types$categorical_vars), "\n")
  cat("Variables flagged for encoding:", length(variable_types$vars_to_encode), "\n")
  
  if (length(variable_types$numeric_vars) == 0 && length(variable_types$categorical_vars) == 0) {
    stop("No valid variables found after variable type diagnosis", call. = FALSE)
  }
  
  # ----- Update Stat object or return list -------------------------------
  if (inherits(object, "Stat")) {
    cat("Updating 'Stat' object...\n")
    object@variable.types <- variable_types
    cat("The 'Stat' object has been updated with the following slots:\n")
    cat("- 'variable.types' slot updated.\n")
    return(object)
  }
  
  return(variable_types)
}


#' Gaze Analysis for Group Comparison
#'
#' This function performs gaze analysis for group comparison using the specified formula,
#' method, and settings. It provides a detailed result, optionally saving it as a Word document.
#' The gaze analysis can be performed based on group columns in the data, and various settings
#' are available to customize the output.
#'
#' @import here
#' @import officer
#' @import flextable
#' @import stats
#' @import methods
#' @import autoReg
#' @param data A data frame containing the data to analyze.
#' @param formula A formula specifying the model to fit. If NULL, the formula is automatically created
#' based on the provided group columns (default is NULL).
#' @param group_cols A vector of column names representing the grouping variables (default is NULL).
#' @param digits The number of digits to display for the result (default is 1).
#' @param show.p A logical value indicating whether to display p-values in the output (default is TRUE).
#' @param gaze_method An integer between 1 and 5 representing the gaze analysis method to use (default is 3).
#' @param save_word A logical value indicating whether to save the results to a Word document (default is TRUE).
#' @param save_dir The directory to save the Word document (default is the "PrognosiX/gaze_baseline" folder).
#'
#' @returns A data frame or matrix with the gaze analysis results. If `save_word` is TRUE, the results
#' are also saved as a Word document in the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#'my_data=stat_obj_test@clean.data
#' # Example 1: Performing gaze analysis with a formula and custom settings
#'result <- gaze_analysis(data = my_data,formula = ~ SWAB + AGE,digits = 2,show.p = TRUE,
#'                     gaze_method = 3,save_word = TRUE, save_dir = "./")
#' # Example 2: Using the default formula based on group columns
#'  result <- gaze_analysis(data = my_data,group_cols = c("SWAB"),
#'                         digits = 1,show.p = FALSE,save_word = TRUE,
#'                         save_dir = "./",gaze_method = 1)
#' }
gaze_analysis <- function(data,
                          formula = NULL,
                          group_cols = NULL,
                          digits = 1,
                          show.p = TRUE,
                          gaze_method = 3,
                          save_word = TRUE,
                          save_dir = NULL) {

  if (!is.data.frame(data)) stop("The input 'data' must be a data frame.")

  if (is.null(formula)) {
    if (!is.null(group_cols)) {
      if (length(group_cols) == 1 && group_cols %in% colnames(data)) {
        formula <- as.formula(paste(group_cols, "~ ."))
        cat("Using formula with one group:", deparse(formula), "\n")
      } else if (length(group_cols) > 1 && all(group_cols %in% colnames(data))) {
        formula <- as.formula(paste(paste(group_cols, collapse = " + "), "~ ."))
        cat("Using formula with multiple groups:", deparse(formula), "\n")
      } else {
        stop("Group columns not found in data.")
      }
    } else {
      formula <- as.formula("~ .")
      cat("Using default formula: ~ .\n")
    }
  } else if (!inherits(formula, "formula")) {
    stop("The input 'formula' must be a valid formula.")
  }

  if (!is.numeric(digits) || digits < 0 || digits != as.integer(digits))
    stop("The input 'digits' must be a non-negative integer.")

  if (!is.logical(show.p))
    stop("The input 'show.p' must be a logical value (TRUE or FALSE).")

  if (!is.numeric(gaze_method) || gaze_method < 1 || gaze_method > 5 || gaze_method != as.integer(gaze_method))
    stop("The input 'gaze_method' must be an integer between 1 and 5.")

  tryCatch({
    cat("Running gaze analysis with method:", gaze_method, "\n")
    result <- gaze(formula, data, digits = digits, show.p = show.p, method = gaze_method)

    cat("Result type:", class(result), "\n")
    if (is.data.frame(result) || is.matrix(result)) {
      result <- myft(result)
      cat("Gaze analysis completed successfully.\n")

      if (save_word) {
        doc <- read_docx()
        doc <- doc %>%
          body_add_flextable(result) %>%
          body_add_par("Gaze Analysis Results", style = "heading 1")

        if (!dir.exists(save_dir)) {
          dir.create(save_dir, recursive = TRUE)
        }

        word_filename <- file.path(save_dir, "gaze_analysis.docx")
        print(doc, target = word_filename)
        cat("Word file saved to:", word_filename, "\n")
      }

      return(result)
    } else {
      stop("The result is not a data frame or matrix.")
    }
  }, error = function(e) {
    stop("An error occurred while performing the gaze analysis: ", e$message)
  })
}

#' Statistical Gaze Analysis
#' Safe gaze analysis for Table 1 (auto‑removes problematic columns)
#'
#' Wrapper around \code{stat_gaze_analysis} that first drops zero‑variance
#' and all‑missing columns from \code{clean.data} to avoid errors in
#' \code{autoReg::gaze}.
#'
#' @param object A \code{Stat} object or data frame.
#' @param formula Optional formula.
#' @param group_col Group column name.
#' @param digits Number of digits for statistics.
#' @param show.p Logical. Show p‑values?
#' @param gaze_method Integer 1–5.
#' @param save_word Logical. Save Word file?
#' @param save_dir Output directory. If \code{NULL} and \code{save_word = TRUE},
#'   uses \code{./tables/}.
#' @return Updated \code{Stat} object or flextable.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # --- Prepare example data ---
#' library(autoReg)
#' library(flextable)
#' 
#' data("mtcars")
#' mtcars$cyl <- factor(mtcars$cyl)
#' mtcars$vs  <- factor(mtcars$vs, labels = c("V-shaped", "Straight"))
#' 
#' # Create a Stat object with "vs" as the grouping variable
#' stat_obj <- CreateStatObject(raw.data = mtcars, group_col = "vs")
#' 
#' # --- Example 1: Basic usage with 2 groups (automatic formula) ---
#' stat_obj <- stat_gaze_analysis(stat_obj,
#'                                digits = 2,
#'                                show.p = TRUE,
#'                                gaze_method = 3)
#' # View the formatted table
#' stat_obj@baseline.table
#' 
#' # --- Example 2: Direct use with a data frame, no Word output ---
#' result <- stat_gaze_analysis(mtcars,
#'                              group_col = "cyl",    # >2 groups
#'                              digits = 1,
#'                              show.p = FALSE,
#'                              gaze_method = 3,
#'                              save_word = FALSE)
#' result
#' 
#' # --- Example 3: Custom formula to compare selected variables ---
#' stat_obj <- stat_gaze_analysis(stat_obj,
#'                                formula = vs ~ mpg + hp + wt,
#'                                digits = 2,
#'                                show.p = TRUE,
#'                                save_word = FALSE)
#' 
#' # --- Example 4: Different gaze method (method 1: t-test / Wilcoxon) ---
#' stat_obj <- stat_gaze_analysis(stat_obj,
#'                                gaze_method = 1,
#'                                digits = 1,
#'                                save_word = FALSE)
#' 
#' # --- Example 5: Save as Word document with custom directory ---
#' stat_obj <- stat_gaze_analysis(stat_obj,
#'                                save_word = TRUE,
#'                                save_dir = "my_results/tables")
#' # The file is saved as my_results/tables/gaze_analysis.docx
#' 
#' # --- Example 6: Suppress p‑values (descriptive statistics only) ---
#' stat_obj <- stat_gaze_analysis(stat_obj,
#'                                show.p = FALSE,
#'                                save_word = FALSE)
#' 
#' # --- Example 7: Multi‑group comparison (>2 groups) ---
#' result_multi <- stat_gaze_analysis(mtcars,
#'                                    group_col = "cyl",
#'                                    digits = 2,
#'                                    show.p = TRUE,
#'                                    save_word = FALSE)
#'}
stat_gaze_analysis <- function(object,
                               formula = NULL,
                               group_col = "group",
                               digits = 1,
                               show.p = TRUE,
                               gaze_method = 3,
                               save_word = TRUE,
                               save_dir = NULL) {
  
  cat("Starting stat_gaze_analysis (safe version)...\n")
  
  # ---- 1. Extract data ----
  if (inherits(object, "Stat")) {
    data <- slot(object, "clean.data")
    if (is.null(data) || nrow(data) == 0) {
      data <- slot(object, "raw.data")
    }
    group_col <- slot(object, "group_col")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }
  
  if (is.null(data) || nrow(data) == 0)
    stop("No valid data found in the input.")
  
  # ---- 2. Sanitize grouping column ----
  if (!is.null(group_col) && !group_col %in% colnames(data))
    group_col <- NULL
  
  # ---- 3. Drop problematic columns ----
  # (a) Remove all‑NA columns
  na_cols <- sapply(data, function(x) all(is.na(x)))
  if (any(na_cols)) {
    cat("Removing all‑NA columns:", paste(names(data)[na_cols], collapse = ", "), "\n")
    data <- data[, !na_cols, drop = FALSE]
  }
  
  # (b) Remove zero‑variance numeric columns
  numeric_cols <- sapply(data, is.numeric)
  if (any(numeric_cols)) {
    zero_var <- sapply(data[, numeric_cols, drop = FALSE], var, na.rm = TRUE) == 0
    if (any(zero_var)) {
      cat("Removing zero‑variance numeric columns:",
          paste(names(data[, numeric_cols])[zero_var], collapse = ", "), "\n")
      data <- data[, !(colnames(data) %in% names(data[, numeric_cols])[zero_var]), drop = FALSE]
    }
  }
  
  # (c) Remove factor columns with only one level (excluding group_col)
  factor_cols <- sapply(data, is.factor)
  if (any(factor_cols)) {
    single_level <- sapply(data[, factor_cols, drop = FALSE],
                           function(x) length(levels(droplevels(x))) < 2)
    # keep group_col if present
    if (!is.null(group_col) && group_col %in% names(single_level))
      single_level[group_col] <- FALSE
    if (any(single_level)) {
      cat("Removing single‑level factor columns:",
          paste(names(data[, factor_cols])[single_level], collapse = ", "), "\n")
      data <- data[, !(colnames(data) %in% names(data[, factor_cols])[single_level]), drop = FALSE]
    }
  }
  
  cat("Data prepared for gaze analysis. Number of rows:", nrow(data),
      "Columns:", ncol(data), "\n")
  
  # ---- 4. Default save_dir for Word output ----
  if (save_word && is.null(save_dir)) {
    if (exists("get_output_dir")) {
      save_dir <- get_output_dir("m1", "baseline_tables")
    } else {
      save_dir <- file.path(".", "tables")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  }
  
  # ---- 5. Call gaze_analysis with cleaned data ----
  result <- gaze_analysis(data,
                          formula = formula,
                          group_cols = group_col,
                          digits = digits,
                          show.p = show.p,
                          gaze_method = gaze_method,
                          save_word = save_word,
                          save_dir = save_dir)
  
  print(result)
  
  # ---- 6. Update Stat object if needed ----
  if (inherits(object, "Stat")) {
    object@baseline.table <- result
    cat("'baseline.table' slot updated.\n")
    return(object)
  }
  
  return(result)
}
