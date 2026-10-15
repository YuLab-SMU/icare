#' Modify column names of a data frame
#'
#' This function modifies the column names of a data frame by performing the following:
#' 1. Converts all column names to lowercase.
#' 2. Replaces spaces with underscores.
#' 3. Removes any non-alphanumeric characters except for underscores.
#' 4. Ensures all column names are unique by appending numbers to duplicated names.
#'
#' @param data A data frame whose column names need to be modified.
#' @param tolower Logical. If \code{TRUE}, column names are converted to
#'   lowercase before other modifications. Default is \code{FALSE}.
#' @returns A data frame with modified column names.
#' @export
#'
#' @examples
#' data <- data.frame("First Name" = c("John", "Jane"), "Age" = c(25, 30))
#' modified_data <- modify_column_names(data)
#' print(modified_data)

modify_column_names <- function(data,tolower=F) {
  if (is.null(names(data))) stop("Column names cannot be NULL!")
  if(tolower){ new_names <- tolower(names(data)) }   
  new_names=names(data)
  new_names <- gsub(" ", "_", new_names)
  new_names <- gsub("[^[:alnum:]_]", "", new_names)
  if (anyDuplicated(new_names)) {
    warning("Duplicate column names found; made unique.")
    new_names <- make.unique(new_names)
  }
  names(data) <- new_names
  return(data)
}


#' Class to store statistical analysis results
#'
#' This class is designed to store and manage various components of statistical analysis.
#' It includes raw and cleaned data, processed results, descriptive statistics, correlation results,
#' variable types, and metadata related to the analysis.
#'
#' @import methods
#' @slot raw.data A data.frame containing the raw input data.
#' @slot clean.data A data.frame containing the cleaned data after preprocessing.
#' @slot info.data A data.frame containing additional information related to the data.
#' @slot scale.data A data.frame containing the scaled version of the data.
#' @slot meta.featurename A character vector containing metadata feature names.
#' @slot group_col The grouping column in the data, which can be of any type.
#' @slot baseline.table A table that stores baseline information, can be of any type.
#' @slot process.info A list containing information about the processing steps.
#' @slot compute.descriptive A list containing results of descriptive statistics computations.
#' @slot corr.result A list containing the results of correlation analysis.
#' @slot var.result A list containing the results of variance analysis.
#' @slot variable.types A list specifying the types of the variables in the dataset.
#'
#' @returns An object of class 'Stat' containing the slots mentioned above.
#' @export
Stat <- setClass(
  Class = "Stat",
  slots = c(
    raw.data = "data.frame",
    clean.data = "data.frame",
    info.data = "data.frame",
    scale.data = "data.frame",
    meta.featurename = "character",
    group_col = "ANY",
    baseline.table = "ANY",
    process.info = "list",
    compute.descriptive = "list",
    corr.result = "list",
    var.result = "list",
    variable.types = "list"
  ),
  prototype = list(
    raw.data = data.frame(),
    clean.data = data.frame(),
    info.data = data.frame(),
    scale.data = data.frame(),
    meta.featurename = character(0),  
    group_col = NULL,
    baseline.table = NULL,
    process.info = list(),
    compute.descriptive = list(),
    corr.result = list(),
    var.result = list(),
    variable.types = list()
  ),
  validity = function(object) {
    # Check data frame slots
    if (!is.data.frame(object@raw.data)) {
      return("Slot 'raw.data' must be a data.frame")
    }
    if (!is.data.frame(object@clean.data)) {
      return("Slot 'clean.data' must be a data.frame")
    }
    if (!is.data.frame(object@info.data)) {
      return("Slot 'info.data' must be a data.frame")
    }
    if (!is.data.frame(object@scale.data)) {
      return("Slot 'scale.data' must be a data.frame")
    }
    
    # Check character slot
    if (!is.character(object@meta.featurename)) {
      return("Slot 'meta.featurename' must be a character vector")
    }
    
    # Check list slots
    if (!is.list(object@process.info)) {
      return("Slot 'process.info' must be a list")
    }
    if (!is.list(object@compute.descriptive)) {
      return("Slot 'compute.descriptive' must be a list")
    }
    if (!is.list(object@corr.result)) {
      return("Slot 'corr.result' must be a list")
    }
    if (!is.list(object@var.result)) {
      return("Slot 'var.result' must be a list")
    }
    if (!is.list(object@variable.types)) {
      return("Slot 'variable.types' must be a list")
    }
    
    # Row alignment. info.data is append-only: it must COVER clean.data, not
    # equal it, because operations that drop samples from clean.data leave
    # info.data untouched and labels are aligned by row name at read time.
    if (nrow(object@clean.data) > 0 && nrow(object@info.data) > 0) {
      if (!all(rownames(object@clean.data) %in% rownames(object@info.data))) {
        return("Row names of 'clean.data' must all be present in 'info.data'")
      }
    }

    # Label location. Note this only runs on new()/validObject(); `@<-` does not
    # re-validate, so .stat_validate() is called explicitly at write sites.
    if (length(object@group_col) == 1 && !is.null(object@group_col)) {
      gc <- as.character(object@group_col)
      for (nm in c("raw.data", "clean.data")) {
        d <- if (nm == "raw.data") object@raw.data else object@clean.data
        if (nrow(d) > 0 && gc %in% colnames(d)) {
          return(sprintf("Label column '%s' must not appear in %s", gc, nm))
        }
      }
      if (nrow(object@info.data) > 0 && !(gc %in% colnames(object@info.data))) {
        return(sprintf("Label column '%s' not found in info.data", gc))
      }
    }

    return(TRUE)
  }
)

#' Create a Stat object for statistical analysis
#'
#' This function creates a `Stat` object, which is used to store and manage various components
#' of statistical analysis, such as raw data, cleaned data, additional metadata, and processing
#' information. It performs basic checks and preparation on the input data before creating the object.
#' 
#' @import methods
#' @param raw.data A data.frame containing the raw features. Defaults to an empty
#'   data frame. Must not contain the group column.
#' @param clean.data A data.frame containing pre-cleaned features. Defaults to an
#'   empty data frame. Must not contain the group column.
#' @param info.data A data.frame of sample metadata. This is the only place the
#'   group/outcome column may live.
#' @param group_col A character string naming the grouping column **in `info.data`**.
#' @param ... Additional arguments passed to methods (not used in this function).
#' @param na.action Character string specifying how to handle NA values.
#'   Options are "allow" (keep NA values, default), "omit" (remove rows with NA),
#'   or "error" (stop if NA values are found).
#' @returns An object of class `Stat`, which contains the processed data and metadata as slots.
#' @export
#'
#' @details
#' **Label location.** `group_col` must be a column of `info.data`. It must NOT be
#' a column of `raw.data` or `clean.data`; passing it there is an error, with no
#' automatic migration. The reason is that a label sitting in the feature matrix
#' is indistinguishable from a feature to any code that selects columns by type:
#' a numeric `0`/`1` label is silently consumed by `PlotPCA()` and by
#' `PlotCorrelationHeatmap()`, which inflates the apparent group separation
#' without any warning. Keeping the label out of the matrix makes that class of
#' bug unrepresentable rather than merely discouraged.
#'
#' `info.data` is append-only. Operations that drop samples from `clean.data`
#' leave `info.data` untouched; labels are aligned by row name at read time.
#'
#' @examples
#' raw <- data.frame(id = 1:10, value = rnorm(10))
#' inf <- data.frame(group = rep(c("A","B"), each = 5))
#' rownames(raw) <- rownames(inf) <- paste0("S", 1:10)
#' stat_new <- CreateStatObject(raw.data = raw, info.data = inf, group_col = "group")
#' print(stat_new)
CreateStatObject <- function(
    raw.data = data.frame(),
    clean.data = data.frame(),
    info.data = data.frame(),
    group_col = "group",
    na.action = c("allow", "omit", "error"),
    ...
) {

  na.action <- match.arg(na.action)
  raw.data <- as.data.frame(raw.data)
  clean.data <- as.data.frame(clean.data)
  info.data <- as.data.frame(info.data)

  # Store original group_col before any column name modifications
  original_group_col <- group_col

  # Handle NA values based on na.action
  handle_na <- function(data, data_name) {
    if (is.null(data) || nrow(data) == 0) return(data)
    if (!any(is.na(data))) return(data)

    switch(na.action,
           "allow" = {
             message("Note: NA values found in ", data_name, " but kept (na.action='allow')")
             return(data)
           },
           "omit" = {
             warning("NA values found in ", data_name, ". Removing ", sum(is.na(data)), " NA values.")
             return(na.omit(data))
           },
           "error" = {
             stop("NA values found in ", data_name, ". Set na.action='omit' to remove or na.action='allow' to keep.")
           }
    )
  }

  raw.data <- handle_na(raw.data, "raw.data")
  clean.data <- handle_na(clean.data, "clean.data")

  # info.data is append-only, so it must NOT be na.omit()'d on its own: its row
  # mask would differ from the feature frames' and the setequal check below would
  # then abort. Row removal is expressed once, at the end, by re-indexing
  # info.data to the feature frame. na.action="error" must still see label NAs,
  # so that case is reported here.
  if (nrow(info.data) > 0 && any(is.na(info.data))) {
    switch(na.action,
           "allow" = message("Note: NA values found in info.data but kept (na.action='allow')"),
           "error" = stop("NA values found in info.data. Set na.action='omit' to remove or na.action='allow' to keep."),
           "omit"  = NULL)
  }

  clean_symbol_values <- function(data) {
    for(col in colnames(data)) {
      if(is.character(data[[col]])) {

        if(any(grepl("^[<>]", data[[col]]))) {

          data[[col]] <- gsub("[<>]", "", data[[col]])

          data[[col]] <- as.numeric(data[[col]])
          warning(paste("Removed >/< symbols from column", col, "and converted to numeric"))
        }
      }
    }
    return(data)
  }

  if (nrow(raw.data) == 0 && nrow(clean.data) == 0) {
    stop("At least one of 'raw.data' or 'clean.data' must be provided and not empty.")
  }


  if(nrow(raw.data) > 0) {
    raw.data <- clean_symbol_values(raw.data)
  }
  if(nrow(clean.data) > 0) {
    clean.data <- clean_symbol_values(clean.data)
  }

  # ---------------------------------------------------------------------------
  # Label location contract (checked BEFORE any column renaming, so that the
  # names reported in the error are the ones the user actually wrote).
  # ---------------------------------------------------------------------------
  if (!is.null(original_group_col)) {
    in_info  <- original_group_col %in% colnames(info.data)
    in_clean <- original_group_col %in% colnames(clean.data)
    in_raw   <- original_group_col %in% colnames(raw.data)
    in_matrix <- in_clean || in_raw

    if (in_matrix) {
      where <- c(if (in_raw) "raw.data", if (in_clean) "clean.data")
      stop(sprintf(
        paste0(
          "group_col '%s' is present in the feature matrix (%s), but icare requires ",
          "the label column to live in info.data.\n\n",
          "Why: a label inside the feature matrix is indistinguishable from a ",
          "feature to any code that selects columns by type. A numeric 0/1 label ",
          "is silently consumed by PlotPCA() and PlotCorrelationHeatmap(), which ",
          "inflates the apparent group separation without warning.\n\n",
          "Fix -- move it, don't copy it:\n",
          "  inf$%s <- <your group labels>\n",
          "  stat_obj <- CreateStatObject(raw.data = mat, info.data = inf,\n",
          "                               group_col = \"%s\", na.action = \"allow\")\n",
          "  # and delete the `mat$%s <- ...` line\n"),
        original_group_col, paste(where, collapse = " and "),
        original_group_col, original_group_col, original_group_col),
        call. = FALSE)
    }

    if (!in_info) {
      if (nrow(info.data) == 0) {
        stop(sprintf(
          paste0(
            "group_col '%s' was requested but info.data is empty.\n\n",
            "The label column must be supplied through info.data:\n",
            "  inf$%s <- <your group labels>          # rownames must match the features\n",
            "  stat_obj <- CreateStatObject(raw.data = mat, info.data = inf, group_col = \"%s\")"),
          original_group_col, original_group_col, original_group_col),
          call. = FALSE)
      }
      stop(sprintf("Column '%s' is not in info.data (columns present: %s). ",
                   original_group_col, paste(colnames(info.data), collapse = ", ")),
           "group_col is resolved against info.data, not against the feature matrix, ",
           "and is matched verbatim -- note that data.frame() and read.csv() rewrite ",
           "spaces and punctuation unless check.names = FALSE.",
           call. = FALSE)
    }
  }

  if (nrow(info.data) > 0) {
    if (nrow(raw.data) > 0 && !setequal(rownames(info.data), rownames(raw.data))) {
      stop("Row names in 'info.data' do not match 'raw.data' (content mismatch).")
    }

    if (nrow(clean.data) > 0 && !setequal(rownames(info.data), rownames(clean.data))) {
      stop("Row names in 'info.data' do not match 'clean.data' (content mismatch).")
    }
  }

  # raw.data and clean.data are filtered independently by handle_na() above, so
  # na.action="omit" can leave them describing different sample sets. Warn: every
  # downstream read aligns labels to one specific frame.
  if (nrow(raw.data) > 0 && nrow(clean.data) > 0 &&
      !identical(rownames(raw.data), rownames(clean.data))) {
    warning("'raw.data' and 'clean.data' cover different samples. ",
            "Labels are aligned against 'raw.data'; make the two frames agree ",
            "to keep sample-level results unambiguous.")
  }

  prepare_data <- function(data, data_name) {
    if (nrow(data) == 0) {
      return(data)
    }

    if (anyDuplicated(rownames(data))) {
      warning(paste("Duplicate row names found in", data_name, "; they have been made unique."))
      rownames(data) <- make.unique(rownames(data))
    }

    if (anyDuplicated(colnames(data))) {
      warning(paste("Duplicate column names found in", data_name, "; they have been made unique."))
      colnames(data) <- make.unique(colnames(data))
    }

    if (is.null(colnames(data))) {
      stop(paste(data_name, "is missing column names."))
    }

    data <- modify_column_names(data)

    cat(paste("Data prepared for", data_name, "\n"))

    return(data)
  }

  # Track the label by POSITION before renaming, so that group_col resolves to
  # the right column no matter what modify_column_names() does (spaces, illegal
  # characters, or duplicate names made unique). Resolving by re-deriving the
  # mangled string would silently drift the moment the renamer changes.
  group_idx <- if (!is.null(original_group_col) && original_group_col %in% colnames(info.data)) {
    which(colnames(info.data) == original_group_col)[1]
  } else {
    NA_integer_
  }

  raw.data <- prepare_data(raw.data, "raw.data")
  clean.data <- prepare_data(clean.data, "clean.data")
  # info.data goes through the SAME renamer as the feature frames. Previously it
  # was skipped, so group_col = "Diag Group" resolved against a renamed feature
  # matrix but not against an unrenamed info.data.
  info.data <- prepare_data(info.data, "info.data")

  if (!is.na(group_idx)) {
    resolved <- colnames(info.data)[group_idx]
    if (!identical(resolved, original_group_col)) {
      message("Note: group_col '", original_group_col, "' was renamed to '",
              resolved, "' in info.data")
    }
    group_col <- resolved
  }

  if (nrow(raw.data) > 0) {
    meta.featurename <- as.character(colnames(raw.data))
  } else if (nrow(clean.data) > 0) {
    meta.featurename <- as.character(colnames(clean.data))
  } else {
    meta.featurename <- character()
  }

  # Align info.data to the feature frame. This is the single place where a row
  # set is narrowed, so raw/clean/scale/info can never disagree by construction:
  # info.data keeps every row it has, and reads align to whatever frame they ask
  # for (see .stat_labels()).
  target <- if (nrow(clean.data) > 0) clean.data else if (nrow(raw.data) > 0) raw.data else NULL

  if (is.null(target)) {
    if (nrow(info.data) > 0) {
      warning("No feature data supplied; info.data was kept as-is.")
    }
  } else if (nrow(info.data) == 0) {
    # Only reachable when group_col is NULL: a label was never requested, so an
    # empty metadata frame is legitimate.
    info.data <- data.frame(row.names = rownames(target))
  } else {
    missing_rows <- setdiff(rownames(target), rownames(info.data))
    if (length(missing_rows) > 0) {
      stop(sprintf("info.data is missing %d row name(s) present in the feature frame ",
                   length(missing_rows)),
           sprintf("(e.g. %s). ", paste(head(missing_rows, 3), collapse = ", ")),
           "Both frames must describe the same samples.", call. = FALSE)
    }
    info.data <- info.data[rownames(target), , drop = FALSE]
    rownames(info.data) <- rownames(target)
  }

  cat("Final info.data prepared.\n")

  # The label must not survive in any feature frame. This is a construction-time
  # invariant, so assert it rather than repair it -- the earlier check already
  # rejected it, and this catches anything a future edit might reintroduce.
  if (!is.null(group_col)) {
    for (nm in c("raw.data", "clean.data")) {
      d <- if (nm == "raw.data") raw.data else clean.data
      if (nrow(d) > 0 && group_col %in% colnames(d)) {
        stop(sprintf("Internal consistency error: '%s' ended up in %s after preparation.",
                     group_col, nm), call. = FALSE)
      }
    }
  }

  Stat <- new(
    Class = 'Stat',
    raw.data = raw.data,
    info.data = info.data,
    clean.data = clean.data,
    meta.featurename = meta.featurename,
    process.info = list(),
    group_col = group_col
  )

  cat("Stat object created.\n")

  return(Stat)
}
