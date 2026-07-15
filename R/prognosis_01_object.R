#' PrognosiX Class
#'
#' @import methods
#' @slot clean.data Data frame.
#' @slot info.data Data frame.
#' @slot time_col Any.
#' @slot status_col Any.
#' @slot survival.data Data frame.
#' @slot baseline.table Any.
#' @slot variable.types List.
#' @slot survival.var List.
#' @slot sub.data Data frame.
#' @slot univariate.analysis List.
#' @slot split.data List.
#' @slot feature.result List.
#' @slot filtered.set List.
#' @slot survival.model Any.
#' @slot best.model List.
#' @slot subgroup.risk List.
#' @export
#' @examples
#' \dontrun{
#'   # Create a PrognosiX object
#'   obj <- new("PrognosiX",
#'              clean.data = data.frame(gene1 = c(1, 2, 3)),
#'              info.data = data.frame(time = c(10, 20, 30), status = c(1, 0, 1)),
#'              time_col = "time",
#'              status_col = "status")
#' }
PrognosiX <- setClass(
  Class = 'PrognosiX',
  slots = c(
    clean.data = 'data.frame',
    info.data = 'data.frame',
    time_col = 'ANY',
    status_col = 'ANY',
    survival.data ='data.frame',
    baseline.table = 'ANY',
    variable.types = 'list',
    survival.var = 'list',
    sub.data = 'data.frame',
    univariate.analysis= 'list',
    split.data = 'list',
    feature.result = 'list',
    filtered.set = 'list',
    survival.model = 'ANY',
    best.model = 'list',
    subgroup.risk = 'list'
  ),
  prototype = list(
    clean.data = data.frame(),
    info.data = data.frame(),
    time_col = NULL,
    status_col = NULL,
    survival.data = data.frame(),
    baseline.table = NULL,
    variable.types = list(),
    survival.var = list(),
    sub.data = data.frame(),
    univariate.analysis= list(),
    split.data = list(),
    filtered.set = list(),
    feature.result = list(),
    survival.model = NULL,
    best.model = list(),
    subgroup.risk = list()
  ),
  validity = function(object) {
    # Check data frame slots
    if (!is.data.frame(object@clean.data)) {
      return("Slot 'clean.data' must be a data.frame")
    }
    if (!is.data.frame(object@info.data)) {
      return("Slot 'info.data' must be a data.frame")
    }
    if (!is.data.frame(object@survival.data)) {
      return("Slot 'survival.data' must be a data.frame")
    }
    if (!is.data.frame(object@sub.data)) {
      return("Slot 'sub.data' must be a data.frame")
    }
    
    # Check character columns
    if (!is.null(object@time_col) && !is.character(object@time_col)) {
      return("Slot 'time_col' must be a character string or NULL")
    }
    if (!is.null(object@status_col) && !is.character(object@status_col)) {
      return("Slot 'status_col' must be a character string or NULL")
    }
    
    # Check list slots
    if (!is.list(object@variable.types)) {
      return("Slot 'variable.types' must be a list")
    }
    if (!is.list(object@survival.var)) {
      return("Slot 'survival.var' must be a list")
    }
    if (!is.list(object@univariate.analysis)) {
      return("Slot 'univariate.analysis' must be a list")
    }
    if (!is.list(object@split.data)) {
      return("Slot 'split.data' must be a list")
    }
    if (!is.list(object@feature.result)) {
      return("Slot 'feature.result' must be a list")
    }
    if (!is.list(object@filtered.set)) {
      return("Slot 'filtered.set' must be a list")
    }
    if (!is.list(object@best.model)) {
      return("Slot 'best.model' must be a list")
    }
    if (!is.list(object@subgroup.risk)) {
      return("Slot 'subgroup.risk' must be a list")
    }
    
    # Check survival data has required columns
    if (nrow(object@survival.data) > 0 && !is.null(object@time_col)) {
      if (!object@time_col %in% colnames(object@survival.data)) {
        return(sprintf("Time column '%s' not found in survival.data", object@time_col))
      }
    }
    if (nrow(object@survival.data) > 0 && !is.null(object@status_col)) {
      if (!object@status_col %in% colnames(object@survival.data)) {
        return(sprintf("Status column '%s' not found in survival.data", object@status_col))
      }
    }
    
    # Check row name consistency
    if (nrow(object@clean.data) > 0 && nrow(object@info.data) > 0) {
      if (!setequal(rownames(object@clean.data), rownames(object@info.data))) {
        return("Row names of 'clean.data' and 'info.data' should match")
      }
    }
    
    return(TRUE)
  }
)


`%||%` <- function(a,b) if(!is.null(a)&&length(a)>0) a else b

# Helper to save plots – renamed as requested
save_plot_sur <- function(p, name, w=8, h=6) {
  if(inherits(p, "ggplot")) ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, w, h)
}


#' Create PrognosiX Object
#'
#' @param clean.data Clean data.
#' @param info.data Info data.
#' @param time_col Time col.
#' @param status_col Status col.
#' @param baseline.table Baseline table.
#' @param variable.types Variable types.
#' @param survival.var Survival var.
#' @param survival.data Survival data.
#' @param sub.data Sub data.
#' @param univariate.analysis Univariate analysis.
#' @param split.data Split data.
#' @param feature.result Feature result.
#' @param filtered.set Filtered set.
#' @param survival.model Survival model.
#' @param best.model Best model.
#' @param subgroup.risk Subgroup risk.
#' @param object Input object.
#' @return A PrognosiX S4 object.
#' @export
#' 
#' @examples
#' \dontrun{
#'   # Create from data frames
#'   clean <- data.frame(gene1 = c(1, 2, 3), gene2 = c(4, 5, 6),
#'                       row.names = c("S1", "S2", "S3"))
#'   info <- data.frame(time = c(10, 20, 30), status = c(1, 0, 1),
#'                      row.names = c("S1", "S2", "S3"))
#'   obj <- CreatePrognosiXObject(clean.data = clean, info.data = info,
#'                                time_col = "time", status_col = "status")
#'
#'   # Create from existing Stat object
#'   # obj <- CreatePrognosiXObject(object = stat_obj, time_col = "time", status_col = "status")
#' }
CreatePrognosiXObject <- function(
    clean.data = NULL,
    info.data = data.frame(),
    time_col = "time",  
    status_col = "status",  
    baseline.table = NULL,
    variable.types = list(),
    survival.var = list(),
    survival.data = data.frame(),
    sub.data = data.frame(),
    univariate.analysis = list(),
    split.data = list(),
    feature.result = list(),
    filtered.set = list(),
    survival.model = NULL,
    best.model = list(),
    subgroup.risk = list(),
    object = NULL
) {
  
  if (is.null(clean.data) && is.null(object)) {
    stop("At least one of 'clean.data' or 'object' must be provided.")
  }
  if (!is.null(object)) {
    cat("Extracting data from provided object...\n")
    
    if (!inherits(object, "Stat") && !inherits(object, "PrognosiX") && !inherits(object, "Subtyping")) {
      stop("The 'object' parameter must be an instance of class 'Stat' or 'PrognosiX'.")
    }
    
    if (inherits(object, "Stat")) {
      clean.data <- ExtractCleanData(object)
      info.data <- ExtractInfoData(object)
      
      if (is.null(clean.data) || nrow(clean.data) == 0) {
        stop("Failed to extract valid clean data from the provided 'Stat' object.")
      }
      
      if (is.null(info.data) || nrow(info.data) == 0) {
        info.data <- data.frame(row.names = rownames(clean.data))
        cat("info.data was created from clean.data with", nrow(info.data), "rows.\n")
      }
      
    } else if (inherits(object, "PrognosiX")) {
      clean.data <- object@clean.data
      info.data <- object@info.data
      sub.data <- object@sub.data
      time_col <- object@time_col
      status_col <- object@status_col
    } else if (inherits(object, "Subtyping")) {
      if (!is.null(slot(object, "clustered.data")) && nrow(slot(object, "clustered.data")) > 0) {
        clean.data <- slot(object, "clustered.data")
        cat("Using 'clustered.data' from 'Subtyping' object as 'clean.data'.\n")
      } else {
        clean.data <- slot(object, "clean.data")
        cat("'clustered.data' not found. Using 'clean.data' from 'Subtyping' object.\n")
      }
      info.data <- slot(object, "info.data")
    }
  }
  
  if (is.null(clean.data) || nrow(clean.data) == 0) {
    stop("'clean.data' must be provided and not empty.")
  }
  if (nrow(info.data) == 0) {
    if (!is.null(time_col) && !is.null(status_col)) {
      info.data <- clean.data[, c(time_col, status_col), drop = FALSE]
      clean.data <- clean.data[, !colnames(clean.data) %in% c(time_col, status_col)]
      rownames(info.data) <- rownames(clean.data)
      colnames(info.data)[colnames(info.data) == time_col]<- "time"
      colnames(info.data)[colnames(info.data) == status_col]<- "status"
      cat("info.data created from clean.data with columns:", time_col, "and", status_col, "\n")
    } else {
      stop("Both 'time_col' and 'status_col' must be provided.")
    }
  }
  
  # Match info.data rows with clean.data rows with safety check
  if (nrow(info.data) > 0) {
    common_rows <- intersect(rownames(clean.data), rownames(info.data))
    missing_rows <- setdiff(rownames(clean.data), rownames(info.data))
    
    if (length(missing_rows) > 0) {
      warning("The following rows in clean.data are not found in info.data (", 
              length(missing_rows), " rows): ", 
              paste(head(missing_rows, 5), collapse = ", "),
              ifelse(length(missing_rows) > 5, " ...", ""))
    }
    
    if (length(common_rows) == 0) {
      warning("No matching rows found between clean.data and info.data. Creating empty info.data.")
      info.data <- data.frame(row.names = rownames(clean.data))
    } else {
      if (length(common_rows) < nrow(clean.data)) {
        warning("Only ", length(common_rows), " of ", nrow(clean.data), 
                " clean.data rows found in info.data. Subsetting clean.data to match.")
        clean.data <- clean.data[common_rows, , drop = FALSE]
      }
      info.data <- info.data[rownames(clean.data), , drop = FALSE]
    }
  }
  
  if (nrow(info.data) != 0) {
    if (!is.null(time_col) && !is.null(status_col)) {
      # If time_col/status_col are in clean.data but not info.data, move them
      if (time_col %in% colnames(clean.data) && !(time_col %in% colnames(info.data))) {
        info.data[[time_col]] <- clean.data[[time_col]]
        clean.data <- clean.data[, !colnames(clean.data) %in% time_col, drop=FALSE]
      }
      if (status_col %in% colnames(clean.data) && !(status_col %in% colnames(info.data))) {
        info.data[[status_col]] <- clean.data[[status_col]]
        clean.data <- clean.data[, !colnames(clean.data) %in% status_col, drop=FALSE]
      }
      colnames(info.data)[colnames(info.data) == time_col] <- "time"
      colnames(info.data)[colnames(info.data) == status_col] <- "status"
    } else {
      exists <- c("time", "status") %in% colnames(info.data)
      if (!all(exists)) {
        stop("Both 'time_col' and 'status_col' must be provided or 'time' and 'status' columns must exist.")
      }
    }
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
    
    return(data)
  }
  
  time_col <- "time"
  status_col <- "status"
  
  clean.data <- prepare_data(clean.data, "clean.data")
  info.data <- prepare_data(info.data, "info.data")
  
  
  if (!identical(rownames(clean.data), rownames(info.data))) {
    info.data <- info.data[rownames(clean.data), ]
    cat("info.data row names synchronized with clean.data.\n")
  }
  
  if (any(is.na(clean.data))) {
    stop("clean.data contains missing values. Please clean the data.")
  }
  
  if (status_col %in% colnames(info.data)) {
    info.data[[status_col]] <- as.numeric(as.character(info.data[[status_col]]))
  } else {
    warning("Column not found: ", status_col)
  }
  
  survival.data <- cbind(clean.data, info.data[, c(time_col, status_col)])
  cat("Removing rows with missing values in time or status columns...\n")  
  survival.data <- survival.data[complete.cases(survival.data[, c(time_col, status_col)]), ]
  
  PrognosiX <- new(
    Class = 'PrognosiX',
    clean.data = clean.data,
    info.data = info.data,
    sub.data = sub.data,
    univariate.analysis = univariate.analysis,
    time_col = time_col,
    status_col = status_col,
    survival.data = survival.data,
    baseline.table = baseline.table,
    variable.types = variable.types,
    survival.var = survival.var,
    split.data = split.data,
    feature.result = feature.result,
    filtered.set = filtered.set,
    survival.model = survival.model,
    best.model = best.model,
    subgroup.risk = subgroup.risk
  )
  
  cat("PrognosiX object created successfully.\n")
  
  return(PrognosiX)
}
