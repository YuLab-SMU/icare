#' Run Univariate Cox Analysis
#'
#' This function has been moved to module_4_pn5_0_cox_univariate.R.
#' Calling this function will result in an error or redirection.
#'
#' @keywords internal
run_univariate_cox_analysis <- function(data,
                                        time_col = "time",
                                        status_col = "status",
                                        selected_vars = NULL,
                                        covariates = NULL,
                                        P_value =0.05) {
  stop("This function is deprecated in module_3. Please use the version in module_4 which is now the canonical implementation.")
}

#' Sur Cox Univariate Analysis
#'
#' @param object SurObj.
#' @param formula Formula.
#' @param status_col Status col.
#' @param time_col Time col.
#' @param selected_vars Selected vars.
#' @param covariates Covariates.
#' @param P_value P value.
#' @param save_dir Save dir.
#' @export
Sur_cox_univariate_analysis <- function(object,
                                        formula = NULL,
                                        status_col = "status",
                                        time_col = "time",
                                        selected_vars = NULL,
                                        covariates = NULL,
                                        P_value =0.05,
                                        save_dir = here('SurObj', "univariate_analysis")) {
  
  if (inherits(object, 'SurObj')) {
    example_data <- slot(object, "survival.data")
    cat("Using survival data from SurObj object...\n")
  } else if (is.data.frame(object)) {
    example_data <- object
  } else {
    stop("Input must be an object of class 'SurObj' or a data frame.")
  }
  
  if (is.null(example_data) || nrow(example_data) == 0)
    stop("No valid data found in the input.")
  
  
  
  example_data[[status_col]] <- as.numeric(example_data[[status_col]])
  
  results <- run_univariate_cox_analysis(
    data = example_data,
    time_col = time_col,
    status_col = status_col,
    selected_vars = selected_vars,
    covariates = covariates,
    P_value = P_value
  )
  
  if (inherits(object, 'SurObj')) {
    cat("Updating 'SurObj' object...\n")
    
    object@univariate.analysis[["all_univariate_results"]] <- results
    
    cat("The 'SurObj' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }
  
  return(results)
}
