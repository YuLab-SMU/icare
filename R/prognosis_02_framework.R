#' Extract Survival Task from PrognosiX or TaskSurv Object
#'
#' This helper function extracts or returns a \code{TaskSurv} object from either
#' a \code{PrognosiX} S4 object or a \code{TaskSurv} object. It serves as a
#' unified interface for downstream functions that require a survival task.
#'
#' @param object An object of class \code{PrognosiX} or \code{TaskSurv}.
#'   If \code{PrognosiX}, the function extracts the survival data, time column,
#'   and event column to create a \code{TaskSurv} object.
#'
#' @return A \code{\link[mlr3proba]{TaskSurv}} object.
#'
#' @seealso \code{\link{surv_create_surv_task}} for creating tasks from data frames
#' @export
#' 
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' # From TaskSurv
#' data("veteran", package = "survival")
#' task <- TaskSurv$new("veteran", backend = veteran, time = "time", event = "status")
#' task_out <- surv_extract_task(task)
#' 
#' # From PrognosiX object (if available)
#' # task_out <- surv_extract_task(prog_obj)
#' }
surv_extract_task <- function(object) {
  .check_prognosis_packages()
  if (inherits(object, 'PrognosiX')) {
    data <- object@survival.data
    time_col <- object@time_col
    event_col <- object@status_col
    return(surv_create_surv_task(data, time_col, event_col))
  } else if (inherits(object, 'TaskSurv')) {
    return(object)
  }
  stop('Input must be a PrognosiX object or TaskSurv')
}

# ==============================================================================

#' Available Survival Learners
#'
#' A character vector of all survival learner IDs currently available in the
#' \code{mlr3} environment. This is generated at package load time by filtering
#' \code{mlr3::mlr_learners$keys()} for learners with the \code{"surv."} prefix.
#'
#' @format A character vector of learner IDs, e.g., \code{"surv.coxph"},
#'   \code{"surv.ranger"}, etc. Returns an empty vector if \code{mlr3} is not installed.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("mlr3", quietly = TRUE)) {
#'   head(surv_keys, 5)}
#' }
surv_keys <- if (requireNamespace("mlr3", quietly = TRUE)) {
  mlr3::mlr_learners$keys()[grep("^surv\\.", mlr3::mlr_learners$keys())]
} else {
  character(0)
}

#' Flexible Search Space Manager for Survival Analysis
#'
#' Returns a predefined hyperparameter search space for a given survival learner.
#' The search space includes sensible ranges for tuning, with dynamic scaling of
#' parameters like \code{mtry} based on the number of features in the task.
#'
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.ranger"}).
#' @param object Optional. A \code{TaskSurv} or \code{PrognosiX} object used to
#'   determine the number of features for dynamic parameter scaling. If \code{NULL},
#'   defaults to 100 features.
#'
#' @return A \code{\link[paradox]{ParamSet}} object defining the hyperparameter
#'   search space. If no predefined space exists for the learner, attempts to
#'   fetch from \code{mlr3tuningspaces}; if that fails, returns an empty \code{ParamSet}.
#'
#' @details
#' The function maintains predefined search spaces for over 30 survival learners,
#' organized into categories:
#' \itemize{
#'   \item \strong{Random Forests \& Ensemble Trees}: \code{surv.ranger}, \code{surv.rfsrc}, \code{surv.aorsf}, etc.
#'   \item \strong{Gradient Boosting}: \code{surv.xgboost.cox}, \code{surv.gbm}, \code{surv.mboost}, etc.
#'   \item \strong{Regularized Regression}: \code{surv.glmnet}, \code{surv.cv_glmnet}, \code{surv.penalized}, etc.
#'   \item \strong{Decision Trees \& SVM}: \code{surv.rpart}, \code{surv.ctree}, \code{surv.svm}, etc.
#'   \item \strong{Splines \& Flexible Models}: \code{surv.flexreg}, \code{surv.flexspline}, etc.
#'   \item \strong{Neural Networks}: \code{surv.survdnn}
#'   \item \strong{Non-parametric}: \code{surv.kaplan}, \code{surv.nelson} (empty ParamSet)
#' }
#'
#' @seealso \code{\link{surv_get_tuning_config}} for tuning strategy
#' @export
#' 
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' 
#' # Get search space for random forest
#' ps <- surv_get_search_space("surv.ranger")
#' print(ps)
#' 
#' # With task for dynamic scaling
#' data("veteran", package = "survival")
#' task <- TaskSurv$new("veteran", backend = veteran, time = "time", event = "status")
#' ps_dynamic <- surv_get_search_space("surv.rfsrc", object = task)
#' }
surv_get_search_space <- function(learner_id, object = NULL) {
  .check_prognosis_packages()
  task <- if (!is.null(object)) surv_extract_task(object) else NULL
  
  # Get the number of features for dynamic scaling of parameters like mtry
  p <- if (!is.null(task)) length(task$feature_names) else 100
  
  search_spaces <- list(
    
    # === 1. Random Forests & Ensemble Trees ===
    "surv.ranger"      = ps(num.trees = p_int(100, 1000), mtry.ratio = p_dbl(0.1, 0.8), min.node.size = p_int(1, 20), splitrule = p_fct(c("logrank", "extratrees"))),
    "surv.rfsrc"       = ps(ntree = p_int(100, 1000), mtry = p_int(1, max(2, floor(sqrt(p)*2))), nodesize = p_int(1, 20), splitrule = p_fct(c("logrank", "random"))),
    "surv.aorsf"       = ps(n_tree = p_int(100, 500), leaf_min_events = p_int(1, 10), split_min_events = p_int(5, 20)),
    "surv.blockforest" = ps(n_trees = p_int(100, 500), block.weights = p_fct(c("proportional", "equal"))),
    "surv.cforest"     = ps(ntree = p_int(100, 500), mtry = p_int(1, max(2, floor(sqrt(p)))), mincriterion = p_dbl(0.5, 0.99)),
    "surv.bart"        = ps(num_trees = p_int(20, 200), k = p_dbl(1, 3), power = p_dbl(1, 3)),
    
    # === 2. Gradient Boosting Machines (GBM) ===
    "surv.xgboost.cox" = ps(nrounds = p_int(50, 500), eta = p_dbl(1e-3, 0.3, logscale = TRUE), max_depth = p_int(2, 8), subsample = p_dbl(0.5, 1)),
    "surv.xgboost.aft" = ps(nrounds = p_int(50, 500), eta = p_dbl(1e-3, 0.3, logscale = TRUE), aft_loss_distribution = p_fct(c("normal", "logistic")), max_depth = p_int(2, 8)),
    "surv.gbm"         = ps(n.trees = p_int(100, 1000), interaction.depth = p_int(1, 5), shrinkage = p_dbl(1e-3, 0.1, logscale = TRUE), n.minobsinnode = p_int(2, 15)),
    "surv.mboost"      = ps(mstop = p_int(50, 500), nu = p_dbl(0.01, 0.2), baselearner = p_fct(c("bbs", "bols", "btree"))),
    "surv.blackboost"  = ps(mstop = p_int(50, 500), maxdepth = p_int(2, 8), nu = p_dbl(0.01, 0.2)),
    "surv.gamboost"    = ps(mstop = p_int(50, 500), nu = p_dbl(0.01, 0.2)),
    "surv.glmboost"    = ps(mstop = p_int(50, 500), nu = p_dbl(0.01, 0.2)),
    "surv.coxboost"    = ps(stepno = p_int(10, 200), penalty = p_dbl(1, 100, logscale = TRUE)),
    "surv.cv_coxboost" = ps(maxstepno = p_int(50, 200), penalty = p_dbl(1, 100, logscale = TRUE)),
    
    # === 3. Regularized & Penalized Regression ===
    "surv.glmnet"      = ps(alpha = p_dbl(0, 1), lambda = p_dbl(1e-4, 1, logscale = TRUE)),
    "surv.cv_glmnet"   = ps(alpha = p_dbl(0, 1)), 
    "surv.penalized"   = ps(lambda1 = p_dbl(0, 20), lambda2 = p_dbl(0, 20)),
    "surv.priority_lasso" = ps(block1.penalization = p_dbl(0, 1), lambda.type = p_fct(c("lambda.min", "lambda.1se"))),
    "surv.cv_ncvsurv"  = ps(penalty = p_fct(c("MCP", "SCAD", "lasso")), alpha = p_dbl(0.1, 1)),
    "surv.coxph"       = ps(ties = p_fct(c("efron", "breslow"))),
    
    # === 4. Decision Trees & Support Vector Machines ===
    "surv.rpart"       = ps(cp = p_dbl(1e-4, 0.1, logscale = TRUE), maxdepth = p_int(1, 30)),
    "surv.ctree"       = ps(mincriterion = p_dbl(0.5, 0.99), minsplit = p_int(2, 30), minbucket = p_int(1, 20)),
    "surv.svm"         = ps(type = p_fct(c("regression", "vanbelle1", "vanbelle2")), kernel = p_fct(c("lin_kernel", "rbf_kernel", "poly_kernel")), mu = p_dbl(0, 1)),
    
    # === 5. Splines & Flexible Parametric Models ===
    "surv.flexreg"     = ps(dist = p_fct(c("weibull", "gengamma", "genf", "gompertz"))),
    "surv.flexspline"  = ps(k = p_int(1, 10), scale = p_fct(c("hazard", "odds", "normal"))),
    "surv.gam.cox"     = ps(select = p_lgl()), 
    
    # === 6. Neural Networks & Baseline Estimators ===
    "surv.survdnn"     = ps(epochs = p_int(10, 100), lr = p_dbl(1e-4, 1e-2, logscale = TRUE), batch_size = p_int(16, 128)),
    "surv.kaplan"      = ps(), # Non-parametric (No Tuning)
    "surv.nelson"      = ps()  # Non-parametric (No Tuning)
  )
  
  # --- Selection Logic ---
  if (learner_id %in% names(search_spaces)) {
    return(search_spaces[[learner_id]])
  } else {
    # Attempt to fetch from mlr3tuningspaces (Expert Default) if not in the list
    t_space <- tryCatch({lts(learner_id)$values}, error = function(e) NULL)
    if (!is.null(t_space)) return(t_space)
    
    message(sprintf("[-] Info: No predefined space for '%s', returning empty ParamSet.", learner_id))
    return(ps())
  }
}

#' Get Recommended Tuning Configuration
#'
#' Provides a recommended tuning strategy (tuner and terminator) based on the
#' complexity of the specified survival learner. Complex models (e.g., random
#' forests, XGBoost) default to random search, while simpler models use grid search.
#'
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param tuning_budget An integer specifying the number of evaluations allowed
#'   during tuning. Default is \code{50}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{tuner}{A \code{\link[mlr3tuning]{Tuner}} object, or \code{NULL}
#'       if no tuning is needed.}
#'     \item{terminator}{A \code{\link[mlr3tuning]{Terminator}} object, or
#'       \code{NULL} if no tuning is needed.}
#'   }
#'
#' @details
#' The function distinguishes between complex learners (requiring more
#' sophisticated tuning strategies) and simple learners:
#' \itemize{
#'   \item \strong{Complex learners}: \code{surv.ranger}, \code{surv.xgboost.*},
#'     \code{surv.gbm}, \code{surv.cforest} → random search with budget.
#'   \item \strong{Simple learners}: For learners with parameters, grid search
#'     with adaptive resolution based on the number of parameters.
#'   \item \strong{No-tuning learners}: \code{surv.kaplan}, \code{surv.nelson}
#'     → returns \code{NULL} for both components.
#' }
#'
#' @examples
#' \dontrun{
#' # Get tuning config for random forest
#' config <- surv_get_tuning_config("surv.ranger", tuning_budget = 100)
#' print(config$tuner$id)
#' 
#' # For Cox PH (simple learner)
#' config_cox <- surv_get_tuning_config("surv.coxph")
#' print(config_cox$tuner$id)
#' }
#'
#' @seealso \code{\link{surv_get_search_space}} for available search spaces
#' @export
surv_get_tuning_config <- function(learner_id, tuning_budget = 50) {
  .check_prognosis_packages()
  # Select appropriate tuning strategy based on algorithm characteristics
  complex_learners <- c("surv.ranger", "surv.xgboost", "surv.gbm", "surv.cforest")
  
  if (learner_id %in% complex_learners) {
    # Complex models use random search or Bayesian optimization
    tuner <- tnr("random_search")
    terminator <- trm("evals", n_evals = tuning_budget)
  } else {
    # Simple models use grid search
    search_space <- surv_get_search_space(learner_id)
    n_params <- length(search_space$params)
    
    if (n_params == 0) {
      tuner <- NULL
      terminator <- NULL
    } else {
      resolution <- ceiling(tuning_budget^(1/n_params))
      tuner <- tnr("grid_search", resolution = resolution)
      terminator <- trm("evals", n_evals = tuning_budget)
    }
  }
  
  list(tuner = tuner, terminator = terminator)
}

# ==============================================================================
# 2. Data Processing Module
# ==============================================================================

#' Create a Survival Analysis Task
#'
#' Creates a \code{TaskSurv} object from a data frame for use with \code{mlr3proba}
#' survival analysis workflows. The function automatically coerces the data to
#' \code{data.table} for optimized performance.
#'
#' @param data A data frame containing the dataset with survival time and event
#'   status columns.
#' @param time_col A character string specifying the name of the column containing
#'   survival times.
#' @param event_col A character string specifying the name of the column containing
#'   event status (typically 1 for event, 0 for censored).
#' @param id A character string specifying the task identifier. Default is
#'   \code{"survival_task"}.
#'
#' @return A \code{\link[mlr3proba]{TaskSurv}} object ready for use with
#'   \code{mlr3} survival learners.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data("veteran", package = "survival")
#' 
#' task <- surv_create_surv_task(
#'   data = veteran,
#'   time_col = "time",
#'   event_col = "status",
#'   id = "veteran_task"
#' )
#' print(task)
#' }
#'
#' @seealso \code{\link{surv_extract_task}} for extracting tasks from other objects
#' @export
surv_create_surv_task <- function(data, time_col, event_col, id = "survival_task") {
  .check_prognosis_packages()
  # Coerce to data.table for optimized performance in mlr3
  data <- as.data.table(data)
  
  TaskSurv$new(
    id = id,
    backend = data,
    time = time_col,
    event = event_col
  )
}

# ==============================================================================
# 3. Model Training & Tuning Module
# ==============================================================================

#' Instantiate and Configure a Survival Learner
#'
#' Creates a survival learner instance with appropriate \code{predict_type} and
#' automatically handles categorical features by adding an encoding pipeline
#' if necessary.
#'
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param task A \code{TaskSurv} object used to determine feature types and
#'   learner capabilities.
#'
#' @return A configured \code{\link[mlr3]{Learner}} object, possibly wrapped in
#'   a \code{PipeOp} pipeline if encoding is required.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Instantiates the learner using \code{lrn(learner_id)}.
#'   \item Sets \code{predict_type} to \code{"distr"} if available, otherwise
#'     \code{"crank"}.
#'   \item Checks if the task contains factor/character features and if the
#'     learner supports them. If not, adds a \code{po("encode")} pipeline.
#' }
#'
#' @keywords internal
#' @noRd
surv_get_learner <- function(learner_id, task) {
  lrn_obj <- lrn(learner_id)
  
  if ("distr" %in% lrn_obj$predict_types) {
    lrn_obj$predict_type <- "distr"
  } else if ("crank" %in% lrn_obj$predict_types) {
    lrn_obj$predict_type <- "crank"
  }
  
  task_ftypes <- task$feature_types$type
  unsupported_factors <- ("factor" %in% task_ftypes || "character" %in% task_ftypes) && 
                         !("factor" %in% lrn_obj$feature_types)
  
  if (unsupported_factors) {
    lrn_obj <- po("encode", method = "treatment") %>>% lrn_obj
    lrn_obj <- as_learner(lrn_obj)
    
    if ("distr" %in% lrn_obj$predict_types) {
      lrn_obj$predict_type <- "distr"
    } else if ("crank" %in% lrn_obj$predict_types) {
      lrn_obj$predict_type <- "crank"
    }
  }
  return(lrn_obj)
}

#' Train and tune a survival learner with hyperparameter optimization
#'
#' @param object A survival task object (e.g., from \code{mlr3} package).
#' @param learner_id Character string identifying the learner (e.g., \code{"surv.coxph"}).
#' @param search_space A \code{\link[paradox]{ParamSet}} defining the hyperparameter
#'   search space. If \code{NULL}, a default is generated.
#' @param resampling A \code{\link[mlr3]{Resampling}} object. If \code{NULL},
#'   defaults to 5-fold cross-validation.
#' @param measure A \code{\link[mlr3]{Measure}} object. If \code{NULL},
#'   defaults to the concordance index (\code{surv.cindex}).
#' @param tuning_budget Integer number of evaluations allowed during tuning.
#' @param tuner A \code{\link[mlr3tuning]{Tuner}} object. If \code{NULL},
#'   random search is used.
#' @param seed Integer seed for reproducibility.
#'
#' @return A list with four components:
#' \describe{
#' \item{learner}{The trained learner object with optimal parameters.}
#' \item{best_params}{List of best hyperparameter values.}
#' \item{tuning_result}{Data frame of all evaluated parameter sets and performance.}
#' \item{cv_performance}{Numeric cross‑validated performance score.}
#' }
#' @details
#' The function extracts a task from \code{object}, instantiates the specified
#' learner, and – if a non‑empty search space is provided – tunes its
#' hyperparameters using the chosen resampling and measure. The best parameters
#' are then applied to the learner, which is retrained on the full task.
#' After training, predictions are validated to ensure the model produces
#' sensible \code{crank} values. If the search space is empty, training
#' proceeds directly without tuning.
#'
#' @examples
#' \dontrun{
#' library(mlr3)
#' library(mlr3proba)
#' task <- tsk("lung")
#' result <- surv_train_and_tune(task, "surv.coxph", tuning_budget = 10)
#' print(result$best_params)
#' }
#'
#' @seealso \code{\link{surv_extract_task}}, \code{\link{surv_get_learner}}
#' @export
surv_train_and_tune <- function(object,
                           learner_id,
                           search_space = NULL,
                           resampling = NULL,
                           measure = NULL,
                           tuning_budget = 50,
                           tuner = NULL,
                           seed = 123) {
  task <- surv_extract_task(object)
  
  set.seed(seed)
  
  # 1. Instantiate the learner and smartly configure predict_type
  learner <- tryCatch({
    surv_get_learner(learner_id, task)
  }, error = function(e) {
    stop(sprintf("Failed to instantiate learner '%s': %s\nConsider running mlr3extralearners::install_learner('%s')", learner_id, e$message, learner_id))
  })
  
  # 2. Set Search Space (Dynamically pass task to scale features)
  if (is.null(search_space)) {
    search_space <- surv_get_search_space(learner_id, object = task)
  }
  
  # Skip tuning if there are no hyperparameters defined
  if (length(search_space$params) == 0) {
    message(sprintf("[-] Learner '%s' requires no tuning. Training directly...", learner_id))
    learner$train(task)
    return(list(
      learner = learner,
      best_params = list(),
      tuning_result = NULL,
      cv_performance = NA
    ))
  }
  
  # 3. Configure default evaluation strategies
  if (is.null(resampling)) resampling <- rsmp("cv", folds = 5)
  if (is.null(measure)) measure <- msr("surv.cindex")
  if (is.null(tuner)) tuner <- tnr("random_search") 
  
  terminator <- trm("evals", n_evals = tuning_budget)
  
  # 4. Create Tuning Instance
  instance <- TuningInstanceSingleCrit$new(
    task = task,
    learner = learner,
    resampling = resampling,
    measure = measure,
    search_space = search_space,
    terminator = terminator
  )
  
  # 5. Execute Tuning
  message(sprintf("[*] Starting tuning for learner '%s' (Budget: %d evals)...", learner_id, tuning_budget))
  tuner$optimize(instance)
  
  # 6. Train the final model on the full dataset using the best hyperparameters
  learner$param_set$values <- instance$result_learner_param_vals
  
  # Try training with error handling
  train_success <- tryCatch({
    learner$train(task)
    TRUE
  }, error = function(e) {
    message(sprintf("[!] Training failed with error: %s", e$message))
    FALSE
  })
  
  if (!train_success) {
    stop(sprintf("Model training failed for learner '%s'. Please check data quality and feature suitability.", learner_id))
  }
  
  # Verify predictions work (mlr3 learners may not have state$trained flag)
  test_pred <- tryCatch({
    pred <- learner$predict(task)
    if (is.null(pred$crank)) {
      stop("Model predictions do not contain crank values")
    }
    if (length(unique(pred$crank)) < 2) {
      warning(sprintf("Model predictions have low variance (unique values: %d). Risk stratification may not be meaningful.", length(unique(pred$crank))))
    }
    TRUE
  }, error = function(e) {
    message(sprintf("[!] Prediction test failed: %s", e$message))
    FALSE
  })
  
  if (!test_pred) {
    stop(sprintf("Model training verification failed for learner '%s'. Predictions are invalid.", learner_id))
  }
  
  # 7. Return Results
  list(
    learner = learner,
    best_params = instance$result_learner_param_vals,
    tuning_result = instance$result,
    cv_performance = instance$result_y[[1]] # Extract single-criterion score
  )
}

# ==============================================================================
# 4. Model Evaluation Module
# ==============================================================================

#' Evaluate Model Performance on a Survival Task
#'
#' Evaluates a trained survival learner on a given task using specified measures.
#' The function automatically selects appropriate measures based on the learner's
#' \code{predict_type}.
#'
#' @param learner A trained \code{\link[mlr3]{Learner}} object.
#' @param object A \code{TaskSurv} or \code{PrognosiX} object. The function
#'   extracts the task using \code{surv_extract_task()}.
#' @param measures A list of \code{\link[mlr3]{Measure}} objects. If \code{NULL},
#'   automatically selects:
#'   \itemize{
#'     \item \code{surv.cindex} (always included)
#'     \item \code{surv.graf} if learner supports \code{"distr"} predictions
#'   }
#'
#' @return A data frame with one row containing the performance metrics for each
#'   requested measure.
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' # Auto-select measures
#' perf <- surv_evaluate_model(learner, task)
#' print(perf)
#' 
#' # Custom measures
#' perf_custom <- surv_evaluate_model(
#'   learner, task,
#'   measures = list(msr("surv.cindex"), msr("surv.logloss"))
#' )
#' }
#'
#' @seealso \code{\link{surv_train_and_tune}}, \code{\link{surv_benchmark_learners}}
#' @export
surv_evaluate_model <- function(learner, object, measures = NULL) {
  task <- surv_extract_task(object)
  
  # Smart measure selection based on learner capabilities
  if (is.null(measures)) {
    if (learner$predict_type == "distr") {
      measures <- list(msr("surv.cindex"), msr("surv.graf")) # C-index & Brier Score
    } else {
      measures <- list(msr("surv.cindex")) # Fallback to C-index only
    }
  }
  
  # Generate predictions on the task (Note: usually done on test data)
  predictions <- learner$predict(task)
  
  # Calculate scores safely
  scores <- sapply(measures, function(m) {
    tryCatch({
      predictions$score(m)
    }, error = function(e) {
      NA_real_
    })
  })
  
  # Format results
  names(scores) <- sapply(measures, function(m) m$id)
  as.data.frame(t(scores))
}

# ==============================================================================
# 5. Batch Benchmark Module
# ==============================================================================

#' Batch Train and Benchmark Multiple Learners with CV Performance
#'
#' Trains and evaluates multiple survival learners using cross-validation.
#' Optionally performs hyperparameter tuning for each learner before benchmarking.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param learner_ids A character vector of learner IDs (e.g., \code{c("surv.coxph", "surv.ranger")}).
#' @param tune A logical value. Should hyperparameter tuning be performed?
#'   Default is \code{TRUE}.
#' @param resampling A \code{\link[mlr3]{Resampling}} strategy. If \code{NULL},
#'   defaults to 5-fold cross-validation.
#' @param measures A list of evaluation measures. If \code{NULL}, uses
#'   \code{surv.cindex} for CV and training evaluation.
#' @param tuning_budget An integer specifying the number of tuning evaluations
#'   when \code{tune = TRUE}. Default is \code{50}.
#'
#' @return A list where each element corresponds to a learner and contains:
#'   \describe{
#'     \item{learner}{The trained learner object.}
#'     \item{best_params}{A list of best hyperparameter values (if tuning was performed).}
#'     \item{cv_performance}{A numeric cross-validated C-index.}
#'     \item{performance}{A data frame of training set metrics.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' veteran$celltype <- as.factor(veteran$celltype)
#' task <- surv_create_surv_task(veteran, "time", "status", "veteran_task")
#' 
#' # Without tuning (fast)
#' results <- surv_benchmark_learners(
#'   object = task,
#'   learner_ids = c("surv.coxph", "surv.ranger"),
#'   tune = FALSE
#' )
#' 
#' # With tuning (slower)
#' \dontrun{
#' results_tuned <- surv_benchmark_learners(
#'   object = task,
#'   learner_ids = c("surv.coxph", "surv.ranger"),
#'   tune = TRUE,
#'   tuning_budget = 20
#' )
#' }
#' }
#'
#' @seealso \code{\link{surv_summarize_benchmark}} for summarizing results
#' @export
surv_benchmark_learners <- function(object,
                                    learner_ids,
                                    tune = TRUE,
                                    resampling = NULL,
                                    measures = NULL,
                                    tuning_budget = 50) {
  task <- surv_extract_task(object)
  
  # Default resampling: 5‑fold CV
  if (is.null(resampling)) {
    resampling <- rsmp("cv", folds = 5)
  }
  
  # Default measure for tuning and CV
  if (is.null(measures)) {
    measures <- list(msr("surv.cindex"))
  }
  cv_measure <- measures[[1]]   # Use first measure for CV
  
  results <- list()
  
  for (learner_id in learner_ids) {
    message(sprintf("\n========== Processing Learner: %s ==========", learner_id))
    
    learner_result <- tryCatch({
      if (tune) {
        # Tuning workflow (includes CV performance from tuning)
        tune_res <- surv_train_and_tune(
          object = task,
          learner_id = learner_id,
          resampling = resampling,
          measure = cv_measure,
          tuning_budget = tuning_budget
        )
        learner <- tune_res$learner
        best_params <- tune_res$best_params
        cv_perf <- tune_res$cv_performance
      } else {
        # No tuning: use a FRESH clone for CV (the trained learner must not be
        # reused for CV, as it has already seen all data -- that would invalidate
        # the CV estimate).
        learner_for_cv <- surv_get_learner(learner_id, task)
        rr <- resample(task, learner_for_cv, resampling, store_models = FALSE)
        cv_perf <- rr$aggregate(cv_measure)
        # Now train the final model on the full dataset for downstream use
        learner <- surv_get_learner(learner_id, task)
        learner$train(task)
        best_params <- list()
      }
      
      # Training set (apparent) performance
      train_perf <- surv_evaluate_model(learner, task, measures)
      
      list(
        learner = learner,
        best_params = best_params,
        cv_performance = cv_perf,
        performance = train_perf
      )
    }, error = function(e) {
      message(sprintf("✗ [%s] Failed: %s", learner_id, e$message))
      NULL
    })
    
    if (!is.null(learner_result)) {
      results[[learner_id]] <- learner_result
      message(sprintf("✓ [%s] Completed successfully (CV C‑index = %.4f)", 
                      learner_id, learner_result$cv_performance))
    }
  }
  
  return(results)
}

#' Summarize Benchmark Results into a Leaderboard
#'
#' Converts the output from \code{surv_benchmark_learners} into a sorted
#' leaderboard data frame for easy comparison of model performance.
#'
#' @param benchmark_results The list output from \code{\link{surv_benchmark_learners}}.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{learner}{Learner ID.}
#'     \item{cv_score}{Cross-validated C-index score.}
#'     \item{...}{Additional performance metrics from training set evaluation.}
#'   }
#'   The data frame is sorted by \code{cv_score} in descending order.
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' 
#' # Run benchmark (without tuning for speed)
#' results <- surv_benchmark_learners(
#'   object = task,
#'   learner_ids = c("surv.coxph", "surv.ranger"),
#'   tune = FALSE
#' )
#' 
#' # Summarize
#' summary_df <- surv_summarize_benchmark(results)
#' print(summary_df)
#' }
#'
#' @seealso \code{\link{surv_benchmark_learners}}
#' @export
surv_summarize_benchmark <- function(benchmark_results) {
  perf_list <- lapply(names(benchmark_results), function(learner_id) {
    res <- benchmark_results[[learner_id]]
    if (!is.null(res)) {
      # Extract cv_score (if missing, use training C‑index as fallback)
      cv_score <- if (!is.null(res$cv_performance) && !is.na(res$cv_performance)) {
        res$cv_performance
      } else if (!is.null(res$performance$surv.cindex)) {
        warning(sprintf("No CV score for %s, using training C‑index as fallback.", learner_id))
        res$performance$surv.cindex
      } else {
        NA_real_
      }
      df <- data.frame(
        learner = learner_id,
        cv_score = cv_score,
        stringsAsFactors = FALSE
      )
      cbind(df, res$performance)
    }
  })
  perf_df <- data.table::rbindlist(perf_list, fill = TRUE)
  if (nrow(perf_df) > 0 && "cv_score" %in% colnames(perf_df)) {
    perf_df <- perf_df[order(-perf_df$cv_score), ]
  }
  return(as.data.frame(perf_df))
}

# ==============================================================================
# 6. Utility Functions Module
# ==============================================================================

#' List Available Survival Learners
#'
#' Returns a character vector of all survival analysis learner IDs currently
#' available in the \code{mlr3} environment.
#'
#' @return A character vector of learner IDs with the \code{"surv."} prefix.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("mlr3", quietly = TRUE)) {
#'   available <- surv_list_available_learners()
#'   print(head(available, 10))
#' }
#' }
#'
#' @seealso \code{\link{surv_keys}} for the global list of survival learners
#' @export
surv_list_available_learners <- function() {
  surv_learners <- mlr_learners$keys()[grep("^surv\\.", mlr_learners$keys())]
  return(surv_learners)
}


# ==============================================================================
# 7. Advanced Visualization & Interpretability Module
# ==============================================================================
#' Plot Risk Stratification Kaplan-Meier Curve (Training or Validation)
#'
#' Generates Kaplan-Meier survival curves stratified by risk groups defined by
#' a model's predicted risk scores (crank values). Supports multiple cutoff
#' determination methods including median, tertile, quartile, and p-value optimization.
#'
#' @param learner A trained \code{mlr3} learner that outputs \code{crank} predictions.
#' @param object A \code{TaskSurv} or \code{PrognosiX} object (can be training or validation task).
#' @param cutoff_method A character string specifying the method for determining
#'   risk group cutoffs. Must be one of \code{"median"}, \code{"tertile"},
#'   \code{"quartile"}, \code{"p_optimize"}, or \code{"custom"}.
#' @param custom_cutoffs A numeric vector of custom cutoffs (required when
#'   \code{cutoff_method = "custom"}). Length 1 → binary split, length 2 → three
#'   groups, length 3 → four groups, etc.
#' @param n_boot An integer specifying the number of bootstrap samples for
#'   \code{p_optimize} method. Default is \code{10}.
#' @param fraction A numeric value specifying the subsample fraction for
#'   \code{p_optimize}. Default is \code{0.1}.
#' @param conf_int A logical value. Should confidence intervals be shown?
#'   Default is \code{FALSE}.
#' @param risk_table A logical value. Should a risk table be shown below the plot?
#'   Default is \code{FALSE}.
#' @param palette_name A character string specifying the Wes Anderson palette name.
#'   Default is \code{"AsteroidCity1"}.
#' @param show_cutoff A logical value. Should cutoffs be displayed in the plot subtitle?
#'   Default is \code{TRUE}.
#' @param title An optional custom plot title. If \code{NULL}, auto-generated.
#'
#' @return A \code{\link[survminer]{ggsurvplot}} object with the Kaplan-Meier plot.
#'   The cutoffs used and group distribution are stored as attributes.
#'
#' @details
#' The cutoff methods work as follows:
#' \itemize{
#'   \item \code{"median"}: Splits at the median risk score (binary groups).
#'   \item \code{"tertile"}: Splits at the 1/3 and 2/3 quantiles (3 groups).
#'   \item \code{"quartile"}: Splits at the 1/4, 1/2, and 3/4 quantiles (4 groups).
#'   \item \code{"p_optimize"}: Finds the cutoff that maximizes log-rank test
#'     significance across bootstrap samples (binary groups).
#'   \item \code{"custom"}: Uses user-provided cutoffs for arbitrary group numbers.
#' }
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' # Median split
#' p <- surv_plot_risk_km(learner, task, cutoff_method = "median")
#' print(p)
#' 
#' # Tertile split (3 groups)
#' p_tertile <- surv_plot_risk_km(learner, task, cutoff_method = "tertile")
#' 
#' # Custom cutoffs
#' p_custom <- surv_plot_risk_km(
#'   learner, task,
#'   cutoff_method = "custom",
#'   custom_cutoffs = c(-1, 0, 1)
#' )
#' }
#'
#' @seealso \code{\link{get_cf}} for extracting cutoffs from the plot object
#' @export
surv_plot_risk_km <- function(learner, object, 
                              cutoff_method = c("median", "tertile", "quartile", "p_optimize"),
                              custom_cutoffs = NULL,
                              n_boot = 10, fraction = 0.1,
                              conf_int = FALSE,
                              risk_table = FALSE,
                              palette_name = "AsteroidCity1",
                              show_cutoff = TRUE,
                              title = NULL) {
  cutoff_method <- match.arg(cutoff_method)
  task <- surv_extract_task(object)
  if (!requireNamespace("survminer", quietly = TRUE)) 
    stop("Please install 'survminer'")
  
  # Predict risk scores (crank)
  if ("crank" %in% learner$predict_types) learner$predict_type <- "crank"
  lp <- learner$predict(task)$crank
  
  if (length(unique(lp)) < 2) 
    stop("Cannot perform risk stratification: all predictions are identical.")
  
  surv_data <- as.data.frame(task$data(cols = task$target_names))
  surv_data$lp <- lp
  
  # Determine cutoffs
  if (cutoff_method == "custom") {
    if (is.null(custom_cutoffs)) 
      stop("cutoff_method = 'custom' requires custom_cutoffs.")
    cutoffs_used <- sort(custom_cutoffs)
    cat(sprintf("[*] Using custom cutoffs: %s\n", 
                paste(round(cutoffs_used, 4), collapse = ", ")))
  } else {
    if (cutoff_method == "median") {
      cutoffs_used <- median(lp, na.rm = TRUE)
    } else if (cutoff_method == "tertile") {
      q <- quantile(lp, probs = c(1/3, 2/3), na.rm = TRUE)
      if (length(unique(q)) < 2) {
        warning("Tertile cutoffs not unique, falling back to median.")
        cutoffs_used <- median(lp)
        cutoff_method <- "median"
      } else {
        cutoffs_used <- as.numeric(q)
      }
    } else if (cutoff_method == "quartile") {
      q <- quantile(lp, probs = c(1/4, 2/4, 3/4), na.rm = TRUE)
      if (length(unique(q)) < 3) {
        warning("Quartile cutoffs not unique, falling back to median.")
        cutoffs_used <- median(lp)
        cutoff_method <- "median"
      } else {
        cutoffs_used <- as.numeric(q)
      }
    } else if (cutoff_method == "p_optimize") {
      # Only meaningful on training set; if used on validation, it will compute based on validation data
      cat(sprintf("[*] Running p_optimize: %d bootstraps with %.1f%% fraction...\n", 
                  n_boot, fraction * 100))
      n_samples <- nrow(surv_data)
      sample_size <- max(10, floor(n_samples * fraction))
      best_cutoffs <- numeric(n_boot)
      for (i in 1:n_boot) {
        idx <- sample(1:n_samples, size = sample_size, replace = TRUE)
        boot_data <- surv_data[idx, ]
        boot_lp <- sort(unique(boot_data$lp))
        if (length(boot_lp) > 2) boot_lp <- boot_lp[-c(1, length(boot_lp))]
        if (length(boot_lp) < 2) {
          best_cutoffs[i] <- median(boot_data$lp, na.rm = TRUE)
          next
        }
        p_vals <- numeric(length(boot_lp))
        for (j in seq_along(boot_lp)) {
          cutoff <- boot_lp[j]
          group <- ifelse(boot_data$lp > cutoff, "High", "Low")
          if (length(unique(group)) == 2 && min(table(group)) >= 3) {
            fit_diff <- tryCatch(
              survival::survdiff(survival::Surv(time, status) ~ group, data = boot_data),
              error = function(e) NULL
            )
            if (!is.null(fit_diff) && !is.na(fit_diff$chisq)) {
              p_vals[j] <- 1 - pchisq(fit_diff$chisq, length(fit_diff$n) - 1)
            } else {
              p_vals[j] <- NA
            }
          } else {
            p_vals[j] <- NA
          }
        }
        if (all(is.na(p_vals))) {
          best_cutoffs[i] <- median(boot_data$lp, na.rm = TRUE)
        } else {
          best_cutoffs[i] <- boot_lp[which.min(p_vals)]
        }
      }
      cutoffs_used <- median(best_cutoffs, na.rm = TRUE)
      cat(sprintf("  -> Optimal threshold: %.4f\n", cutoffs_used))
    }
  }
  
  # Assign risk groups (unified logic)
  risk_group <- .assign_risk_groups(lp, cutoffs_used)
  surv_data$risk_group <- risk_group
  surv_data$lp <- NULL
  
  # Fit KM curves
  fit <- survival::survfit(survival::Surv(time, status) ~ risk_group, data = surv_data)
  actual_levels <- levels(droplevels(risk_group))
  n_groups <- length(actual_levels)
  
  # Colors
  if (requireNamespace("wesanderson", quietly = TRUE)) {
    color_palette <- wesanderson::wes_palette(palette_name, n = max(n_groups, 3), 
                                              type = "continuous")[1:n_groups]
  } else {
    color_palette <- RColorBrewer::brewer.pal(min(n_groups, 8), "Set1")
  }
  
  # Title and subtitle
  if (is.null(title)) {
    title <- paste("Risk Stratification:", learner$id)
  }
  subtitle_text <- paste("Method:", cutoff_method)
  if (show_cutoff && !is.null(cutoffs_used)) {
    cutoff_disp <- paste(round(cutoffs_used, 3), collapse = ", ")
    subtitle_text <- paste0(subtitle_text, " | Cutoffs: ", cutoff_disp)
  }
  
  p <- survminer::ggsurvplot(
    fit, data = surv_data,
    pval = (n_groups > 1), pval.method = (n_groups > 1),
    conf.int = conf_int, risk.table = risk_table,
    title = title, subtitle = subtitle_text,
    palette = color_palette,
    ggtheme = ggprism::theme_prism() + 
      ggplot2::theme(legend.title = ggplot2::element_blank()),
    legend = "right",
    risk.table.y.text.col = TRUE, risk.table.y.text = TRUE,
    risk.table.fontsize = 3.5
  )
  attr(p, "cutoffs_used") <- cutoffs_used
  attr(p, "group_distribution") <- table(risk_group)
  return(p)
}

# Internal helper for risk group assignment (used by above)
.assign_risk_groups <- function(lp, cutoffs) {
  cutoffs <- sort(cutoffs)
  n_cuts <- length(cutoffs)
  if (n_cuts == 1) {
    factor(ifelse(lp > cutoffs, "High Risk", "Low Risk"),
           levels = c("Low Risk", "High Risk"))
  } else if (n_cuts == 2) {
    cut(lp, breaks = c(-Inf, cutoffs[1], cutoffs[2], Inf),
        labels = c("Low Risk", "Medium Risk", "High Risk"))
  } else {
    labels <- c("Low Risk", paste("Risk Group", 2:(n_cuts)), "High Risk")
    if (length(labels) != n_cuts + 1) 
      labels <- paste("Group", 1:(n_cuts+1))
    cut(lp, breaks = c(-Inf, cutoffs, Inf), labels = labels)
  }
}


#' Extract Cutoffs from a Risk Stratification Plot
#'
#' Retrieves the cutoffs used to create risk groups in a Kaplan-Meier plot
#' generated by \code{surv_plot_risk_km}.
#'
#' @param km_plot A \code{ggsurvplot} object returned by \code{\link{surv_plot_risk_km}}.
#'
#' @return A numeric vector of cutoffs used in the plot, or \code{NULL} if not available.
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' p <- surv_plot_risk_km(learner, task, cutoff_method = "median")
#' cutoffs <- get_cf(p)
#' print(cutoffs)
#' }
#'
#' @seealso \code{\link{surv_plot_risk_km}}
#' @export
get_cf <- function(km_plot) {
  attr(km_plot, "cutoffs_used")
}


#' Generate Clinical Nomogram for Survival Model
#'
#' Creates a nomogram for visualizing and predicting survival probabilities
#' at specified time points using a fitted Cox proportional hazards model.
#' The nomogram is generated using the \code{rms} package.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param selected_features A character vector of feature names to include in the
#'   nomogram. Defaults to the top 5 features from the task.
#' @param time_points A numeric vector of time points at which to predict survival
#'   probabilities (e.g., \code{c(3, 5)}).
#' @param time_unit An optional character string specifying the time unit label
#'   (e.g., \code{"days"}, \code{"months"}, \code{"years"}). If \code{NULL},
#'   no unit is shown. Default is \code{NULL}.
#'
#' @return A \code{nomogram} object (invisibly) and prints the nomogram plot.
#'   The object can be used for further customization.
#'
#' @note
#' This function uses \code{rms::datadist()} for model fitting. The \code{datadist}
#' is set locally to avoid polluting the global environment.
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' 
#' # Generate nomogram with top 5 features
#' nom <- surv_generate_nomogram(
#'   object = task,
#'   selected_features = c("age", "karno", "diagtime", "celltype", "prior"),
#'   time_points = c(3, 5),
#'   time_unit = "months"
#' )
#' }
#'
#' @seealso \code{\link[rms]{nomogram}}, \code{\link[rms]{cph}}
#' @export
surv_generate_nomogram <- function(object, selected_features = NULL, time_points = c(3, 5), time_unit = NULL) {
  task <- surv_extract_task(object)
  if (!requireNamespace("rms", quietly = TRUE)) stop("Please install 'rms'")
  
  data <- as.data.frame(task$data())
  
  # Default to top 5 features if none provided
  if (is.null(selected_features)) {
    selected_features <- head(task$feature_names, 5)
  }
  
  # rms::datadist is strictly required for rms::cph nomograms
  dd <- rms::datadist(data)
  options(datadist = "dd")
  
  # Refit the model using rms::cph (mlr3 uses survival::coxph internally)
  formula_str <- paste("survival::Surv(time, status) ~", paste(selected_features, collapse = " + "))
  cph_fit <- rms::cph(as.formula(formula_str), data = data, x = TRUE, y = TRUE, surv = TRUE)
  
  # Create survival prediction functions for specific time points
  surv_obj <- rms::Survival(cph_fit)
  surv_funcs <- lapply(time_points, function(t) {
    function(x) surv_obj(t, x)
  })
  
  # Build funlabel based on time_unit
  if (is.null(time_unit)) {
    funlabel <- paste0(time_points, " Survival")
  } else {
    # Capitalize first letter of unit
    unit_label <- paste0(toupper(substring(time_unit, 1, 1)), substring(time_unit, 2))
    funlabel <- paste0(time_points, "-", unit_label, " Survival")
  }
  
  # Build nomogram
  nom <- rms::nomogram(
    cph_fit, 
    fun = surv_funcs,
    funlabel = funlabel,
    lp = FALSE # Hide linear predictor to save space
  )
  
  plot(nom)
  return(nom)
}

#' SurvSHAP(t) Explanations for Survival Models — Production Version
#'
#' Computes time-dependent SurvSHAP(t) explanations for survival predictions.
#' Supports both local (individual patient) and global (population-level)
#' interpretations using Kernel SHAP.
#'
#' @param learner A trained \code{mlr3} \code{LearnerSurv} object.
#' @param task A \code{TaskSurv} object.
#' @param type A character string specifying the explanation type. Must be one of
#'   \code{"local"} (explain specific patients) or \code{"global"} (population-level
#'   importance averaged over many patients).
#' @param n_explain An integer specifying the number of observations to explain.
#'   For \code{type = "local"}: number of patients to explain. For
#'   \code{type = "global"}: number of observations to aggregate over.
#'   Default is \code{NULL} (auto-set to 20 for local, 50 for global).
#' @param n_background An integer specifying the background data size for Kernel SHAP.
#'   Default is \code{50L}.
#' @param n_timepoints An integer specifying how many evaluation time points to use.
#'   If \code{NULL}, uses all available time points. Default is \code{NULL}.
#' @param n_top_features An integer specifying the top features to show in plots.
#'   Default is \code{6L}.
#' @param bar_color A character string specifying the hex color for bar plots.
#'   Default is \code{"#2980b9"}.
#' @param seed An integer seed for reproducibility. Default is \code{123L}.
#' @param verbose A logical value. Print progress messages. Default is \code{TRUE}.
#'
#' @return A list with six components:
#'   \describe{
#'     \item{shap_long}{A tidy data frame with columns \code{feature}, \code{time},
#'       \code{shap_value}, and \code{observation}.}
#'     \item{explainer}{A \code{survex} explainer object.}
#'     \item{eval_times}{The evaluation time points actually used.}
#'     \item{plots}{A list with \code{bar_plot} (feature importance) and
#'       \code{line_plot} (SHAP dynamics over time).}
#'     \item{original_features}{A data frame with the original feature values.}
#'   }
#'
#' @details
#' The function uses the \code{survex} package for efficient SHAP computation.
#' Key features:
#' \itemize{
#'   \item \strong{Time-dependent}: SHAP values are computed at multiple time points.
#'   \item \strong{Sampling acceleration}: Sub-samples time points and observations for speed.
#'   \item \strong{Native mlr3 support}: Handles factor/character columns without errors.
#' }
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' # Global explanation (population-level)
#' shap_result <- surv_explain_shap(
#'   learner = learner,
#'   task = task,
#'   type = "global",
#'   n_explain = 20,
#'   n_background = 20,
#'   n_timepoints = 5
#' )
#' 
#' # View bar plot
#' print(shap_result$plots$bar_plot)
#' 
#' # Local explanation (single patient)
#' shap_local <- surv_explain_shap(
#'   learner = learner,
#'   task = task,
#'   type = "local",
#'   n_explain = 1,
#'   n_background = 20
#' )
#' }
#'
#' @seealso \code{\link{surv_plot_shap_beeswarm}} for visualizing SHAP values
#' @export
surv_explain_shap <- function(
    learner,
    task,
    type               = c("local", "global"),
    n_explain          = NULL,
    n_background       = 50L,
    n_timepoints       = NULL,
    n_top_features     = 6L,
    bar_color          = "#2980b9",
    seed               = 123L,
    verbose            = TRUE) {
  
  type <- match.arg(type)
  set.seed(seed)
  
  # ---- Input validation ----
  if (!inherits(learner, "LearnerSurv"))
    stop("learner must be mlr3::LearnerSurv. Got: ", class(learner)[1])
  if (!inherits(task, "TaskSurv"))
    stop("task must be mlr3::TaskSurv. Got: ", class(task)[1])
  
  if (n_background < 10L)
    warning("n_background < 10: estimates may be very noisy. Recommended >= 50.")
  
  # ---- Set defaults ----
  if (is.null(n_explain)) {
    n_explain <- if (type == "local") 20L else 50L
  }
  n_explain    <- min(as.integer(n_explain), nrow(task$data()))
  n_background <- min(as.integer(n_background), nrow(task$data()))
  
  if (verbose) {
    cat(sprintf("[SurvSHAP] mode=%s, n_explain=%d, n_background=%d\n",
                type, n_explain, n_background))
  }
  
  # ---- Build explainer ----
  data     <- as.data.frame(task$data())
  features <- data[, task$feature_names, drop = FALSE]
  
  target   <- survival::Surv(data[[task$target_names[1L]]],
                             data[[task$target_names[2L]]])
  
  predict_surv_fn <- function(model, newdata, times) {
    old_type <- model$predict_type
    model$predict_type <- "distr"
    on.exit(model$predict_type <- old_type, add = TRUE)
    
    time_col <- task$target_names[1L]
    status_col <- task$target_names[2L]
    
    tmp_df <- newdata
    tmp_df[[time_col]] <- 1
    tmp_df[[status_col]] <- 0
    
    tmp_task <- mlr3proba::TaskSurv$new("tmp", as.data.frame(tmp_df),
                                        time = time_col, event = status_col)
    pred <- model$predict(tmp_task)
    t(as.matrix(pred$distr$survival(times)))
  }
  
  explainer <- survex::explain_survival(
    learner, data = features, y = target,
    predict_survival_function = predict_surv_fn,
    label = learner$id, verbose = FALSE
  )
  
  # ---- Time-point subsampling for speed ----
  eval_times_full <- explainer$times
  if (!is.null(n_timepoints) && length(eval_times_full) > n_timepoints) {
    idx_subsample <- seq(1L, length(eval_times_full),
                         length.out = as.integer(n_timepoints))
    eval_times_use <- eval_times_full[idx_subsample]
    if (verbose) {
      cat(sprintf("[SurvSHAP] Subsampling times: %d → %d points (%.1fx speedup)\n",
                  length(eval_times_full), length(eval_times_use),
                  length(eval_times_full) / length(eval_times_use)))
    }
    explainer$times <- eval_times_use
  }
  
  # ---- Sample observations to explain ----
  obs_idx <- sample(nrow(features), size = n_explain, replace = FALSE)
  new_obs <- features[obs_idx, , drop = FALSE]
  
  # ---- Compute SHAP (either local or global mode) ----
  if (type == "global") {
    if (verbose) cat("[SurvSHAP] Computing GLOBAL SurvSHAP(t)...\n")
    shap_obj <- tryCatch(
      survex::model_survshap(
        explainer,
        new_observation    = new_obs,
        N                  = n_background,
        calculation_method = "kernelshap",
        aggregation_method = "mean_absolute"
      ),
      error = function(e) {
        stop("model_survshap failed. Diagnostic:\n  ", conditionMessage(e))
      }
    )
  } else {
    if (verbose) cat("[SurvSHAP] Computing LOCAL SurvSHAP(t)...\n")
    shap_obj <- tryCatch(
      survex::predict_parts(
        explainer,
        new_observation = new_obs,
        type            = "survshap",
        output_type     = "survival",
        N               = n_background
      ),
      error = function(e) {
        stop("predict_parts failed. Diagnostic:\n  ", conditionMessage(e))
      }
    )
  }
  
  # ---- Convert to tidy long format ----
  shap_long <- .survshap_to_long(shap_obj, obs_idx, features)
  
  if (is.null(shap_long) || nrow(shap_long) == 0L)
    stop("No SHAP rows returned. Check task/learner compatibility.")
  
  # ---- Generate plots ----
  plots <- .generate_shap_plots(shap_long, n_top_features, bar_color,
                                learner$id, type)
  
  if (verbose) {
    cat(sprintf("[SurvSHAP] Complete. %d observations × %d features × %d times.\n",
                length(unique(shap_long$observation)),
                length(unique(shap_long$feature)),
                length(unique(shap_long$time[!is.na(shap_long$time)]))))
  }
  
  list(
    shap_long         = shap_long,
    explainer         = explainer,
    eval_times        = if (!is.null(n_timepoints)) eval_times_use else eval_times_full,
    plots             = plots,
    original_features = features
  )
}


# =========================================================================
# Beeswarm / Violin Plot 
# =========================================================================

#' Survival SHAP Beeswarm/Violin Summary Plot
#'
#' This function creates a beeswarm or violin summary plot for SHAP (SHapley 
#' Additive exPlanations) values from survival models. It aggregates SHAP values 
#' across observations and optionally at a specific time point, displaying the 
#' most important features based on mean absolute SHAP values.
#'
#' @param shap_result A list object returned by a survival SHAP calculation function,
#'   typically containing `shap_long` (long-format SHAP values) and optionally
#'   `original_features` (original feature values for coloring).
#' @param time_point A numeric value specifying the time point at which to aggregate
#'   SHAP values. If `NULL` or all SHAP values have missing time, SHAP values are
#'   averaged across all time points. Default is `NULL`.
#' @param top_n An integer specifying the number of top features to display based
#'   on mean absolute SHAP values. Default is `8L`.
#' @param method A character string specifying the plot type. Must be one of
#'   `"beeswarm"` (default) or `"violin"`.
#' @param color_low A character string specifying the color for low feature values
#'   in the color gradient. Default is `"#2c7bb6"` (blue).
#' @param color_high A character string specifying the color for high feature values
#'   in the color gradient. Default is `"#d7191c"` (red).
#' @param title A character string for the plot title. If `NULL`, a default title
#'   is generated with time point and sample size information. Default is `NULL`.
#'
#' @return A \code{ggplot} object representing the SHAP summary plot. Points are
#'   colored by feature values if `original_features` is available in `shap_result`.
#'
#' @details
#' The function first aggregates SHAP values from the long-format input. If a
#' `time_point` is specified, SHAP values are filtered to the closest available
#' evaluation time before aggregation; otherwise, values are averaged across all
#' time points. Features are ranked by mean absolute SHAP value, and the top
#' \code{top_n} features are displayed.
#'
#' Two visualization methods are supported:
#' \itemize{
#'   \item \code{"beeswarm"}: Points are distributed along the y-axis to avoid
#'     overlap, with optional coloring by feature values.
#'   \item \code{"violin"}: Violin plots show the distribution of SHAP values,
#'     with quasi-random points overlaid for individual observations.
#' }
#'
#' When `original_features` is provided in `shap_result`, points are colored by
#' the corresponding feature values. Numeric features use a continuous color
#' gradient, while categorical features use discrete colors. The function uses
#' the \pkg{ggbeeswarm} package for beeswarm and quasirandom geometries.
#'
#' @note
#' The `shap_result` object must contain a `shap_long` data frame with columns
#' `observation`, `feature`, `shap_value`, and optionally `time`. If
#' `original_features` is provided, it should be a data frame with observations
#' as row names and features as columns.
#'
#' This function requires the following packages: \pkg{ggplot2}, \pkg{ggbeeswarm},
#' \pkg{ggprism}, \pkg{dplyr}, and \pkg{tidyr}.
#'
#' @importFrom dplyr group_by summarise filter arrange desc slice rename left_join
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_vline theme element_text
#'   labs scale_color_gradient geom_violin
#' @importFrom ggbeeswarm geom_beeswarm geom_quasirandom
#' @importFrom ggprism theme_prism
#' @export
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(dplyr)
#' library(ggplot2)
#'
#' # Example SHAP result structure (simulated)
#' set.seed(123)
#' shap_long <- data.frame(
#'   observation = rep(1:50, each = 5),
#'   feature = rep(c("age", "sex", "bmi", "stage", "treatment"), 50),
#'   shap_value = rnorm(250, 0, 1),
#'   time = rep(c(12, 24, 36), length.out = 250)
#' )
#'
#' original_features <- data.frame(
#'   age = rnorm(50, 65, 10),
#'   sex = factor(sample(c("M", "F"), 50, replace = TRUE)),
#'   bmi = rnorm(50, 28, 5),
#'   stage = factor(sample(1:4, 50, replace = TRUE)),
#'   treatment = factor(sample(c("A", "B"), 50, replace = TRUE))
#' )
#' rownames(original_features) <- 1:50
#'
#' shap_result <- list(
#'   shap_long = shap_long,
#'   original_features = original_features
#' )
#'
#' # Basic beeswarm plot (averaged across time)
#' surv_plot_shap_beeswarm(shap_result, top_n = 5)
#'
#' # Plot at specific time point with custom colors
#' surv_plot_shap_beeswarm(
#'   shap_result,
#'   time_point = 24,
#'   top_n = 6,
#'   color_low = "darkblue",
#'   color_high = "darkred"
#' )
#'
#' # Violin plot method
#' surv_plot_shap_beeswarm(
#'   shap_result,
#'   method = "violin",
#'   top_n = 4,
#'   title = "SHAP Summary - Violin Plot"
#' )
#'
#' # Without feature coloring
#' shap_result_no_color <- list(shap_long = shap_long)
#' surv_plot_shap_beeswarm(shap_result_no_color, top_n = 5)
#' }
surv_plot_shap_beeswarm <- function(shap_result,
                                    time_point = NULL,
                                    top_n      = 8L,
                                    method     = c("beeswarm", "violin"),
                                    color_low  = "#2c7bb6",
                                    color_high = "#d7191c",
                                    title      = NULL) {
  
  method    <- match.arg(method)
  shap_long <- shap_result$shap_long
  
  if (is.null(shap_long) || nrow(shap_long) == 0L)
    stop("shap_long is empty.")
  
  n_obs <- length(unique(shap_long$observation))
  
  if (is.null(time_point) || all(is.na(shap_long$time))) {
    shap_agg   <- shap_long %>%
      dplyr::group_by(observation, feature) %>%
      dplyr::summarise(shap = mean(shap_value, na.rm = TRUE), .groups = "drop")
    time_label <- ""
  } else {
    eval_times <- sort(unique(shap_long$time[!is.na(shap_long$time)]))
    closest_t  <- eval_times[which.min(abs(eval_times - time_point))]
    shap_agg   <- shap_long %>%
      dplyr::filter(abs(time - closest_t) < 1e-9) %>%
      dplyr::rename(shap = shap_value)
    time_label <- sprintf(" (t = %.1f)", closest_t)
  }
  
  feat_imp <- shap_agg %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_abs = mean(abs(shap), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice(seq_len(min(as.integer(top_n), dplyr::n())))
  
  top_features <- feat_imp$feature
  plot_df      <- shap_agg %>% dplyr::filter(feature %in% top_features)
  plot_df$feature <- factor(plot_df$feature, levels = rev(top_features))
  
  color_by_value <- FALSE
  is_numeric_scale <- FALSE
  
  if (!is.null(shap_result$original_features)) {
    var_vals <- shap_result$original_features
    
    var_vals[] <- lapply(var_vals, as.character)
    
    var_long <- tidyr::pivot_longer(
      cbind(observation = rownames(var_vals), var_vals),
      cols = -observation, names_to = "feature", values_to = "feature_value"
    )
    plot_df <- dplyr::left_join(plot_df, var_long, by = c("observation", "feature"))
    color_by_value <- !all(is.na(plot_df$feature_value))

    plot_df$feature_value_num <- suppressWarnings(as.numeric(as.character(plot_df$feature_value)))
    is_numeric_scale <- !all(is.na(plot_df$feature_value_num))
  }
  
  if (!requireNamespace("ggbeeswarm", quietly = TRUE))
    stop("Package 'ggbeeswarm' required.")
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = shap, y = feature))
  
  if (method == "beeswarm") {
    if (color_by_value) {
      if (is_numeric_scale) {
        p <- p + ggbeeswarm::geom_beeswarm(ggplot2::aes(color = feature_value_num), size = 2, cex = 2)
      } else {
        p <- p + ggbeeswarm::geom_beeswarm(ggplot2::aes(color = factor(feature_value)), size = 2, cex = 2)
      }
    } else {
      p <- p + ggbeeswarm::geom_beeswarm(ggplot2::aes(color = feature), size = 2, cex = 2)
    }
  } else {
    fill_aes <- if (color_by_value) {
      if (is_numeric_scale) ggplot2::aes(color = feature_value_num) else ggplot2::aes(color = factor(feature_value))
    } else {
      ggplot2::aes(color = feature)
    }
    p <- p + ggplot2::geom_violin(fill = "gray90", alpha = 0.5) + ggbeeswarm::geom_quasirandom(fill_aes, size = 2, width = 0.3)
  }
  
  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (color_by_value) {
    if (is_numeric_scale) {
      p <- p + ggplot2::scale_color_gradient(low = color_low, high = color_high, name = "Feature Value")
    } else {
      p <- p + ggplot2::labs(color = "Category")
    }
  }
  
  if (is.null(title))
    title <- sprintf("SurvSHAP Summary%s (%d patients)", time_label, n_obs)
  p + ggplot2::labs(title = title, x = "SHAP value", y = NULL)
}
# ==============================================================================
# 9. Clinical Reporting & Subgroup Analysis
# ==============================================================================

#' Generate Subgroup Forest Plot
#'
#' Creates a forest plot showing hazard ratios for risk score in different
#' subgroups defined by categorical variables. This helps assess whether the
#' model's prognostic effect is consistent across patient subgroups.
#'
#' @param learner A trained \code{mlr3} \code{LearnerSurv} object.
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param subgroup_vars A character vector of categorical column names to analyze.
#' @param prog Optional. A \code{PrognosiX} S4 object. If \code{NULL}, the function
#'   attempts to detect it from parent environments.
#' @param base_size A numeric value specifying the base font size for ggplot2.
#'   Default is \code{14}.
#' @param palette_name A character string specifying the Wes Anderson palette name.
#'   Default is \code{"Darjeeling1"}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{plot}{A \code{ggplot} object showing the forest plot.}
#'     \item{data}{A data frame containing the hazard ratio estimates for each subgroup.}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the task and data from the input object.
#'   \item For each subgroup variable, fits a Cox model with risk score as predictor.
#'   \item Computes hazard ratios and 95\% confidence intervals.
#'   \item Creates a forest plot with subgroups ordered by HR magnitude.
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' veteran$celltype <- as.factor(veteran$celltype)
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' # Forest plot for celltype subgroups
#' result <- surv_plot_subgroup_forest(
#'   learner = learner,
#'   object = task,
#'   subgroup_vars = "celltype"
#' )
#' print(result$plot)
#' }
surv_plot_subgroup_forest <- function(learner, object, subgroup_vars, prog = NULL,
                                      base_size = 14, palette_name = "Darjeeling1") {
  
  # =========================================================================
  # Phase 1: Polymorphic Input Parsing & Smart Prog Object Recovery
  # =========================================================================
  if (inherits(object, "TaskSurv")) {
    task <- object
    row_filter <- task$row_ids     # Extract row indices of current Task to isolate train/val sets
    cohort_name <- task$id
    
    # [Smart Dual-Track Recovery] Prioritize explicitly passed prog; if NULL, search parent environments
    if (is.null(prog)) {
      prog <- get0("prog", envir = parent.frame(), inherits = TRUE)
      if (is.null(prog)) prog <- get0("prog", envir = .GlobalEnv)
    }
    
    # Ultimate Defense: If prog is completely missing and features are not within the Task, throw an error
    missing_in_task <- setdiff(subgroup_vars, task$backend$cols)
    if (length(missing_in_task) > 0L && is.null(prog)) {
      stop(sprintf(
        "Error: Subgroup variables [%s] not found in mlr3 Task backend.\nNo 'prog' object detected in current environments. Please provide explicitly: prog = your_object",
        paste(missing_in_task, collapse = ", ")
      ))
    }
    
  } else if (inherits(object, "PrognosiX")) {
    prog <- object
    cohort_name <- "Full Cohort"
    
    # Extract Full Task if available
    task <- tryCatch({
      surv_extract_task(prog)
    }, error = function(e) NULL)
    
    if (!is.null(task)) {
      row_filter <- task$row_ids
    } else {
      # Fallback: Use all row positions if task extraction fails
      row_filter <- 1:nrow(prog@survival.data)
    }
  } else {
    stop("Error: 'object' must be an instance of 'TaskSurv' or 'PrognosiX'.")
  }
  
  # Dynamically retrieve target endpoint column names for the current survival analysis
  if (!is.null(task)) {
    time_col   <- task$target_names[1L]
    status_col <- task$target_names[2L]
  } else {
    time_col   <- prog@time_col %||% "time"
    status_col <- prog@status_col %||% "status"
  }
  
  # =========================================================================
  # Phase 2: Absolutely Safe Row Alignment & Data Slicing
  # =========================================================================
  # FIX: Use Base R 'inherits' instead of 'hasSlot' to ensure cross-platform safety
  if (!is.null(prog) && inherits(prog, "PrognosiX") && nrow(prog@survival.data) > 0) {
    raw_source <- as.data.frame(prog@survival.data)
    
    # Validate the existence of features in the data source
    available_vars <- intersect(subgroup_vars, colnames(raw_source))
    missing_vars   <- setdiff(subgroup_vars, colnames(raw_source))
    
    if (length(missing_vars) > 0L) {
      message(paste("Warning: The following variables were not found in prog and will be ignored:", 
                    paste(missing_vars, collapse = ", ")))
    }
    
    # Perform row slicing based on row_filter to ensure absolute alignment with current Task
    data <- raw_source[row_filter, unique(c(time_col, status_col, available_vars)), drop = FALSE]
  } else {
    # Fallback: Directly fetch data from Task Backend
    cols_to_fetch <- unique(c(time_col, status_col, subgroup_vars))
    valid_cols <- intersect(cols_to_fetch, task$backend$cols)
    data <- as.data.frame(task$backend$data(rows = row_filter, cols = valid_cols))
  }
  
  # Mount prediction risk scores (Crank) — physical row order corresponds strictly 1:1
  if (!is.null(task)) {
    predictions <- learner$predict(task)
    data$risk_score <- predictions$crank
  } else {
    stop("Error: Failed to construct mlr3 Task, unable to calculate risk scores.")
  }
  
  # Extract valid variables eventually used for the loop
  valid_vars <- intersect(subgroup_vars, colnames(data))
  if (length(valid_vars) == 0L) {
    stop("Error: No valid subgroup variables found in the parsed dataset.")
  }
  
  results <- list()
  
  # =========================================================================
  # Phase 3: Robust Subgroup Cox Regression Loop
  # =========================================================================
  for (var in valid_vars) {
    if (all(is.na(data[[var]]))) next
    
    levels_var <- unique(na.omit(data[[var]]))
    for (lev in levels_var) {
      subset_data <- data[data[[var]] == lev & !is.na(data[[var]]), ]
      
      n_total  <- nrow(subset_data)
      n_events <- sum(subset_data[[status_col]], na.rm = TRUE)
      
      # Robustness Filtering: If a subgroup has too few events (< 5), skip to prevent singularity
      if (n_events < 5L) {
        message(sprintf("Info: Skipping %s - %s due to insufficient events (n=%d, events=%d).", 
                        var, lev, n_total, n_events))
        next 
      }
      
      fit <- tryCatch({
        survival::coxph(
          survival::Surv(subset_data[[time_col]], subset_data[[status_col]]) ~ risk_score, 
          data = subset_data
        )
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        hr <- exp(stats::coef(fit))
        ci <- exp(confint(fit))
        
        if (is.na(hr) || is.infinite(hr) || any(is.na(ci))) next
        
        results[[paste0(var, "_", lev)]] <- data.frame(
          Variable = var,
          Subgroup = as.character(lev),
          N        = n_total,
          HR       = as.numeric(hr),
          Lower    = as.numeric(ci[1L]),
          Upper    = as.numeric(ci[2L])
        )
      }
    }
  }
  
  # =========================================================================
  # Phase 4: Forest Plot Rendering
  # =========================================================================
  if (length(results) == 0L) {
    warning("Warning: No subgroups met the statistical criteria. Execution aborted.")
    return(NULL)
  }
  
  res_df <- do.call(rbind, results)
  rownames(res_df) <- NULL
  
  res_df$Label <- paste0(res_df$Variable, " - ", res_df$Subgroup, " (n=", res_df$N, ")")
  res_df$Label <- factor(res_df$Label, levels = rev(res_df$Label))
  
  # Palette availability and safety check
  if (requireNamespace("wesanderson", quietly = TRUE)) {
    all_pals <- names(wesanderson::wes_palettes)
    active_pal <- if (palette_name %in% all_pals) palette_name else "Darjeeling1"
    colors <- wesanderson::wes_palette(active_pal, n = max(nrow(res_df), 2L), type = "continuous")
  } else {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(min(nrow(res_df), 8L), "Set1"))(nrow(res_df))
  }
  
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = HR, y = Label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper), height = 0.2, size = 0.8) +
    ggplot2::geom_point(ggplot2::aes(color = Label), size = 3) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = paste("Subgroup Analysis of Risk Score (", cohort_name, ")"),
      x = "Hazard Ratio (95% CI)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = base_size * 0.8)
    ) +
    ggplot2::scale_x_log10()
  
  print(p)
  return(list(plot = p, data = res_df))
}
# ==============================================================================
# 9. Stability & Sensitivity Analysis Module
# ==============================================================================

#' Robust Feature Stability Analysis for Sparse Survival Models
#'
#' Evaluates the stability of feature selection under repeated subsampling using
#' \code{glmnet::cv.glmnet} with Cox regression. Provides the Jaccard stability
#' index and selection frequencies for each feature.
#'
#' @param object A \code{PrognosiX} object, data frame, or \code{TaskSurv}.
#' @param time_col A character string specifying the name of the time column
#'   (required if \code{object} is a data frame). Default is \code{"time"}.
#' @param status_col A character string specifying the name of the status column
#'   (required if \code{object} is a data frame). Default is \code{"status"}.
#' @param n_repeat An integer specifying the number of subsampling iterations.
#'   Default is \code{30}.
#' @param train_ratio A numeric value specifying the proportion of data to sample
#'   each iteration. Default is \code{0.8}.
#' @param alpha Elastic net mixing parameter: 1 = LASSO, 0 = Ridge, 0.5 = elastic net.
#'   Default is \code{1}.
#' @param palette_name A character string specifying the name of the Wes Anderson
#'   palette. Default is \code{"AsteroidCity1"}.
#' @param seed An integer seed for reproducibility. Default is \code{2025}.
#' @param verbose A logical value. Print progress messages. Default is \code{TRUE}.
#'
#' @return A list with four components:
#'   \describe{
#'     \item{stability_index}{Jaccard stability index (mean pairwise Jaccard similarity).}
#'     \item{frequencies}{A data frame of features and their selection frequencies.}
#'     \item{plot}{A \code{ggplot} object of the top 15 features.}
#'     \item{success_rate}{The proportion of successful iterations.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data("veteran", package = "survival")
#' 
#' # Stability analysis on data frame
#' stab <- surv_analyze_feature_stability(
#'   object = veteran,
#'   time_col = "time",
#'   status_col = "status",
#'   n_repeat = 20,
#'   alpha = 1
#' )
#' 
#' print(paste("Stability Index:", round(stab$stability_index, 3)))
#' print(stab$plot)
#' }
#'
#' @seealso \code{\link{surv_feature_selection_multi}} for multi-method selection
#' @export
surv_analyze_feature_stability <- function(object,
                                           time_col = "time",
                                           status_col = "status",
                                           n_repeat = 30,
                                           train_ratio = 0.8,
                                           alpha = 1,
                                           palette_name = "AsteroidCity1",
                                           seed = 2025,
                                           verbose = TRUE) {
  
  # Extract data
  if (inherits(object, "PrognosiX")) {
    data <- object@survival.data
    time_col <- object@time_col
    status_col <- object@status_col
  } else if (is.data.frame(object)) {
    data <- object
  } else if (inherits(object, "TaskSurv")) {
    data <- as.data.frame(object$data())
    time_col <- object$target_names[1]
    status_col <- object$target_names[2]
  } else {
    stop("object must be a PrognosiX object, data frame, or TaskSurv")
  }
  
  # Required packages
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required for stability analysis.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")
  
  set.seed(seed)
  total_samples <- nrow(data)
  feature_cols <- setdiff(colnames(data), c(time_col, status_col))
  
  # One-hot encode factor features (robust, avoids mlr3 encoding issues)
  encode_data <- function(df) {
    df_encoded <- df[, c(time_col, status_col), drop = FALSE]
    for (col in feature_cols) {
      x <- df[[col]]
      if (is.factor(x) || is.character(x)) {
        x <- as.factor(x)
        dummies <- model.matrix(~ x - 1)
        colnames(dummies) <- paste0(col, "_", levels(x))
        df_encoded <- cbind(df_encoded, dummies)
      } else {
        df_encoded[[col]] <- x
      }
    }
    df_encoded <- df_encoded[, !duplicated(colnames(df_encoded))]
    colnames(df_encoded) <- make.names(colnames(df_encoded))
    return(df_encoded)
  }
  
  data_encoded <- encode_data(data)
  encoded_features <- setdiff(colnames(data_encoded), c(time_col, status_col))
  
  if (length(encoded_features) == 0) {
    stop("No features available after encoding. Check input data.")
  }
  
  if (verbose) {
    message(sprintf("\n[*] Starting Feature Stability Analysis (%d iterations)...", n_repeat))
    message(sprintf("    Subsample ratio: %.0f%% | Elastic net alpha = %.1f", train_ratio * 100, alpha))
    message(sprintf("    Original features: %d | Encoded features: %d", 
                    length(feature_cols), length(encoded_features)))
  }
  
  selected_sets <- list()
  success_count <- 0
  
  # Progress bar if interactive
  if (verbose && interactive() && requireNamespace("utils", quietly = TRUE)) {
    pb <- utils::txtProgressBar(min = 0, max = n_repeat, style = 3)
  } else {
    pb <- NULL
  }
  
  for (i in seq_len(n_repeat)) {
    idx <- sample(total_samples, size = floor(total_samples * train_ratio))
    sub_data <- data_encoded[idx, , drop = FALSE]
    
    y <- survival::Surv(sub_data[[time_col]], sub_data[[status_col]])
    x <- as.matrix(sub_data[, encoded_features, drop = FALSE])
    
    fit <- tryCatch({
      glmnet::cv.glmnet(x = x, y = y, family = "cox", alpha = alpha, nfolds = 5)
    }, error = function(e) {
      if (verbose && i <= 3) warning(sprintf("Iter %d failed: %s", i, e$message))
      NULL
    })
    
    if (is.null(fit)) {
      selected_sets[[i]] <- character(0)
      next
    }
    
    coefs <- as.matrix(stats::coef(fit, s = "lambda.min"))
    selected <- rownames(coefs)[coefs[, 1] != 0]
    selected <- setdiff(selected, "(Intercept)")
    
    if (length(selected) == 0) selected <- "none"
    selected_sets[[i]] <- selected
    success_count <- success_count + 1
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
  }
  
  if (!is.null(pb)) close(pb)
  
  if (verbose) {
    message(sprintf("\n    Successful iterations: %d / %d", success_count, n_repeat))
  }
  
  if (success_count == 0) {
    warning("No successful iterations. Returning empty result.")
    return(list(stability_index = NA, frequencies = data.frame(), plot = NULL, success_rate = 0))
  }
  
  valid_sets <- selected_sets[sapply(selected_sets, length) > 0 & 
                                sapply(selected_sets, function(x) !all(x == "none"))]
  
  if (length(valid_sets) == 0) {
    warning("No features selected in any iteration. Try increasing sample size or changing alpha.")
    return(list(stability_index = NA, frequencies = data.frame(), plot = NULL, success_rate = 0))
  }
  
  # Selection frequencies
  all_selected <- unlist(valid_sets)
  freq_tab <- sort(table(all_selected), decreasing = TRUE)
  freq_df <- data.frame(
    Feature = names(freq_tab),
    Frequency = as.numeric(freq_tab) / length(valid_sets)
  )
  
  # Jaccard stability index
  stab_index <- NA
  if (length(valid_sets) > 1) {
    jaccard_vals <- c()
    for (i in seq_len(length(valid_sets) - 1)) {
      for (j in (i + 1):length(valid_sets)) {
        inter <- length(intersect(valid_sets[[i]], valid_sets[[j]]))
        union <- length(union(valid_sets[[i]], valid_sets[[j]]))
        if (union > 0) jaccard_vals <- c(jaccard_vals, inter / union)
      }
    }
    stab_index <- mean(jaccard_vals, na.rm = TRUE)
  }
  
  # Plot top 15 features with optional wesanderson color
  plot_df <- head(freq_df, 15)
  p <- NULL
  if (nrow(plot_df) > 0) {
    # Determine bar color
    if (requireNamespace("wesanderson", quietly = TRUE)) {
      bar_color <- wesanderson::wes_palette(palette_name, n = 1, type = "continuous")
    } else {
      bar_color <- "#2980b9"  # default blue
      if (verbose) message("    Note: Install 'wesanderson' for more colorful plots.")
    }
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = stats::reorder(Feature, Frequency), y = Frequency)) +
      ggplot2::geom_bar(stat = "identity", fill = bar_color, width = 0.6) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::labs(
        x = "Feature",
        y = "Selection Frequency",
        title = "Feature Selection Stability (LASSO/Elastic Net)",
        subtitle = sprintf(
          "Subsampling ratio: %.0f%% | Iterations: %d | Jaccard Index: %.3f",
          train_ratio * 100, length(valid_sets), stab_index
        )
      ) +
      ggprism::theme_prism(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(p)
  }
  
  if (verbose) {
    message(sprintf("  ✓ Stability index (Jaccard): %.3f", stab_index))
  }
  
  list(
    stability_index = stab_index,
    frequencies = freq_df,
    plot = p,
    success_rate = success_count / n_repeat
  )
}


#' Analyze Model Performance Sensitivity
#'
#' Evaluates how model performance (C-index) changes under varying conditions
#' such as sample size or censoring rate. This helps assess the robustness of
#' the model to data perturbations.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param analysis_type A character string specifying the type of sensitivity
#'   analysis. Must be one of \code{"sample_size"} or \code{"censoring"}.
#' @param param_values A numeric vector of parameter values to test. For
#'   \code{"sample_size"}, these are proportions (e.g., \code{c(0.3, 0.5, 0.7, 0.9, 1.0)}).
#'   For \code{"censoring"}, these are additional censoring proportions
#'   (e.g., \code{c(0.1, 0.2, 0.3, 0.5)}). If \code{NULL}, uses reasonable defaults.
#' @param palette_name A character string specifying the Wes Anderson palette name.
#'   Default is \code{"AsteroidCity1"}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{results}{A data frame with columns: \code{Parameter}, \code{C_Index}, and \code{SE}.}
#'     \item{plot}{A \code{ggplot} object showing the sensitivity trajectory.}
#'   }
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' 
#' # Sample size sensitivity
#' sens <- surv_analyze_model_sensitivity(
#'   object = task,
#'   learner_id = "surv.coxph",
#'   analysis_type = "sample_size",
#'   param_values = c(0.3, 0.5, 0.7, 0.9)
#' )
#' print(sens$plot)
#' }
surv_analyze_model_sensitivity <- function(object, learner_id, analysis_type = c("sample_size", "censoring"), param_values = NULL, palette_name = "AsteroidCity1") {
  task <- surv_extract_task(object)
  
  analysis_type <- match.arg(analysis_type)
  original_data <- as.data.table(task$data())
  
  message(sprintf("\n[*] Starting Sensitivity Analysis: %s...", analysis_type))
  
  results <- data.frame(Parameter = numeric(), C_Index = numeric(), SE = numeric())
  
  # Default parameters if not provided
  if (is.null(param_values)) {
    if (analysis_type == "sample_size") param_values <- c(0.3, 0.5, 0.7, 0.9, 1.0)
    # For censoring: values are ADDITIONAL censoring PROPORTIONS (0 = no extra censoring,
    # 0.5 = censor an additional 50% of currently observed events at random times).
    # This is NOT the rate parameter of an exponential distribution.
    if (analysis_type == "censoring") param_values <- c(0.1, 0.2, 0.3, 0.5)
  }
  
  for (val in param_values) {
    
    cv_scores <- numeric(5L)
    
    for (fold in seq_len(5L)) {
      temp_data <- data.table::copy(original_data)
      
      if (analysis_type == "sample_size") {
        keep_idx  <- sample(nrow(temp_data), size = max(20L, floor(nrow(temp_data) * val)))
        temp_data <- temp_data[keep_idx, ]
      } else if (analysis_type == "censoring") {
        # val is the ADDITIONAL CENSORING PROPORTION: the fraction of currently
        # observed events (status == 1) that will be randomly censored.
        # We randomly select val*100% of event rows and replace their observed
        # event time with a uniform random time in [0, observed_time).
        # This correctly controls the additional censoring fraction and avoids
        # the misinterpretation of 'rate' in rexp().
        time_col   <- task$target_names[1L]
        status_col <- task$target_names[2L]
        event_idx  <- which(temp_data[[status_col]] == 1)
        n_to_censor <- floor(length(event_idx) * val)
        if (n_to_censor > 0L) {
          censor_idx <- sample(event_idx, size = n_to_censor, replace = FALSE)
          # Censoring time uniform in (0, observed_time] for each selected subject
          new_times  <- stats::runif(n_to_censor, min = 0,
                                     max = temp_data[[time_col]][censor_idx])
          new_times  <- pmax(new_times, .Machine$double.eps)  # avoid zero times
          temp_data[[time_col]][censor_idx]   <- new_times
          temp_data[[status_col]][censor_idx] <- 0L
        }
      }
      
      # 2. Create Temporary Task & Evaluate
      temp_task <- surv_create_surv_task(temp_data, task$target_names[1], task$target_names[2], id = "temp_sens")
      
      res <- tryCatch({
        rr <- resample(temp_task, lrn(learner_id), rsmp("cv", folds = 3), store_models = FALSE)
        rr$aggregate(msr("surv.cindex"))
      }, error = function(e) NA)
      
      cv_scores[fold] <- res
    }
    
    # 3. Store Results (Mean and Standard Error)
    results <- rbind(results, data.frame(
      Parameter = val,
      C_Index = mean(cv_scores, na.rm = TRUE),
      SE = stats::sd(cv_scores, na.rm = TRUE) / sqrt(sum(!is.na(cv_scores)))
    ))
  }
  
  # Generate Plot
  x_label <- ifelse(analysis_type == "sample_size",
                    "Proportion of Total Sample Size",
                    "Additional Censoring Proportion (fraction of events re-censored)")
  
  # Get colors from wesanderson palette
  if (requireNamespace("wesanderson", quietly = TRUE)) {
    colors <- wesanderson::wes_palette(palette_name, n = 3, type = "continuous")
    line_color <- colors[1]
    point_color <- colors[2]
    error_color <- colors[3]
  } else {
    line_color <- "#E74C3C"
    point_color <- "#C0392B"
    error_color <- "#7F8C8D"
  }
  
  p <- ggplot2::ggplot(results, ggplot2::aes(x = Parameter, y = C_Index)) +
    ggplot2::geom_line(color = line_color, size = 1) +
    ggplot2::geom_point(color = point_color, size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = C_Index - SE, ymax = C_Index + SE), width = 0.05, color = error_color) +
    ggplot2::scale_y_continuous(limits = c(0.45, 1)) +
    ggplot2::labs(
      x = x_label, 
      y = "Cross-Validated C-Index (± SE)", 
      title = "Model Sensitivity Analysis",
      subtitle = sprintf("Model: %s | Perturbation: %s", learner_id, analysis_type)
    ) +
    ggprism::theme_prism(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  print(p)
  message("  ✓ Sensitivity analysis complete.")
  
  list(results = results, plot = p)
}

#' Feature Ablation Sensitivity Analysis
#'
#' Evaluates the impact of removing individual features on model performance.
#' For each feature, the model is retrained without that feature and the change
#' in C-index is recorded.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param features_to_test A character vector of feature names to test. If \code{NULL},
#'   tests all features in the task. Default is \code{NULL}.
#'
#' @return A list with three components:
#'   \describe{
#'     \item{results}{A data frame with columns: \code{Feature_Removed}, \code{New_CIndex},
#'       and \code{Performance_Drop}.}
#'     \item{plot}{A \code{ggplot} object showing the performance drop for each feature.}
#'     \item{baseline}{The baseline C-index with all features.}
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Calculates baseline performance using all features.
#'   \item For each feature, creates a task without that feature and evaluates performance.
#'   \item Reports the drop in C-index for each feature.
#'   \item Creates a bar plot sorted by impact magnitude.
#' }
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' 
#' # Ablation analysis on all features
#' ablation <- surv_analyze_feature_ablation(
#'   object = task,
#'   learner_id = "surv.coxph"
#' )
#' print(ablation$plot)
#' 
#' # Test specific features
#' ablation_subset <- surv_analyze_feature_ablation(
#'   object = task,
#'   learner_id = "surv.coxph",
#'   features_to_test = c("age", "karno", "celltype")
#' )
#' }
#'
#' @export
surv_analyze_feature_ablation <- function(object, learner_id, features_to_test = NULL) {
  task <- surv_extract_task(object)
  
  if (is.null(features_to_test)) {
    features_to_test <- task$feature_names
  }
  
  message(sprintf("\n[*] Starting Feature Ablation Analysis for %d features...", length(features_to_test)))
  
  # 1. Calculate Baseline Performance (using all features)
  learner <- surv_get_learner(learner_id, task)
  
  # Use a small 3-fold CV to get a stable baseline
  resampling <- rsmp("cv", folds = 3)
  baseline_rr <- resample(task, learner, resampling, store_models = FALSE)
  baseline_cindex <- baseline_rr$aggregate(msr("surv.cindex"))
  
  results <- list()
  
  # 2. Iterate and Remove Each Feature
  for (feat in features_to_test) {
    # Create a task without the specific feature
    remaining_feats <- setdiff(task$feature_names, feat)
    temp_task <- task$clone()$select(remaining_feats)
    
    # Evaluate model performance without this feature
    temp_rr <- tryCatch({
      resample(temp_task, learner, resampling, store_models = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(temp_rr)) {
      new_cindex <- temp_rr$aggregate(msr("surv.cindex"))
      drop <- baseline_cindex - new_cindex
      
      results[[feat]] <- data.frame(
        Feature_Removed = feat,
        New_CIndex = new_cindex,
        Performance_Drop = drop
      )
    }
  }
  
  res_df <- do.call(rbind, results)
  res_df <- res_df[order(-res_df$Performance_Drop), ] # Sort by impact
  
  # 3. Generate Visualization
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = stats::reorder(Feature_Removed, Performance_Drop), y = Performance_Drop)) +
    ggplot2::geom_bar(stat = "identity", fill = ifelse(res_df$Performance_Drop > 0, "#E67E22", "#95A5A6")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Feature Removed",
      y = "Drop in C-Index (Baseline - New)",
      title = "Feature Ablation Sensitivity Analysis",
      subtitle = sprintf("Baseline C-Index (Full Model): %.4f", baseline_cindex)
    ) +
    ggprism::theme_prism()
  
  print(p)
  message("  ✓ Feature ablation analysis complete.")
  
  return(list(results = res_df, plot = p, baseline = baseline_cindex))
}

# ==============================================================================
# 11. Clinical Validation: Calibration & Time-dependent AUC
# ==============================================================================
#' Plot Calibration Curve with Quantitative Metrics
#'
#' Calibration is assessed by grouping predicted survival probabilities into
#' quantile bins and comparing the bin mean predicted probability to the
#' Kaplan-Meier observed survival probability in that bin.
#'
#' IMPORTANT: When evaluated on the TRAINING SET the curve reflects "apparent"
#' (optimistic) calibration. For valid calibration use an independent test set
#' or cross-validated out-of-fold predictions.
#'
#' ICI follows Austin & Steyerberg (2019): loess of observed ~ predicted fitted
#' on bin centres (span = 1.0), then |loess(x) - x| integrated over the
#' prediction range. Falls back to trapezoidal rule when fewer than 4 bins.
#'https://doi.org/10.1002/sim.8281
#'
#' @param learner Trained mlr3 learner with predict_type "distr"
#' @param object TaskSurv or PrognosiX object
#' @param time_point Numeric. Time horizon for calibration.
#' @param n_bins Integer. Number of bins (default 10; minimum 5 recommended).
#' @param apparent Logical. If TRUE (default) labels the plot as APPARENT
#'   (training-set) calibration and prints a reminder.
#' @param print_metrics Logical. Print metrics to console (default TRUE).
#' @param show_ici Logical. Include ICI in plot subtitle (default TRUE).
#' @return ggplot object; calibration metrics attached as attribute "calibration_metrics".
#' @export
#' 
#' @examples
#' \dontrun{
#' # Requires trained learner with distr predict_type and task
#' # cal_plot <- surv_plot_calibration(learner, task, time_point = 365)
#' }
surv_plot_calibration <- function(learner, object, time_point, n_bins = 10,
                                  apparent = TRUE,
                                  print_metrics = TRUE, show_ici = TRUE) {
  # ---- 1. Environment Protection: Prevent altering external R6 learner state ----
  if (!("distr" %in% learner$predict_types)) {
    warning("Learner does not support 'distr' predict_type. Calibration skipped.")
    return(NULL)
  }
  old_predict_type <- learner$predict_type
  on.exit({ learner$predict_type <- old_predict_type }, add = TRUE)
  
  # ---- 2. Data Preparation and Cleaning ----
  cal_df <- .prepare_cal_data(learner = learner, object = object, 
                              time_point = time_point, n_bins = n_bins)
  if (is.null(cal_df)) return(NULL)
  
  if (apparent) {
    message(paste(
      "[Calibration] Evaluating on the TRAINING SET.",
      "This is APPARENT (optimistic) calibration.",
      "For valid calibration use an independent test set or CV out-of-fold predictions."
    ))
  }
  
  # ---- 3. Metric Computations ----
  metrics <- .compute_cal_metrics(cal_df)
  
  if (print_metrics) {
    cat("\n========== Calibration Metrics at t =", time_point, "==========\n")
    if (apparent) cat("** APPARENT (training-set) calibration -- interpret with caution **\n")
    cat(sprintf("Calibration slope    (ideal 1): %.4f\n", metrics$slope))
    cat(sprintf("Calibration intercept (ideal 0): %.4f\n", metrics$intercept))
    cat(sprintf("R-squared:                       %.4f\n", metrics$r_squared))
    cat(sprintf("Mean Absolute Error (MAE):       %.4f\n", metrics$mae))
    cat(sprintf("Integrated Calibration Index:    %s\n",
                ifelse(is.na(metrics$ici), "NA", sprintf("%.4f", metrics$ici))))
    cat(sprintf("E50 (median abs error):          %.4f\n", metrics$e50))
    cat(sprintf("E90 (90th pct abs error):        %.4f\n", metrics$e90))
    cat("====================================================\n")
  }
  
  # ---- 4. Plot Configurations ----
  slope_label <- round(metrics$slope, 3)
  ici_label   <- ifelse(is.na(metrics$ici), "NA", round(metrics$ici, 3))
  plot_title  <- ifelse(apparent, 
                        paste("Apparent Calibration at t =", time_point, "(training set)"),
                        paste("Calibration Curve at t =", time_point))
  sub_title   <- ifelse(show_ici, 
                        paste(learner$id, "| Slope =", slope_label, "| ICI =", ici_label),
                        learner$id)
  
  p <- ggplot2::ggplot(cal_df, ggplot2::aes(x = predicted, y = observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "red", linewidth = 0.8) +
    ggplot2::geom_line(color = "#2980b9", linewidth = 0.9) +
    ggplot2::geom_point(size = 3, color = "#2980b9") +
    ggplot2::labs(x = "Predicted Survival Probability",
                  y = "Observed Survival Probability (KM)",
                  title = plot_title, subtitle = sub_title) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) # Core fix: Prevent early data clipping from ruining the path line
  
  # Graceful fallback if ggprism is not available
  if (requireNamespace("ggprism", quietly = TRUE)) {
    p <- p + ggprism::theme_prism()
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  
  attr(p, "calibration_metrics") <- metrics
  attr(p, "calibration_data")    <- cal_df
  return(p)
}


#' Plot Comparison Calibration (Apparent Training vs Validation)
#'
#' Overlays calibration curves for the training set (apparent/optimistic) and
#' a held-out validation set on the same plot.
#'
#' @param learner Trained mlr3 learner (must support "distr" predict_type)
#' @param train_task Training TaskSurv object
#' @param val_task Validation TaskSurv object
#' @param time_point Numeric. Evaluation time point for calibration.
#' @param n_bins Integer. Number of bins (default 10; minimum 5).
#' @param print_metrics Logical. Print metrics for both datasets (default TRUE).
#' @return A ggplot object showing calibration curves for training and validation.
#' @export
#' 
#' @examples
#' \dontrun{
#' # Requires trained learner, training task and validation task
#' # comp_plot <- surv_plot_comparison_calibration(learner, train_task, val_task, time_point = 365)
#' }
surv_plot_comparison_calibration <- function(learner, train_task, val_task,
                                             time_point, n_bins = 10,
                                             print_metrics = TRUE) {
  # ---- 1. Environment Protection: Avoid polluting the external R6 learner state ----
  if (!("distr" %in% learner$predict_types)) {
    warning("Learner does not support 'distr' predict_type. Calibration skipped.")
    return(NULL)
  }
  old_predict_type <- learner$predict_type
  on.exit({ learner$predict_type <- old_predict_type }, add = TRUE)
  
  # ---- 2. Extract underlying data using shared internal function to eliminate redundancy ----
  cal_train <- .prepare_cal_data(learner, train_task, time_point, n_bins)
  cal_val   <- .prepare_cal_data(learner, val_task, time_point, n_bins)
  
  if (is.null(cal_train) || is.null(cal_val)) {
    warning("Insufficient bins or missing columns for calibration comparison.")
    return(NULL)
  }
  
  m_train <- .compute_cal_metrics(cal_train)
  m_val   <- .compute_cal_metrics(cal_val)
  
  if (print_metrics) {
    cat("\n========== Calibration Metrics at t =", time_point, "==========\n")
    cat("** TRAINING (apparent -- optimistic for flexible learners): **\n")
    cat(sprintf("  Slope=%.4f | Intercept=%.4f | R2=%.4f | MAE=%.4f | ICI=%.4f\n",
                m_train$slope, m_train$intercept, m_train$r_squared, m_train$mae, m_train$ici))
    cat("** VALIDATION (unbiased estimate): **\n")
    cat(sprintf("  Slope=%.4f | Intercept=%.4f | R2=%.4f | MAE=%.4f | ICI=%.4f\n",
                m_val$slope, m_val$intercept, m_val$r_squared, m_val$mae, m_val$ici))
    cat("  Note: Training calibration is always expected to be more optimistic.\n")
    cat("=================================================================\n")
  }
  
  cal_train$Dataset <- "Training (apparent)"
  cal_val$Dataset   <- "Validation"
  cal_all <- rbind(cal_train, cal_val)
  
  ici_tr  <- ifelse(is.na(m_train$ici), "NA", round(m_train$ici, 3))
  ici_val <- ifelse(is.na(m_val$ici),   "NA", round(m_val$ici,   3))
  
  p <- ggplot2::ggplot(cal_all, ggplot2::aes(x = predicted, y = observed, color = Dataset)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "gray40", linewidth = 0.8) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = c("Training (apparent)" = "#2980b9", "Validation" = "#e74c3c")) +
    ggplot2::labs(x = "Predicted Survival Probability",
                  y = "Observed Survival Probability (KM)",
                  title = paste("Calibration Comparison at t =", time_point),
                  subtitle = paste0(learner$id, " | Train ICI = ", ici_tr, " (apparent) | Val ICI = ", ici_val)) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) # Core fix: Avoid line breakage due to strict axis clipping
  
  if (requireNamespace("ggprism", quietly = TRUE)) {
    p <- p + ggprism::theme_prism()
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  
  return(p)
}

#' Plot Continuous Time-dependent AUC (Dynamic AUC)
#' @param learner Trained mlr3 learner
#' @param task TaskSurv object
#' @return ggplot object
#' @export
surv_plot_time_dependent_auc <- function(learner, object) {
  task <- surv_extract_task(object)
  if (!requireNamespace("risksetROC", quietly = TRUE)) stop("Please install 'risksetROC'")
  
  data <- as.data.frame(task$data())
  pred <- learner$predict(task)
  
  # Extract Time, Status, and Marker (linear predictor or crank)
  # risksetROC needs 'marker' where higher value means higher risk
  stime <- data[[task$target_names[1]]]
  status <- data[[task$target_names[2]]]
  marker <- pred$crank 
  
  # Get unique event times for plotting (as in your snippet)
  utimes <- sort(unique(stime[status == 1]))
  
  message("[*] Calculating Time-dependent AUC across all event times...")
  
  # Use risksetROC to calculate AUC at each time point
  out <- risksetROC::risksetAUC(
    Stime = stime,
    status = status,
    marker = marker,
    tmax = max(stime) * 0.95, # Avoid instability at the very end
    plot = FALSE
  )
  
  # Prepare data for ggplot
  auc_df <- data.frame(
    times = utimes,
    tAUC = out$AUC[match(utimes, out$utimes)]
  )
  # Remove NAs
  auc_df <- na.omit(auc_df)
  
  # Create the plot similar to your provided ggplot snippet
  p <- ggplot2::ggplot(auc_df, ggplot2::aes(x = times, y = tAUC)) +
    ggplot2::geom_step(direction = "vh", color = "#2C3E50", size = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0.5, ymax = tAUC), fill = "#3498DB", alpha = 0.1) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ggplot2::labs(
      x = "Evaluation Time Points",
      y = "AUC",
      title = "Time-Dependent AUC (Dynamic Accuracy)",
      subtitle = sprintf("Model: %s | Integrated AUC: %.3f", learner$id, mean(auc_df$tAUC))
    ) +
    ggplot2::ylim(0.4, 1.0) +
    ggprism::theme_prism()
  
  print(p)
  return(auc_df)
}

#' Feature Selection Pipeline: Univariate to Lasso
#' @param task TaskSurv object
#' @param p_threshold P-value threshold for univariate filtering (default 0.05)
#' @return A list containing the filtered task and the univariate results table
#' @export
surv_filter_features_clinical <- function(object, p_threshold = 0.05) {
  task <- surv_extract_task(object)
  
  message("[*] Starting Feature Selection: Univariate Cox Filtering...")
  
  data <- as.data.frame(task$data())
  features <- task$feature_names
  target_time <- task$target_names[1]
  target_status <- task$target_names[2]
  
  unv_results <- list()
  
  for (feat in features) {
    formula_str <- sprintf("survival::Surv(%s, %s) ~ %s", target_time, target_status, feat)
    fit <- tryCatch({
      survival::coxph(as.formula(formula_str), data = data)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      s <- summary(fit)
      p_val <- s$coefficients[1, "Pr(>|z|)"]
      unv_results[[feat]] <- data.frame(
        Feature = feat,
        HR = s$conf.int[1, "exp(coef)"],
        P_Value = p_val,
        Lower_CI = s$conf.int[1, "lower .95"],
        Upper_CI = s$conf.int[1, "upper .95"]
      )
    }
  }
  
  unv_df <- do.call(rbind, unv_results)
  # Filter significant features
  significant_feats <- unv_df$Feature[unv_df$P_Value < p_threshold]
  
  message(sprintf("  [✓] Univariate filter complete: %d -> %d features", length(features), length(significant_feats)))
  
  # Return a task with only significant features
  new_task <- task$clone()$select(significant_feats)
  
  # Visualization: P-value Lollipop Plot
  p <- ggplot2::ggplot(unv_df, ggplot2::aes(x = stats::reorder(Feature, -P_Value), y = -log10(P_Value))) +
    ggplot2::geom_segment(ggplot2::aes(xend = Feature, yend = 0), color = "grey") +
    ggplot2::geom_point(ggplot2::aes(color = P_Value < p_threshold), size = 3) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "red") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Univariate Feature Filtering", y = "-log10(P-value)", x = "Features") +
    ggprism::theme_prism()
  
  print(p)
  
  return(list(task = new_task, table = unv_df, plot = p))
}


#' Benchmark Multiple Survival Algorithms
#'
#' Compares the performance of multiple survival learners using cross-validation.
#' The function automatically creates a benchmarking design, runs the comparison,
#' and generates a boxplot visualization of the C-index scores across folds.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object. The function
#'   extracts the task using \code{surv_extract_task()}.
#' @param learners_list A list of \code{Learner} objects to benchmark. If \code{NULL},
#'   defaults to four commonly used survival learners:
#'   \itemize{
#'     \item \code{surv.coxph} (Cox Proportional Hazards)
#'     \item \code{surv.cv_glmnet} (LASSO)
#'     \item \code{surv.ranger} (Random Forest)
#'     \item \code{surv.xgboost} (XGBoost)
#'   }
#' @param resampling A \code{\link[mlr3]{Resampling}} object. If \code{NULL},
#'   defaults to 5-fold cross-validation.
#'
#' @return A list with three components:
#'   \describe{
#'     \item{bmr}{A \code{\link[mlr3benchmark]{BenchmarkResult}} object
#'       containing all benchmark results.}
#'     \item{table}{A data frame of aggregated performance metrics (C-index)
#'       for each learner.}
#'     \item{plot}{A \code{ggplot} object showing the performance distribution
#'       across CV folds for each learner.}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the survival task from \code{object}.
#'   \item If no learners are provided, instantiates default learners.
#'   \item Ensures all learners use the same \code{predict_type} (preferring
#'     \code{"distr"} if available) for fair comparison.
#'   \item Creates a benchmark design and runs the benchmark.
#'   \item Aggregates performance using the concordance index (\code{surv.cindex}).
#'   \item Generates a boxplot comparing the performance distribution.
#' }
#'
#' @note
#' The default learners require additional packages:
#' \itemize{
#'   \item \code{surv.coxph}: built-in to \code{mlr3proba}
#'   \item \code{surv.cv_glmnet}: requires \code{mlr3learners} and \code{glmnet}
#'   \item \code{surv.ranger}: requires \code{mlr3learners} and \code{ranger}
#'   \item \code{surv.xgboost}: requires \code{mlr3extralearners} and \code{xgboost}
#' }
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' 
#' # Run benchmark with default learners
#' bm <- surv_run_algorithm_benchmark(object = task)
#' 
#' # View performance table
#' print(bm$table)
#' 
#' # Custom learners
#' library(mlr3learners)
#' custom_learners <- list(
#'   lrn("surv.coxph", id = "CoxPH"),
#'   lrn("surv.ranger", id = "RF", num.trees = 100)
#' )
#' 
#' bm_custom <- surv_run_algorithm_benchmark(
#'   object = task,
#'   learners_list = custom_learners,
#'   resampling = rsmp("cv", folds = 3)
#' )
#' }
#'
#' @seealso
#' \code{\link{surv_benchmark_learners}} for more detailed benchmarking with tuning support,
#' \code{\link{surv_summarize_benchmark}} for summarizing results
#'
#' @export
surv_run_algorithm_benchmark <- function(object, learners_list = NULL, resampling = NULL) {
  task <- surv_extract_task(object)
  
  if (is.null(learners_list)) {
    learners_list <- list(
      lrn("surv.coxph", id = "CoxPH"),
      lrn("surv.cv_glmnet", id = "Lasso"),
      lrn("surv.ranger", id = "RandomForest"),
      lrn("surv.xgboost", id = "XGBoost")
    )
  }
  
  # Ensure all learners use the same predict_type for fair comparison
  for (l in learners_list) {
    if ("distr" %in% l$predict_types) l$predict_type <- "distr"
  }
  
  if (is.null(resampling)) resampling <- rsmp("cv", folds = 5)
  
  message("[*] Running Benchmark: Comparing Algorithms via 5-fold CV...")
  
  design <- benchmark_grid(task, learners_list, resampling)
  bmr <- benchmark(design)
  
  # Measure performance
  measures <- list(msr("surv.cindex"))
  perf_tab <- bmr$aggregate(measures)
  
  # Visualization: Boxplot of performance across folds
  p <- autoplot(bmr, measure = msr("surv.cindex")) +
    ggplot2::labs(title = "Algorithm Performance Comparison", subtitle = "5-Fold Cross-Validation C-Index") +
    ggprism::theme_prism() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  print(p)
  
  return(list(bmr = bmr, table = perf_tab, plot = p))
}

# ==============================================================================
# 12. Pipeline Helper Functions
# ==============================================================================

#' Check Data Quality for Survival Analysis
#' @param data Data frame to check
#' @param time_col Time column name
#' @param event_col Event/status column name
#' @return Invisible NULL, stops on error
#' @keywords internal
check_data_quality <- function(data, time_col, event_col) {
  # Check if data is a data frame
  if (!is.data.frame(data)) {
    stop("train_data must be a data frame.")
  }
  
  # Check required columns exist
  if (!time_col %in% colnames(data)) {
    stop(sprintf("Time column '%s' not found in data.", time_col))
  }
  if (!event_col %in% colnames(data)) {
    stop(sprintf("Event column '%s' not found in data.", event_col))
  }
  
  # Check for missing values in target columns
  time_missing <- sum(is.na(data[[time_col]]))
  event_missing <- sum(is.na(data[[event_col]]))
  
  if (time_missing > 0) {
    warning(sprintf("Time column has %d missing values.", time_missing))
  }
  if (event_missing > 0) {
    warning(sprintf("Event column has %d missing values.", event_missing))
  }
  
  # Check time is numeric
  if (!is.numeric(data[[time_col]])) {
    stop("Time column must be numeric.")
  }
  
  # Check event is binary (0/1)
  event_vals <- unique(na.omit(data[[event_col]]))
  if (!all(event_vals %in% c(0, 1))) {
    warning("Event column should be binary (0/1). Found values: ", paste(event_vals, collapse = ", "))
  }
  
  # Check for non-finite times
  non_finite <- sum(!is.finite(data[[time_col]]))
  if (non_finite > 0) {
    stop(sprintf("Time column has %d non-finite values (Inf, -Inf, or NaN).", non_finite))
  }
  
  # Check for negative times
  neg_times <- sum(data[[time_col]] < 0, na.rm = TRUE)
  if (neg_times > 0) {
    warning(sprintf("Time column has %d negative values.", neg_times))
  }
  
  # Check censoring rate
  event_rate <- mean(data[[event_col]] == 1, na.rm = TRUE)
  if (event_rate < 0.05) {
    warning(sprintf("Very low event rate (%.1f%%). Model may be unstable.", event_rate * 100))
  }
  if (event_rate > 0.95) {
    warning(sprintf("Very high event rate (%.1f%%). Check data coding.", event_rate * 100))
  }
  
  message(sprintf("[✓] Data quality check passed. N=%d, Events=%.1f%%", 
                  nrow(data), event_rate * 100))
  invisible(NULL)
}

#' Create Step Directory for Pipeline
#' @param base_dir Base output directory
#' @param step_num Step number
#' @param step_name Step name
#' @return Full path to created directory
#' @keywords internal
create_step_dir <- function(base_dir, step_num, step_name) {
  dir_name <- sprintf("Step%02d_%s", step_num, step_name)
  full_path <- file.path(base_dir, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
  }
  message(sprintf("  → Step %d: %s", step_num, step_name))
  return(full_path)
}


# ==============================================================================
# 13. Multi-Dataset Validation & Comparison Module
# ==============================================================================

#' Evaluate Model on a New Dataset (Validation/Test)
#' @param learner Trained mlr3 learner
#' @param test_data New data frame (validation/test set)
#' @param task_ref Reference task to ensure column consistency
#' @return Prediction object
#' @export
surv_predict_on_validation <- function(learner, test_data, task_ref) {
  # Ensure the test data has the same structure/encoding as training data
  test_task <- surv_create_surv_task(
    data = test_data, 
    time_col = task_ref$target_names[1], 
    event_col = task_ref$target_names[2],
    id = "validation_set"
  )
  
  # Perform prediction
  pred <- learner$predict(test_task)
  return(list(task = test_task, prediction = pred))
}

#' Plot Comparison Time-dependent AUC (Train vs Validation)
#' @param learner Trained learner
#' @param train_task Training task
#' @param val_task Validation task
#' @return ggplot object
#' @export
surv_plot_comparison_auc <- function(learner, train_task, val_task) {
  if (!requireNamespace("risksetROC", quietly = TRUE)) stop("Please install 'risksetROC'")
  
  # Calculate for Training
  pred_train <- learner$predict(train_task)
  auc_train <- risksetROC::risksetAUC(
    Stime = train_task$data()[[train_task$target_names[1]]],
    status = train_task$data()[[train_task$target_names[2]]],
    marker = pred_train$crank,
    tmax = max(train_task$data()[[train_task$target_names[1]]]) * 0.9,
    plot = FALSE
  )
  
  # Calculate for Validation
  pred_val <- learner$predict(val_task)
  auc_val <- risksetROC::risksetAUC(
    Stime = val_task$data()[[val_task$target_names[1]]],
    status = val_task$data()[[val_task$target_names[2]]],
    marker = pred_val$crank,
    tmax = max(val_task$data()[[val_task$target_names[1]]]) * 0.9,
    plot = FALSE
  )
  
  # Combine Data
  df_train <- data.frame(Time = auc_train$utimes, AUC = auc_train$AUC, Group = "Training Set")
  df_val <- data.frame(Time = auc_val$utimes, AUC = auc_val$AUC, Group = "Validation Set")
  df_all <- rbind(df_train, df_val)
  
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = Time, y = AUC, color = Group)) +
    ggplot2::geom_step(size = 1) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
    ggplot2::scale_color_manual(values = c("Training Set" = "#2980b9", "Validation Set" = "#e74c3c")) +
    ggplot2::labs(title = "Time-Dependent AUC Comparison", subtitle = learner$id) +
    ggprism::theme_prism()
  
  print(p)
  return(p)
}

#' Multi‑strategy feature selection for survival analysis
#'
#' Applies multiple feature selection methods to a survival task and returns a
#' consensus set of features. Supports 12+ algorithms including univariate Cox,
#' penalized Cox (LASSO, Ridge, Elastic Net), random forest importance,
#' XGBoost importance, VIMP, stepwise, stability selection, mRMR, and Boruta.
#'
#' @param object A \code{PrognosiX} object or a \code{TaskSurv} (mlr3 survival task).
#' @param methods Character vector of method names to apply. Available methods:
#'   \describe{
#'     \item{\code{"uni_cox"}}{Univariate Cox regression; keeps features with p < \code{p_threshold}.}
#'     \item{\code{"lasso"}}{LASSO penalized Cox (alpha = 1, lambda.min).}
#'     \item{\code{"ridge"}}{Ridge penalized Cox (alpha = 0, lambda.min).}
#'     \item{\code{"enet"}}{Elastic net (alpha = 0.5, lambda.min).}
#'     \item{\code{"rf_imp"}}{Random forest (ranger) permutation importance; keeps top \code{top_ratio} features.}
#'     \item{\code{"rfsrc_imp"}}{Random survival forest (randomForestSRC) importance; keeps top \code{top_ratio}.}
#'     \item{\code{"xgb_imp"}}{XGBoost gain importance; keeps top \code{top_ratio}.}
#'     \item{\code{"vimp"}}{VIMP variable importance from \code{randomForestSRC::vimp} (recommended).}
#'     \item{\code{"boruta"}}{Boruta wrapper algorithm (requires \code{Boruta} package; disabled by default).}
#'     \item{\code{"stepwise"}}{Stepwise Cox regression (both directions, AIC). Only for low‑dimensional data (p < 30).}
#'     \item{\code{"stab_sel"}}{Stability selection using \code{c060::stabpath} with Lasso.}
#'     \item{\code{"mrmr"}}{Minimum Redundancy Maximum Relevance using Cox risk score as proxy (approximate).}
#'   }
#' @param p_threshold Numeric. P‑value threshold for univariate Cox (default 0.05).
#' @param top_ratio Numeric. For importance‑based methods (RF, XGBoost, VIMP), keep this proportion of top features (default 0.5).
#' @param combine Character. How to combine results from different methods:
#'   \itemize{
#'     \item \code{"union"} – take union of all selected feature sets.
#'     \item \code{"intersection"} – take intersection (common features).
#'     \item \code{"freq"} – keep features selected by at least \code{freq_cutoff} methods.
#'   }
#' @param freq_cutoff Integer. Minimum number of methods that must select a feature when \code{combine = "freq"} (default 2).
#' @param verbose Logical. Print progress messages (default TRUE).
#' @param use_boruta Logical. Enable Boruta (can be slow and may fail on survival data). Default FALSE.
#'
#' @return A list with three components:
#'   \describe{
#'     \item{\code{selected}}{Character vector of finally selected feature names.}
#'     \item{\code{method_table}}{Data frame with rows = all features, columns = each method, indicating selection status (TRUE/FALSE).}
#'     \item{\code{method_results}}{List of raw outputs from each method (e.g., fitted models, importance vectors) for further inspection.}
#'   }
#'
#' @importFrom stats coef as.formula
#' @importFrom utils head
#' @importFrom survival Surv coxph
#' @importFrom MASS stepAIC
#' @importFrom mlr3 lrn
#' @importFrom mlr3proba TaskSurv
#' @importFrom randomForestSRC vimp
#' @importFrom c060 stabpath
#' @export
#' @examples
#' \dontrun{
#' # Load veteran data and create PrognosiX object
#' data("veteran", package = "survival")
#' stat <- CreateStatObject(raw.data = veteran, clean.data = veteran, group_col = "status")
#' prog <- Stat_to_PrognosiX(stat, time_col = "time", status_col = "status")
#'
#' # Run feature selection with 4 robust methods
#' fs <- surv_feature_selection_multi(
#'   object = prog,
#'   methods = c("uni_cox", "lasso", "rf_imp", "vimp"),
#'   p_threshold = 0.1,
#'   top_ratio = 0.5,
#'   combine = "freq",
#'   freq_cutoff = 2
#' )
#'
#' # View selected features
#' print(fs$selected)
#' # View which methods selected each feature
#' head(fs$method_table)
#' }
surv_feature_selection_multi <- function(object,
                                         methods = c("uni_cox", "lasso", "rf_imp", "vimp"),
                                         p_threshold = 0.05,
                                         top_ratio = 0.5,
                                         combine = c("union", "intersection", "freq"),
                                         freq_cutoff = 2,
                                         verbose = TRUE,
                                         use_boruta = FALSE) {
  
  combine <- match.arg(combine)
  if (!requireNamespace("mlr3", quietly = TRUE))
    stop("Package 'mlr3' required.")
  
  # Extract task and data
  if (inherits(object, "PrognosiX")) {
    task <- surv_extract_task(object)
  } else if (inherits(object, "TaskSurv")) {
    task <- object
  } else {
    stop("object must be a PrognosiX or TaskSurv object.")
  }
  
  data <- as.data.frame(task$data())
  time_var <- task$target_names[1]
  status_var <- task$target_names[2]
  all_features <- task$feature_names
  
  if (verbose) cat("\n[Multi‑Feature Selection] Methods:", paste(methods, collapse=", "), "\n")
  
  # Helper: keep top proportion of features based on importance vector
  keep_top <- function(imp_vec, ratio) {
    if (length(imp_vec) == 0) return(character(0))
    n_keep <- max(1, floor(length(imp_vec) * ratio))
    names(sort(imp_vec, decreasing = TRUE)[1:n_keep])
  }
  
  selection_list <- list()
  raw_results <- list()
  
  # ----------------------------------------------------------------------------
  # 1. Univariate Cox
  # ----------------------------------------------------------------------------
  if ("uni_cox" %in% methods) {
    if (verbose) cat("  - Running univariate Cox (p <", p_threshold, ")...\n")
    uni <- tryCatch({
      surv_filter_features_clinical(object, p_threshold = p_threshold)
    }, error = function(e) NULL)
    if (!is.null(uni)) {
      selection_list$uni_cox <- uni$task$feature_names
      raw_results$uni_cox <- uni$table
    } else {
      selection_list$uni_cox <- character(0)
    }
  }
  
  # ----------------------------------------------------------------------------
  # 2-4. glmnet family: LASSO, Ridge, Elastic Net
  # ----------------------------------------------------------------------------
  if ("lasso" %in% methods || "ridge" %in% methods || "enet" %in% methods) {
    if (verbose) cat("  - Running glmnet family (LASSO, Ridge, Elastic Net)...\n")
    
    # Helper to extract non-zero coefficients from glmnet model
    extract_glmnet_features <- function(lrn_obj, alpha_val) {
      selected <- character(0)
      if (!is.null(lrn_obj)) {
        coef_mat <- tryCatch(as.matrix(stats::coef(lrn_obj$model, s = "lambda.min")), error = function(e) NULL)
        if (!is.null(coef_mat)) {
          selected <- rownames(coef_mat)[abs(coef_mat[,1]) > 1e-6]
          selected <- setdiff(selected, "(Intercept)")
        }
      }
      return(selected)
    }
    
    if ("lasso" %in% methods) {
      lasso_lrn <- tryCatch(lrn("surv.cv_glmnet", alpha = 1, s = "lambda.min")$train(task), error = function(e) NULL)
      selection_list$lasso <- extract_glmnet_features(lasso_lrn, 1)
      raw_results$lasso <- lasso_lrn
    }
    
    if ("ridge" %in% methods) {
      ridge_lrn <- tryCatch(lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")$train(task), error = function(e) NULL)
      selection_list$ridge <- extract_glmnet_features(ridge_lrn, 0)
      raw_results$ridge <- ridge_lrn
    }
    
    if ("enet" %in% methods) {
      enet_lrn <- tryCatch(lrn("surv.cv_glmnet", alpha = 0.5, s = "lambda.min")$train(task), error = function(e) NULL)
      selection_list$enet <- extract_glmnet_features(enet_lrn, 0.5)
      raw_results$enet <- enet_lrn
    }
  }
  
  # ----------------------------------------------------------------------------
  # 5. Random Forest (ranger) permutation importance
  # ----------------------------------------------------------------------------
  if ("rf_imp" %in% methods) {
    if (verbose) cat("  - Running Random Forest (ranger) importance...\n")
    rf_lrn <- tryCatch({
      lrn("surv.ranger", importance = "permutation")$train(task)
    }, error = function(e) NULL)
    selected <- character(0)
    if (!is.null(rf_lrn) && !is.null(rf_lrn$importance())) {
      imp <- rf_lrn$importance()
      selected <- keep_top(imp, top_ratio)
    }
    selection_list$rf_imp <- selected
    raw_results$rf_imp <- if(!is.null(rf_lrn)) rf_lrn$importance() else NULL
  }
  
  # ----------------------------------------------------------------------------
  # 6. Random Survival Forest (rfsrc) VIMP importance
  # ----------------------------------------------------------------------------
  if ("rfsrc_imp" %in% methods) {
    if (verbose) cat("  - Running Random Survival Forest (rfsrc) importance...\n")
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
      warning("Package 'randomForestSRC' not installed. Skipping rfsrc_imp.")
    } else {
      rfsrc_lrn <- tryCatch({
        lrn("surv.rfsrc", importance = "permute")$train(task)
      }, error = function(e) NULL)
      selected <- character(0)
      if (!is.null(rfsrc_lrn) && !is.null(rfsrc_lrn$importance())) {
        imp <- rfsrc_lrn$importance()
        selected <- keep_top(imp, top_ratio)
      }
      selection_list$rfsrc_imp <- selected
      raw_results$rfsrc_imp <- if(!is.null(rfsrc_lrn)) rfsrc_lrn$importance() else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 7. XGBoost gain importance
  # ----------------------------------------------------------------------------
  if ("xgb_imp" %in% methods) {
    if (verbose) cat("  - Running XGBoost importance (gain)...\n")
    xgb_lrn <- tryCatch({
      lrn("surv.xgboost.cox", nrounds = 50)$train(task)
    }, error = function(e) NULL)
    selected <- character(0)
    if (!is.null(xgb_lrn) && !is.null(xgb_lrn$importance())) {
      imp_df <- xgb_lrn$importance()
      if (nrow(imp_df) > 0) {
        imp_vec <- setNames(imp_df$Gain, imp_df$Feature)
        selected <- keep_top(imp_vec, top_ratio)
      }
    }
    selection_list$xgb_imp <- selected
    raw_results$xgb_imp <- if(!is.null(xgb_lrn)) xgb_lrn$importance() else NULL
  }
  
  # ----------------------------------------------------------------------------
  # 8. VIMP (Variable Importance) from randomForestSRC — Standalone algorithm
  # ----------------------------------------------------------------------------
  if ("vimp" %in% methods) {
    if (verbose) cat("  - Running VIMP variable importance...\n")
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
      warning("Package 'randomForestSRC' not installed. Skipping vimp.")
    } else {
      vimp_fit <- tryCatch({
        randomForestSRC::vimp(
          formula = as.formula(paste("Surv(", time_var, ",", status_var, ") ~ .")),
          data = data[, c(all_features, time_var, status_var)],
          importance = "permute"
        )
      }, error = function(e) NULL)
      selected <- character(0)
      if (!is.null(vimp_fit)) {
        imp <- vimp_fit$importance
        selected <- keep_top(imp, top_ratio)
      }
      selection_list$vimp <- selected
      raw_results$vimp <- if(!is.null(vimp_fit)) vimp_fit else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 9. Boruta (optional, may fail on survival data)
  # ----------------------------------------------------------------------------
  if (use_boruta && "boruta" %in% methods) {
    if (verbose) cat("  - Running Boruta (optional, may be slow)...\n")
    if (!requireNamespace("Boruta", quietly = TRUE)) {
      warning("Package 'Boruta' not installed. Skipping boruta.")
    } else {
      selected <- character(0)
      tryCatch({
        cox_lrn <- lrn("surv.coxph")$train(task)
        risk <- cox_lrn$predict(task)$crank
        boruta_data <- data[, all_features, drop = FALSE]
        boruta_data$.risk <- risk
        boruta_fit <- Boruta::Boruta(.risk ~ ., data = boruta_data, doTrace = 0)
        selected <- names(boruta_fit$finalDecision[boruta_fit$finalDecision == "Confirmed"])
      }, error = function(e) {
        warning("Boruta failed: ", e$message)
      })
      selection_list$boruta <- selected
      raw_results$boruta <- if(exists("boruta_fit")) boruta_fit else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 10. Stepwise Cox (low-dimensional only)
  # ----------------------------------------------------------------------------
  if ("stepwise" %in% methods) {
    if (verbose) cat("  - Running stepwise Cox (direction = both, AIC)...\n")
    selected <- character(0)
    tryCatch({
      # Prevent high-dimensional explosion: take top 20 features from uni_cox filter first
      pre_feats <- if (length(selection_list$uni_cox) > 0) {
        head(selection_list$uni_cox, min(20, length(selection_list$uni_cox)))
      } else {
        head(all_features, min(20, length(all_features)))
      }
      if (length(pre_feats) > 0) {
        full_form <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~", 
                                      paste(pre_feats, collapse = " + ")))
        cox_full <- survival::coxph(full_form, data = data)
        step_mod <- MASS::stepAIC(cox_full, direction = "both", trace = 0)
        selected <- names(stats::coef(step_mod))
      }
    }, error = function(e) {
      warning("Stepwise failed: ", e$message)
    })
    selection_list$stepwise <- selected
    raw_results$stepwise <- if(exists("step_mod")) step_mod else NULL
  }
  
  # ----------------------------------------------------------------------------
  # 11. Stability Selection (via c060::stabpath)
  # ----------------------------------------------------------------------------
  if ("stab_sel" %in% methods) {
    if (verbose) cat("  - Running stability selection (c060::stabpath)...\n")
    selected <- character(0)
    if (!requireNamespace("c060", quietly = TRUE)) {
      warning("Package 'c060' not installed. Skipping stab_sel.")
    } else {
      tryCatch({
        x_mat <- as.matrix(data[, all_features, drop = FALSE])
        y_mat <- survival::Surv(data[[time_var]], data[[status_var]])
        stab_path <- c060::stabpath(y = y_mat, x = x_mat, steps = 50)
        selection_prob <- stab_path$stabpath
        if (!is.null(selection_prob) && ncol(selection_prob) > 0) {
          prob_avg <- colMeans(selection_prob, na.rm = TRUE)
          selected <- names(prob_avg)[prob_avg > 0.6]
        }
      }, error = function(e) {
        warning("Stability selection failed: ", e$message)
      })
      selection_list$stab_sel <- selected
      raw_results$stab_sel <- if(exists("stab_path")) stab_path else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 12. mRMR (using Cox risk score as proxy)
  # ----------------------------------------------------------------------------
  if ("mrmr" %in% methods) {
    if (verbose) cat("  - Running mRMR (praznik with Cox risk proxy)...\n")
    selected <- character(0)
    if (!requireNamespace("praznik", quietly = TRUE)) {
      warning("Package 'praznik' not installed. Skipping mrmr.")
    } else {
      tryCatch({
        cox_lrn <- lrn("surv.coxph")$train(task)
        risk <- cox_lrn$predict(task)$crank
        mrmr_res <- praznik::MRMR(X = data[, all_features, drop = FALSE], 
                                  Y = risk, k = min(20, length(all_features)))
        selected <- mrmr_res$selection
      }, error = function(e) {
        warning("mRMR failed: ", e$message)
      })
      selection_list$mrmr <- selected
      raw_results$mrmr <- if(exists("mrmr_res")) mrmr_res else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 13. Combine results
  # ----------------------------------------------------------------------------
  method_names <- names(selection_list)
  if (length(method_names) == 0) {
    warning("No method succeeded. Returning all features.")
    selected_final <- all_features
  } else {
    method_table <- data.frame(Feature = all_features)
    for (m in method_names) {
      method_table[[m]] <- all_features %in% selection_list[[m]]
    }
    
    if (combine == "union") {
      selected_final <- unique(unlist(selection_list))
    } else if (combine == "intersection") {
      selected_final <- Reduce(intersect, selection_list)
    } else if (combine == "freq") {
      counts <- rowSums(method_table[, -1])
      selected_final <- all_features[counts >= freq_cutoff]
    }
    
    if (verbose) {
      cat(sprintf("\n[Combined] %s selection: %d features out of %d\n",
                  toupper(combine), length(selected_final), length(all_features)))
    }
  }
  
  return(list(
    selected = selected_final,
    method_table = method_table,
    method_results = raw_results
  ))
}

#' Decision Curve Analysis for One or More Survival Models
#'
#' Computes and plots decision curves at a specified time point for one or more
#' survival models using standard Kaplan-Meier corrections via the 'dcurves' package.
#' Supports highly flexible aesthetic configurations for colors and linetypes.
#'
#' @param learners A **named list** of trained `mlr3` survival learners.
#'   Each element must be a learner that supports `"distr"` predictions.
#'   Example: `list("Ranger" = learner_ranger)`.
#' @param object A `TaskSurv` or `PrognosiX` object containing the validation data.
#' @param eval_time Numeric. The time point at which to evaluate event probabilities.
#' @param thresholds Numeric vector. Risk thresholds (probabilities of event) at
#'   which net benefit is calculated. Defaults to `seq(0.01, 0.99, length.out = 50)`.
#' @param colors Character vector or single string. Can be:
#'   - A built-in palette keyword: `"default"`, `"clinical"`, `"vibrant"`, or `"jama"`.
#'   - A fully/partially named vector, e.g. `c(Ranger = "#2c7fb8")`. Missing reference
#'     strategies (`TreatAll`, `TreatNone`) will be filled automatically with high-contrast distinct colors.
#' @param linetypes Named character vector of line types for the strategies.
#'   Can be partially named (e.g., `c(Ranger = "solid")`); missing reference styles
#'   will default to distinct `"dashed"` and `"dotted"` profiles automatically.
#' @param include_reference Logical. Should the "Treat All" and "Treat None" curves be added? Default is `TRUE`.
#' @param ylim Numeric vector of length 2. Y-axis limits for net benefit. Defaults to `c(-0.05, NA)`.
#' @param title Character. Custom plot title.
#' @param subtitle Character. Custom plot subtitle.
#' @param print_stats Logical. Should summary statistics be printed to the console? Default `TRUE`.
#' @param clin_range Numeric vector of length 2. Default is `c(0.05, 0.5)`.
#'
#' @return A list with three components: `plot`, `table`, and `summary`.
#'
#' @importFrom survival Surv
#' @importFrom dcurves dca
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_line labs scale_color_manual scale_linetype_manual coord_cartesian theme_minimal
#' @importFrom tidyr pivot_wider
#' @export
#' @examples
#' \dontrun{
#' # Requires trained learners with distr predict_type and task
#' # dca_result <- plot_dca_survival(list("Cox" = cox_learner), task, eval_time = 365)
#' }
plot_dca_survival <- function(learners,
                                object,
                                eval_time,
                                thresholds = seq(0.01, 0.99, length.out = 50),
                                colors = NULL,
                                linetypes = NULL,
                                include_reference = TRUE,
                                ylim = c(-0.05, NA),
                                title = NULL,
                                subtitle = NULL,
                                print_stats = TRUE,
                                clin_range = c(0.05, 0.5)) {
    
    # ---- 1. Input validation and Task extraction  ------------------------------
    if (!is.list(learners) || is.null(names(learners))) {
      stop("'learners' must be a named list (e.g., list(Model1 = lrn1, Model2 = lrn2))")
    }
    if (length(learners) == 0) stop("At least one learner is required.")
    if (!is.numeric(eval_time) || length(eval_time) != 1) stop("'eval_time' must be a single number.")
    if (!is.numeric(clin_range) || length(clin_range) != 2 || clin_range[1] >= clin_range[2]) {
      stop("'clin_range' must be a numeric vector of length 2 with min < max.")
    }
    
    if (inherits(object, "TaskSurv")) {
      task <- object
    } else if (inherits(object, "PrognosiX")) {
      task <- surv_extract_task(object)
    } else {
      stop("object must be a TaskSurv or PrognosiX object.")
    }
    
    data_df <- as.data.frame(task$data())
    time_var <- task$target_names[1]
    status_var <- task$target_names[2]
    
    # ---- 2. Safely extract predicted event probability from mlr3 learners [1 - S(t)] ---------
    dca_data <- data_df[, c(time_var, status_var), drop = FALSE]
    
    for (nm in names(learners)) {
      lrn <- learners[[nm]]
      if (!"distr" %in% lrn$predict_types) {
        stop(sprintf("Learner '%s' does not support 'distr' predictions.", lrn$id))
      }
      lrn$predict_type <- "distr"
      pred <- lrn$predict(task)
      
      surv_prob <- as.numeric(pred$distr$survival(eval_time))
      dca_data[[nm]] <- 1 - surv_prob
    }
    
    # ---- 3. Standard survival DCA calculation using dcurves -----------------------
    formula_str <- sprintf("survival::Surv(%s, %s) ~ %s", 
                           time_var, status_var, 
                           paste(names(learners), collapse = " + "))
    
    dca_obj <- dcurves::dca(
      formula = as.formula(formula_str),
      data = dca_data,
      time = eval_time,
      thresholds = thresholds
    )
    
    nb_table <- as.data.frame(tibble::as_tibble(dca_obj))
    
    colnames(nb_table)[colnames(nb_table) == "variable"]    <- "Strategy"
    colnames(nb_table)[colnames(nb_table) == "threshold"]   <- "Threshold"
    colnames(nb_table)[colnames(nb_table) == "net_benefit"] <- "NetBenefit"
    
    if (!include_reference) {
      nb_table <- nb_table[!nb_table$Strategy %in% c("all", "none"), ]
    } else {
      nb_table$Strategy <- ifelse(nb_table$Strategy == "all", "TreatAll",
                                  ifelse(nb_table$Strategy == "none", "TreatNone", nb_table$Strategy))
    }
    
    all_strategies <- unique(nb_table$Strategy)
    model_names <- names(learners)
    
    # ---- 4. Automated aesthetic settings (high-contrast colors and line types) ----------------------------
    default_ref_colors <- c("TreatAll" = "#E69F00", "TreatNone" = "#000000") 
    default_ref_lts    <- c("TreatAll" = "dashed", "TreatNone" = "dotted")
    
    if (is.null(colors)) {
      hues <- seq(15, 375, length.out = length(model_names) + 1)
      mod_cols <- hcl(hues[1:length(model_names)], l = 55, c = 90)
      names(mod_cols) <- model_names
      colors <- c(mod_cols, default_ref_colors)
    } else if (is.character(colors) && length(colors) == 1) {
      pal_choice <- tolower(colors)
      if (pal_choice == "clinical") {
        mod_cols <- c("#00A087FF", "#3C5488FF", "#4DBBD5FF")[1:length(model_names)]
        names(mod_cols) <- model_names
        colors <- c(mod_cols, "TreatAll" = "#E64B35FF", "TreatNone" = "#111111")
      } else if (pal_choice == "vibrant") {
        mod_cols <- c("#0073C2FF", "#CD534CFF", "#7AA6C2FF")[1:length(model_names)]
        names(mod_cols) <- model_names
        colors <- c(mod_cols, "TreatAll" = "#EFC000FF", "TreatNone" = "#000000")
      } else if (pal_choice == "jama") {
        mod_cols <- c("#374E55FF", "#DF8F44FF", "#00A1D5FF")[1:length(model_names)]
        names(mod_cols) <- model_names
        colors <- c(mod_cols, "TreatAll" = "#B24745FF", "TreatNone" = "#79AF97FF")
      } else {
        hues <- seq(15, 375, length.out = length(model_names) + 1)
        mod_cols <- hcl(hues[1:length(model_names)], l = 55, c = 90)
        names(mod_cols) <- model_names
        colors <- c(mod_cols, default_ref_colors)
      }
    } else {
      missing_refs <- setdiff(c("TreatAll", "TreatNone"), names(colors))
      if (length(missing_refs) > 0) {
        colors <- c(colors, default_ref_colors[missing_refs])
      }
    }
    
    if (is.null(linetypes)) {
      lt <- rep("solid", length(all_strategies))
      names(lt) <- all_strategies
      if ("TreatAll" %in% all_strategies) lt["TreatAll"] <- "dashed"
      if ("TreatNone" %in% all_strategies) lt["TreatNone"] <- "dotted"
      linetypes <- lt
    } else {
      missing_lts <- setdiff(c("TreatAll", "TreatNone"), names(linetypes))
      if (length(missing_lts) > 0) {
        linetypes <- c(linetypes, default_ref_lts[missing_lts])
      }
      unmapped_models <- setdiff(model_names, names(linetypes))
      if (length(unmapped_models) > 0) {
        extra_lts <- rep("solid", length(unmapped_models))
        names(extra_lts) <- unmapped_models
        linetypes <- c(linetypes, extra_lts)
      }
    }
    
    # ---- 5. Build ggplot (using coord_cartesian to avoid data clipping warnings) -----------------
    if (is.null(title)) title <- "Decision Curve Analysis"
    if (is.null(subtitle)) subtitle <- paste("Time =", eval_time, "| Validation set")
    
    p <- ggplot2::ggplot(nb_table, ggplot2::aes(x = Threshold, y = NetBenefit, 
                                                color = Strategy, linetype = Strategy)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_linetype_manual(values = linetypes) +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = "Risk Threshold (Probability of Event)",
                    y = "Net Benefit",
                    color = "Strategy", linetype = "Strategy") +
      ggplot2::coord_cartesian(ylim = ylim) # Core fix: avoid warnings and preserve full curve continuity
    
    if (requireNamespace("ggprism", quietly = TRUE)) {
      p <- p + ggprism::theme_prism()
    } else {
      p <- p + ggplot2::theme_minimal()
    }
    
    # ---- 6. Metric evaluation within clinically relevant range (smartly adapts to user-defined risk threshold) ---------------------
    max_avail_thr <- max(nb_table$Threshold, na.rm = TRUE)
    if (clin_range[2] > max_avail_thr) {
      clin_range[2] <- max_avail_thr # Auto-align evaluation upper bound to user input max threshold, prevent out-of-bounds
    }
    
    summary_stats <- data.frame()
    
    for (nm in model_names) {
      model_nb <- nb_table[nb_table$Strategy == nm, ]
      idx <- which(model_nb$Threshold >= clin_range[1] & model_nb$Threshold <= clin_range[2])
      if (length(idx) == 0) idx <- seq_len(nrow(model_nb))
      
      nb_clin <- model_nb$NetBenefit[idx]
      thr_clin <- model_nb$Threshold[idx]
      
      max_nb   <- max(nb_clin, na.rm = TRUE)
      best_thr <- thr_clin[which.max(nb_clin)]
      avg_nb   <- mean(nb_clin, na.rm = TRUE)
      
      treat_all_clin <- nb_table[nb_table$Strategy == "TreatAll" & nb_table$Threshold %in% thr_clin, "NetBenefit"]
      avg_gain <- if(length(treat_all_clin) == length(nb_clin)) mean(nb_clin - treat_all_clin, na.rm = TRUE) else NA
      
      summary_stats <- rbind(summary_stats, data.frame(
        Model = nm,
        ClinRange_Min = clin_range[1],
        ClinRange_Max = clin_range[2],
        Max_NetBenefit = round(max_nb, 4),
        Threshold_at_Max = round(best_thr, 4),
        Avg_NetBenefit = round(avg_nb, 4),
        Avg_NetBenefit_Gain = round(avg_gain, 4)  # Removed invalid and non-standard AUC_NetBenefit metric
      ))
    }
    
    if (print_stats) {
      cat("\n========== DCA Summary (Clinical Range: [", clin_range[1], ", ", clin_range[2], "]) ==========\n", sep = "")
      cat("  Calculated using standard Kaplan-Meier survival adjustments via 'dcurves'.\n")
      print(summary_stats, row.names = FALSE)
      cat("========================================================================\n")
    }
    
    # ---- 7. Matrix form output (Core fix: remove unique attribute columns, ensure wide table perfect alignment) -------
    nb_table_core <- nb_table[, c("Threshold", "Strategy", "NetBenefit")]
    dca_wide <- tidyr::pivot_wider(nb_table_core, names_from = Strategy, values_from = NetBenefit)
    
    invisible(list(
      plot = p,
      table = as.data.frame(dca_wide),
      summary = summary_stats
    ))
  }
                           
#' List all available feature selection methods in surv_feature_selection_multi
#'
#' @param verbose If TRUE, prints the table to console. If FALSE, returns the data frame.
#' @return A data frame with columns: Method, Description, RequiredPackages, Recommendation.
#' @export
#'
#' @examples
#' list_surv_feature_methods()
list_surv_feature_methods <- function(verbose = TRUE) {
  methods_df <- data.frame(
    Method = c(
      "uni_cox", "lasso", "ridge", "enet", "rf_imp", "rfsrc_imp",
      "xgb_imp", "vimp", "boruta", "stepwise", "stab_sel", "mrmr"
    ),
    Description = c(
      "Univariate Cox regression (p < threshold)",
      "LASSO penalized Cox (lambda.min)",
      "Ridge penalized Cox",
      "Elastic net (alpha = 0.5, cross-validated)",
      "Random forest (ranger) permutation importance, keep top ratio",
      "Random survival forest (randomForestSRC) importance, keep top ratio",
      "XGBoost gain importance, keep top ratio",
      "VIMP variable importance from randomForestSRC (robust, recommended)",
      "Boruta wrapper algorithm (default OFF, may fail on survival data)",
      "Stepwise Cox regression (both directions, AIC) – low-dimensional only",
      "Stability selection with Lasso via c060::stabpath",
      "Minimum Redundancy Maximum Relevance (Cox risk proxy) – approximate"
    ),
    RequiredPackages = c(
      "survival (built-in)",
      "glmnet (via mlr3learners)",
      "glmnet",
      "glmnet",
      "ranger (mlr3learners)",
      "randomForestSRC",
      "xgboost (mlr3extralearners)",
      "randomForestSRC",
      "Boruta, randomForest (optional)",
      "MASS",
      "c060",
      "praznik"
    ),
    Recommendation = c(
      "★ ★ ★ ★ ★ (must-have)",
      "★ ★ ★ ★ ★ (top choice)",
      "★ ★ ★ ★ ☆ (high collinearity)",
      "★ ★ ★ ★ ★ (often best glmnet)",
      "★ ★ ★ ★ ☆ (nonlinear effects)",
      "★ ★ ★ ★ ☆ (survival-specialized)",
      "★ ★ ★ ★ ☆ (handles missing data)",
      "★ ★ ★ ★ ★ (stable, official RF method)",
      "★ ★ ☆ ☆ ☆ (use with caution, OFF by default)",
      "★ ★ ☆ ☆ ☆ (only for low dimension, p < 30)",
      "★ ★ ★ ★ ☆ (robust for high-dim)",
      "★ ★ ★ ☆ ☆ (approximate, informative only)"
    ),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("\n============================================================\n")
    cat("Survival Feature Selection Methods (updated)\n")
    cat("============================================================\n\n")
    print(methods_df, row.names = FALSE)
    cat("\n[Note] Use combine = 'union', 'intersection', or 'freq'.\n")
    cat("Recommended production set: c('uni_cox', 'lasso', 'vimp', 'rf_imp')\n")
  }
  
  invisible(methods_df)
}

# Short alias
print_surv_feature_methods <- function() {
  list_surv_feature_methods(verbose = TRUE)
}

# =========================================================================
# Internal Helpers
# =========================================================================
#' Convert SurvSHAP Result to Long Format
#'
#' Internal helper function to convert the output from \code{survex} SHAP
#' calculations into a tidy long-format data frame.
#'
#' @param shap_obj The SHAP object returned by \code{survex}.
#' @param obs_idx The indices of observations explained.
#' @param features The original feature data frame.
#'
#' @return A data frame with columns: \code{time}, \code{feature}, \code{shap_value},
#'   and \code{observation}.
#'
#' @keywords internal
#' @noRd
.survshap_to_long <- function(shap_obj, obs_idx, features) {
  result <- shap_obj$result
  times  <- shap_obj$eval_times
  var_vals <- features[obs_idx, , drop = FALSE]
  obs_names <- rownames(var_vals)
  if (is.null(obs_names)) obs_names <- as.character(obs_idx)
  
  if (is.data.frame(result)) {
    df_long <- .df_to_long(result, times, obs_names[1L])
    return(df_long)
  }
  
  if (is.list(result)) {
    n <- length(result)
    long_list <- lapply(seq_len(n), function(i) {
      .df_to_long(as.data.frame(result[[i]]), times, obs_names[i])
    })
    return(do.call(rbind, long_list))
  }
  
  stop("Cannot parse SurvSHAP result. Unexpected structure: ", class(result)[1])
}

.df_to_long <- function(df, times, obs_label) {
  df <- df[, !names(df) %in% c("times", "time", "_times_", "id", "B"),
           drop = FALSE]
  
  n_rows <- nrow(df)
  t_vec <- if (!is.null(times) && length(times) == n_rows) {
    times
  } else {
    rn <- rownames(df)
    t_parsed <- suppressWarnings(as.numeric(sub("^t=", "", rn)))
    if (all(!is.na(t_parsed))) t_parsed else seq_len(n_rows)
  }
  
  df_long <- tidyr::pivot_longer(
    cbind(data.frame(.time = t_vec, stringsAsFactors = FALSE), df),
    cols = -".time", names_to = "feature", values_to = "shap_value"
  )
  names(df_long)[names(df_long) == ".time"] <- "time"
  df_long$observation <- obs_label
  as.data.frame(df_long)
}

.generate_shap_plots <- function(shap_long, n_top, bar_col, model_id, type) {
  feat_imp <- shap_long %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice(seq_len(min(n_top, dplyr::n())))
  
  feat_imp$feature <- factor(feat_imp$feature, levels = rev(feat_imp$feature))
  
  n_obs_used <- length(unique(shap_long$observation))
  title_text <- sprintf("SurvSHAP Importance (%s mode, n=%d): %s",
                        toupper(type), n_obs_used, model_id)
  
  p_bar <- ggplot2::ggplot(feat_imp, ggplot2::aes(x = feature, y = mean_abs)) +
    ggplot2::geom_bar(stat = "identity", fill = bar_col, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title_text, x = "Feature", y = "Mean |SHAP|") +
    ggprism::theme_prism()
  
  if (!all(is.na(shap_long$time))) {
    top_feats  <- feat_imp$feature
    line_data  <- shap_long %>%
      dplyr::filter(feature %in% top_feats) %>%
      dplyr::group_by(feature, time) %>%
      dplyr::summarise(mean_shap = mean(shap_value, na.rm = TRUE),
                       .groups = "drop")
    
    p_line <- ggplot2::ggplot(line_data,
                              ggplot2::aes(x = time, y = mean_shap, color = feature)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::labs(title = "SHAP Dynamics over Time",
                    x = "Time", y = "Average SHAP Value") +
      ggprism::theme_prism()
  } else {
    p_line <- NULL
  }
  
  list(bar_plot = p_bar, line_plot = p_line)
}

.prognosis_optional_packages <- c(
  "mlr3", "mlr3proba", "mlr3tuning", "mlr3learners",
  "mlr3extralearners", "survival", "tidyverse", "paradox", "data.table"
)

.check_prognosis_packages <- function() {
  missing <- .prognosis_optional_packages[
    !vapply(.prognosis_optional_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0) {
    stop(
      "Missing packages required for PrognosiX framework: ",
      paste(missing, collapse = ", "),
      ". Install them before using prognosis-related functions."
    )
  }
  invisible(TRUE)
}


# Unified pipeline to extract predictions and compute Kaplan-Meier observed rates per bin
.prepare_cal_data <- function(learner, object, time_point, n_bins) {
  if (inherits(object, "TaskSurv")) {
    task <- object
  } else if (inherits(object, "PrognosiX")) {
    task <- surv_extract_task(object)
  } else {
    stop("object must be a TaskSurv or PrognosiX object")
  }
  
  n_bins <- max(n_bins, 5L)
  learner$predict_type <- "distr"
  pred <- learner$predict(task)
  
  surv_prob <- .extract_surv_prob(pred$distr, time_point, task)
  if (is.null(surv_prob)) return(NULL)
  surv_prob <- as.numeric(surv_prob)
  
  if (length(surv_prob) != task$nrow) {
    warning("Length of survival probabilities does not match task rows. Calibration skipped.")
    return(NULL)
  }
  
  data_df <- as.data.frame(task$data())
  time    <- data_df[[task$target_names[1L]]]
  status  <- data_df[[task$target_names[2L]]]
  
  df <- data.frame(pred = surv_prob, time = time, status = status)
  df <- df[order(df$pred), ]
  
  breaks <- unique(quantile(df$pred, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE))
  if (length(breaks) < 3L) {
    warning("Not enough unique predicted probabilities to form bins. Calibration skipped.")
    return(NULL)
  }
  df$bin <- cut(df$pred, breaks = breaks, include.lowest = TRUE)
  
  obs_surv <- sapply(split(df, df$bin), function(bin_data) {
    if (nrow(bin_data) < 2L) return(NA_real_)
    km <- survival::survfit(survival::Surv(time, status) ~ 1, data = bin_data)
    sp <- summary(km, times = time_point, extend = TRUE)$surv
    if (length(sp) == 0L) NA_real_ else sp[[1L]]
  })
  
  bin_centers <- tapply(df$pred, df$bin, mean, na.rm = TRUE)
  cal_df <- na.omit(data.frame(predicted = as.numeric(bin_centers),
                               observed  = as.numeric(obs_surv)))
  
  if (nrow(cal_df) < 2L) {
    warning("Not enough valid bins for calibration (need at least 2).")
    return(NULL)
  }
  return(cal_df)
}

# Unified interface to calculate calibration metrics
.compute_cal_metrics <- function(cal_df) {
  metrics <- list()
  lm_fit            <- lm(observed ~ predicted, data = cal_df)
  metrics$slope     <- unname(stats::coef(lm_fit)[2L])
  metrics$intercept <- unname(stats::coef(lm_fit)[1L])
  metrics$r_squared <- summary(lm_fit)$r.squared
  
  errors            <- abs(cal_df$observed - cal_df$predicted)
  metrics$mae       <- mean(errors, na.rm = TRUE)
  metrics$ici       <- .compute_ici(cal_df)
  metrics$e50       <- unname(quantile(errors, 0.5, na.rm = TRUE))
  metrics$e90       <- unname(quantile(errors, 0.9, na.rm = TRUE))
  return(metrics)
}

# Extractor for predicted survival probabilities from mlr3 distribution objects
.extract_surv_prob <- function(distr, time_point, task) {
  if (inherits(distr, "Matdist") || is.environment(distr)) {
    if (!is.null(distr$survival)) {
      sp <- distr$survival(time_point)
      return(if (is.matrix(sp)) as.numeric(sp) else sp)
    } else if (!is.null(distr$cdf)) {
      sp <- 1 - distr$cdf(time_point)
      return(if (is.matrix(sp)) as.numeric(sp) else sp)
    }
    warning("No survival or cdf method in distr object.")
    return(NULL)
  } else if (is.matrix(distr) || is.array(distr)) {
    times <- attr(distr, "times")
    if (is.null(times)) {
      data_tmp <- as.data.frame(task$data())
      times <- sort(unique(data_tmp[[task$target_names[1L]]][data_tmp[[task$target_names[2L]]] == 1]))
    }
    t_idx <- which.min(abs(times - time_point))
    sp <- if (length(dim(distr)) == 2L) distr[, t_idx] else distr[, t_idx, 1L]
    return(as.numeric(sp))
  }
  warning("Unknown distr type")
  NULL
}

# Integrated Calibration Index (ICI) computation engine
.compute_ici <- function(cal_df) {
  if (nrow(cal_df) < 2L) return(NA_real_)
  tryCatch({
    if (nrow(cal_df) >= 4L) {
      lo    <- loess(observed ~ predicted, data = cal_df,
                     span = 1.0, degree = 1L, surface = "direct")
      x_seq <- seq(min(cal_df$predicted), max(cal_df$predicted), length.out = 200L)
      y_hat <- predict(lo, newdata = data.frame(predicted = x_seq))
      valid <- !is.na(y_hat)
      if (sum(valid) < 2L) return(NA_real_)
      x_v <- x_seq[valid]; d_v <- abs(y_hat[valid] - x_v)
    } else {
      x_v <- cal_df$predicted; d_v <- abs(cal_df$observed - x_v)
    }
    rng <- diff(range(x_v))
    if (rng < 1e-6) return(NA_real_)
    sum(diff(x_v) * (d_v[-length(d_v)] + d_v[-1L]) / 2L) / rng
  }, error = function(e) NA_real_)
}

