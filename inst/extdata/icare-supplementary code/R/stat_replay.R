# =============================================================================
# icare_stat_replay_helpers.R
# Used only by Part C (Module 1 / Stat object). Content is taken verbatim from
# the uploaded test-leakage.R:
#   build_stat()         fits a Stat object on the given samples
#                        (imputation / outlier handling / encoding / scaling)
#   replay_stat_data()   applies the fitted parameters to new samples as-is
#                        (no re-estimation)
#   transform_rows()     alias of replay_stat_data
# Note: these functions were written to work around a problem in the installed
#       icare's outlier replay path. Once process_new_data() is fixed in the
#       package, merge this logic into the package instead of keeping it in a
#       script.
# =============================================================================
build_stat <- function(raw_data, info_data) {
  if (!is.data.frame(raw_data) || !is.data.frame(info_data)) stop("Inputs must be data frames")
  if (!"group" %in% names(info_data)) stop("info_data must contain group")
  if (nrow(raw_data) != nrow(info_data)) stop("raw_data and info_data row counts differ")
  if (!setequal(rownames(raw_data), rownames(info_data))) stop("Sample IDs do not match")
  if (anyDuplicated(names(raw_data))) stop("Duplicate input feature names are not supported")
  info_data <- info_data[rownames(raw_data), , drop = FALSE]
  original_names <- names(raw_data)
  raw_data <- icare::modify_column_names(raw_data)
  so <- icare::CreateStatObject(raw.data = raw_data, info.data = info_data,
                               group_col = "group", na.action = "allow")
  so <- icare::stat_diagnose_variable_type(so)
  so <- icare::stat_convert_variables(so, save_data = FALSE)
  template <- so@clean.data
  so <- icare::stat_miss_processed(so, impute_method = "median_mode", save_data = FALSE)
  # These replacement values are fitted BEFORE training outliers are replaced,
  # exactly as handle_outliers(..., impute_value='median') does during fitting.
  before_outliers <- so@clean.data
  numeric_cols <- names(before_outliers)[vapply(before_outliers, is.numeric, logical(1))]
  outlier_fills <- vapply(before_outliers[numeric_cols], stats::median,
                          numeric(1), na.rm = TRUE)
  so <- icare::stat_detect_and_mark_outliers(so, method = "iqr", save_data = FALSE)
  so <- icare::stat_handle_outliers(so, method = "impute", impute_value = "median", save_data = FALSE)
  before_encode <- so@clean.data
  encode_vars <- icare::diagnose_variable_type(before_encode, group_col = "group")$vars_to_encode
  encodings <- lapply(before_encode[encode_vars], function(x) unique(x[!is.na(x)]))
  so <- icare::stat_onehot_encode(so, save_data = FALSE)
  so <- icare::stat_normalize_process(so, method = "auto", save_data = FALSE)
  so@process.info$replay_state <- list(
    version = 1L, original_names = original_names, columns = names(template),
    template = template[0, , drop = FALSE], outlier_fills = outlier_fills,
    encodings = encodings, output_columns = names(so@scale.data))
  so
}

replay_stat_data <- function(stat_object, new_data) {
  if (!inherits(stat_object, "Stat") || !is.data.frame(new_data))
    stop("Expected a Stat object and a new-data data frame")
  state <- stat_object@process.info$replay_state
  if (is.null(state) || !identical(state$version, 1L))
    stop("No replay state: rebuild this object using the repaired build_stat()")
  if (anyDuplicated(names(new_data))) stop("Duplicate new-data feature names")
  d <- new_data
  if (!setequal(names(d), state$columns)) {
    mapping <- stats::setNames(state$columns, state$original_names)
    mapped <- unname(mapping[names(d)])
    if (anyNA(mapped) || anyDuplicated(mapped) || !setequal(mapped, state$columns))
      stop("New-data columns do not match the training feature schema")
    names(d) <- mapped
  }
  d <- d[, state$columns, drop = FALSE]
  # Freeze variable types at training time rather than re-diagnosing a batch.
  for (col in names(state$template)) {
    target <- state$template[[col]]
    if (is.factor(target)) {
      unknown <- !is.na(d[[col]]) & !as.character(d[[col]]) %in% levels(target)
      if (any(unknown)) stop("Unseen category in ", col)
      d[[col]] <- factor(d[[col]], levels = levels(target), ordered = is.ordered(target))
    } else if (is.numeric(target)) {
      x <- if (is.numeric(d[[col]])) d[[col]] else
        suppressWarnings(as.numeric(as.character(d[[col]])))
      if (any(!is.na(d[[col]]) & !is.finite(x))) stop("Non-finite or nonnumeric input in ", col)
      d[[col]] <- x
    }
  }
  info <- stat_object@process.info
  fills <- info$missing_info$imputation_info$imputation_values
  for (col in names(fills)) {
    if (col %in% names(d) && anyNA(d[[col]])) {
      fill <- fills[[col]]$used_value
      if (length(fill) != 1L || is.na(fill)) stop("No fitted missing-value replacement for ", col)
      d[[col]][is.na(d[[col]])] <- fill
    }
  }
  detection <- info$outlier_detection
  if (!is.null(detection)) {
    for (col in intersect(names(detection$stats), names(state$outlier_fills))) {
      x <- d[[col]]; pars <- detection$stats[[col]]
      is_outlier <- switch(detection$method,
        iqr = x < pars$Q1 - detection$threshold * pars$IQR |
              x > pars$Q3 + detection$threshold * pars$IQR,
        zscore = abs((x - pars$mean) / pars$sd) > detection$threshold,
        stop("Unsupported fitted outlier detector: ", detection$method))
      if (any(is_outlier, na.rm = TRUE)) {
        fill <- state$outlier_fills[[col]]
        if (!is.finite(fill)) stop("No finite training outlier replacement for ", col)
        d[[col]][which(is_outlier)] <- fill
      }
    }
  }
  # Use the same training levels, order and column names as one_hot_encode().
  for (col in names(state$encodings)) {
    for (value in state$encodings[[col]])
      d[[paste(col, value, sep = "_")]] <- as.integer(d[[col]] == value)
    d[[col]] <- NULL
  }
  if (!setequal(names(d), state$output_columns)) stop("Encoded feature schema differs from training")
  d <- d[, state$output_columns, drop = FALSE]
  norm_cols <- names(info$normalization)
  if (length(norm_cols) && (is.null(info$normalization_params) ||
      !all(norm_cols %in% names(info$normalization_params)) ||
      any(vapply(info$normalization_params[norm_cols], is.null, logical(1)))))
    stop("Stored normalization parameters are incomplete; refit the training object")
  # Delegate only the stored normalization transforms to icare. This local copy
  # cannot re-enter its broken outlier path or mutate the caller's object.
  normalizer <- stat_object
  normalizer@process.info <- info[c("normalization", "normalization_params")]
  out <- as.data.frame(icare::process_new_data(normalizer, d, save_data = FALSE))
  if (!identical(rownames(out), rownames(new_data)) || nrow(out) != nrow(new_data))
    stop("Replay changed the input row count or order")
  if (!identical(names(out), state$output_columns)) stop("Replay changed the output feature schema")
  out
}

transform_rows <- function(so, rows) replay_stat_data(so, rows)
