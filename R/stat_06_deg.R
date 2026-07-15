
#' Compute log2 group means for logFC calculation
#'
#' Internal helper. Given a data frame and group column, returns log2 of per-group
#' means with pseudocount added only when zero means exist.
#'
#' @param mat   Data frame containing numeric columns and group column.
#' @param group_col  Name of the grouping column.
#' @return A tibble: n_groups rows, columns = numeric features, values = log2(mean).
#' @keywords internal
.compute_log2_means <- function(mat, group_col) {
  means <- mat %>%
    group_by_at(group_col) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::select(-dplyr::all_of(group_col))
  
  raw_vals <- unlist(means)
  if (any(raw_vals == 0, na.rm = TRUE)) {
    positive_vals <- raw_vals[raw_vals > 0]
    pseudocount   <- if (length(positive_vals) > 0) min(positive_vals, na.rm = TRUE) * 0.01 else 1e-6
    message("Zero group means detected — adding pseudocount (", round(pseudocount, 8), ") before log2 transform.")
    means <- log2(means + pseudocount)
  } else {
    means <- log2(means)
  }
  return(means)
}

#' Perform Batch Wilcoxon Test for Multiple Variables (Binary Classification)
#'
#' Runs pairwise Wilcoxon tests between two groups for every numeric variable.
#' The test itself always operates on `mat` (whatever data is passed in).
#' logFC is computed from `logfc_mat`, which is the appropriate data for ratio
#' calculation.
#'
#' @param mat        Data frame used for the Wilcoxon test (can be clean, scaled, or log-transformed).
#' @param group_col  Column name defining the two groups.
#' @param logfc_mat  Data frame used exclusively for logFC calculation.
#' @param logfc_type One of `"log2ratio"` (default) or `"diff"`.
#' @param p_threshold Numeric; only features with p < p_threshold are retained. Defaults to 0.05.
#' @param save_dir   Directory to save results.
#' @param save_data  Logical; whether to write CSV output.
#' @param csv_filename  Output filename.
#'
#' @returns Data frame with W, p, mean_x/y, median_x/y, sd_x/y, p.adjust, logFC, change.
#' @export
#'
#' @examples
#' \dontrun{
#' deg_res <- batch_Wilcoxon(stat_obj_test@clean.data, group_col = "SWAB", p_threshold = 0.05)
#' }
batch_Wilcoxon <- function(mat,
                           group_col    = "group",
                           logfc_mat    = NULL,
                           logfc_type   = "log2ratio",
                           p_threshold  = 0.05,
                           save_dir     = NULL,
                           save_data    = FALSE,
                           csv_filename = "last_test_sig.csv") {
  
  # Validate that exactly two groups exist
  group_levels <- unique(na.omit(as.character(mat[[group_col]])))
  if (length(group_levels) != 2) {
    stop("batch_Wilcoxon requires exactly 2 groups in '", group_col,
         "'. Found: ", paste(group_levels, collapse = ", "))
  }
  
  # ── 1. Wilcoxon test ──────────────────────────────────────────────────────
  test.fun <- function(dat, col) {
    index <- unique(dat[[group_col]])
    sigs  <- wilcox.test(
      dat[dat[[group_col]] == index[1], col],
      dat[dat[[group_col]] == index[2], col]
    )
    data.frame(
      W        = sigs$statistic,
      p        = sigs$p.value,
      mean_x   = mean(dat[dat[[group_col]] == index[1], col]),
      mean_y   = mean(dat[dat[[group_col]] == index[2], col]),
      median_x = median(dat[dat[[group_col]] == index[1], col]),
      median_y = median(dat[dat[[group_col]] == index[2], col])
    )
  }
  
  mat[[group_col]] <- as.factor(as.character(mat[[group_col]]))
  numeric_cols <- sapply(mat, is.numeric)
  numeric_cols[group_col] <- FALSE
  mat_num <- mat[, c(names(numeric_cols)[numeric_cols], group_col)]
  
  feat_cols <- colnames(mat_num)[colnames(mat_num) != group_col]
  tests     <- do.call(rbind, lapply(feat_cols, function(x) test.fun(mat_num, x)))
  rownames(tests) <- feat_cols
  
  test_sig          <- tests[tests$p < p_threshold, , drop = FALSE]
  test_sig$p.adjust <- p.adjust(test_sig$p, method = "bonferroni")
  test_sig          <- test_sig[order(test_sig$p), ]
  
  # ── 2. SD (from test mat) ─────────────────────────────────────────────────
  sd_file <- mat_num %>%
    group_by_at(group_col) %>%
    dplyr::summarise_all(sd) %>%
    t()
  colnames(sd_file)   <- sd_file[1, ]
  sd_file             <- as.data.frame(sd_file[-1, ])
  sd_file$id          <- rownames(sd_file)
  colnames(sd_file)[1:2] <- paste0("sd_", colnames(sd_file)[1:2])
  
  # ── 3. logFC ──────────────────────────────────────────────────────────────
  ref_mat <- if (!is.null(logfc_mat) && nrow(logfc_mat) > 0) logfc_mat else mat_num
  
  if (is.null(logfc_mat) && logfc_type == "log2ratio") {
    raw_check <- unlist(ref_mat[, colnames(ref_mat) != group_col])
    if (any(raw_check < 0, na.rm = TRUE)) {
      warning("logfc_mat not supplied and mat contains negative values. ",
              "log2ratio is not meaningful for negative data. ",
              "Falling back to 'diff' (mean difference).")
      logfc_type <- "diff"
    }
  }
  
  ref_mat[[group_col]] <- as.factor(as.character(ref_mat[[group_col]]))
  
  if (logfc_type == "log2ratio") {
    log2_means <- .compute_log2_means(ref_mat, group_col)
    shared <- intersect(rownames(test_sig), colnames(log2_means))
    logFC  <- log2_means[2, shared, drop = FALSE] - log2_means[1, shared, drop = FALSE]
  } else {
    raw_means <- ref_mat %>%
      group_by_at(group_col) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::select(-dplyr::all_of(group_col))
    shared <- intersect(rownames(test_sig), colnames(raw_means))
    logFC  <- raw_means[2, shared, drop = FALSE] - raw_means[1, shared, drop = FALSE]
  }
  
  logFC     <- as.data.frame(t(logFC))
  logFC$id  <- rownames(logFC)
  colnames(logFC)[1] <- "logFC"
  
  # ── 4. Merge & annotate ───────────────────────────────────────────────────
  test_sig$id   <- rownames(test_sig)
  last_test_sig <- merge(test_sig, sd_file,  by = "id")
  last_test_sig <- merge(last_test_sig, logFC, by = "id")
  last_test_sig <- last_test_sig[order(last_test_sig$p), ]
  
  last_test_sig$change <- as.factor(
    ifelse(last_test_sig$p < p_threshold,
           ifelse(last_test_sig$logFC >  0.5, "Up",
                  ifelse(last_test_sig$logFC < -0.5, "Down", "Stable")),
           "Stable")
  )
  
  if (save_data) {
    if (is.null(save_dir)) {
      stop("'save_dir' cannot be NULL when 'save_data' is TRUE. Please provide a valid directory path.")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    full_path <- file.path(save_dir, csv_filename)
    write.csv(last_test_sig, file = full_path, row.names = FALSE)
    cat("DEG results saved to:", full_path, "\n")
  }
  
  return(last_test_sig)
}


#' Internal function: Wilcoxon test for two specific groups
#'
#' Helper function that performs Wilcoxon test between two groups.
#' Used internally by multi-class DEG functions.
#'
#' @param mat Data frame with numeric features and group column.
#' @param group_col Name of the grouping column (must have exactly 2 levels).
#' @param logfc_mat Optional data frame for logFC calculation.
#' @param logfc_type One of "log2ratio" or "diff".
#'
#' @return Data frame with test statistics and logFC.
#' @keywords internal
.wilcoxon_two_groups <- function(mat, 
                                 group_col = "group",
                                 logfc_mat = NULL,
                                 logfc_type = "log2ratio") {

  if (!is.factor(mat[[group_col]])) {
    mat[[group_col]] <- factor(mat[[group_col]], levels = unique(mat[[group_col]]))
  }
  
  group_levels <- levels(mat[[group_col]])
  if (length(group_levels) != 2) {
    stop("Internal .wilcoxon_two_groups requires exactly 2 groups, found: ",
         paste(group_levels, collapse = ", "))
  }
  test.fun <- function(dat, col) {
    index <- levels(dat[[group_col]])
    sigs  <- wilcox.test(
      dat[dat[[group_col]] == index[1], col],
      dat[dat[[group_col]] == index[2], col]
    )
    data.frame(
      W        = sigs$statistic,
      p        = sigs$p.value,
      mean_x   = mean(dat[dat[[group_col]] == index[1], col]),
      mean_y   = mean(dat[dat[[group_col]] == index[2], col]),
      median_x = median(dat[dat[[group_col]] == index[1], col]),
      median_y = median(dat[dat[[group_col]] == index[2], col]),
      sd_x     = sd(dat[dat[[group_col]] == index[1], col]),
      sd_y     = sd(dat[dat[[group_col]] == index[2], col])
    )
  }

  numeric_cols <- sapply(mat, is.numeric)
  numeric_cols[group_col] <- FALSE
  mat_num <- mat[, c(names(numeric_cols)[numeric_cols], group_col)]
  
  feat_cols <- colnames(mat_num)[colnames(mat_num) != group_col]
  tests     <- do.call(rbind, lapply(feat_cols, function(x) test.fun(mat_num, x)))
  rownames(tests) <- feat_cols

  tests$p.adjust <- p.adjust(tests$p, method = "bonferroni")
 
  ref_mat <- if (!is.null(logfc_mat) && nrow(logfc_mat) > 0) logfc_mat else mat_num

  if (!is.factor(ref_mat[[group_col]])) {
    ref_mat[[group_col]] <- factor(ref_mat[[group_col]], levels = unique(ref_mat[[group_col]]))
  } else {
    if (!identical(levels(ref_mat[[group_col]]), group_levels)) {
      ref_mat[[group_col]] <- factor(ref_mat[[group_col]], levels = group_levels)
    }
  }

  if (logfc_type == "log2ratio") {
    raw_check <- unlist(ref_mat[, colnames(ref_mat) != group_col])
    if (any(raw_check < 0, na.rm = TRUE)) {
      logfc_type <- "diff"
    }
  }
  
  if (logfc_type == "log2ratio") {
    log2_means <- .compute_log2_means(ref_mat, group_col)
    shared <- intersect(rownames(tests), colnames(log2_means))
    logFC  <- log2_means[2, shared, drop = FALSE] - log2_means[1, shared, drop = FALSE]
  } else {
    raw_means <- ref_mat %>%
      group_by_at(group_col) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::select(-dplyr::all_of(group_col))
    shared <- intersect(rownames(tests), colnames(raw_means))
    logFC  <- raw_means[2, shared, drop = FALSE] - raw_means[1, shared, drop = FALSE]
  }
  
  logFC     <- as.data.frame(t(logFC))
  logFC$id  <- rownames(logFC)
  colnames(logFC)[1] <- "logFC"
  
  tests$id <- rownames(tests)
  result <- merge(tests, logFC, by = "id")
  result <- result[order(result$p), ]
  
  return(result)
}

#' Find All Markers for Multi-class Classification
#'
#' Performs DEG analysis for multi-class data by comparing each class against all others.
#' This function internally converts multi-class problem to binary comparisons (class_i vs. others).
#'
#' @param mat Data frame with numeric features and a group column.
#' @param group_col Name of the grouping column (can have 2+ groups).
#' @param logfc_mat Optional data frame for logFC calculation.
#' @param logfc_type One of "log2ratio" (default) or "diff".
#' @param p_threshold P-value threshold for feature significance. Default 0.05.
#' @param only.pos Logical; if TRUE, only keep features with logFC > 0. Default FALSE.
#' @param save_dir Directory to save results.
#' @param save_data Logical; whether to write CSV output.
#' @param csv_filename Output filename.
#'
#' @return Data frame combining results for all groups, with added "target_group" column.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- batch_Wilcoxon_MultiClass(mat = iris,group_col = "Species",only.pos = FALSE)
#' }
batch_Wilcoxon_MultiClass <- function(mat,
                                      group_col = "group",
                                      logfc_mat = NULL,
                                      logfc_type = "log2ratio",
                                      p_threshold = 0.05,
                                      only.pos = FALSE,
                                      save_dir = NULL,
                                      save_data = FALSE,
                                      csv_filename = "deg_multiclass_all.csv") {
  
  # Get unique groups
  group_levels <- unique(na.omit(as.character(mat[[group_col]])))
  if (length(group_levels) < 2) {
    stop("batch_Wilcoxon_MultiClass requires at least 2 groups. Found: ",
         paste(group_levels, collapse = ", "))
  }
  
  message("Multi-class DEG analysis: ", length(group_levels), " groups detected")
  
  results_list <- list()
  
  # For each group, compare it against all others
  for (i in seq_along(group_levels)) {
    target_group <- group_levels[i]
    message("Processing group: ", target_group, " (", i, "/", length(group_levels), ")")
    
    # Create binary comparison: target_group vs. others
    temp_mat <- mat
    temp_mat[[group_col]] <- as.character(temp_mat[[group_col]])
    temp_mat[[group_col]] <- ifelse(
      temp_mat[[group_col]] == target_group,
      target_group,
      "others"
    )
    temp_mat[[group_col]] <- factor(temp_mat[[group_col]], levels = c("others", target_group))
    
    # Prepare logfc_mat if provided
    temp_logfc_mat <- NULL
    if (!is.null(logfc_mat)) {
      temp_logfc_mat <- logfc_mat
      temp_logfc_mat[[group_col]] <- as.character(temp_logfc_mat[[group_col]])
      temp_logfc_mat[[group_col]] <- ifelse(
        temp_logfc_mat[[group_col]] == target_group,
        target_group,
        "others"
      )
      temp_logfc_mat[[group_col]] <- factor(temp_logfc_mat[[group_col]], 
                                            levels = c("others", target_group))
    }
    
    # Run Wilcoxon test
    deg_result <- .wilcoxon_two_groups(
      temp_mat,
      group_col = group_col,
      logfc_mat = temp_logfc_mat,
      logfc_type = logfc_type
    )
    
    # Filter by p-threshold
    deg_result <- deg_result[deg_result$p < p_threshold, , drop = FALSE]
    
    if (only.pos) {
      deg_result <- deg_result[deg_result$logFC > 0, , drop = FALSE]
    }
    
    # Add annotation
    deg_result$change <- as.factor(
      ifelse(deg_result$logFC > 0.5, "Up",
             ifelse(deg_result$logFC < -0.5, "Down", "Stable"))
    )
    
    deg_result$target_group <- target_group
    
    results_list[[i]] <- deg_result
  }
  
  # Combine all results
  all_results <- do.call(rbind, results_list)
  rownames(all_results) <- NULL
  
  # Sort by group and p-value
  all_results <- all_results[order(all_results$target_group, all_results$p), ]
  
  if (save_data) {
    if (is.null(save_dir)) {
      stop("'save_dir' cannot be NULL when 'save_data' is TRUE. Please provide a valid directory path.")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    full_path <- file.path(save_dir, csv_filename)
    write.csv(all_results, file = full_path, row.names = FALSE)
    cat("Multi-class DEG results saved to:", full_path, "\n")
  }
  
  return(all_results)
}


#' Find Sequential Markers for Ordered Multi-class Classification
#'
#' Performs pairwise DEG analysis between consecutive classes (ordered comparisons).
#' For k classes, performs k-1 comparisons: class_1 vs. class_2, class_2 vs. class_3, etc.
#'
#' @param mat Data frame with numeric features and a group column.
#' @param group_col Name of the grouping column (must have ordered factor levels).
#' @param logfc_mat Optional data frame for logFC calculation.
#' @param logfc_type One of "log2ratio" (default) or "diff".
#' @param p_threshold P-value threshold for feature significance. Default 0.05.
#' @param only.pos Logical; if TRUE, only keep features with logFC > 0. Default TRUE.
#' @param save_dir Directory to save results.
#' @param save_data Logical; whether to write CSV output.
#' @param csv_filename Output filename.
#'
#' @return Data frame combining results for all consecutive pairs.
#' @export
#'
#' @examples
#' \dontrun{
#' iris$Species=factor(iris$Species,levels = unique(iris$Species))
#' result <- batch_Wilcoxon_OrderedMultiClass(mat = iris,group_col = "Species",only.pos = FALSE)
#' }
batch_Wilcoxon_OrderedMultiClass <- function(mat,
                                             group_col = "group",
                                             logfc_mat = NULL,
                                             logfc_type = "log2ratio",
                                             p_threshold = 0.05,
                                             only.pos = TRUE,
                                             save_dir = NULL,
                                             save_data = FALSE,
                                             csv_filename = "deg_ordered_sequential.csv") {
  
  # Get unique groups in order
  group_levels <- levels(as.factor(as.character(mat[[group_col]])))
  if (length(group_levels) < 2) {
    stop("batch_Wilcoxon_OrderedMultiClass requires at least 2 groups. Found: ",
         paste(group_levels, collapse = ", "))
  }
  
  message("Ordered multi-class DEG analysis: ", length(group_levels), " groups detected")
  
  results_list <- list()
  
  # For consecutive pairs
  for (i in 1:(length(group_levels) - 1)) {
    group1 <- group_levels[i]
    group2 <- group_levels[i + 1]
    message("Comparing: ", group1, " vs. ", group2, " (", i, "/", length(group_levels) - 1, ")")
    
    # Keep only these two groups
    temp_mat <- mat[as.character(mat[[group_col]]) %in% c(group1, group2), ]
    temp_mat[[group_col]] <- factor(temp_mat[[group_col]], levels = c(group1, group2))
    
    # Prepare logfc_mat if provided
    temp_logfc_mat <- NULL
    if (!is.null(logfc_mat)) {
      temp_logfc_mat <- logfc_mat[as.character(logfc_mat[[group_col]]) %in% c(group1, group2), ]
      temp_logfc_mat[[group_col]] <- factor(temp_logfc_mat[[group_col]], levels = c(group1, group2))
    }
    
    # Run Wilcoxon test
    deg_result <- .wilcoxon_two_groups(
      temp_mat,
      group_col = group_col,
      logfc_mat = temp_logfc_mat,
      logfc_type = logfc_type
    )
    
    # Filter by p-threshold
    deg_result <- deg_result[deg_result$p < p_threshold, , drop = FALSE]
    
    if (only.pos) {
      deg_result <- deg_result[deg_result$logFC > 0, , drop = FALSE]
    }
    
    # Add annotation
    deg_result$change <- as.factor(
      ifelse(deg_result$logFC > 0.5, "Up",
             ifelse(deg_result$logFC < -0.5, "Down", "Stable"))
    )
    
    deg_result$comparison <- paste0(group1, " vs. ", group2)
    deg_result$target_group <- group2
    
    results_list[[i]] <- deg_result
  }
  
  # Combine all results
  all_results <- do.call(rbind, results_list)
  rownames(all_results) <- NULL
  
  # Sort by comparison order and p-value
  all_results <- all_results[order(all_results$comparison, all_results$p), ]
  
  if (save_data) {
    if (is.null(save_dir)) {
      stop("'save_dir' cannot be NULL when 'save_data' is TRUE. Please provide a valid directory path.")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    full_path <- file.path(save_dir, csv_filename)
    write.csv(all_results, file = full_path, row.names = FALSE)
    cat("Ordered multi-class DEG results saved to:", full_path, "\n")
  }
  
  return(all_results)
}

#' Unified Interface for DEG Analysis (Binary and Multi-class)
#'
#' Automatically detects number of groups and routes to appropriate function.
#' Provides a single entry point for both binary and multi-class DEG analysis.
#'
#' @param mat Data frame with numeric features and a group column.
#' @param group_col Name of the grouping column.
#' @param analysis_type One of:
#'   - "auto" (default): automatically detect based on number of groups
#'   - "binary": force binary classification (2 groups)
#'   - "multiclass": one-vs-rest comparison (2+ groups)
#'   - "ordered": sequential pairwise comparison (2+ groups with order)
#' @param logfc_mat Optional data frame for logFC calculation.
#' @param logfc_type One of "log2ratio" (default) or "diff".
#' @param p_threshold P-value threshold. Default 0.05.
#' @param only.pos Logical; filter to positive logFC. Default depends on analysis_type.
#' @param save_dir Directory to save results.
#' @param save_data Logical; whether to write CSV. Default TRUE.
#' @param csv_filename Output filename. If NULL, auto-generated based on analysis_type.
#'
#' @return Data frame with DEG results. For multi-class: includes "target_group" column.
#'         For ordered: includes "comparison" column.
#' @export
#'
#' @examples
#' \dontrun{
#' # Binary analysis (2 groups) - automatic
#' result_binary <- batch_Wilcoxon_Unified(mtcars,group_col = 'vs')
#' # Multi-class analysis (3+ groups) - automatic
#' result_multi <- batch_Wilcoxon_Unified(iris,group_col ="Species")
#' # Ordered analysis (stages)
#' iris$Species=factor(iris$Species,levels = unique(iris$Species))
#' result_ordered <- batch_Wilcoxon_Unified(
#' iris,group_col = "Species",analysis_type = "ordered")
#' }
batch_Wilcoxon_Unified <- function(mat,
                                   group_col = "group",
                                   analysis_type = "auto",
                                   logfc_mat = NULL,
                                   logfc_type = "log2ratio",
                                   p_threshold = 0.05,
                                   only.pos = NULL,
                                   save_dir = NULL,
                                   save_data = FALSE,
                                   csv_filename = NULL) {
  
  # Detect number of groups
  group_levels <- unique(na.omit(as.character(mat[[group_col]])))
  n_groups <- length(group_levels)
  
  message("Detected ", n_groups, " group(s): ", paste(group_levels, collapse = ", "))
  
  # Auto-detect analysis type if not specified
  if (analysis_type == "auto") {
    if (n_groups == 2) {
      analysis_type <- "binary"
    } else if (n_groups > 2) {
      analysis_type <- "multiclass"
    } else {
      stop("Cannot perform DEG analysis with fewer than 2 groups.")
    }
    message("Auto-detected analysis type: ", analysis_type)
  }
  
  # Validate analysis type
  if (!analysis_type %in% c("binary", "multiclass", "ordered")) {
    stop("analysis_type must be one of: 'auto', 'binary', 'multiclass', 'ordered'")
  }
  
  # Check group count compatibility
  if (analysis_type == "binary" && n_groups != 2) {
    stop("Binary analysis requires exactly 2 groups, but found ", n_groups)
  }
  if (analysis_type %in% c("multiclass", "ordered") && n_groups < 2) {
    stop(analysis_type, " analysis requires at least 2 groups, but found ", n_groups)
  }
  
  # Set only.pos default if not provided
  if (is.null(only.pos)) {
    only.pos <- if (analysis_type == "ordered") TRUE else FALSE
  }
  
  # Auto-generate csv_filename if not provided
  if (is.null(csv_filename)) {
    csv_filename <- switch(analysis_type,
                           binary = "deg_binary_result.csv",
                           multiclass = "deg_multiclass_allvrest.csv",
                           ordered = "deg_multiclass_ordered.csv")
  }
  
  # Route to appropriate function
  result <- switch(analysis_type,
                   binary = batch_Wilcoxon(
                     mat = mat,
                     group_col = group_col,
                     logfc_mat = logfc_mat,
                     logfc_type = logfc_type,
                     p_threshold = p_threshold,
                     save_dir = save_dir,
                     save_data = save_data,
                     csv_filename = csv_filename
                   ),
                   multiclass = batch_Wilcoxon_MultiClass(
                     mat = mat,
                     group_col = group_col,
                     logfc_mat = logfc_mat,
                     logfc_type = logfc_type,
                     p_threshold = p_threshold,
                     only.pos = only.pos,
                     save_dir = save_dir,
                     save_data = save_data,
                     csv_filename = csv_filename
                   ),
                   ordered = batch_Wilcoxon_OrderedMultiClass(
                     mat = mat,
                     group_col = group_col,
                     logfc_mat = logfc_mat,
                     logfc_type = logfc_type,
                     p_threshold = p_threshold,
                     only.pos = only.pos,
                     save_dir = save_dir,
                     save_data = save_data,
                     csv_filename = csv_filename
                   )
  )
  
  message("DEG analysis complete. Analysis type: ", analysis_type)
  
  return(result)
}

# ══════════════════════════════════════════════════════════════════════════════
# stat_var_feature 
# ══════════════════════════════════════════════════════════════════════════════

#' Perform Feature Selection Using Batch Wilcoxon Test
#'
#' Entry point for DEG analysis. Applies the following priority logic:
#'
#' **Wilcoxon test data** (for ranking/p-value):
#'   - Always uses `clean.data` if available (original untransformed values
#'     give the most interpretable Wilcoxon statistics).
#'   - Falls back to `scale.data` if `clean.data` is empty or user sets
#'     `data_type = "scale"` explicitly (e.g. when the input data itself is
#'     already log-transformed or otherwise pre-processed before loading).
#'
#' **logFC calculation data**:
#'   - `clean.data` (raw positive values)  → `logfc_type = "log2ratio"` (default)
#'   - `scale.data` with log-space method  → `logfc_type = "diff"` (mean difference = log ratio)
#'   - `scale.data` with ratio-intact method (`scale`, `min_max`, `max_abs`) → `logfc_type = "log2ratio"`
#'   - `scale.data` with centered method   → logFC still from `clean.data` if available
#'
#' @param object    A `Stat` object or a plain data frame.
#' @param group_col Group column name (auto-read from Stat slot when object is Stat).
#' @param data_type `"auto"` (default): prefer clean, fall back to scale.
#'                  `"clean"`: force use of clean.data for both test and logFC.
#'                  `"scale"`: force use of scale.data for the test; logFC source
#'                  is determined by the scale method.
#' @param p_threshold Passed through to `batch_Wilcoxon`; features with p >= p_threshold
#'                    are excluded. Defaults to 0.05.
#' @param save_dir  Output directory.
#' @param save_data Logical; save CSV output.
#' @param csv_filename  Output filename.
#'
#' @returns Updated `Stat` object (if input is Stat) or a data frame of DEG results.
#' @export
#'
#' @examples
#' \dontrun{
#' # Recommended: let the function choose automatically
#' stat_obj <- stat_var_feature(stat_obj_test, p_threshold = 0.05)
#' }
stat_var_feature <- function(object,
                             group_col    = "group",
                             data_type    = "auto",
                             p_threshold  = 0.05,
                             save_dir     = NULL,
                             save_data    = FALSE,
                             csv_filename = "last_test_sig.csv") {
  # ---- Auto‑generate default save directory ----
  if (save_data && is.null(save_dir)) {
    if (exists("get_output_dir")) {
      save_dir <- get_output_dir("m1", "deg_results")
    } else {
      save_dir <- file.path(".", "results", "deg")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    cat("DEG results will be saved to:", save_dir, "\n")
  }
  # ── 1. Extract data from Stat object or data frame ───────────────────────
  if (inherits(object, "Stat")) {
    
    group_col <- slot(object, "group_col")
    if (length(group_col) == 0 || is.null(group_col)) {
      stop("group_col slot is empty in the Stat object. Please set it when creating the object.")
    }
    
    clean_dat <- ExtractCleanData(object)
    scale_dat <- tryCatch(ExtractScaleData(object), error = function(e) NULL)
    
    has_clean <- !is.null(clean_dat) && nrow(clean_dat) > 0
    has_scale <- !is.null(scale_dat) && nrow(scale_dat) > 0
    
    # Determine scale method recorded in process.info (may be NULL)
    norm_info <- tryCatch(slot(object, "process.info")[["normalization"]], error = function(e) NULL)
    scale_method <- if (!is.null(norm_info) && length(norm_info) > 0) {
      methods_used <- unique(unlist(norm_info))
      if (length(methods_used) == 1) {
        methods_used
      } else {
        most_common <- names(sort(table(unlist(norm_info)), decreasing = TRUE))[1]
        warning("Mixed normalization methods across columns: ",
                paste(methods_used, collapse = ", "),
                ". Using most common ('", most_common, "') to determine logFC strategy.")
        most_common
      }
    } else {
      "none"
    }
    
    # ── Decision table ───────────────────────────────────────────────────
    #
    #  data_type = "auto"  → prefer clean for test; logFC always from clean
    #  data_type = "clean" → same as auto but explicit
    #  data_type = "scale" → user knows their input is already processed;
    #                         test on scale data, logFC strategy from scale_method
    #
    log_space_methods    <- c("log", "box_cox", "yeo_johnson")
    ratio_intact_methods <- c("none", "min_max", "max_abs", "scale")
    centered_methods     <- c("z_score", "center")
    
    if (data_type %in% c("auto", "clean")) {
      
      if (!has_clean) {
        if (!has_scale) stop("No valid data found in clean.data or scale.data.")
        message("clean.data is empty — falling back to scale.data for test. ",
                "If your raw data is already normalised, set data_type = 'scale'.")
        test_dat   <- scale_dat
        logfc_dat  <- scale_dat
        lfc_type   <- if (scale_method %in% log_space_methods) "diff" else "log2ratio"
      } else {
        # Ideal path: clean data for test, clean data for logFC
        test_dat   <- clean_dat
        logfc_dat  <- clean_dat
        lfc_type   <- "log2ratio"
        if (has_scale) {
          message("Using clean.data for Wilcoxon test and logFC calculation (recommended).")
        }
      }
      
    } else if (data_type == "scale") {
      
      if (!has_scale) stop("data_type = 'scale' requested but scale.data is empty.")
      test_dat <- scale_dat
      
      if (scale_method %in% log_space_methods) {
        # scale data is in log space: mean diff = log ratio
        logfc_dat <- scale_dat
        lfc_type  <- "diff"
        message("scale_method='", scale_method, "': using scale.data for logFC (mean diff = log ratio).")
        
      } else if (scale_method %in% centered_methods) {
        # centering broke ratio; recover logFC from clean data if available
        if (has_clean) {
          logfc_dat <- clean_dat
          lfc_type  <- "log2ratio"
          message("scale_method='", scale_method, "': Wilcoxon on scale.data, logFC from clean.data.")
        } else {
          logfc_dat <- scale_dat
          lfc_type  <- "log2ratio"   # will auto-warn inside batch_Wilcoxon if negatives detected
          warning("scale_method='", scale_method, "' and clean.data unavailable. ",
                  "logFC computed from scale.data — values may not be true log2 fold changes.")
        }
        
      } else {
        # ratio-intact (scale x/sd, min_max, max_abs, none): log2ratio directly
        logfc_dat <- scale_dat
        lfc_type  <- "log2ratio"
        message("scale_method='", scale_method, "': ratio relationship intact, using log2ratio for logFC.")
      }
      
    } else {
      stop("data_type must be one of: 'auto', 'clean', 'scale'.")
    }
    
  } else if (is.data.frame(object)) {
    # Plain data frame: take as-is, user is responsible for the data state
    test_dat  <- object
    logfc_dat <- object
    lfc_type  <- "log2ratio"
    # Auto-detect: if negatives present, switch to diff and warn
    num_vals <- unlist(object[, sapply(object, is.numeric)])
    if (any(num_vals < 0, na.rm = TRUE)) {
      warning("Input data frame contains negative values. ",
              "Switching logfc_type to 'diff' (mean difference). ",
              "If data is in log space this is correct; if it is z-scored, ",
              "consider passing a Stat object with clean.data for proper log2FC.")
      lfc_type <- "diff"
    }
  } else {
    stop("Input must be a 'Stat' object or a data frame.")
  }
  
  # ── 2. Validate group column ──────────────────────────────────────────────
  if (!group_col %in% colnames(test_dat)) {
    stop("Group column '", group_col, "' not found in the test data.")
  }
  group_values <- test_dat[[group_col]]
  if (any(is.na(group_values)) || length(unique(na.omit(group_values))) < 2) {
    stop("Group column '", group_col, "' must have at least two distinct non-missing values.")
  }
  
  # ── 3. Run batch Wilcoxon ─────────────────────────────────────────────────
  test_dat  <- data.frame(lapply(test_dat,  function(x) if (is.numeric(x)) as.numeric(x) else x))
  logfc_dat <- data.frame(lapply(logfc_dat, function(x) if (is.numeric(x)) as.numeric(x) else x))
  
  cat("Starting batch Wilcoxon test...\n")
  cat("  Test data  : ", nrow(test_dat), "samples x",
      sum(sapply(test_dat, is.numeric)), "numeric features\n")
  cat("  logFC data : ", if (identical(test_dat, logfc_dat)) "same as test data"
      else paste(nrow(logfc_dat), "samples"), "\n")
  cat("  logFC type :", lfc_type, "\n")
  
  last_test_sig <- batch_Wilcoxon(
    mat          = test_dat,
    group_col    = group_col,
    logfc_mat    = logfc_dat,
    logfc_type   = lfc_type,
    p_threshold  = p_threshold,
    save_dir     = save_dir,
    save_data    = save_data,
    csv_filename = csv_filename
  )
  
  # ── 4. Store result ───────────────────────────────────────────────────────
  if (inherits(object, "Stat")) {
    object@var.result <- list(last_test_sig = last_test_sig)
    cat("Feature selection completed. Significant features:", nrow(last_test_sig), "\n")
    cat("'var.result' slot updated.\n")
    return(object)
  } else {
    cat("Significant features:", nrow(last_test_sig), "\n")
    return(last_test_sig)
  }
}


#' Extract Last Significant Test Results
#'
#' This function extracts the last significant test results stored in the `var.result` slot of a `Stat` object. If
#' the `last_test_sig` result is not available, it returns `NULL`.
#'
#' @param object An object of class `Stat`. The function attempts to extract the `last_test_sig` from the `var.result` slot.
#'
#' @returns Returns the last significant test results stored in the `last_test_sig` slot of the `var.result` list.
#' If the `last_test_sig` is not found, it returns `NULL`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' last_sig <- ExtractLastTestSig(stat_obj)
#' print(last_sig)
#' }
ExtractLastTestSig <- function(object) {
  last_test_sig <- tryCatch(object@var.result[["last_test_sig"]],
                            error = function(e) NULL)
  return(last_test_sig)
}

#' Plot DEG Radar Chart
#'
#' This function generates a radar chart visualizing the differential expression (DEG) results from a data frame.
#' It uses a set of custom parameters to customize the appearance, such as colors, size, and grouping. The chart can be
#' saved to a specified directory as a PDF file.
#'
#'
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr pivot_longer
#' @import wesanderson
#' @import here
#' @param df A data frame containing the DEG results with columns for ID, means, log-fold change, p-values,
#' adjusted p-values, and other necessary data.
#' @param palette_name A string indicating the name of the color palette to use for the chart .
#' @param x_col The column name to use as the x-axis (default: "id").
#' @param y_cols A vector of column names to be plotted on the y-axis (default: c("mean_x", "mean_y")).
#' @param size_col The column name for the size of the points (default: "logFC").
#' @param color_col The column name for the color of the points, typically for p-values (default: "logp").
#' @param fill_col The column name for the fill color, which defines the grouping (default: "change").
#' @param p_adjust_col The column name for the adjusted p-values (default: "p.adjust").
#' @param plot_width The width of the plot (default: 5).
#' @param plot_height The height of the plot (default: 5).
#' @param save_dir The directory where the radar chart image will be saved (default: "here('StatObject', 'deg_info')").
#' @param base_size The base font size for the plot (default: 14).
#' @param title The title of the radar chart (default: "Radar Chart Title").
#'
#' @returns The ggplot object representing the radar chart.
#' @export
#'
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' last_sig <- ExtractLastTestSig(stat_obj)
#' print(last_sig)
#' plot_deg_radarchart(last_sig, palette_name = "Zissou1", title = "Differential Expression Radar")
#' }
plot_deg_radarchart <- function(df,
                                palette_name = "Zissou1",
                                x_col = "id",
                                y_cols = c("mean_x", "mean_y"),
                                size_col = "logFC",
                                color_col = "logp",
                                fill_col = "change",
                                p_adjust_col = "p.adjust",
                                plot_width = 5,
                                plot_height = 5,
                                save_dir = NULL,
                                base_size = 14,
                                title = "Radar Chart Title") {
  
  colors <- wes_palette(n = 3, name = palette_name, type = "continuous")
  
  colors <- as.list(colors)
  
  
  primary_color <- colors[[1]]
  secondary_color <- colors[[2]]
  
  
  df <- df %>%
    mutate(
      logp = -log10(get(p_adjust_col)),
      group = factor(get(fill_col), levels = c("Down", "Stable"))
    )
  
  df_long <- df %>%
    pivot_longer(cols = all_of(y_cols), names_to = "var", values_to = "value")
  
  plot <- ggplot(df_long, aes(x = get(x_col), y = value, fill = group)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.7) +
    geom_point(aes(size = abs(get(size_col)), color = get(color_col)), position = position_dodge(width = 0.9)) +
    coord_polar(start = 0) +
    scale_fill_manual(values = c("Down" = primary_color, "Stable" = secondary_color)) +
    scale_color_gradient(low = primary_color, high = secondary_color) +
    scale_size_continuous(range = c(2, 10)) +
    ggprism::theme_prism(base_size = base_size) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      legend.position = "right",
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
    ) +
    labs(fill = "Change", color = "-log10(p.adjust)", size = "LogFC", title = title) +
    geom_hline(yintercept = 0, color = "gray60", linetype = "dashed")
  
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(save_dir, "radar_chart.pdf"), plot = plot, width = plot_width, height = plot_height,
           device = "pdf")
    cat("Radar chart saved to: ", file.path(save_dir, "radar_chart.pdf"))
  }
  
  return(plot)
}

#' Generate Radar Chart for Variable Features (with auto‑save option)
#'
#' @param object Stat object or data frame.
#' @param group_col Group column name. If NULL, auto‑detect from Stat object.
#' @param palette_name Colour palette.
#' @param plot_width,plot_height Plot dimensions (inches).
#' @param save_dir Output directory. If NULL and save_plots=TRUE, will use
#'   "./figures/deg_info" or get_output_dir() if available.
#' @param save_plots Logical. Save the plot? Default TRUE.
#' @param base_size Base font size.
#' @param title Plot title.
#'
#' @return Updated Stat object or ggplot.
#' @export
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' VarFeature_radarchart(stat_obj, palette_name = "Zissou1", title = "Differential Expression Radar")
#' }
VarFeature_radarchart <- function(object,
                                  group_col = NULL,
                                  palette_name = "Zissou1",
                                  plot_width = 5,
                                  plot_height = 5,
                                  save_dir = NULL,
                                  save_plots = TRUE,
                                  base_size = 14,
                                  title = "Variable Feature Radar Chart") {
  cat("Generating variable feature radar chart...\n")
  
  # ---- Handle default save_dir ----
  if (is.null(save_dir) && save_plots) {
    if (exists("get_output_dir")) {
      save_dir <- get_output_dir("Figures", "deg_info")
    } else {
      save_dir <- file.path(".", "figures", "deg_info")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  } else if (!save_plots) {
    save_dir <- NULL   # force no saving
  }
  
  # ---- Extract data and group_col ----
  if (inherits(object, "Stat")) {
    mat <- ExtractLastTestSig(object)
    if (is.null(group_col)) {
      group_col <- slot(object, "group_col")
    }
    if (length(group_col) == 0) group_col <- NULL
  } else if (is.data.frame(object)) {
    mat <- object
    if (is.null(group_col)) group_col <- "group"
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }
  
  if (is.null(mat) || nrow(mat) == 0) {
    stop("No valid data found in the input.")
  }
  
  # ---- Generate radar chart ----
  plot <- plot_deg_radarchart(
    mat,
    palette_name = palette_name,
    plot_width = plot_width,
    plot_height = plot_height,
    save_dir = save_dir,        # now correctly NULL if save_plots = FALSE
    base_size = base_size,
    title = title
  )
  
  print(plot)
  
  if (inherits(object, "Stat")) {
    object@var.result[["VarFeaturePlot"]] <- plot
    return(object)
  }
  return(plot)
}


#' Generate Volcano Plot for Differential Expression
#'
#' This function generates a volcano plot based on the log-fold change and adjusted p-values of differential expression results.
#' The plot highlights significant and non-significant genes using different colors and includes labels for the significant genes.
#' The volcano plot is saved as an image in the specified directory.
#' @importFrom dplyr mutate %>%
#' @param last_test_sig A data frame containing the results of differential expression analysis. It should have columns
#'                      for log-fold change and adjusted p-values.
#' @param logFC_col The name of the column representing log-fold change values (default: "logFC").
#' @param p_adjust_col The name of the column representing adjusted p-values (default: "p.adjust").
#' @param title The title of the volcano plot (default: "Volcano Plot").
#' @param palette_name A string indicating the color palette to use for the plot.
#' @param save_dir The directory to save the volcano plot image (default: "here('StatObject', 'deg_info')").
#' @param plot_width The width of the plot in inches (default: 5).
#' @param plot_height The height of the plot in inches (default: 5).
#' @param base_size The base font size for the plot (default: 14).
#'
#' @returns A `ggplot` object representing the volcano plot.
#' @export
#'
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' last_sig <- ExtractLastTestSig(stat_obj)
#' print(last_sig)
#' plot_deg_volcano(last_test_sig = last_sig, logFC_col = "logFC", p_adjust_col = "p.adjust")
#' plot_deg_volcano(last_test_sig = last_sig, title = "Custom Volcano Plot", palette_name = "Zissou1")
#' }
plot_deg_volcano <- function(last_test_sig,
                             logFC_col = "logFC",
                             p_adjust_col = "p.adjust",
                             title = "Volcano Plot",
                             palette_name = "Zissou1",
                             save_dir = NULL,
                             plot_width = 5,
                             plot_height = 5,
                             base_size = 14) {
  
  if (is.null(last_test_sig) || nrow(last_test_sig) == 0) {
    stop("No valid data found in last_test_sig.")
  }
  
  last_test_sig <- last_test_sig %>%
    mutate(log_p = -log10(get(p_adjust_col)))
  
  
  colors <- wes_palette(n = 3, name = palette_name, type = "continuous")
  
  colors <- as.list(colors)
  significant_color <- colors[[1]]
  not_significant_color <- colors[[2]]
  
  volcano_plot <- ggplot(last_test_sig, aes_string(x = logFC_col, y = "log_p")) +
    geom_point(aes(color = ifelse(get(p_adjust_col) < 0.05, "Significant", "Not Significant")),
               alpha = 0.6, size = 3) +
    scale_color_manual(values = c("Significant" = significant_color, "Not Significant" = not_significant_color)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") +
    labs(title = title, x = "Log Fold Change", y = "-Log10 Adjusted P-value") +
    ggprism::theme_prism(base_size = base_size) +
    theme(
      legend.title = element_blank()
    ) +
    geom_text(aes(label = id), vjust = -1, hjust = 0.5, size = 3, check_overlap = TRUE, color = "black")
  
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("Created directory: ", save_dir, "\n")
    }
    ggsave(filename = file.path(save_dir, "volcano_plot.pdf"), plot = volcano_plot,
           width = plot_width, height = plot_height,
           device = "pdf")
    cat("Volcano plot saved to: ", file.path(save_dir, "volcano_plot.pdf"), "\n")
  }
  
  cat("Volcano plot saved to: ", file.path(save_dir, "volcano_plot.pdf"), "\n")
  
  return(volcano_plot)
}


#' Generate Volcano Plot for Variable Features
#'
#' This function generates a volcano plot based on the log-fold change and adjusted p-values of variable features from
#' a 'Stat' object or a data frame. It then saves the volcano plot as an image in the specified directory.
#'
#' @param object An object of class 'Stat' or a data frame containing the results of differential expression analysis.
#'               It should include columns for log-fold change and adjusted p-values.
#' @param logFC_col The name of the column representing log-fold change values (default: "logFC").
#' @param p_adjust_col The name of the column representing adjusted p-values (default: "p.adjust").
#' @param title The title of the volcano plot (default: "Volcano Plot").
#' @param palette_name A string indicating the color palette to use for the plot.
#' @param plot_width The width of the plot in inches (default: 5).
#' @param plot_height The height of the plot in inches (default: 5).
#' @param save_dir The directory to save the volcano plot image (default: "here('StatObject', 'deg_info')").
#' @param base_size The base font size for the plot (default: 14).
#'
#' @returns A `ggplot` object representing the volcano plot. If the input is a 'Stat' object, it also updates the
#'          'Stat' object with the volcano plot.
#' @export
#'
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' VarFeature_volcano(object = stat_obj, logFC_col = "logFC", p_adjust_col = "p.adjust")
#' VarFeature_volcano(object = stat_obj, title = "Custom Volcano Plot", palette_name = "Zissou1")
#' }
VarFeature_volcano <- function(object,
                               logFC_col = "logFC",
                               p_adjust_col = "p.adjust",
                               title = "Volcano Plot",
                               palette_name = "Zissou1",
                               plot_width = 5,
                               plot_height = 5,
                               save_dir = NULL,
                               base_size = 14) {
  
  if (inherits(object, "Stat")) {
    last_test_sig <- ExtractLastTestSig(object)
  } else if (is.data.frame(object)) {
    last_test_sig <- object
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }
  
  if (is.null(last_test_sig) || nrow(last_test_sig) == 0) {
    stop("No valid data found in last_test_sig.")
  }
  
  volcano_plot <- plot_deg_volcano(last_test_sig,
                                   logFC_col = logFC_col,
                                   p_adjust_col = p_adjust_col,
                                   title = title,
                                   palette_name = palette_name,
                                   plot_width = plot_width,
                                   plot_height = plot_height,
                                   base_size = base_size)
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(save_dir, "volcano_plot.pdf"), plot = volcano_plot,
           width = plot_width, height = plot_height, device = "pdf")
  }
  
  if (inherits(object, "Stat")) {
    object@var.result[["VolcanoPlot"]] <- volcano_plot
    return(object)
  }
  
  return(volcano_plot)
}


#' Generate box Plot for Differentially Expressed Genes
#'
#' This function generates a violin plot for the top differentially expressed genes based on their log-fold change
#' and statistical significance. It also calculates the Wilcoxon test p-values for group comparisons and annotates
#' the plot with significance stars. The plot can be saved to a specified directory.
#' @importFrom dplyr %>% group_by arrange desc top_n mutate case_when
#' @importFrom dplyr summarise left_join
#' @importFrom tidyr pivot_longer
#'
#' @import stats
#' @param last_test_sig A data frame containing the results of differential expression analysis with columns
#'                      such as 'change' (e.g., 'Stable', 'Upregulated', 'Downregulated') and 'logFC' (log-fold change).
#' @param data A data frame containing the expression data with rows as samples and columns as features.
#' @param control The label of the control group (default: "health").
#' @param case The label of the case group (default: "cancer").
#' @param top_n The number of top significant features to display (default: 5).
#' @param palette_name The name of the color palette for the plot.
#' @param name_identity The identity for the type of analysis, default is "deg".
#' @param save_plots Logical value to indicate whether to save the plot (default: TRUE).
#' @param save_dir The directory to save the plot (default: "here('StatObject', 'deg_info')").
#' @param plot_width The width of the plot in inches (default: 5).
#' @param plot_height The height of the plot in inches (default: 5).
#' @param base_size The base font size for the plot (default: 14).
#' @param title The title of the plot (default: "Violin Plot").
#' @param group_col The column name used to group the data (default: 'group').
#' @param pseudocount A small positive value added before log10 transform to avoid
#'                    log(0) / NaN. Default: 1 (i.e. log10(value + 1)).
#'
#' @returns A ggplot object representing the violin plot.
#' @export
#'
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' last_sig <- ExtractLastTestSig(stat_obj)
#' print(last_sig)
#' plot_deg_boxplot(last_test_sig = last_sig, data = stat_obj@clean.data,group='SWAB',save_dir = "./")
#' plot_deg_boxplot(last_test_sig = last_sig,data = stat_obj@clean.data,group='SWAB',save_dir = "./",top_n = 3)
#' }
plot_deg_boxplot <- function(last_test_sig,
                                data,
                                control = 'health',
                                case = 'cancer',
                                top_n = 5,
                                palette_name = "Royal1",
                                name_identity = 'deg',
                                save_plots = TRUE,
                                save_dir = NULL,
                                plot_width = 5,
                                plot_height = 5,
                                base_size = 14,
                                title = "Violin Plot",
                                group_col = 'group',
                                pseudocount = 1) {
  
  cat("Data columns:", colnames(data), "\n")
  
  cat("Filtering significant features...\n")
  left <- last_test_sig[last_test_sig$change != 'Stable', ]
  
  for_label <- left %>%
    group_by(change) %>%
    arrange(desc(abs(logFC)), .by_group = TRUE) %>%
    top_n(n = top_n, wt = abs(logFC)) %>%
    dplyr::filter(!is.na(id))
  
  cat("Filtered features (for_label):\n")
  print(for_label)
  
  if (nrow(for_label) == 0) {
    cat("No significant features found. Exiting the function.\n")
    return(NULL)
  }
  
  n <- nrow(for_label)
  n <- min(n, top_n * 2)
  
  selected_ids <- for_label$id[seq_len(n)]
  selected_columns <- c(selected_ids, group_col)
  cat("Selected columns:", selected_columns, "\n")
  missing_columns <- setdiff(selected_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("Missing columns in data:", paste(missing_columns, collapse = ", ")))
  }
  
  box_test <- data[, selected_columns]
  box_test <- pivot_longer(box_test, cols = all_of(selected_ids), values_to = 'value', names_to = 'id')
  box_test[[group_col]] <- as.factor(box_test[[group_col]])
  
  cat("Calculating significance differences...\n")
  p_values <- box_test %>%
    group_by(id) %>%
    summarise(p_value = wilcox.test(value ~ .data[[group_col]])$p.value, .groups = 'drop') %>%
    mutate(significance = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  y_positions <- box_test %>%
    group_by(id) %>%
    summarise(y_position = max(log10(abs(value) + pseudocount)) + 0.1)
  
  p_values <- p_values %>%
    left_join(y_positions, by = "id")
  
  box_test <- left_join(box_test, p_values)
  index <- box_test$id[order(box_test$y_position, decreasing = FALSE)]
  box_test$id <- factor(box_test$id, levels = unique(index))
  
  cat("Drawing violin plot...\n")
  p <- ggplot(data = box_test, aes(x = id, y = log10(value + pseudocount), fill = .data[[group_col]])) +
    geom_boxplot(outlier.alpha = 0) +
    geom_text(data = box_test, aes(x = id, y = y_position, label = significance),
              size = 3, color = "black") +
    scale_fill_manual(values = wes_palette(palette_name)) +
    ggprism::theme_prism(base_size = base_size) +
    theme(
      legend.key.size = unit(0.4, 'cm'),
      legend.text = element_text(size = base_size * 0.5),
      legend.title = element_text(size = base_size * 0.6),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    labs(title = title) +
    labs(x = "Features", y = paste0("Log10(Value + ", pseudocount, ")"))
  
  if (save_plots) {
    if (is.null(save_dir)) {
      stop("'save_dir' cannot be NULL when 'save_plots' is TRUE. Please provide a valid directory path.")
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    output_path <- file.path(save_dir, "violinplot.pdf")
    ggsave(filename = output_path, plot = p, height = plot_height, width = plot_width, device = "pdf")
    
    cat("Plot saved to:", output_path, "\n")
  }
  
  return(p)
}

#' Generate Boxplot Plot for Differential Features (fixed)
#'
#' @param object Stat object or data frame.
#' @param control Control group label.
#' @param case Case group label.
#' @param top_n Number of top features to show.
#' @param palette_name Colour palette.
#' @param name_identity Analysis identifier.
#' @param data_type "clean" or "scale".
#' @param save_dir Output directory. If NULL, auto‑creates one under
#'   `./figures/deg_info/`.
#' @param save_plots Logical. Save the plot? Default TRUE.
#' @param plot_width,plot_height Plot dimensions (inches).
#' @param base_size Base font size.
#' @param group_col Group column name. If NULL, read from Stat object.
#' @param pseudocount Pseudocount for log10 transform.
#'
#' @return Updated Stat object or ggplot.
#' @export
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' VarFeature_boxplot(object = stat_obj)
#' }
VarFeature_boxplot <- function(object,
                                  control = 'health',
                                  case = 'cancer',
                                  top_n = 5,
                                  palette_name = "Royal1",
                                  name_identity = 'deg',
                                  data_type = "clean",
                                  save_dir = NULL,
                                  save_plots = TRUE,
                                  plot_width = 5,
                                  plot_height = 5,
                                  base_size = 10,
                                  group_col = NULL,
                                  pseudocount = 1) {
  # ---- Handle default save_dir ----
  if (is.null(save_dir) && save_plots) {
    if (exists("get_output_dir")) {
      save_dir <- get_output_dir("Figures", "deg_info")
    } else {
      save_dir <- file.path(".", "figures", "deg_info")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  }
  
  # ---- Extract data and group_col ----
  if (inherits(object, "Stat")) {
    last_test_sig <- ExtractLastTestSig(object)
    
    # If group_col not provided, read from Stat slot
    if (is.null(group_col)) {
      group_col <- slot(object, "group_col")
    }
    if (length(group_col) == 0) group_col <- NULL
    
    cat("Extracting data from 'Stat' object...\n")
    data <- if (data_type == "clean") ExtractCleanData(object) else ExtractScaleData(object)
    
    if (!is.null(group_col) && !group_col %in% colnames(data)) {
      warning("group_col '", group_col, "' not found in data. Setting to NULL.")
      group_col <- NULL
    }
  } else if (is.data.frame(object)) {
    last_test_sig <- object
    data          <- object
    if (is.null(group_col)) group_col <- "group"
  } else {
    stop("Input must be an object of class 'Stat' or a data frame.")
  }
  
  if (is.null(last_test_sig) || nrow(last_test_sig) == 0)
    stop("No valid data found in last_test_sig.")
  
  # ---- Call the core plotting function ----
  boxplot <- plot_deg_boxplot(
    last_test_sig = last_test_sig,
    data = data,
    control = control,
    case = case,
    top_n = top_n,
    name_identity = name_identity,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    palette_name = palette_name,
    group_col = group_col,
    pseudocount = pseudocount
  )
  
  print(boxplot)
  
  if (inherits(object, "Stat")) {
    object@var.result[["boxplot"]] <- boxplot
    return(object)
  }
  
  return(boxplot)
}

####----here
#' Plot ROC Curve for Differential Expression Gene Analysis
#'
#' This function generates a ROC curve for the top features identified in differential expression analysis.
#' It filters the significant features, calculates the ROC curve, and visualizes the AUC (Area Under the Curve)
#' for each feature. The plot can be saved to a specified directory.
#'
#' @importFrom dplyr %>% group_by arrange desc top_n
#' @importFrom pROC roc ggroc
#' @param deg_test A data frame containing the results of differential expression analysis.
#'                 It should have a column `change` indicating whether a feature is "Stable" or changed.
#' @param mat_test A data frame containing the expression values of features across samples.
#'                 It should have columns representing features (e.g., genes) and rows representing samples.
#' @param group_col The column in `mat_test` that indicates the grouping of samples (e.g., control vs case).
#'                  Default is `'group'`.
#' @param control The label used to indicate the control group in the data. Default is `'health'`.
#' @param case The label used to indicate the case group in the data. Default is `'lung cancer'`.
#' @param top_n The number of top features to consider based on log fold change. Default is 5.
#' @param palette_name The name of the palette to be used for the ROC curves.
#' @param name_identity The name or identifier used for the features in the `deg_test` data frame. Default is `'deg'`.
#' @param save_plots A logical value indicating whether the plot should be saved. Default is `TRUE`.
#' @param save_dir The directory path where the plot will be saved. Default is `get_output_dir("StatObject", "deg_info")`.
#' @param plot_width The width of the saved plot in cm. Default is 5.
#' @param plot_height The height of the saved plot in cm. Default is 5.
#' @param base_size The base font size for the plot. Default is 14.
#' @param title The title of the plot. Default is `'ROC Curve'`.
#' @param normalize_auc Logical, if `TRUE`, for any feature with AUC < 0.5,
#'   the predictor values are negated (multiplied by -1) and the ROC curve is
#'   recalculated so that the new AUC becomes 1 - original AUC, ensuring all
#'   curves lie above the diagonal. Default is `TRUE`.
#' @returns A `ggplot` object representing the ROC curve for the top features.
#'          If `save_plots` is `TRUE`, the plot is saved to the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' last_sig <- ExtractLastTestSig(stat_obj)
#' print(last_sig)
#' plot_deg_Roc_plot(deg_test = last_sig, mat_test = stat_obj@clean.data,group_col = 'SWAB',save_dir = "./")
#' }
plot_deg_Roc_plot <- function(deg_test,
                              mat_test,
                              group_col = 'group',
                              control = 'health',
                              case = 'lung cancer',
                              top_n = 5,
                              palette_name = "Royal1",
                              name_identity = 'deg',
                              save_plots = TRUE,
                              save_dir = NULL,
                              plot_width = 5,
                              plot_height = 5,
                              base_size = 10,
                              title = 'ROC Curve',
                              normalize_auc = TRUE) {  
  left <- deg_test[deg_test$change != 'Stable', ]
  for_label <- left %>%
    group_by(change) %>%
    arrange(desc(abs(logFC)), .by_group = TRUE) %>%
    top_n(n = top_n, wt = abs(logFC)) %>%
    dplyr::filter(!is.na(id))
  
  cat("Filtered features (for_label):\n")
  print(for_label)
  
  if (nrow(for_label) == 0) {
    cat("No significant features found. Exiting the function.\n")
    return(NULL)
  }
  
  if (!all(for_label$id %in% colnames(mat_test))) {
    stop("Some features in for_label$id are not found in mat_test.")
  }
  
  roc_tem <- list()
  auc_tem <- numeric()
  
  for (j in 1:min(nrow(for_label), top_n * 2)) {
    id <- for_label$id[j]
    cat("Processing feature:", id, "\n")
    
    if (!is.numeric(mat_test[[id]])) {
      stop(paste("Feature", id, "is not numeric."))
    }
    
    pred <- mat_test[[id]]
    roc_obj <- roc(mat_test[[group_col]], pred)
    
    if (normalize_auc && roc_obj$auc < 0.5) {
      cat("  AUC < 0.5, reversing predictor for", id, "\n")
      roc_obj <- roc(mat_test[[group_col]], -pred)
    }
    
    roc_tem[[j]] <- roc_obj
    auc_tem[j] <- round(roc_obj$auc, 2)
  }
  
  n_curves <- length(roc_tem)
  base_colors <- tryCatch(
    as.vector(wes_palette(n = min(n_curves, 5), name = palette_name, type = "discrete")),
    error = function(e) as.vector(wes_palette(name = palette_name))
  )
  if (n_curves > length(base_colors)) {
    colors <- colorRampPalette(base_colors)(n_curves)
  } else {
    colors <- base_colors[seq_len(n_curves)]
  }
  
  p <- ggroc(roc_tem) +
    ggprism::theme_prism(base_size = base_size) +
    theme(legend.position = 'none') +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
    scale_color_manual(values = colors) +
    labs(title = title)
  
  for (j in 1:length(auc_tem)) {
    p <- p + annotate("text", x = 0.25, y = 0.01 + (j - 1) * 0.05,
                      label = paste(for_label$id[j], "AUC:", auc_tem[j]),
                      size = 5, color = colors[j])
  }
  
  if (save_plots) {
    if (is.null(save_dir)) {
      stop("'save_dir' cannot be NULL when 'save_plots' is TRUE. Please provide a valid directory path.")
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    output_path <- file.path(save_dir, "roc_plot.pdf")
    ggsave(filename = output_path, plot = p, height = plot_height, width = plot_width, device = "pdf")
    cat("Saved plots to: ", output_path, "\n")
  }
  
  return(p)
}


#' Generate ROC Curve for Variable Features 
#'
#' @param object Stat object, or named list with $last_test_sig and $data.
#' @param group_col Group column name. If NULL, auto‑detect from Stat object.
#' @param control Control group label.
#' @param case Case group label.
#' @param top_n Number of top features to show.
#' @param palette_name Colour palette.
#' @param name_identity Analysis identifier.
#' @param data_type "clean" or "scale".
#' @param save_dir Output directory. If NULL and save_plots=TRUE, will use
#'   "./figures/deg_info" or get_output_dir() if available.
#' @param save_plots Logical. Save the plot? Default TRUE.
#' @param plot_width,plot_height Plot dimensions (inches).
#' @param base_size Base font size.
#' @param normalize_auc Logical, if `TRUE`, for any feature with AUC < 0.5,
#'   the predictor values are negated (multiplied by -1) and the ROC curve is
#'   recalculated so that the new AUC becomes 1 - original AUC, ensuring all
#'   curves lie above the diagonal. Default is `TRUE`.
#' @return Updated Stat object or ggplot.
#' @export
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' VarFeature_ROC(stat_obj)
#' }
VarFeature_ROC <- function(object,
                           group_col = NULL,
                           control = 'health',
                           case = 'cancer',
                           top_n = 5,
                           palette_name = "Royal1",
                           name_identity = 'deg',
                           data_type = "clean",
                           save_dir = NULL,
                           save_plots = TRUE,
                           plot_width = 5,
                           plot_height = 5,
                           base_size = 10,
                           normalize_auc = TRUE) {
  
  # ---- Handle default save_dir ----
  if (is.null(save_dir) && save_plots) {
    if (exists("get_output_dir")) {
      save_dir <- get_output_dir("Figures", "deg_info")
    } else {
      save_dir <- file.path(".", "figures", "deg_info")
    }
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  }
  
  # ---- Extract data and group_col ----
  if (inherits(object, "Stat")) {
    last_test_sig <- ExtractLastTestSig(object)
    
    if (is.null(group_col)) {
      group_col <- slot(object, "group_col")
    }
    if (length(group_col) == 0) group_col <- NULL
    
    data <- if (data_type == "clean") {
      ExtractCleanData(object)
    } else {
      ExtractScaleData(object)
    }
    
    if (!is.null(group_col) && !group_col %in% colnames(data)) {
      stop(paste("Group column", group_col, "not found in data. Available columns:",
                 paste(colnames(data), collapse = ", ")))
    }
    
  } else if (is.list(object) && !is.data.frame(object)) {
    if (is.null(object$last_test_sig) || is.null(object$data)) {
      stop("When passing a list to VarFeature_ROC, it must contain both ",
           "$last_test_sig (DEG result data frame) and $data (expression matrix).")
    }
    last_test_sig <- object$last_test_sig
    data          <- object$data
    
    if (is.null(group_col)) group_col <- "group"
    if (!group_col %in% colnames(data)) {
      stop(paste("Group column", group_col, "not found in data frame."))
    }
    
  } else {
    stop("Input must be a 'Stat' object or a named list with $last_test_sig and $data. ",
         "A plain DEG result data frame is not accepted here because the ROC ",
         "calculation requires the raw expression matrix as well.")
  }
  
  if (is.null(last_test_sig) || nrow(last_test_sig) == 0) {
    stop("No valid data found in last_test_sig.")
  }
  
  # ---- Generate ROC plot ----
  roc_plot <- plot_deg_Roc_plot(
    deg_test = last_test_sig,
    mat_test = data,
    group_col = group_col,
    control = control,
    case = case,
    top_n = top_n,
    palette_name = palette_name,
    name_identity = name_identity,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    normalize_auc = normalize_auc
  )
  print(roc_plot)
  if (inherits(object, "Stat")) {
    object@var.result[["Rocplot"]] <- roc_plot
    return(object)
  }
  
  return(roc_plot)
}
