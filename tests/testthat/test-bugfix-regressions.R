# =============================================================================
# Regression tests for bugs fixed alongside the label-location refactor
# =============================================================================

make_scaled_stat <- function(seed = 11, n = 60) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n))
  feat$f1 <- feat$f1 + rep(c(1.2, 0), each = n / 2)
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = factor(rep(c("HC", "Case"), each = n / 2)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  suppressMessages(suppressWarnings({
    obj <- stat_miss_processed(obj, impute_method = "median_mode", save_data = FALSE)
    obj <- stat_normalize_process(obj, method = "auto", save_data = FALSE)
    obj
  }))
}

test_that("stat_compute_descriptive groups by a label that is not named 'group'", {
  # The worker's own default is the literal string "group", and the wrapper used
  # to discard the slot value, so any other label name silently produced
  # ungrouped descriptives.
  set.seed(2)
  n <- 40
  feat <- data.frame(value = rnorm(n))
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(diagnosis = rep(c("HC", "Case"), each = n / 2),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "diagnosis", na.action = "allow")

  obj <- suppressMessages(stat_compute_descriptive(obj))
  expect_false(is.null(obj@compute.descriptive$Group_Counts))
  expect_setequal(names(obj@compute.descriptive$Group_Counts), c("HC", "Case"))
})

test_that("compute_descriptive_stats() tolerates group_col = NULL", {
  # `NULL %in% colnames(d)` is logical(0), so `if (logical(0))` aborted.
  set.seed(3)
  d <- data.frame(a = rnorm(20), b = rnorm(20))
  expect_error(
    suppressMessages(compute_descriptive_stats(d, group_col = NULL)),
    NA
  )
})

test_that("batch_Wilcoxon_MultiClass accepts a label-free logfc_mat", {
  # stat_var_feature() hands over clean.data, which is label-free by contract;
  # subscripting the label out of it used to raise "subscript out of bounds".
  obj <- make_scaled_stat()
  m <- CreateModelObject(object = obj, group_col = "group")
  feat_only <- m@clean.df[, c("f1", "f2", "f3"), drop = FALSE]

  res <- expect_error(
    batch_Wilcoxon_MultiClass(m@clean.df, group_col = "group", logfc_mat = feat_only),
    NA
  )
  expect_true(nrow(res) >= 1)
})

test_that("one_hot_encode keeps clean.data, info.data and survival.data in step", {
  obj <- make_scaled_stat()
  info <- obj@info.data
  info$time <- rexp(nrow(info), 0.05)
  info$status <- rbinom(nrow(info), 1, 0.6)
  obj@info.data <- info
  px <- CreatePrognosiXObject(object = obj, group_col = "group",
                              time_col = "time", status_col = "status")

  info <- px@info.data
  info$time[1:6] <- NA
  px@info.data <- info

  px <- suppressMessages(one_hot_encode(px, group_col = "group", save_data = FALSE))

  expect_equal(nrow(px@clean.data), nrow(px@info.data))
  expect_equal(nrow(px@info.data), nrow(px@survival.data))
  expect_identical(rownames(px@clean.data), rownames(px@info.data))
  expect_silent(validObject(px))
})

test_that("DownsampleObject matches groups by row name, not by position", {
  obj <- make_scaled_stat(seed = 21, n = 60)
  # Reverse info.data row order: a positional match would pair the wrong labels.
  obj@info.data <- obj@info.data[rev(rownames(obj@info.data)), , drop = FALSE]

  d <- suppressMessages(DownsampleObject(obj, n = 5, per_group = TRUE,
                                         group_col = "group", seed = 1))
  labs <- .stat_labels(d)
  expect_equal(nrow(d@clean.data), 10)
  counts <- as.integer(table(labs[[1]]))
  expect_equal(counts, c(5L, 5L))
})

test_that("the survival univariate scan never tests the label as a feature", {
  set.seed(6)
  n <- 100
  feat <- data.frame(g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n), g4 = rnorm(n))
  rownames(feat) <- paste0("P", seq_len(n))
  info <- data.frame(group = sample(c(0, 1), n, TRUE),
                     time = rexp(n, 0.04),
                     status = rbinom(n, 1, 0.7),
                     age = round(rnorm(n, 62, 9)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  obj <- suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                              save_data = FALSE))
  px <- CreatePrognosiXObject(object = obj, group_col = "group",
                              time_col = "time", status_col = "status")

  res <- suppressMessages(surv_filter_features_clinical(px, include_clinical = TRUE))
  expect_false("group" %in% res$table$Feature)
})

test_that("CreatePrognosiXObject resolves normalised time/status names", {
  # time_col is normalised by CreateStatObject, and the constructor then looked
  # the ORIGINAL name up in info.data -- an opaque "replacement has 0 rows".
  set.seed(8)
  n <- 40
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n))
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = rep(c("HC", "Case"), each = n / 2),
                     row.names = rownames(feat))
  info[["Time To Event"]] <- rexp(n, 0.05)
  info[["Event Status"]] <- rbinom(n, 1, 0.6)

  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  obj <- suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                              save_data = FALSE))

  norm_time <- setdiff(colnames(obj@info.data), c("group", "Event_Status"))
  px <- expect_error(
    CreatePrognosiXObject(object = obj, group_col = "group",
                          time_col = norm_time, status_col = "Event_Status"),
    NA
  )
  expect_true(norm_time %in% colnames(px@survival.data))
})
