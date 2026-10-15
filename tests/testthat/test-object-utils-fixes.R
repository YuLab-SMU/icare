# =============================================================================
# Object-utils fixes: metadata preservation across merge/rename/inspect
# =============================================================================

make_stat <- function(seed = 8, n = 40) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n))
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(diagnosis = factor(rep(c("HC", "Case"), each = n / 2)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, clean.data = feat, info.data = info,
                          group_col = "diagnosis", na.action = "allow")
  suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                       save_data = FALSE))
}

make_tm <- function(seed = 2, n = 40, rows = NULL) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n))
  rownames(feat) <- if (is.null(rows)) paste0("S", seq_len(n)) else rows
  info <- data.frame(diagnosis = factor(rep(c("HC", "Case"), each = n / 2)),
                     age = round(rnorm(n, 60, 8)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, clean.data = feat, info.data = info,
                          group_col = "diagnosis", na.action = "allow")
  suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                       save_data = FALSE))
  CreateModelObject(object = obj, group_col = "diagnosis")
}

test_that("MergeObjects keeps info.data and the label name for Train_Model", {
  # CreateModelObject(data = ...) rebuilds from an empty info.data, silently
  # dropping the merged metadata.
  m1 <- make_tm(seed = 2)
  m2 <- make_tm(seed = 3, rows = paste0("T", 1:40))

  merged <- MergeObjects(list(m1, m2))
  expect_s4_class(merged, "Train_Model")
  expect_equal(merged@group_col, "diagnosis")
  expect_equal(nrow(merged@info.data),
               nrow(m1@info.data) + nrow(m2@info.data))
  expect_true("diagnosis" %in% colnames(merged@info.data))
  expect_true("age" %in% colnames(merged@info.data))
  expect_silent(validObject(merged))
})

make_px <- function(seed = 5, n = 40, rows = NULL) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n))
  rownames(feat) <- if (is.null(rows)) paste0("S", seq_len(n)) else rows
  info <- data.frame(diagnosis = factor(rep(c("HC", "Case"), each = n / 2)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, clean.data = feat, info.data = info,
                          group_col = "diagnosis", na.action = "allow")
  obj <- suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                              save_data = FALSE))
  info <- obj@info.data
  info$time <- rexp(n, 0.05); info$status <- rbinom(n, 1, 0.6)
  obj@info.data <- info
  CreatePrognosiXObject(object = obj, group_col = "diagnosis",
                        time_col = "time", status_col = "status")
}

test_that("MergeObjects carries group_col for PrognosiX", {
  px  <- make_px(seed = 5)
  px2 <- make_px(seed = 6, rows = paste0("T", 1:40))

  merged <- MergeObjects(list(px, px2))
  expect_equal(merged@group_col, "diagnosis")   # used to be NULL
})

test_that("RenameFeatures updates group_col when the outcome is renamed", {
  m <- make_tm(seed = 6)
  renamed <- RenameFeatures(m, "diagnosis", "label")
  expect_equal(renamed@group_col, "label")
  expect_true("label" %in% colnames(renamed@clean.df))
  expect_false("diagnosis" %in% colnames(renamed@clean.df))
  expect_silent(validObject(renamed))

  # A collision with an existing feature is rejected rather than duplicated.
  expect_error(RenameFeatures(m, "f1", "diagnosis"), "already exist")
})

test_that("InspectObject surfaces the new slots", {
  m <- make_tm(seed = 7)
  out <- capture.output(InspectObject(m))
  expect_true(any(grepl("info.data", out)))

  obj <- make_stat(seed = 8)
  info <- obj@info.data
  info$time <- rexp(nrow(info), 0.05); info$status <- rbinom(nrow(info), 1, 0.6)
  obj@info.data <- info
  px <- CreatePrognosiXObject(object = obj, group_col = "diagnosis",
                              time_col = "time", status_col = "status")
  out2 <- capture.output(InspectObject(px))
  expect_true(any(grepl("group_col", out2)))
})

test_that("one_hot_encode keeps clean/survival/info/sub frames in step", {
  obj <- make_stat(seed = 10)
  info <- obj@info.data
  info$time <- rexp(nrow(info), 0.05); info$status <- rbinom(nrow(info), 1, 0.6)
  obj@info.data <- info
  px <- CreatePrognosiXObject(object = obj, group_col = "diagnosis",
                              time_col = "time", status_col = "status")
  px@sub.data <- data.frame(note = rep("x", nrow(px@clean.data)),
                            row.names = rownames(px@clean.data))

  # Introduce missing survival endpoints so the drop path fires.
  info <- px@info.data
  info$time[1:5] <- NA
  px@info.data <- info

  px2 <- suppressMessages(one_hot_encode(px, group_col = "diagnosis", save_data = FALSE))
  n_row <- nrow(px2@clean.data)
  expect_equal(nrow(px2@info.data), n_row)
  expect_equal(nrow(px2@survival.data), n_row)
  expect_equal(nrow(px2@sub.data), n_row)      # used to keep all rows
  expect_identical(rownames(px2@clean.data), rownames(px2@sub.data))
})
