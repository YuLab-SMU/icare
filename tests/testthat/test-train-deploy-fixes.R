# =============================================================================
# Train/deploy chain fixes
# =============================================================================

make_stat <- function(seed = 9, n = 40) {
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

make_trained <- function(seed = 2, n = 60) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n), f4 = rnorm(n))
  lab  <- rep(c("HC", "Case"), each = n / 2)
  feat$f1 <- feat$f1 + ifelse(lab == "HC", 1.4, 0)
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = factor(lab), row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  obj <- suppressMessages(suppressWarnings({
    o <- stat_miss_processed(obj, impute_method = "median_mode", save_data = FALSE)
    o <- stat_normalize_process(o, method = "auto", save_data = FALSE)
    o
  }))
  m <- CreateModelObject(object = obj, group_col = "group")
  idx <- rownames(m@clean.df)
  set.seed(seed)
  tr <- sample(idx, floor(0.7 * length(idx)))
  m@filtered.set <- list(training = m@clean.df[tr, , drop = FALSE],
                         testing  = m@clean.df[setdiff(idx, tr), , drop = FALSE])
  m@split.data  <- m@filtered.set
  ModelTrainAnalysis(m, methods = c("glm"), metric_selection = "ROC",
                     seed = seed, save_plots = FALSE)
}

test_that("ModelTrainAnalysis records preprocess_columns on the main path", {
  # deploy_clinlab_app() reads this field; it used to be set only by
  # PreprocessModelData, so the deployment fallback chain was dead code.
  m <- make_trained()
  cols <- m@process.info$preprocess_columns
  expect_false(is.null(cols))
  expect_true(length(cols) > 0)
  expect_false("group" %in% cols)          # the outcome must not be a feature
  expect_true(all(cols %in% colnames(m@clean.df)))

  d <- ModelDeployment(m)
  expect_identical(d$preprocess_columns, cols)
})

test_that("one_hot_encode() never one-hot-encodes the Train_Model outcome", {
  # The default group_col = "group" did not match an outcome named anything
  # else, so the label was treated as an ordinary categorical and encoded away,
  # producing an invalid Train_Model.
  obj <- make_trained()
  m <- obj

  out <- suppressMessages(one_hot_encode(m, save_data = FALSE))
  expect_true("group" %in% colnames(out@clean.df))
  expect_silent(validObject(out))

  # Same guarantee for an outcome with a non-default name.
  obj2 <- make_stat(seed = 9)
  m2 <- CreateModelObject(object = obj2, group_col = "diagnosis")
  out2 <- suppressMessages(one_hot_encode(m2, save_data = FALSE))
  expect_true("diagnosis" %in% colnames(out2@clean.df))
  expect_silent(validObject(out2))
})

test_that(".extract_xy rejects an outcome-less Train_Model with a clear error", {
  # cd[[NULL]] produced base R's "attempt to select less than one element in
  # get1index"; group_col = NULL is the class prototype, so this is reachable.
  m <- CreateModelObject(data = mtcars)
  expect_error(icare:::.extract_xy(m), "no group_col")
})
