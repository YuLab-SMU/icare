# =============================================================================
# logFC strategy must match the actual data content
# =============================================================================
# Regression for: stat_var_feature() inferred the logFC strategy from the
# recorded normalisation method and hardcoded "log2ratio" on the clean-data
# path. When clean.data contained negative values (log-ratio or already
# centred input), log2() of a negative group mean is undefined, so every
# affected feature came back with logFC = NaN and only a bare
# "NaNs produced" warning. The batch_Wilcoxon() guard that would have caught
# this was dead code on the main path, because stat_var_feature() always
# supplies logfc_mat.

make_stat <- function(seed = 1, n = 60, shift = c(1.2, 0)) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n))
  feat$f1 <- feat$f1 + rep(shift, each = n / 2)
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = factor(rep(c("HC", "Case"), each = n / 2)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                       save_data = FALSE))
}

capture_warn <- function(expr) {
  w <- NULL
  withCallingHandlers(
    expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }
  )
  w
}

test_that("negative clean.data falls back to diff instead of producing NaN logFC", {
  obj <- make_stat()
  expect_true(min(obj@clean.data$f1) < 0)   # the condition that used to break

  obj <- suppressMessages(stat_var_feature(obj, p_threshold = 0.05))
  res <- ExtractLastTestSig(obj)

  expect_false(anyNA(res$logFC))
  # diff IS the group mean difference, so it must match an independent check.
  d <- obj@clean.data
  mu <- tapply(d$f1, icare:::.stat_labels(obj)[[1]], mean)
  expect_equal(res$logFC[res$id == "f1"], unname(mu["HC"] - mu["Case"]),
               tolerance = 1e-8)
  expect_identical(attr(res, "logfc_type"), "diff")
  expect_identical(attr(res, "logfc_data_src"), "clean.data")
})

test_that("an actionable warning replaces the bare 'NaNs produced'", {
  obj <- make_stat()
  warnings <- capture_warn(suppressMessages(stat_var_feature(obj, p_threshold = 0.05)))
  expect_true(any(grepl("negative values", warnings)))
  expect_false(any(warnings == "NaNs produced"))
})

test_that("all-positive data keeps log2ratio", {
  obj <- make_stat(seed = 2, shift = c(6, 4))   # shifts the whole distribution up
  obj@clean.data <- abs(obj@clean.data) + 1

  res <- ExtractLastTestSig(suppressMessages(stat_var_feature(obj, p_threshold = 0.05)))
  expect_identical(attr(res, "logfc_type"), "log2ratio")
  expect_false(anyNA(res$logFC))
})

test_that(".compute_log2_means names the offending feature instead of going NaN", {
  set.seed(3)
  d <- data.frame(f1 = c(-1, -2, 1, 2), f2 = c(1, 2, 3, 4),
                  group = rep(c("A", "B"), each = 2))
  expect_error(
    icare:::.compute_log2_means(d, group_col = "group"),
    "group means are negative"
  )
})

test_that("an empty significant-feature set returns a well-formed empty result", {
  # Regression: with no feature passing the threshold, `shared` became empty,
  # logFC lost its `id` column and merge() aborted with an opaque
  # "'by' must specify a uniquely valid column".
  set.seed(4)
  n <- 40
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n))   # no real group difference
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = rep(c("HC", "Case"), each = n / 2),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  obj <- suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                              save_data = FALSE))

  obj <- expect_warning(
    suppressMessages(stat_var_feature(obj, p_threshold = 1e-12)),
    "No features passed"
  )
  res <- ExtractLastTestSig(obj)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
})

test_that("a row-dropped object stays aligned and still runs the full pipeline", {
  # info.data is append-only: dropping rows from clean.data must leave labels
  # aligned by name and the object usable downstream.
  set.seed(5)
  n <- 40
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n))
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = rep(c("HC", "Case"), each = n / 2),
                     row.names = rownames(feat))
  feat$f1[1:8] <- NA; feat$f2[9:16] <- NA

  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  obj <- suppressMessages(stat_miss_remove(obj, miss_threshold = 25, save_data = FALSE))

  expect_true(all(rownames(obj@clean.data) %in% rownames(obj@info.data)))
  labels <- icare:::.stat_labels(obj)
  expect_equal(nrow(labels), nrow(obj@clean.data))
  expect_silent(validObject(obj))

  obj <- suppressMessages(suppressWarnings({
    x <- stat_miss_processed(obj, impute_method = "median_mode", save_data = FALSE)
    x <- stat_normalize_process(x, method = "auto", save_data = FALSE)
    x <- stat_var_feature(x, p_threshold = 1e-12)   # force the empty path
    x
  }))
  expect_s3_class(ExtractLastTestSig(obj), "data.frame")
})
