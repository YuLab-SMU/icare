# =============================================================================
# Label TYPE must not change the analysis
# =============================================================================
# Before the label moved out of the feature matrix, whether a label behaved as a
# label or as a feature depended on its R type: a numeric 0/1 column was picked
# up by every "select numeric columns" step, while "0"/"1" was not. These tests
# pin the opposite property: the type carries no meaning.

build_stat <- function(labels, seed = 4, n = 60) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n), f4 = rnorm(n))
  # Give f1 a real group difference and leave the others null, so the DEG result
  # is a meaningful assertion rather than noise.
  feat$f1 <- feat$f1 + ifelse(as.character(labels) == as.character(labels)[1], 1.4, 0)
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = labels, row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  suppressMessages(suppressWarnings({
    obj <- stat_diagnose_variable_type(obj)
    obj <- stat_convert_variables(obj, save_data = FALSE)
    obj <- stat_miss_processed(obj, impute_method = "median_mode", save_data = FALSE)
    obj <- stat_normalize_process(obj, method = "auto", save_data = FALSE)
    obj
  }))
}

test_that("a character label survives the cleaning chain instead of becoming NA", {
  # stat_convert_variables() used to drop group_col, so the label was treated as
  # an ordinary variable and as.numeric()'d -- character labels became all-NA.
  obj <- build_stat(rep(c("HC", "Case"), each = 30))
  labels <- obj@clean.data
  expect_false("group" %in% colnames(labels))
  expect_false(anyNA(.stat_labels(obj)[[1]]))
  expect_setequal(unique(.stat_labels(obj)[[1]]), c("HC", "Case"))
})

test_that("numeric, character and factor labels give the same DEG ranking", {
  lab_chr <- rep(c("0", "1"), each = 30)
  lab_num <- rep(c(0, 1), each = 30)
  lab_fac <- factor(lab_chr, levels = c("0", "1"))

  res <- lapply(list(num = lab_num, chr = lab_chr, fac = lab_fac), function(l) {
    obj <- build_stat(l)
    # log2ratio of z-scored values legitimately warns about NaN; that is a
    # separate pre-existing wart and not what this test is about.
    obj <- suppressWarnings(suppressMessages(stat_var_feature(obj, p_threshold = 0.05)))
    r <- ExtractLastTestSig(obj)
    r[order(r$id), c("id", "p")]
  })

  expect_equal(res$num$id, res$chr$id)
  expect_equal(res$num$id, res$fac$id)
  expect_equal(res$num$p, res$chr$p, tolerance = 1e-12)
  expect_equal(res$num$p, res$fac$p, tolerance = 1e-12)
  expect_true("f1" %in% res$num$id)
})

test_that("the label column never appears in a numeric feature list", {
  for (l in list(rep(c(0, 1), each = 30), rep(c("HC", "Case"), each = 30))) {
    obj <- build_stat(l)
    for (frame in list(obj@raw.data, obj@clean.data, obj@scale.data)) {
      if (nrow(frame) == 0) next
      numeric_cols <- names(frame)[vapply(frame, is.numeric, logical(1))]
      expect_false("group" %in% numeric_cols)
    }
  }
})
