# =============================================================================
# Stratified splitting must be a function of the LABELS, not their storage type
# =============================================================================
# Regression from the 1.1.0 -> user-report incident: the same cohort produced
# different train/test partitions (421/180 vs 422/179) depending on whether the
# outcome was stored as numeric 0/1 or as a factor with readable labels, and
# the resulting NMF rank choice flipped from K=2 to K=5.
#
# The mechanism: caret::createDataPartition() branches on is.numeric(y) and
# bins numeric input with cut(y, unique(quantile(...))). A 0/1 label collapses
# to ONE bin, so the "stratified" split silently degenerated into a single
# random sample -- stratification was requested and not delivered. A factor
# splits by level, which is what stratify_by actually means.
#
# These tests pin both halves: (a) the partition must not depend on the
# storage class, and (b) stratification must really be stratified.

.make_subtyping <- function(labels, seed = 3, n = 60) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n))
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = labels, row.names = rownames(feat))
  stat_obj <- CreateStatObject(raw.data = feat, clean.data = feat,
                               info.data = info,
                               group_col = "group", na.action = "allow")
  ConvertObject(stat_obj, to = "Subtyping")
}

test_that("the split partition does not depend on the label storage type", {
  n <- 60
  binary <- rep(c(0, 1), each = n / 2)

  # The same underlying labels, three storage classes a user might hand over.
  variants <- list(
    numeric  = binary,
    character = as.character(binary),
    factor   = factor(binary, levels = c(0, 1), labels = c("HC", "Case"))
  )

  partitions <- lapply(names(variants), function(nm) {
    sub <- .make_subtyping(variants[[nm]])
    set.seed(123)
    sp <- SplitSubtypingObject(sub, p = 0.7, stratify_by = "group")
    rownames(sp$train@clean.data)
  })
  names(partitions) <- names(variants)

  # Essence: identical sample sets, not merely identical counts. The bug this
  # replaces produced the same 422/179 count for character and factor but a
  # DIFFERENT 421/180 partition for numeric -- so a count check would pass.
  expect_identical(partitions$numeric, partitions$factor)
  expect_identical(partitions$character, partitions$factor)
})

test_that("stratify_by really stratifies (the old numeric path did not)", {
  n <- 60
  binary <- rep(c(0, 1), each = n / 2)
  sub <- .make_subtyping(binary)

  set.seed(123)
  sp <- SplitSubtypingObject(sub, p = 0.7, stratify_by = "group")

  train_tab <- table(sp$train@info.data$group)
  test_tab  <- table(sp$test@info.data$group)

  # 70% of each level, rounded -- not "one random sample of the whole cohort".
  expect_equal(as.integer(train_tab[["0"]]), ceiling(n / 2 * 0.7))
  expect_equal(as.integer(train_tab[["1"]]), ceiling(n / 2 * 0.7))
  expect_equal(as.integer(test_tab[["0"]]), n / 2 - ceiling(n / 2 * 0.7))
  expect_equal(as.integer(test_tab[["1"]]), n / 2 - ceiling(n / 2 * 0.7))

  # Both levels must actually appear on both sides of the split.
  expect_true(all(train_tab > 0))
  expect_true(all(test_tab > 0))
})

test_that("a single-level stratifier falls back to an unstratified split", {
  n <- 60
  sub <- .make_subtyping(rep(1L, n))   # one level only

  expect_warning(
    sp <- SplitSubtypingObject(sub, p = 0.7, stratify_by = "group"),
    "distinct value"
  )
  # The stratifier is dropped, so the sample size follows the unstratified
  # path. (Assert the invariant that matters -- the split still covers the
  # whole cohort exactly once -- rather than caret's exact rounding, which
  # differs between its numeric and factor branches.)
  expect_gt(nrow(sp$train@clean.data), 0)
  expect_equal(nrow(sp$train@clean.data) + nrow(sp$test@clean.data), n)
  expect_equal(
    length(intersect(rownames(sp$train@clean.data), rownames(sp$test@clean.data))),
    0
  )
})

test_that("NAs in the stratifier do not silently produce a degenerate split", {
  n <- 60
  labels <- rep(c("HC", "Case"), each = n / 2)
  labels[1:3] <- NA
  sub <- .make_subtyping(labels)

  # na.omit drops the NA rows before counting levels, so the split still
  # stratifies over the levels that exist rather than crashing or silently
  # falling back.
  set.seed(123)
  sp <- SplitSubtypingObject(sub, p = 0.7, stratify_by = "group")
  expect_true(nrow(sp$train@clean.data) > 0)
  expect_equal(nrow(sp$train@clean.data) + nrow(sp$test@clean.data), n)
})

# -----------------------------------------------------------------------------
# Splitting an ALREADY-normalised object must keep the scaling
# -----------------------------------------------------------------------------
# Regression: CreateSubtypingObject() was called without `scale.data`, so
# every split silently reset that slot to 0x0. A split performed AFTER
# normalisation then handed downstream clustering an empty scale.data, and
# (because .extract_data() falls back to clean.data when scale.data is empty)
# the analysis quietly ran on UNSCALED features with no message at all.

test_that("splitting a normalised object keeps scale.data and group_col", {
  n <- 60
  binary <- rep(c(0, 1), each = n / 2)
  sub <- .make_subtyping(binary)

  # Normalise first, so scale.data is non-empty at split time.
  sub <- Sub_normalize_process(sub, normalize_method = "min_max", group_col = "group")
  expect_true(nrow(sub@scale.data) > 0)

  sp <- SplitSubtypingObject(sub, p = 0.7, stratify_by = "group")

  # Essence: the scaling SURVIVES the split, row-aligned to each side.
  expect_true(nrow(sp$train@scale.data) > 0 && nrow(sp$test@scale.data) > 0)
  expect_identical(rownames(sp$train@scale.data), rownames(sp$train@clean.data))
  expect_identical(rownames(sp$test@scale.data),  rownames(sp$test@clean.data))

  # And the values are the ORIGINAL scaling of those samples, not recomputed.
  expect_equal(as.matrix(sp$train@scale.data),
               as.matrix(sub@scale.data[rownames(sp$train@clean.data), , drop = FALSE]))
  expect_equal(as.matrix(sp$test@scale.data),
               as.matrix(sub@scale.data[rownames(sp$test@clean.data), , drop = FALSE]))

  # group_col must survive too, or the label becomes unreachable downstream.
  expect_identical(sp$train@group_col, "group")
  expect_identical(sp$test@group_col,  "group")
})

test_that("a split with no scale.data still works (the normal order)", {
  n <- 60
  sub <- .make_subtyping(rep(c("HC", "Case"), each = n / 2))
  expect_equal(nrow(sub@scale.data), 0)   # not normalised yet

  sp <- SplitSubtypingObject(sub, p = 0.7, stratify_by = "group")
  expect_equal(nrow(sp$train@scale.data), 0)
  expect_equal(nrow(sp$train@clean.data) + nrow(sp$test@clean.data), n)
})

test_that("a scale.data that does not cover clean.data is an error, not a silent gap", {
  n <- 60
  sub <- .make_subtyping(rep(c(0, 1), each = n / 2))
  sub <- Sub_normalize_process(sub, normalize_method = "min_max", group_col = "group")

  # Drop one sample from scale.data only: a positional/naive subset would
  # silently produce a short frame and misalign every later step.
  broken <- sub
  broken@scale.data <- sub@scale.data[-1, , drop = FALSE]
  expect_error(
    SplitSubtypingObject(broken, p = 0.7, stratify_by = "group"),
    "does not cover"
  )
})

test_that("a Date stratifier is split within levels, not per distinct date", {
  n <- 60
  # A Date column is NOT numeric, so without coercion it reaches caret's
  # factor branch as dates and caret stratifies on every distinct value.
  dates <- as.Date("2020-01-01") + rep(0:1, each = n / 2)
  sub <- .make_subtyping(dates)

  set.seed(123)
  sp <- SplitSubtypingObject(sub, p = 0.7, stratify_by = "group")
  # Two date levels -> both sides keep the 70/30 ratio, and no level is empty.
  train_tab <- table(sp$train@info.data$group)
  expect_length(train_tab, 2)
  expect_true(all(train_tab > 0))
  expect_equal(as.integer(train_tab[1]), ceiling(n / 2 * 0.7))
})

test_that("a single-sample object is a clear error, not a caret message", {
  n <- 60
  sub <- .make_subtyping(rep(c(0, 1), each = n / 2))
  # Narrow to one sample: too small to split, and caret's own message
  # ("y must have at least 2 data points") names neither icare nor the cause.
  tiny <- CreateSubtypingObject(
    clean.data = sub@clean.data[1, , drop = FALSE],
    info.data  = sub@info.data[1, , drop = FALSE],
    group_col  = "group"
  )
  expect_error(SplitSubtypingObject(tiny, p = 0.7, stratify_by = "group"),
               "at least 2")
})

