# =============================================================================
# Label location contract
# =============================================================================
# The label must live in info.data and nowhere else. These tests pin that
# contract, because the failure mode it prevents (a numeric label silently
# treated as a feature) produces plausible-looking but wrong figures rather
# than an error.

make_feat <- function(n = 30, seed = 1) {
  set.seed(seed)
  d <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n))
  rownames(d) <- paste0("S", seq_len(n))
  d
}

make_info <- function(feat, labels) {
  data.frame(group = labels, row.names = rownames(feat))
}

test_that("a label supplied via info.data is accepted and stays out of the matrix", {
  feat <- make_feat()
  info <- make_info(feat, rep(c("HC", "Case"), each = 15))

  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")

  expect_s4_class(obj, "Stat")
  expect_identical(obj@group_col, "group")
  expect_false("group" %in% colnames(obj@raw.data))
  expect_false("group" %in% colnames(obj@clean.data))
  expect_true("group" %in% colnames(obj@info.data))
  expect_false("group" %in% obj@meta.featurename)
})

test_that("a label inside the feature matrix is a hard error with no migration", {
  feat <- make_feat()
  lab <- rep(c("HC", "Case"), each = 15)
  info <- data.frame(row.names = rownames(feat))

  expect_error(
    CreateStatObject(raw.data = cbind(feat, group = lab), info.data = info,
                     group_col = "group", na.action = "allow"),
    "present in the feature matrix"
  )
  # clean.data is checked too, not just raw.data
  expect_error(
    CreateStatObject(clean.data = cbind(feat, group = lab), info.data = info,
                     group_col = "group", na.action = "allow"),
    "present in the feature matrix"
  )
})

test_that("an unknown group_col reports the columns that do exist", {
  feat <- make_feat()
  info <- make_info(feat, rep(c("HC", "Case"), each = 15))

  expect_error(
    CreateStatObject(raw.data = feat, info.data = info, group_col = "nope"),
    "not in info.data"
  )
})

test_that("group_col = NULL is allowed: not every analysis compares groups", {
  feat <- make_feat()
  obj <- CreateStatObject(raw.data = feat, group_col = NULL, na.action = "allow")

  expect_null(obj@group_col)
  # info.data is still synthesised with the sample row names, but carries no columns.
  expect_equal(ncol(obj@info.data), 0)
  expect_equal(nrow(obj@info.data), nrow(feat))
})

test_that("group_col resolves through column-name normalisation", {
  feat <- make_feat()
  info <- data.frame(row.names = rownames(feat))
  info[["Diag Group"]] <- rep(c("HC", "Case"), each = 15)

  obj <- suppressMessages(
    CreateStatObject(raw.data = feat, info.data = info,
                     group_col = "Diag Group", na.action = "allow")
  )

  expect_false(identical(obj@group_col, "Diag Group"))
  expect_true(obj@group_col %in% colnames(obj@info.data))
  expect_setequal(unique(.stat_labels(obj)[[1]]), c("HC", "Case"))
})

test_that("the S4 validity method rejects a label that leaks back into clean.data", {
  feat <- make_feat()
  info <- make_info(feat, rep(c("HC", "Case"), each = 15))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")

  obj@clean.data <- cbind(obj@raw.data, group = info$group)
  expect_error(validObject(obj), "must not appear in clean.data")
})

test_that("info.data may keep rows that clean.data has dropped", {
  # info.data is append-only: dropping samples from clean.data must not
  # invalidate the object, and labels must still align by row name.
  feat <- make_feat()
  info <- make_info(feat, rep(c("HC", "Case"), each = 15))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  obj@clean.data <- obj@raw.data[1:10, , drop = FALSE]

  expect_silent(validObject(obj))
  labels <- .stat_labels(obj)
  expect_equal(nrow(labels), 10)
  expect_identical(rownames(labels), rownames(obj@clean.data))
})

test_that("legacy objects that keep the label in the matrix are migrated on read", {
  feat <- make_feat()
  lab <- rep(c("HC", "Case"), each = 15)
  info <- data.frame(row.names = rownames(feat))
  legacy <- CreateStatObject(raw.data = feat, info.data = info,
                             group_col = NULL, na.action = "allow")
  # Recreate the pre-1.1 layout by hand.
  legacy@group_col <- "group"
  legacy@raw.data <- cbind(legacy@raw.data, group = lab)
  legacy@clean.data <- legacy@raw.data

  migrated <- suppressMessages(.stat_adopt_legacy(legacy))
  expect_false("group" %in% colnames(migrated@raw.data))
  expect_false("group" %in% colnames(migrated@clean.data))
  expect_true("group" %in% colnames(migrated@info.data))
  expect_silent(validObject(migrated))
})

test_that(".stat_features() actively rejects a frame that carries the label", {
  # The assertion never fires on a correctly-built object, because the contract
  # keeps the label out of the matrix in the first place. Exercise it directly
  # so a future regression that turns the assertion into a silent setdiff()
  # cannot pass unnoticed.
  feat <- make_feat()
  info <- make_info(feat, rep(c("HC", "Case"), each = 15))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  obj@clean.data <- cbind(obj@raw.data, group = info$group)   # force the violation

  expect_error(.stat_features(obj, "clean"), "Internal consistency error")
})

test_that(".stat_features() on a well-formed object returns the frame untouched", {
  feat <- make_feat()
  info <- make_info(feat, rep(c("HC", "Case"), each = 15))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  expect_identical(.stat_features(obj, "clean"), obj@clean.data)
})

test_that("a missing label is a hard error, never a silent NULL", {
  # .stat_labels() is the single read path for labels; if it ever degrades to
  # NULL the caller would silently produce an ungrouped analysis -- the exact
  # failure mode the contract exists to remove.
  feat <- make_feat()
  info <- make_info(feat, rep(c("HC", "Case"), each = 15))

  # group_col = NULL: an object that never had a label must not degrade to NULL.
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = NULL, na.action = "allow")
  expect_error(.stat_labels(obj), "no group_col")

  # A cleared slot is the same failure mode.
  obj2 <- CreateStatObject(raw.data = feat, info.data = info,
                           group_col = "group", na.action = "allow")
  obj2@group_col <- NULL
  expect_error(.stat_labels(obj2), "no group_col")
})

test_that("a label column absent from info.data is a hard error", {
  feat <- make_feat()
  obj <- CreateStatObject(raw.data = feat,
                          info.data = data.frame(other = 1:nrow(feat),
                                                 row.names = rownames(feat)),
                          group_col = NULL, na.action = "allow")
  obj@group_col <- "group"   # name that does not exist in info.data
  expect_error(.stat_labels(obj), "not in info.data")
})
