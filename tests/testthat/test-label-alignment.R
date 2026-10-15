# =============================================================================
# Label alignment must be BY ROW NAME, pinned at the value level
# =============================================================================
# Mutation testing showed two protections were pinned only by side effects:
# reverting DownsampleObject to positional grouping passed green, and a
# positional .stat_labels() was caught only by an unrelated error-message
# regex. These tests assert the essence, not the side effect.

make_stat <- function(seed = 7, n = 40) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n))
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(diagnosis = factor(rep(c("HC", "Case"), each = n / 2)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "diagnosis", na.action = "allow")
  suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                       save_data = FALSE))
}

test_that(".stat_labels() aligns BY ROW NAME even when info.data is reordered", {
  obj <- make_stat()
  scrambled <- obj
  scrambled@info.data <- obj@info.data[sample(rownames(obj@info.data)), , drop = FALSE]

  rows <- rownames(scrambled@clean.data)
  got <- icare:::.stat_labels(scrambled, rows = rows)[[1]]
  want <- scrambled@info.data[rows, "diagnosis"]
  expect_identical(as.character(got), as.character(want))
  # And it is NOT the positional reading of the scrambled frame.
  expect_false(identical(as.character(got),
                         as.character(scrambled@info.data[["diagnosis"]])[seq_along(rows)]))
})

test_that("DownsampleObject keeps per-group identity, not just counts", {
  obj <- make_stat(seed = 21)
  # Reverse info.data: a positional implementation pairs the wrong labels.
  obj@info.data <- obj@info.data[rev(rownames(obj@info.data)), , drop = FALSE]
  before_named <- obj@info.data[rownames(obj@clean.data), "diagnosis", drop = FALSE][[1]]
  names(before_named) <- rownames(obj@clean.data)

  d <- suppressMessages(DownsampleObject(obj, n = 5, per_group = TRUE,
                                         group_col = "diagnosis", seed = 1))
  expect_equal(nrow(d@clean.data), 10)

  # Essence: every kept sample's label equals the label that sample had BEFORE
  # downsampling -- not merely that the counts balance. (Extract by row name:
  # a bare df[[col]] drops the names, which is exactly the bug class in play.)
  kept <- rownames(d@clean.data)
  expect_identical(as.character(icare:::.stat_labels(d)[[1]]),
                   as.character(before_named[kept]))
  expect_equal(as.integer(table(icare:::.stat_labels(d)[[1]])), c(5L, 5L))
})

test_that("the Stat -> Train_Model round trip preserves label VALUES", {
  obj <- make_stat()
  m <- CreateModelObject(object = obj, group_col = "diagnosis")

  rows <- rownames(m@clean.df)
  expect_identical(as.character(m@clean.df[["diagnosis"]]),
                   as.character(obj@info.data[rows, "diagnosis"]))

  # And the way back, including with reordered info.data in between.
  back <- ConvertObject(m, to = "Stat")
  expect_identical(as.character(back@info.data[rows, "diagnosis"]),
                   as.character(obj@info.data[rows, "diagnosis"]))
})

test_that("PrognosiX validity rejects a label inside clean.data", {
  obj <- make_stat()
  info <- obj@info.data
  info$time <- rexp(nrow(info), 0.05); info$status <- rbinom(nrow(info), 1, 0.6)
  obj@info.data <- info
  px <- CreatePrognosiXObject(object = obj, group_col = "diagnosis",
                              time_col = "time", status_col = "status")

  leak <- px
  leak@clean.data <- cbind(leak@clean.data, diagnosis = obj@info.data[rownames(leak@clean.data), "diagnosis"])
  expect_error(validObject(leak), "must not appear in clean.data")

  # Same for a group_col that does not exist in info.data.
  dangling <- px
  dangling@group_col <- "nope"
  expect_error(validObject(dangling), "not found in info.data")
})
