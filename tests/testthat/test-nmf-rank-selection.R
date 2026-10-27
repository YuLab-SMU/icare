# =============================================================================
# NMF rank selection must not read evidence into a flat curve
# =============================================================================
# Regression from the same incident as test-split-stratify.R: moving ONE sample
# between train and test moved the cophenetic curve enough to flip the chosen
# rank from K=2 to K=5 (0.970 -> 0.914 at K=2, while K=4/K=5 sat at 0.94).
# A single-sample change must not redefine the subtypes.
#
# These tests pin the selection contract: report the margin, and when the sweep
# cannot separate the candidates, take the smallest tied rank rather than
# whatever `which.max` happens to return first.

.ms <- function(ranks, values) {
  data.frame(rank = ranks, cophenetic = values, stringsAsFactors = FALSE)
}

test_that("a clear winner is taken as-is and its margin is reported", {
  r <- icare:::.nmf_select_rank(.ms(2:6, c(0.9704, 0.9251, 0.8745, 0.8359, 0.7387)))
  expect_equal(r$best_rank, 2)
  expect_false(r$ambiguous)
  # 0.9704 - 0.9251 = 0.0453, comfortably above the default 0.02 tolerance.
  expect_equal(r$margin, 0.9704 - 0.9251, tolerance = 1e-6)
  expect_equal(r$tied, 2)
})

test_that("a near-flat curve resolves to the SMALLEST tied rank", {
  # The real offending curve: argmax says K=5, but K=4 is within 0.005.
  r <- icare:::.nmf_select_rank(.ms(2:6, c(0.9145, 0.9015, 0.9429, 0.9478, 0.7776)))
  expect_true(r$ambiguous)
  expect_equal(r$best_rank, 4)          # not 5, which which.max would pick
  expect_equal(r$tied, c(4, 5))
  expect_equal(r$metric_value, 0.9429) # the value AT the chosen rank
  # In a tie, `margin` is the cost of choosing parsimony over the argmax:
  # negative by construction, and the argmax is still reported.
  expect_equal(r$margin, 0.9429 - 0.9478, tolerance = 1e-6)
  expect_equal(r$best_overall, 5)
  expect_equal(r$best_overall_value, 0.9478)
})

test_that("an exact tie is resolved deterministically, not by sweep order", {
  r <- icare:::.nmf_select_rank(.ms(2:4, c(0.9, 0.9, 0.5)))
  expect_equal(r$best_rank, 2)
  expect_true(r$ambiguous)
  expect_equal(r$tied, c(2, 3))
  expect_equal(r$metric_value, 0.9)
  expect_equal(r$margin, 0)   # nothing is given up: the values are equal
})

test_that("tolerance is honoured: widening it widens the tie set", {
  ms <- .ms(2:6, c(0.9704, 0.9251, 0.8745, 0.8359, 0.7387))
  expect_false(icare:::.nmf_select_rank(ms, tie_tol = 0.02)$ambiguous)
  # 0.9704 - 0.9251 = 0.0453, so a 0.05 tolerance pulls K=3 into the tie.
  wide <- icare:::.nmf_select_rank(ms, tie_tol = 0.05)
  expect_true(wide$ambiguous)
  expect_equal(wide$tied, c(2, 3))
  expect_equal(wide$best_rank, 2)   # smallest still wins
})

test_that("NA and single-rank sweeps do not break selection", {
  # NA must be skipped, not propagated into `which.max` (which returns NA).
  r <- icare:::.nmf_select_rank(.ms(2:3, c(0.9, NA)))
  expect_equal(r$best_rank, 2)
  expect_equal(r$margin, Inf)     # nothing to compare against
  expect_false(r$ambiguous)

  single <- icare:::.nmf_select_rank(.ms(3, 0.9))
  expect_equal(single$best_rank, 3)
  expect_equal(single$margin, Inf)
})

test_that("an unusable metric fails loudly instead of guessing", {
  expect_error(icare:::.nmf_select_rank(.ms(2:4, c(NA, NA, NA))),
               "not usable")
  expect_error(icare:::.nmf_select_rank(data.frame(rank = 2:4, other = 1:3)),
               "must contain 'rank' and 'cophenetic'")
})

test_that("a negative or non-scalar tie_tol is rejected, not silently ignored", {
  # Regression: tie_tol < 0 made `vals >= best_val - tie_tol` false for every
  # row, so `tied` came back EMPTY while best_rank still held the argmax --
  # a self-contradictory result that a caller would read as "no tie".
  ms <- .ms(2:6, c(0.9145, 0.9015, 0.9429, 0.9478, 0.7776))
  expect_error(icare:::.nmf_select_rank(ms, tie_tol = -1), "non-negative")
  expect_error(icare:::.nmf_select_rank(ms, tie_tol = NA), "non-negative")
  expect_error(icare:::.nmf_select_rank(ms, tie_tol = c(0.01, 0.02)),
               "non-negative")

  # tie_tol = 0 is legal and means "exact ties only".
  r <- icare:::.nmf_select_rank(.ms(2:4, c(0.9, 0.9, 0.5)), tie_tol = 0)
  expect_equal(r$best_rank, 2)
  expect_equal(r$tied, c(2, 3))
})

test_that("a duplicated rank row cannot leak a vector into metric_value", {
  # NMF's own summary() never repeats a rank, but a caller-supplied frame
  # might; the scalar contract must hold regardless.
  ms <- data.frame(rank = c(2, 2, 3), cophenetic = c(0.90, 0.90, 0.50),
                   stringsAsFactors = FALSE)
  r <- icare:::.nmf_select_rank(ms)
  expect_length(r$metric_value, 1)
  expect_equal(r$metric_value, 0.9)
  expect_equal(r$tied, 2)          # unique(): a RANK, not two rows
  expect_false(r$ambiguous)
})

test_that("rank order in the sweep does not affect the outcome", {
  # sort() on the tie set must make "smallest rank wins" independent of the
  # row order the sweep happened to be produced in.
  asc  <- icare:::.nmf_select_rank(.ms(2:5, c(0.90, 0.91, 0.905, 0.5)))
  desc <- icare:::.nmf_select_rank(.ms(5:2, c(0.50, 0.905, 0.91, 0.90)))
  expect_equal(asc$best_rank, desc$best_rank)
  expect_equal(asc$tied, desc$tied)
  expect_equal(asc$best_rank, 2)
})
