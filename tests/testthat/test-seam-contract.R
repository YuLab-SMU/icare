# =============================================================================
# Cross-module seam contract
# =============================================================================
# The label must survive every class conversion with its VALUES intact, aligned
# by row NAME (not position), and feature sets derived from one frame must stay
# intersectable with the original colnames.

make_stat <- function(seed = 2, n = 40) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n), f4 = rnorm(n))
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(diagnosis = factor(rep(c("HC", "Case"), each = n / 2)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "diagnosis", na.action = "allow")
  suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                       save_data = FALSE))
}

test_that("PrognosiX -> Subtyping does not crash on the missing scale.data slot", {
  # BLOCKER regression: the PrognosiX branch read @scale.data unconditionally.
  obj <- make_stat()
  info <- obj@info.data
  info$time <- rexp(nrow(info), 0.05); info$status <- rbinom(nrow(info), 1, 0.6)
  obj@info.data <- info
  px <- CreatePrognosiXObject(object = obj, group_col = "diagnosis",
                              time_col = "time", status_col = "status")
  sub <- expect_error(ConvertObject(px, to = "Subtyping"), NA)
  expect_true("diagnosis" %in% colnames(sub@info.data))
})

test_that("Subtyping/PrognosiX -> Train_Model keeps the resolved label name", {
  # The branches used to hardcode "group", discarding the name .to_stat()
  # resolved; a label named anything else then failed downstream.
  obj <- make_stat()
  sub <- ConvertObject(obj, to = "Subtyping")
  expect_equal(sub@group_col, "diagnosis")

  tm <- ConvertObject(sub, to = "Train_Model")
  expect_equal(tm@group_col, "diagnosis")
  expect_true("diagnosis" %in% colnames(tm@clean.df))

  info <- obj@info.data
  info$time <- rexp(nrow(info), 0.05); info$status <- rbinom(nrow(info), 1, 0.6)
  obj@info.data <- info
  px <- CreatePrognosiXObject(object = obj, group_col = "diagnosis",
                              time_col = "time", status_col = "status")
  tm2 <- ConvertObject(px, to = "Train_Model")
  expect_equal(tm2@group_col, "diagnosis")
})

test_that("feature names survive feature elimination unmangled", {
  # make.names() in .extract_xy() renamed "1q21" to "X1q21"; the recorded
  # feature sets then intersected against the ORIGINAL colnames to nothing,
  # and ApplyFeatureSelection() silently left an outcome-only model.
  set.seed(4)
  n <- 40
  feat <- data.frame(`1q21` = rnorm(n), check.names = FALSE)
  feat[["2p22"]] <- rnorm(n)
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(diagnosis = factor(rep(c("HC", "Case"), each = n / 2)),
                     row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "diagnosis", na.action = "allow")
  obj <- suppressMessages(stat_miss_processed(obj, impute_method = "median_mode",
                                              save_data = FALSE))
  m <- CreateModelObject(object = obj, group_col = "diagnosis")

  x <- icare:::.extract_xy(m)$x
  expect_setequal(colnames(x), c("1q21", "2p22"))   # no make.names mangling

  af <- icare:::ApplyFeatureSelection(m, c("1q21", "2p22"))
  expect_true("1q21" %in% colnames(af@clean.df))
  expect_true("2p22" %in% colnames(af@clean.df))
})

test_that("SubsetObject(features=) on Train_Model keeps the outcome column", {
  obj <- make_stat()
  m <- CreateModelObject(object = obj, group_col = "diagnosis")
  m2 <- SubsetObject(m, features = c("f1", "f2"))
  expect_true("diagnosis" %in% colnames(m2@clean.df))
  expect_silent(validObject(m2))
})

test_that("conversion refuses to fabricate NA labels when info.data under-covers", {
  # .object_split_label() used to fill uncovered rows with NA labels.
  m <- CreateModelObject(object = make_stat(), group_col = "diagnosis")
  m@info.data <- m@info.data[1:10, , drop = FALSE]   # simulate under-coverage

  expect_error(icare:::ConvertObject(m, to = "Stat"), "does not cover")
})

test_that("SplitSubtypingObject stratifies by row name, not position", {
  obj <- make_stat(seed = 3, n = 50)
  info <- obj@info.data
  info[["diagnosis"]] <- factor(c(rep("HC", 40), rep("Case", 10)))  # 40/10
  obj@info.data <- info
  s <- ConvertObject(obj, to = "Subtyping")

  split <- SplitSubtypingObject(s, p = 0.5, stratify_by = "diagnosis", seed = 1)
  expect_equal(as.integer(table(split$train@info.data$diagnosis)), c(5L, 20L))

  # Scrambling info.data row order must not defeat stratification...
  s2 <- s
  s2@info.data <- s2@info.data[sample(rownames(s2@info.data)), , drop = FALSE]
  split2 <- SplitSubtypingObject(s2, p = 0.5, stratify_by = "diagnosis", seed = 1)
  expect_equal(as.integer(table(split2$train@info.data$diagnosis)), c(5L, 20L))
  # ...and the label VALUES must still match the original object, per row name.
  common <- intersect(rownames(split2$train@info.data), rownames(info))
  expect_identical(as.character(split2$train@info.data[common, "diagnosis"]),
                   as.character(info[common, "diagnosis"]))
})

test_that("an unknown stratify_by warns instead of silently not stratifying", {
  s <- ConvertObject(make_stat(), to = "Subtyping")
  expect_warning(SplitSubtypingObject(s, stratify_by = "nope"),
                 "not a column of info.data")
})

test_that("Subtyping -> PrognosiX keeps the label name so survival analysis excludes it", {
  # prognosis_02_framework.R reads @group_col to drop the label from the
  # candidate-feature scan; a NULL there makes the label a Cox feature.
  obj <- make_stat()
  info <- obj@info.data
  info$time <- rexp(nrow(info), 0.05); info$status <- rbinom(nrow(info), 1, 0.6)
  obj@info.data <- info
  sub <- ConvertObject(obj, to = "Subtyping")
  px  <- ConvertObject(sub, to = "PrognosiX", time_col = "time", status_col = "status")
  expect_equal(px@group_col, "diagnosis")
})
