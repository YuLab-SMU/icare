# =============================================================================
# The label must never reach an unsupervised computation
# =============================================================================
# Regression for the original defect: a numeric 0/1 label stored inside the
# feature matrix was selected by "give me the numeric columns" and fed to
# prcomp(), manufacturing apparent group separation.

make_case <- function(labels, seed = 11, n = 80) {
  set.seed(seed)
  feat <- data.frame(f1 = rnorm(n), f2 = rnorm(n), f3 = rnorm(n),
                     f4 = rnorm(n), f5 = rnorm(n))
  feat$f1 <- feat$f1 + ifelse(as.character(labels) == as.character(labels)[1], 1.3, 0)
  rownames(feat) <- paste0("S", seq_len(n))
  info <- data.frame(group = labels, row.names = rownames(feat))
  obj <- CreateStatObject(raw.data = feat, info.data = info,
                          group_col = "group", na.action = "allow")
  suppressMessages(suppressWarnings({
    obj <- stat_miss_processed(obj, impute_method = "median_mode", save_data = FALSE)
    obj
  }))
}

# PC1 group separation when the label is numeric (the leaking case) vs when it
# is not. The old code produced ~5x larger separation for the numeric case.
pca_separation <- function(obj) {
  b <- ggplot2::ggplot_build(PlotPCA(obj, color_by = "group", save_plot = FALSE))$data[[1]]
  diff(tapply(b$x, b$colour, mean))
}

test_that("PCA input columns never include the label, whatever its type", {
  for (l in list(rep(c(0, 1), each = 40),
                 rep(c("0", "1"), each = 40),
                 factor(rep(c("HC", "Case"), each = 40)))) {
    obj <- make_case(l)
    cols_in_pca <- names(obj@clean.data)[vapply(obj@clean.data, is.numeric, logical(1))]
    expect_false("group" %in% cols_in_pca)
  }
})

test_that("PCA separation does not depend on the label's type", {
  num <- pca_separation(make_case(rep(c(0, 1), each = 40)))
  chr <- pca_separation(make_case(rep(c("0", "1"), each = 40)))
  fac <- pca_separation(make_case(factor(rep(c("0", "1"), each = 40))))

  # Identical feature frames, so the separation must be identical.
  expect_equal(num, chr, tolerance = 1e-8)
  expect_equal(num, fac, tolerance = 1e-8)
})

test_that("the correlation heatmap does not correlate the label with features", {
  obj <- make_case(rep(c(0, 1), each = 40))
  expect_error(PlotCorrelationHeatmap(obj, save_plot = FALSE), NA)
  expect_false("group" %in% names(obj@clean.data)[vapply(obj@clean.data, is.numeric, logical(1))])
})
