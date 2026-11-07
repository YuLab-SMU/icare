# Run from the extracted starter-kit root, after adding your own model.rds.
# Rscript make_example.R classification|subtyping|survival
args <- commandArgs(trailingOnly = TRUE)
kind <- if (length(args)) args[[1]] else stop("Supply an application type.")
stopifnot(kind %in% c("classification", "subtyping", "survival"))
suppressPackageStartupMessages(library(icare))
app_dir <- file.path("apps", kind)
m <- readRDS(file.path(app_dir, "model.rds"))
if (kind == "classification") {
  fields <- setdiff(names(m@train.models[[1]]$trainingData), c(".outcome", m@group_col))
  ref <- m@split.data$training[, fields, drop = FALSE]
} else if (kind == "subtyping") {
  ref <- m@clean.data
} else {
  ref <- m@clean.data[, m@best.model$features, drop = FALSE]
}
# Artificial marginal-quantile profiles: not observed patients or validation data.
values <- lapply(ref, function(x) {
  if (is.numeric(x)) {
    if (!any(is.finite(x))) stop("A numeric predictor has no finite reference values.")
    as.numeric(quantile(x[is.finite(x)], c(.25, .5, .75)))
  } else {
    choices <- if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
    if (!length(choices)) stop("A categorical predictor has no reference levels.")
    rep(choices[[1]], 3)
  }
})
example <- data.frame(SampleID = sprintf("DEMO_%02d", 1:3), values, check.names = FALSE)
dir.create(file.path(app_dir, "www"), recursive = TRUE, showWarnings = FALSE)
write.csv(example, file.path(app_dir, "www/example.csv"), row.names = FALSE)
message("Review the units, category values and plausibility of www/example.csv before use.")
