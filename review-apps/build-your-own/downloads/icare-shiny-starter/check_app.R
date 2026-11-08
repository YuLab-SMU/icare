# Rscript check_app.R classification|subtyping|survival [subtyping_method]
# Checks the supplied templates with your own fitted objects and artificial example CSV.
args <- commandArgs(trailingOnly = TRUE)
kind <- if (length(args)) args[[1]] else stop("Supply an application type.")
stopifnot(kind %in% c("classification", "subtyping", "survival"))
app_dir <- normalizePath(file.path("apps", kind), mustWork = TRUE)
e <- new.env(parent = globalenv())
local({
  old <- setwd(app_dir); on.exit(setwd(old))
  sys.source("app.R", envir = e)
})
stopifnot(inherits(e$app, "shiny.appobj"))
sample_file <- file.path(app_dir, "www/example.csv")
example <- read.csv(sample_file, row.names = 1, check.names = FALSE)
stopifnot(nrow(example) > 0)
if (kind == "classification") {
  for (field in names(example)) {
    ref <- e$deployment$ref_data[[field]]
    if (is.factor(ref)) example[[field]] <- factor(example[[field]], levels = levels(ref))
  }
  method <- e$deployment$model_list[[1]]
  p <- e$deployment$predict_fn(example, selected_model = method)
  stopifnot(nrow(p) == nrow(example), all(is.finite(as.matrix(p))),
            all(abs(rowSums(p) - 1) < 1e-6))
  shiny::testServer(e$app$serverFuncSource(), {
    vals <- lapply(as.list(example[1, , drop = FALSE]), function(x) if(is.factor(x)) as.character(x) else x)
    names(vals) <- paste0("in_", names(vals))
    do.call(session$setInputs, c(vals, list(model_choice = method, threshold = .5, go = 1)))
    stopifnot(length(report_data()$probabilities) >= 2, nzchar(output$res_ui$html), nzchar(output$plot))
    stopifnot(length(jsonlite::fromJSON(output$download_json)$probabilities) >= 2)
  })
} else if (kind == "subtyping") {
  method <- if (length(args) >= 2) args[[2]] else "nmf"
  stopifnot(method %in% c("nmf", "lpa", "kmeans"))
  shiny::testServer(e$app$serverFuncSource(), {
    vals <- setNames(as.list(example[1, , drop = FALSE]), paste0("in_", names(example)))
    do.call(session$setInputs, c(vals, list(input_mode = "Manual", sid = "DEMO_01", meth = method, go = 1)))
    stopifnot(nrow(results()@info.data) == 1)
    session$setInputs(input_mode = "Batch (CSV)", meth = method,
      up = data.frame(datapath = sample_file), go = 2)
    stopifnot(nrow(results()@info.data) == nrow(example),
      all(!is.na(results()@info.data[[paste0("cluster_", method)]])),
      nrow(read.csv(output$dl_res)) == nrow(example))
  })
} else {
  shiny::testServer(e$app$serverFuncSource(), {
    vals <- setNames(as.list(example[1, , drop = FALSE]), paste0("f_", names(example)))
    do.call(session$setInputs, c(vals, list(input_mode = "Manual", sid = "DEMO_01",
      cutoff_m = "median", show_sc = TRUE, go = 1)))
    stopifnot(nrow(preds()) == 1, is.finite(preds()$risk_score),
      nzchar(output$risk_ui$html), nzchar(output$score_ui$html), nzchar(output$res_tbl))
    for (method in c("median", "tertile", "custom")) {
      session$setInputs(input_mode = "Batch Upload (CSV)", cutoff_m = method,
        custom_thresholds = "10,20", show_sc = TRUE,
        up_file = data.frame(datapath = sample_file), go = 2)
      stopifnot(nrow(preds()) == nrow(example), all(is.finite(preds()$risk_score)),
        all(!is.na(preds()$risk_group)), nzchar(output$res_tbl),
        nrow(read.csv(output$dl_res)) == nrow(example))
    }
  })
}
message("PASS: ", kind, " prediction, Shiny server interaction and download checks.")
message("Also test all methods you intend to expose, then repeat in the hosted browser.")
