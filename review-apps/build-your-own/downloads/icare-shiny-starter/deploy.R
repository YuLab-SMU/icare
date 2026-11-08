# Run from the starter-kit root after registering your own account locally.
# Rscript deploy.R survival YOUR_ACCOUNT my-survival-app
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript deploy.R TYPE ACCOUNT APP_NAME")
kind <- args[[1]]; account <- args[[2]]; app_name <- args[[3]]
stopifnot(kind %in% c("classification", "subtyping", "survival"))
app_dir <- file.path("apps", kind)
files <- c("app.R", "model.rds", "www/example.csv")
if (kind == "classification") files <- c(files, "preprocess.rds")
stopifnot(all(file.exists(file.path(app_dir, files))))
registered <- rsconnect::accounts()
if (!any(registered$name == account & registered$server == "shinyapps.io"))
  stop("Register your own shinyapps.io account in your local R console first.")
deps <- rsconnect::appDependencies(app_dir, appFiles = files, appMode = "shiny")
if (kind == "survival") {
  required <- c("mlr3", "mlr3proba", "mlr3tuning", "mlr3tuningspaces",
    "mlr3learners", "mlr3extralearners", "tidyverse", "paradox", "data.table", "survival", "ranger")
  missing <- setdiff(required, deps$Package)
  if (length(missing)) stop("Deployment manifest misses: ", paste(missing, collapse = ", "))
}
dir.create("reports", showWarnings = FALSE)
write.csv(deps, file.path("reports", paste0(kind, "-dependencies.csv")), row.names = FALSE)
capture.output(sessionInfo(), file = file.path("reports", paste0(kind, "-session.txt")))
rsconnect::deployApp(appDir = app_dir, appFiles = files, appName = app_name,
  account = account, server = "shinyapps.io", appMode = "shiny",
  appVisibility = "public", launch.browser = FALSE, forceUpdate = FALSE)
message("Now open the hosted URL and verify a manual prediction, supported CSV upload and download.")
