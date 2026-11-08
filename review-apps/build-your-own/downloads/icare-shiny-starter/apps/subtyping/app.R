library(shiny)
library(icare)
library(shinydashboard)
library(DT)
if (!requireNamespace("NMF", quietly = TRUE)) stop("The NMF package is required.")
library(nnls)
library(mclust)
options(shiny.maxRequestSize = 10 * 1024^2)

model <- readRDS("model.rds")
stopifnot(inherits(model, "Subtyping"))
manager <- icare::New_Sub_Manager(model)
features <- names(model@clean.data)
var_dict <- data.frame(Feature = features,
                      Description = gsub("_", " ", features),
                      Units = rep("Replace with the actual unit", length(features)))
app <- icare::launch_sub_deploy_app(
  manager, title = "ICARE | Subtyping", var_dict = var_dict,
  project_info = list(
    abstract = paste("Describe your fitted subtyping methods and study here.",
                     "Use manual entry or upload a CSV whose first column contains unique sample IDs.",
                     "Example CSV profiles are artificial and are provided for interface testing."),
    citation = "Replace with your study citation."
  )
)
app
