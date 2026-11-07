library(shiny)
library(icare)
library(shinydashboard)
library(DT)
library(mlr3)
library(mlr3proba)
library(mlr3extralearners)
# icare checks these optional packages dynamically when prediction is requested.
# Declare them explicitly so rsconnect includes them in the cloud bundle.
library(mlr3tuning)
library(mlr3tuningspaces)
library(paradox)
library(mlr3learners)
library(tidyverse)
library(data.table)
library(ranger)
library(survival)
options(shiny.maxRequestSize = 10 * 1024^2)

model <- readRDS("model.rds")
stopifnot(inherits(model, "PrognosiX"))
manager <- icare::New_Prog_Manager(model)
# Exercise the actual prediction path in the deployed runtime before serving UI.
local({
  example <- read.csv("www/example.csv", row.names = 1, check.names = FALSE)
  result <- manager$prog_predict(example, cutoff_method = "median", return_scores = TRUE)
  stopifnot(nrow(result) == nrow(example), all(is.finite(result$risk_score)),
            all(!is.na(result$risk_group)))
  message("ICARE survival startup prediction PASS: ", nrow(result),
          " artificial profiles; scores=", paste(result$risk_score, collapse = ","))
})
features <- model@best.model$features
icare::use_app_theme_grey()
icare::set_prog_app_text(
  title = "ICARE | Survival Risk",
  citation_text = "Replace with your study citation."
)
app <- icare::launch_prog_deploy_app(
  manager, title = "ICARE | Survival Risk",
  var_dict = data.frame(Feature = features, Description = gsub("_", " ", features),
                        Units = rep("Replace with the actual unit", length(features))),
  project_info = list(
    abstract = paste("Describe your fitted survival model and study here.",
                     "Enter the required predictors or upload a batch CSV to obtain risk scores and risk groups.",
                     "Example CSV profiles are artificial and are provided for interface testing."),
    citation = "Replace with your study citation."
  )
)
app
