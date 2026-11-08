library(shiny)
library(icare)
# Explicit imports also let rsconnect identify packages used by saved models.
library(caret)
library(C50)
library(glmnet)
library(kernlab)
library(klaR)
library(randomForest)
library(rpart)
library(MASS)
library(ipred)
library(bslib)
library(plotly)
library(jsonlite)
options(shiny.maxRequestSize = 10 * 1024^2)

model <- readRDS("model.rds")
preproc <- readRDS("preprocess.rds")
features <- setdiff(names(model@train.models[[1]]$trainingData), c(".outcome", model@group_col))
stopifnot(inherits(model, "Train_Model"), length(features) > 0)
deployment <- icare::ModelDeployment(
  model, preproc = preproc,
  class_labels = levels(model@split.data$training[[model@group_col]]),
  model_description = paste(
    "Describe your classification model and intended research use here.",
    "Enter the required predictors, select a model and click INITIATE ANALYSIS.",
    "The probability is a model output for research demonstration."
  )
)
# Include categorical predictors as well as numeric predictors.
# This template uses the caret training schema from the first fitted model.
deployment$preprocess_columns <- features
app <- icare::deploy_clinlab_app(deployment, title = "ICARE | Classification")
app
