#' @title Universal Prediction Function for Clinlabomics
#' @description A standardized prediction interface that handles feature alignment, 
#' preprocessing, and model dispatching for both single models and ensembles.
#'
#' @param object An S4 model object containing trained models and results.
#' @param newdata A data frame containing the features for prediction.
#' @param preproc A caret preProcess object or NULL.
#' @param threshold Numerical value (0-1) for classification cutoff. Default is 0.5.
#' @param class_labels Character vector for outcome labels. Default is c("Neg", "Pos").
#' @param selected_model Character string identifying which model to use (e.g., "Ensemble Stacking", "rf").
#' 
#' @return A data frame containing probability estimates for each class.
#' @import caret
#' @export
#' 
#' @examples
#' \dontrun{
#' # Example usage:
#' preds <- .predict_universal_v4(
#'   object = my_model_obj, 
#'   newdata = test_df, 
#'   selected_model = "Ensemble Stacking"
#' )
#' }
.predict_universal_v4 <- function(object, newdata, preproc = NULL, threshold = 0.5, 
                                  class_labels = NULL, selected_model = "auto") {
  
  newdata <- as.data.frame(newdata)
  if (!is.null(preproc)) newdata <- predict(preproc, newdata)
  
  res <- object@best.model.result
  all_models <- object@train.models 
  
  # Logic to select the specific model from the container
  if (selected_model == "Ensemble Stacking" && !is.null(res$ensemble)) {
    env <- environment(res$ensemble$predict_fn)
    meta_model <- env$meta_model
    model_list <- env$model_list
    reliable_models <- env$reliable_models
    
    # Extract base model probabilities for the positive class
    base_preds <- as.data.frame(lapply(reliable_models, function(nm) {
      as.numeric(predict(model_list[[nm]], newdata, type = "prob")[, 2])
    }))
    colnames(base_preds) <- reliable_models
    
    # Meta-learner inference
    p_val <- if (inherits(meta_model, "glm")) {
      predict(meta_model, base_preds, type = "response")
    } else {
      predict(meta_model, base_preds, type = "prob")[, 2]
    }
    prob_matrix <- as.data.frame(matrix(c(1 - as.numeric(p_val), as.numeric(p_val)), ncol = 2))
    
  } else if (selected_model %in% names(all_models)) {
    # Single algorithm prediction
    prob_matrix <- as.data.frame(predict(all_models[[selected_model]], newdata, type = "prob"))
    
  } else {
    # Default to the best single model identified during training
    prob_matrix <- as.data.frame(predict(res$model, newdata, type = "prob"))
  }
  
  colnames(prob_matrix) <- if (!is.null(class_labels)) class_labels else c("Neg", "Pos")
  return(prob_matrix)
}

#' @title Model Deployment Constructor
#' @description Wraps the model object and preprocessing logic into a deployment-ready list.
#'
#' @param object The trained model object.
#' @param preproc Preprocessing object (e.g., from caret::preProcess).
#' @param class_labels Labels for the classification outcome.
#' @param model_description Short text describing the model's purpose.
#'
#' @return An object of class 'ModelDeployment'.
#' @export
ModelDeployment <- function(object, preproc = NULL, 
                            class_labels = c("Normal", "High Risk"),
                            model_description = "Platelet-based diagnostic ensemble.") {
  
  available_models <- names(object@train.models)
  if (!is.null(object@best.model.result$ensemble)) {
    available_models <- c("Ensemble Stacking", available_models)
  }
  
  deployment <- list(
    object = object,
    preproc = preproc,
    ref_data = as.data.frame(object@split.data$training),
    group_col = object@group_col,
    class_labels = class_labels,
    model_desc = model_description,
    model_list = available_models,
    predict_fn = function(newdata, threshold, model_choice) {
      .predict_universal_v4(object, newdata, preproc, threshold, class_labels, model_choice)
    }
  )
  class(deployment) <- "ModelDeployment"
  return(deployment)
}

#' @title Clinlabomics Terminal: Light-Tech Edition
#' @description Futuristic UI with high-contrast elements on a light grey base for clinical model deployment.
#' 
#' @param deployment A ModelDeployment object.
#' @param title Character string for the application window title.
#'
#' @import shiny bslib
#' @export
deploy_clinlab_app <- function(deployment, title = "Clinical Omics Terminal") {
  
  train_data <- deployment$ref_data
  feat_cols <- setdiff(colnames(train_data), c(".outcome", "group", "Group", deployment$group_col))
  
  # Base Tech Theme Configuration
  tech_theme <- bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#00dfc0", # Electric Cyan
    bg = "#f7f7f7",      # Light Grey Background
    fg = "#1a1d23",      # Dark Text
    base_font = font_google("Inter")
  )
  
  ui <- page_navbar(
    title = title, theme = tech_theme, bg = "#1a1d23", # Dark Tech Header
    
    # Custom CSS for Cyber-Clinical Aesthetics
    header = tags$style(HTML("
      body { background-color: #f7f7f7 !important; }
      .sidebar { background-color: #1a1d23 !important; color: white !important; }
      .sidebar h6 { color: #00dfc0 !important; font-weight: 800; letter-spacing: 1px; }
      .card { border-radius: 12px; border: 1px solid #e0e0e0; box-shadow: 0 4px 15px rgba(0,0,0,0.05); }
      .card-header { background-color: #ffffff !important; border-bottom: 1px solid #eee; font-weight: 700; color: #1a1d23; }
      .btn-primary { background-color: #00dfc0 !important; border: none; color: #1a1d23 !important; font-weight: 800; }
      .irs-bar { background: #00dfc0 !important; border-top: 1px solid #00dfc0 !important; border-bottom: 1px solid #00dfc0 !important; }
      .irs-from, .irs-to, .irs-single { background: #00dfc0 !important; color: #1a1d23 !important; }
    ")),
    
    nav_panel("Diagnostic Terminal",
              layout_sidebar(
                sidebar = sidebar(
                  title = "CORE CONFIG", width = 350,
                  h6("ALGORITHM SELECT"),
                  selectInput("model_choice", NULL, choices = deployment$model_list),
                  
                  h6("SENSITIVITY THRESHOLD"),
                  sliderInput("threshold", NULL, min = 0, max = 1, value = 0.5, step = 0.01),
                  
                  hr(style = "border-top: 1px solid #444;"),
                  h6("BIOMARKER INPUTS"),
                  lapply(feat_cols, function(f) {
                    val <- train_data[[f]]
                    if (is.numeric(val)) {
                      numericInput(paste0("in_", f), f, value = signif(median(val, na.rm=T), 4))
                    } else {
                      selectInput(paste0("in_", f), f, choices = levels(as.factor(val)))
                    }
                  }),
                  actionButton("go", "INITIATE ANALYSIS", class = "btn-primary w-100"),
                  downloadButton("download_json", "GENERATE DATA PACK", class = "btn-outline-light w-100 mt-2")
                ),
                
                layout_column_wrap(
                  width = 1,
                  card(
                    card_header("ANALYTICAL ENGINE STATUS"),
                    layout_column_wrap(
                      width = 1/2,
                      uiOutput("res_ui"),
                      plotly::plotlyOutput("plot", height = "300px")
                    )
                  ),
                  card(
                    card_header("METADATA & SPECIFICATIONS"),
                    p(deployment$model_desc, style = "font-style: italic;"),
                    markdown(paste0(
                      "**Target:** ", deployment$group_col, "  \n",
                      "**System:** Clinlabomics-X v2.1  \n",
                      "**Status:** Synchronized with Training Set"
                    ))
                  )
                )
              )
    )
  )
  
  server <- function(input, output) {
    # Reactive computation of results
    report_data <- eventReactive(input$go, {
      df_input <- as.data.frame(lapply(feat_cols, function(f) {
        val <- input[[paste0("in_", f)]]
        if(is.numeric(train_data[[f]])) as.numeric(val) else val
      }))
      colnames(df_input) <- feat_cols
      probs <- deployment$predict_fn(df_input, input$threshold, input$model_choice)
      list(inputs = as.list(df_input[1,]), probabilities = as.list(probs[1,]), 
           threshold = input$threshold, model_used = input$model_choice, time = Sys.time())
    })
    
    output$res_ui <- renderUI({
      req(res <- report_data())
      pos_prob <- res$probabilities[[2]]
      is_high <- pos_prob >= res$threshold
      label <- if(is_high) names(res$probabilities)[2] else names(res$probabilities)[1]
      
      # Tech colors: Neon Red for high risk, Electric Cyan for low risk
      text_color <- if(is_high) "#ff2d55" else "#00dfc0" 
      
      div(style = "text-align:center; padding: 40px;",
          h4("PREDICTION RESULT", style = "color: #888; font-size: 0.9rem; letter-spacing: 2px;"),
          h1(label, style = paste0("color: ", text_color, "; font-weight: 900; font-size: 3.5rem; text-shadow: 0 0 10px ", text_color, "44;")),
          h5(sprintf("CONFIDENCE SCORE: %.2f%%", pos_prob * 100), style = "color: #555; margin-top: 10px;")
      )
    })
    
    output$plot <- plotly::renderPlotly({
      req(res <- report_data())
      p_vals <- unlist(res$probabilities)
      plotly::plot_ly(x = names(p_vals), y = as.numeric(p_vals), type = "bar", 
                      marker = list(color = c('#e0e0e0', '#00dfc0'))) %>%
        plotly::layout(
          yaxis = list(title = "Probability", range = c(0, 1), gridcolor = "#f0f0f0"),
          xaxis = list(title = ""),
          paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)'
        )
    })
    
    output$download_json <- downloadHandler(
      filename = function() { paste0("Clinlab_Export_", format(Sys.time(), "%Y%m%d"), ".json") },
      content = function(file) { writeLines(jsonlite::toJSON(report_data(), auto_unbox = TRUE, pretty = TRUE), file) }
    )
  }
  
  shinyApp(ui, server)
}
