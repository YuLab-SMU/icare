# =============================================================================
# clinical_analysis.R
# Comprehensive Clinical Analysis Module
# Subgroups · Confounders · Thresholds · Decision Curves
# =============================================================================
# Supports: Train_Model, caret train, ensemble, fine‑tuned models

# 0. Package check -----------------------------------------------------------
#' @keywords internal
.check_clinical_pkgs <- function() {
  required <- c("ggplot2", "dplyr", "tidyr", "wesanderson", "ggprism",
                "pROC", "caret", "reshape2", "scales", "ggrepel",
                "nricens", "corrplot")
  missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Please install them.")
  }
  invisible(TRUE)
}

# ── Internal prediction helper ─────────────────────────────────────────────
#' @keywords internal
.predict_probs <- function(model_obj, newdata, model_name = NULL,
                           positive_class = NULL) {
  if (inherits(model_obj, "Train_Model")) {
    if (is.null(model_name)) {
      best_model <- model_obj@best.model.result$model
      if (is.null(best_model)) best_model <- model_obj@train.models[[1]]
    } else if (model_name == "ensemble") {
      ens <- model_obj@best.model.result$ensemble
      if (is.null(ens)) stop("No ensemble found.")
      return(ens$predict_fn(newdata))
    } else {
      best_model <- model_obj@train.models[[model_name]]
    }
    prob_mat <- predict(best_model, newdata, type = "prob")
    
    # ---- 自动处理正类列 ----
    if (is.matrix(prob_mat) || is.data.frame(prob_mat)) {
      if (ncol(prob_mat) == 1) {
        # 单列向量直接返回
        return(prob_mat[, 1])
      }
      # 如果未指定正类，默认取第二列
      if (is.null(positive_class)) {
        positive_class <- colnames(prob_mat)[2]
        warning("positive_class not specified; using second column: ", positive_class)
        return(prob_mat[, positive_class])
      }
      # 如果指定了正类但列名不存在，降级使用第二列
      if (!positive_class %in% colnames(prob_mat)) {
        warning("Column '", positive_class, "' not found. Using second column.")
        return(prob_mat[, 2])
      }
      return(prob_mat[, positive_class])
    } else {
      # 已经是数值向量
      return(prob_mat)
    }
    
  } else if (inherits(model_obj, "train")) {
    prob_mat <- predict(model_obj, newdata, type = "prob")
    # 同样的逻辑
    if (is.matrix(prob_mat) || is.data.frame(prob_mat)) {
      if (ncol(prob_mat) == 1) return(prob_mat[, 1])
      if (is.null(positive_class)) {
        positive_class <- colnames(prob_mat)[2]
        warning("positive_class not specified; using second column: ", positive_class)
        return(prob_mat[, positive_class])
      }
      if (!positive_class %in% colnames(prob_mat)) {
        warning("Column '", positive_class, "' not found. Using second column.")
        return(prob_mat[, 2])
      }
      return(prob_mat[, positive_class])
    } else {
      return(prob_mat)
    }
  } else {
    stop("model_obj must be a Train_Model or caret train object.")
  }
}
# ── 1. Attach clinical data to model object ────────────────────────────────
#' Attach Clinical Data to Train_Model Object
#'
#' Stores the entire clinical data frame inside \code{@process.info$clinical_data}.
#' Subsequent clinical analysis functions will use this by default.
#'
#' @param object        A \code{Train_Model} object.
#' @param clinical_data A data frame with clinical variables. Rownames must
#'   match the samples in the model's testing set.
#' @return The updated \code{Train_Model} object.
#' @export
AttachClinicalData <- function(object, clinical_data) {
  if (!inherits(object, "Train_Model"))
    stop("object must be a Train_Model.")
  if (!is.data.frame(clinical_data))
    stop("clinical_data must be a data frame.")
  object@process.info$clinical_data <- clinical_data
  cat("Clinical data attached to the model object.\n")
  return(object)
}

# ── 2. Correlation matrix plot ─────────────────────────────────────────────
#' Prediction‑Clinical Correlation Matrix (corrplot Style)
#'
#' Converts categorical variables to numeric, combines them with model
#' predictions, and draws a Spearman correlation matrix using \code{corrplot}.
#' Significant correlations are annotated with stars.
#'
#' @param model_obj,clinical_data,newdata,model_name See \code{\link{ClinicalAnalysis}}.
#' @param save_plot,save_dir,palette_name  Output options.
#' @param ...   Further arguments passed to \code{corrplot::corrplot}.
#' @return Invisibly returns the correlation matrix.
#' @export
PlotClinicalCorrelation <- function(model_obj,
                                    clinical_data = NULL,
                                    newdata       = NULL,
                                    model_name    = NULL,
                                    save_plot     = FALSE,
                                    save_dir      = NULL,
                                    palette_name  = "Royal1",
                                    ...) {
  .check_clinical_pkgs()
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  if (is.null(newdata) && inherits(model_obj, "Train_Model"))
    newdata <- model_obj@filtered.set$testing
  stopifnot(!is.null(newdata), nrow(newdata) == nrow(clinical_data))
  
  probs <- .predict_probs(model_obj, newdata, model_name)
  
  # Convert categorical variables to numeric
  df <- clinical_data
  cat_vars <- names(df)[!sapply(df, is.numeric)]
  for (v in cat_vars) {
    df[[v]] <- as.numeric(as.factor(df[[v]]))
  }
  
  df$Prediction <- probs
  
  cor_res <- corrplot::cor.mtest(df, conf.level = 0.95)
  M <- cor(df, method = "spearman", use = "complete.obs")
  
  cols <- colorRampPalette(
    c(wesanderson::wes_palette(palette_name, 2, type = "discrete")[1],
      "white",
      wesanderson::wes_palette(palette_name, 2, type = "discrete")[2])
  )(200)
  
  draw_corr <- function() {
    corrplot::corrplot(
      M,
      method = "square",
      type = "lower",
      tl.col = "black",
      tl.cex = 0.8,, diag = FALSE,
      p.mat = cor_res$p,
      sig.level = c(0.001, 0.01, 0.05),
      insig = "label_sig",
      pch.cex = 1,
      pch.col = "grey20",
      col = cols,
      title = "Prediction vs Clinical Variables",
      mar = c(0, 0, 2, 0),
      ...
    )
  }
  
  draw_corr()
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    pdf(file.path(save_dir, "clinical_correlation.pdf"), width = 8, height = 8)
    draw_corr()
    dev.off()
    cat("Plot saved to:", file.path(save_dir, "clinical_correlation.pdf"), "\n")
  }
  
  invisible(M)
}

# ── 3. Subgroup forest plot ────────────────────────────────────────────────
#' Subgroup Performance Forest Plot
#'
#' Computes AUC with 95\% CI within each subgroup defined by a categorical
#' clinical variable.
#'
#' @param model_obj,clinical_data,subgroup_var,newdata,model_name As above.
#' @param save_plot,save_dir,palette_name  Output options.
#' @return A ggplot object.
#' @export
PlotSubgroupForest <- function(model_obj,
                               clinical_data = NULL,
                               subgroup_var,
                               newdata      = NULL,
                               model_name   = NULL,
                               save_plot    = FALSE,
                               save_dir     = NULL,
                               palette_name = "Darjeeling1") {
  .check_clinical_pkgs()
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  if (is.null(newdata) && inherits(model_obj, "Train_Model"))
    newdata <- model_obj@filtered.set$testing
  probs <- .predict_probs(model_obj, newdata, model_name)
  
  true <- factor(newdata[[model_obj@group_col]])
  positive <- levels(true)[2]
  
  groups <- unique(clinical_data[[subgroup_var]])
  res_list <- lapply(groups, function(g) {
    idx <- which(clinical_data[[subgroup_var]] == g)
    if (length(idx) < 5) return(NULL)
    roc_obj <- pROC::roc(true[idx], probs[idx], levels = c(levels(true)[1], positive),
                         direction = "auto", quiet = TRUE)
    ci <- pROC::ci.auc(roc_obj)
    data.frame(Subgroup = g, AUC = as.numeric(roc_obj$auc),
               Lower = ci[1], Upper = ci[3], N = length(idx))
  })
  forest_df <- do.call(rbind, res_list)
  forest_df <- forest_df[order(-forest_df$AUC), ]
  forest_df$Subgroup <- factor(forest_df$Subgroup, levels = forest_df$Subgroup)
  
  cols <- wesanderson::wes_palette(palette_name, nrow(forest_df), type = "discrete")
  
  p <- ggplot2::ggplot(forest_df, ggplot2::aes(x = AUC, y = Subgroup, color = Subgroup)) +
    ggplot2::geom_point(size = 3, shape = 18) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper), height = 0.2) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::labs(title = paste("AUC by", subgroup_var),
                  x = "AUC (95% CI)", y = NULL) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("subgroup_", subgroup_var, ".pdf")),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  return(p)
}

# ── 4. Confounder adjustment forest plot ───────────────────────────────────
#' Multivariable Logistic Adjustment Bar Plot (Significance Fill)
#'
#' Fits \code{outcome ~ prediction + clinical_vars} using the binary outcome
#' column from \code{newdata} and plots \code{-log10(p-value)} for each
#' variable. Bar fill represents the significance level (darker = more
#' significant). A dashed reference line marks \code{p = 0.05}.
#'
#' @param model_obj,clinical_data,outcome_var,newdata,model_name As above.
#' @param save_plot,save_dir  Output options.
#' @return A ggplot object.
#' @export
PlotConfounderForest <- function(model_obj,
                                 clinical_data = NULL,
                                 outcome_var,
                                 newdata      = NULL,
                                 model_name   = NULL,
                                 save_plot    = FALSE,
                                 save_dir     = NULL) {
  .check_clinical_pkgs()
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  if (is.null(newdata) && inherits(model_obj, "Train_Model"))
    newdata <- model_obj@filtered.set$testing
  stopifnot(!is.null(newdata), nrow(newdata) == nrow(clinical_data))
  
  probs <- .predict_probs(model_obj, newdata, model_name)
  
  outcome_vec <- newdata[[outcome_var]]
  if (is.null(outcome_vec))
    stop("Outcome column '", outcome_var, "' not found in newdata.")
  
  df <- cbind(clinical_data, prediction = probs)
  df$outcome <- factor(outcome_vec)
  df$outcome_bin <- as.numeric(df$outcome == levels(df$outcome)[2])
  
  covars <- setdiff(names(clinical_data), outcome_var)
  frm <- as.formula(paste("outcome_bin ~ prediction +", paste(covars, collapse = " + ")))
  fit <- glm(frm, data = df, family = binomial)
  
  tem <- summary(fit)$coefficients
  tem <- as.data.frame(tem)
  tem$`-log10P` <- -log10(tem$`Pr(>|z|)`)
  tem <- tem[rownames(tem) != "(Intercept)", ]
  tem$Variable <- rownames(tem)
  
  # Order by -log10P (ascending, for horizontal bar)
  tem$Variable <- factor(tem$Variable, levels = tem$Variable[order(tem$`-log10P`)])
  
  p <- ggplot2::ggplot(tem, ggplot2::aes(x = .data[["-log10P"]], y = .data[["Variable"]],
                                         fill = .data[["-log10P"]])) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradientn(
      colours = c("#f9ddda", "#eda8bd", "#ce78b3", "#9955a8", "#573b88"),
      name = expression(-log[10]("P"))
    ) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                        color = "black", linewidth = 1) +
    ggplot2::annotate("text", x = -log10(0.05) + 0.02, y = 1,
                      label = "p = 0.05", hjust = 0, size = 3.5, color = "black") +
    ggplot2::labs(title = "Confounder Adjustment (Logistic Regression)",
                  x = "-log10(P-value)", y = NULL) +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "confounder_bar.pdf"),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  return(p)
}

# ── 5. Threshold analysis ─────────────────────────────────────────────────
#' Calculate Multi‑Threshold Metrics with Custom Targets
#'
#' Computes a full table of Accuracy, PPV, NPV across all unique thresholds,
#' then finds:
#'   - Youden index (Se + Sp − 1)
#'   - Thresholds that reach a target Sensitivity, Specificity, PPV, or NPV
#'   - Threshold that maximises Accuracy (if requested)
#'
#' @param model_obj,newdata,model_name  As elsewhere.
#' @param target_se      Target Sensitivity (e.g., 0.9). Default NULL.
#' @param target_sp      Target Specificity (e.g., 0.9). Default NULL.
#' @param target_ppv     Target PPV (e.g., 0.95). Default NULL.
#' @param target_npv     Target NPV (e.g., 0.95). Default NULL.
#' @param target_acc     Logical; if TRUE, find threshold maximising Accuracy.
#' @return A list with \code{thresholds}, \code{metrics_df}, \code{probabilities},
#'   \code{true}, \code{positive}, \code{negative}.
#' @export
CalculateThresholds <- function(model_obj,
                                newdata      = NULL,
                                model_name   = NULL,
                                target_se    = NULL,
                                target_sp    = NULL,
                                target_ppv   = NULL,
                                target_npv   = NULL,
                                target_acc   = TRUE) {
  .check_clinical_pkgs()
  if (is.null(newdata) && inherits(model_obj, "Train_Model"))
    newdata <- model_obj@filtered.set$testing
  stopifnot(!is.null(newdata))
  
  probs <- .predict_probs(model_obj, newdata, model_name)
  true  <- factor(newdata[[model_obj@group_col]])
  positive <- levels(true)[2]
  negative <- levels(true)[1]
  
  roc_obj <- pROC::roc(true, probs, levels = c(negative, positive),
                       direction = "auto", quiet = TRUE)
  coords_all <- pROC::coords(roc_obj, "all", ret = c("threshold", "se", "sp"))
  
  youden <- coords_all[which.max(coords_all$se + coords_all$sp - 1), ]
  
  uniq_thr <- sort(unique(round(probs, 4)), decreasing = TRUE)
  met_list <- lapply(uniq_thr, function(t) {
    pred_class <- factor(ifelse(probs > t, positive, negative),
                         levels = c(negative, positive))
    cm <- caret::confusionMatrix(pred_class, true, positive = positive)
    data.frame(Threshold = t,
               Sensitivity = cm$byClass["Sensitivity"],
               Specificity = cm$byClass["Specificity"],
               Accuracy    = cm$overall["Accuracy"],
               PPV         = cm$byClass["Pos Pred Value"],
               NPV         = cm$byClass["Neg Pred Value"],
               F1          = cm$byClass["F1"],
               Precision   = cm$byClass["Precision"],
               Recall      = cm$byClass["Recall"])
  })
  metrics_df <- do.call(rbind, met_list)
  
  .find_thresh <- function(metric_col, target, larger = FALSE) {
    if (larger) {
      idx <- which(metrics_df[[metric_col]] >= target)
      if (length(idx) == 0) return(NULL)
      metrics_df[idx[which.min(metrics_df$Threshold[idx])], ]
    } else {
      dist <- abs(metrics_df[[metric_col]] - target)
      metrics_df[which.min(dist), ]
    }
  }
  
  thresholds <- c(Youden = youden$threshold)
  if (!is.null(target_se)) {
    se_row <- .find_thresh("Sensitivity", target_se, larger = TRUE)
    if (!is.null(se_row)) thresholds <- c(thresholds, Se_Target = se_row$Threshold)
  }
  if (!is.null(target_sp)) {
    sp_row <- .find_thresh("Specificity", target_sp, larger = TRUE)
    if (!is.null(sp_row)) thresholds <- c(thresholds, Sp_Target = sp_row$Threshold)
  }
  if (!is.null(target_ppv)) {
    ppv_row <- .find_thresh("PPV", target_ppv, larger = TRUE)
    if (!is.null(ppv_row)) thresholds <- c(thresholds, PPV_Target = ppv_row$Threshold)
  }
  if (!is.null(target_npv)) {
    npv_row <- .find_thresh("NPV", target_npv, larger = TRUE)
    if (!is.null(npv_row)) thresholds <- c(thresholds, NPV_Target = npv_row$Threshold)
  }
  if (target_acc) {
    acc_max <- metrics_df[which.max(metrics_df$Accuracy), ]
    thresholds <- c(thresholds, MaxAcc = acc_max$Threshold)
  }
  
  list(thresholds    = thresholds,
       metrics_df    = metrics_df,
       probabilities = probs,
       true          = true,
       positive      = positive,
       negative      = negative)
}
#' Calculate Multi‑Threshold Metrics from External Probabilities
#'
#' This function works exactly like \code{CalculateThresholds}, but accepts
#' directly a vector of predicted probabilities and the corresponding true
#' labels, without needing a model object.  It is useful when the prediction
#' scores come from an external source (e.g., literature, another software).
#'
#' @param probs        Numeric vector of predicted probabilities (range 0‑1).
#' @param true         Factor vector of true binary labels.
#' @param positive     Character string specifying the positive class
#'   (e.g., \code{"yes"}).  Must be one of the levels of \code{true}.
#' @param target_se, target_sp, target_ppv, target_npv, target_acc
#'   Same as in \code{CalculateThresholds}.
#' @return A list with the same structure as \code{CalculateThresholds},
#'   suitable for all downstream plotting and analysis functions.
#' @export
CalculateThresholdsFromProbs <- function(probs,
                                         true,
                                         positive,
                                         target_se  = NULL,
                                         target_sp  = NULL,
                                         target_ppv = NULL,
                                         target_npv = NULL,
                                         target_acc = TRUE) {
  .check_clinical_pkgs()
  if (!is.factor(true)) true <- factor(true)
  if (!positive %in% levels(true))
    stop("'positive' must be one of the levels of 'true'.")
  negative <- setdiff(levels(true), positive)[1]
  
  roc_obj <- pROC::roc(true, probs, levels = c(negative, positive),
                       direction = "auto", quiet = TRUE)
  coords_all <- pROC::coords(roc_obj, "all", ret = c("threshold", "se", "sp"))
  
  youden <- coords_all[which.max(coords_all$se + coords_all$sp - 1), ]
  
  # Build metrics table (same as original)
  uniq_thr <- sort(unique(round(probs, 4)), decreasing = TRUE)
  met_list <- lapply(uniq_thr, function(t) {
    pred_class <- factor(ifelse(probs > t, positive, negative),
                         levels = c(negative, positive))
    cm <- caret::confusionMatrix(pred_class, true, positive = positive)
    data.frame(Threshold = t,
               Sensitivity = cm$byClass["Sensitivity"],
               Specificity = cm$byClass["Specificity"],
               Accuracy    = cm$overall["Accuracy"],
               PPV         = cm$byClass["Pos Pred Value"],
               NPV         = cm$byClass["Neg Pred Value"],
               F1          = cm$byClass["F1"],
               Precision   = cm$byClass["Precision"],
               Recall      = cm$byClass["Recall"])
  })
  metrics_df <- do.call(rbind, met_list)
  
  .find_thresh <- function(metric_col, target, larger = FALSE) {
    if (larger) {
      idx <- which(metrics_df[[metric_col]] >= target)
      if (length(idx) == 0) return(NULL)
      metrics_df[idx[which.min(metrics_df$Threshold[idx])], ]
    } else {
      dist <- abs(metrics_df[[metric_col]] - target)
      metrics_df[which.min(dist), ]
    }
  }
  
  thresholds <- c(Youden = youden$threshold)
  if (!is.null(target_se)) {
    se_row <- .find_thresh("Sensitivity", target_se, larger = TRUE)
    if (!is.null(se_row)) thresholds <- c(thresholds, Se_Target = se_row$Threshold)
  }
  if (!is.null(target_sp)) {
    sp_row <- .find_thresh("Specificity", target_sp, larger = TRUE)
    if (!is.null(sp_row)) thresholds <- c(thresholds, Sp_Target = sp_row$Threshold)
  }
  if (!is.null(target_ppv)) {
    ppv_row <- .find_thresh("PPV", target_ppv, larger = TRUE)
    if (!is.null(ppv_row)) thresholds <- c(thresholds, PPV_Target = ppv_row$Threshold)
  }
  if (!is.null(target_npv)) {
    npv_row <- .find_thresh("NPV", target_npv, larger = TRUE)
    if (!is.null(npv_row)) thresholds <- c(thresholds, NPV_Target = npv_row$Threshold)
  }
  if (target_acc) {
    acc_max <- metrics_df[which.max(metrics_df$Accuracy), ]
    thresholds <- c(thresholds, MaxAcc = acc_max$Threshold)
  }
  
  list(thresholds    = thresholds,
       metrics_df    = metrics_df,
       probabilities = probs,
       true          = true,
       positive      = positive,
       negative      = negative)
}
# ── 6. Threshold accuracy / PPV / NPV curve ────────────────────────────────
#' Accuracy/PPV/NPV vs Threshold Plot with Custom Threshold Markers
#'
#' @param thresh_result Output from \code{CalculateThresholds}.
#' @param save_plot,save_dir  Output options.
#' @return A ggplot object.
#' @export
PlotThresholdAccuracy <- function(thresh_result,
                                  save_plot = FALSE,
                                  save_dir  = NULL) {
  df <- thresh_result$metrics_df
  df_long <- tidyr::pivot_longer(df, cols = c("Accuracy", "PPV", "NPV"),
                                 names_to = "Metric", values_to = "Value")
  thr <- thresh_result$thresholds
  
  label_df <- do.call(rbind, lapply(names(thr), function(nm) {
    tval <- thr[nm]
    idx <- which.min(abs(df$Threshold - tval))
    data.frame(Threshold = tval,
               Metric    = c("Accuracy", "PPV", "NPV"),
               Value     = c(df$Accuracy[idx], df$PPV[idx], df$NPV[idx]),
               Label     = nm)
  }))
  
  cols <- c("Accuracy" = "#800000", "PPV" = "#767676", "NPV" = "#cc8214")
  lty  <- c("Accuracy" = "solid",  "PPV" = "dotted",   "NPV" = "dashed")
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Threshold, y = Value,
                                             color = Metric, linetype = Metric)) +
    ggplot2::geom_line(size = 1.1) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_linetype_manual(values = lty) +
    ggplot2::geom_point(data = label_df,
                        ggplot2::aes(x = Threshold, y = Value),
                        size = 3, shape = 8, color = "red") +
    ggrepel::geom_text_repel(data = label_df,
                             ggplot2::aes(label = paste0(Label, "\n(", round(Value, 2), ")")),
                             color = "red", size = 3.5) +
    ggplot2::labs(title = "Accuracy / PPV / NPV vs Threshold",
                  x = "Threshold", y = "Value") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = "bottom")
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "threshold_accuracy.pdf"),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  return(p)
}

# ── 7. Decision zone density plot ─────────────────────────────────────────
#' Decision Density with Threshold Zones
#'
#' @param thresh_result Output from \code{CalculateThresholds}.
#' @param lower_threshold Default: \code{Se_Target} if exists, else Youden.
#' @param upper_threshold Default: \code{Sp_Target} if exists, else Youden.
#' @param save_plot,save_dir  Output options.
#' @return A ggplot object.
#' @export
PlotThresholdDensity <- function(thresh_result,
                                 lower_threshold = NULL,
                                 upper_threshold = NULL,
                                 save_plot = FALSE,
                                 save_dir  = NULL) {
  probs <- thresh_result$probabilities
  true  <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  
  if (is.null(lower_threshold)) {
    lower_threshold <- thresh_result$thresholds["Se_Target"]
    if (is.na(lower_threshold)) lower_threshold <- thresh_result$thresholds["Youden"]
  }
  if (is.null(upper_threshold)) {
    upper_threshold <- thresh_result$thresholds["Sp_Target"]
    if (is.na(upper_threshold)) upper_threshold <- thresh_result$thresholds["Youden"]
  }
  
  df <- data.frame(prob = probs, group = true)
  low_count  <- sum(df$prob <= lower_threshold)
  mid_count  <- sum(df$prob > lower_threshold & df$prob <= upper_threshold)
  high_count <- nrow(df) - low_count - mid_count
  low_npv    <- round(mean(df$group[df$prob <= lower_threshold] == negative) * 100, 1)
  high_ppv   <- round(mean(df$group[df$prob > upper_threshold] == positive) * 100, 1)
  
  cols <- c("#969696", "#fed9a6")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = prob, fill = group)) +
    ggplot2::geom_density(alpha = 0.6, colour = NA) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = lower_threshold,
                      ymin = -Inf, ymax = Inf,
                      fill = "#ffffb3", alpha = 0.15) +
    ggplot2::annotate("rect", xmin = upper_threshold, xmax = Inf,
                      ymin = -Inf, ymax = Inf,
                      fill = "#8dd3c7", alpha = 0.15) +
    ggplot2::geom_vline(xintercept = c(lower_threshold, upper_threshold),
                        linetype = "dashed", color = "grey40") +
    ggplot2::annotate("text", x = lower_threshold/2,
                      y = max(density(df$prob)$y) * 0.9,
                      label = paste0("Low: ", low_count, " (", round(low_count/nrow(df)*100,1),
                                     "%)\nNPV: ", low_npv, "%"), size = 3.5) +
    ggplot2::annotate("text", x = (upper_threshold + 1)/2,
                      y = max(density(df$prob)$y) * 0.9,
                      label = paste0("High: ", high_count, " (", round(high_count/nrow(df)*100,1),
                                     "%)\nPPV: ", high_ppv, "%"), size = 3.5) +
    ggplot2::labs(title = "Prediction Density with Decision Zones",
                  x = "Predicted Probability", y = "Density") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = "top")
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "threshold_density.pdf"),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  return(p)
}

# ── 8. Waterfall plot ──────────────────────────────────────────────────────
#' Waterfall Plot for Threshold Classification
#'
#' @param thresh_result Output from \code{CalculateThresholds}.
#' @param which_threshold Name of the threshold in \code{thresh_result$thresholds}.
#' @param save_plot,save_dir,colors  Output options.
#' @return A ggplot object.
#' @export
PlotThresholdWaterfall <- function(thresh_result,
                                   which_threshold = "Youden",
                                   save_plot = FALSE,
                                   save_dir  = NULL,
                                   colors = c("#f1a340", "#998ec3")) {
  probs <- thresh_result$probabilities
  true  <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  thr <- thresh_result$thresholds[which_threshold]
  
  df <- data.frame(id = seq_along(probs), prob = probs, truth = true)
  df$dif <- df$prob - thr
  df <- df[order(df$dif), ]
  df$predict <- ifelse(df$prob > thr, positive, negative)
  df$correct <- ifelse(df$predict == df$truth, "Correct", "Wrong")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(id, dif), y = dif, fill = correct)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(title = paste("Waterfall Plot —", which_threshold, "Threshold"),
                  x = "Samples (sorted)", y = "Difference from threshold") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("waterfall_", which_threshold, ".pdf")),
                    plot = p, width = 8, height = 5, dpi = 300)
  }
  return(p)
}

# ── 9. Confusion matrix for a chosen threshold ────────────────────────────
#' Confusion Matrix Heatmap (Customizable Colors)
#'
#' @param thresh_result Output from \code{CalculateThresholds}.
#' @param which_threshold Name of the threshold.
#' @param save_plot Logical. Save plot?
#' @param save_dir Output directory.
#' @param fill_colors Character vector of colors for gradient.
#'   Default \code{c("#d8b365", "#f5f5f5", "#5ab4ac")}.
#' @return A ggplot object.
#' @export
PlotThresholdConfusion <- function(thresh_result,
                                   which_threshold = "Youden",
                                   save_plot = FALSE,
                                   save_dir = NULL,
                                   fill_colors = c("#d8b365", "#f5f5f5", "#5ab4ac")) {
  probs <- thresh_result$probabilities
  true  <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  thr <- thresh_result$thresholds[which_threshold]
  
  pred_class <- factor(ifelse(probs > thr, positive, negative),
                       levels = c(negative, positive))
  cm <- caret::confusionMatrix(pred_class, true, positive = positive)
  tab <- as.data.frame(cm$table)
  colnames(tab) <- c("Predicted", "Actual", "Freq")
  tab <- tab %>%
    dplyr::group_by(Actual) %>%
    dplyr::mutate(Pct = round(Freq / sum(Freq) * 100, 1)) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = Actual, y = Predicted, fill = Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Freq, "\n(", Pct, "%)")),
                       size = 5, fontface = "bold") +
    ggplot2::scale_fill_gradientn(colours = fill_colors) +
    ggplot2::labs(title = paste("Confusion Matrix —", which_threshold, "Threshold"),
                  x = "Actual", y = "Predicted", fill = "Count") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, paste0("confusion_", which_threshold, ".pdf")),
                    plot = p, width = 5, height = 4.5, dpi = 300)
  }
  return(p)
}

# ── 10. Multi‑threshold ROC comparison ─────────────────────────────────────
#' Compare Threshold‑Based Classifiers with Original Score (Final, Fixed)
#'
#' @param thresh_result Output from `CalculateThresholds`.
#' @param compare_model A second `thresh_result` for another model (optional).
#' @param compare_label Label for the comparison model.
#' @param save_plot,save_dir  Output options.
#' @return A ggplot object.
#' @export
PlotThresholdROC <- function(thresh_result,
                             compare_model = NULL,
                             compare_label = "Clinician",
                             save_plot = FALSE,
                             save_dir  = NULL) {
  probs <- thresh_result$probabilities
  true  <- thresh_result$true
  positive <- thresh_result$positive
  negative <- thresh_result$negative
  
  roc_main <- pROC::roc(true, probs, levels = c(negative, positive),
                        direction = "auto", quiet = TRUE)
  auc_main <- round(as.numeric(pROC::auc(roc_main)), 3)
  
  thr_names <- names(thresh_result$thresholds)
  # Extract coordinates safely
  sens_sp <- lapply(thr_names, function(nm) {
    tval <- thresh_result$thresholds[nm]
    co <- tryCatch({
      res <- pROC::coords(roc_main, tval, ret = c("se", "sp"), best.method = "closest")
      se <- if (is.list(res)) as.numeric(res$sensitivity[1]) else as.numeric(res[1])
      sp <- if (is.list(res)) as.numeric(res$specificity[1]) else as.numeric(res[2])
      data.frame(Threshold = nm, Sensitivity = se, Specificity = sp)
    }, error = function(e) NULL)
    if (is.null(co) || any(is.na(co[, c("Sensitivity", "Specificity")]))) return(NULL)
    co
  })
  thr_df <- do.call(rbind, sens_sp)
  if (is.null(thr_df) || nrow(thr_df) == 0) stop("No valid threshold coordinates could be extracted.")
  thr_df$FPR <- 1 - thr_df$Specificity
  # Rename Sensitivity to TPR for convenience (optional, we'll just use Sensitivity)
  # We'll use Sensitivity directly in geom_point
  
  roc_df <- data.frame(FPR = 1 - roc_main$specificities,
                       TPR = roc_main$sensitivities)
  n_thr <- nrow(thr_df)
  cols <- wesanderson::wes_palette("Darjeeling1", max(4, n_thr + 1), type = "discrete")
  
  p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = cols[1], size = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    # Explicit mapping for points: x = FPR, y = Sensitivity (= TPR)
    ggplot2::geom_point(data = thr_df, size = 4,
                        mapping = ggplot2::aes(x = FPR, y = Sensitivity),
                        color = cols[2:(n_thr+1)]) +
    ggrepel::geom_text_repel(data = thr_df,
                             mapping = ggplot2::aes(x = FPR, y = Sensitivity, label = Threshold),
                             size = 3.5) +
    ggplot2::annotate("text", x = 0.75, y = 0.25,
                      label = paste0("AUC = ", auc_main), size = 4, color = cols[1]) +
    ggplot2::labs(title = "ROC with Threshold Operating Points",
                  x = "1 – Specificity", y = "Sensitivity") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(compare_model)) {
    roc_comp <- pROC::roc(compare_model$true, compare_model$probabilities,
                          levels = c(negative, positive), direction = "auto", quiet = TRUE)
    auc_comp <- round(as.numeric(pROC::auc(roc_comp)), 3)
    comp_df <- data.frame(FPR = 1 - roc_comp$specificities,
                          TPR = roc_comp$sensitivities)
    p <- p + ggplot2::geom_line(data = comp_df, color = cols[n_thr+2], size = 1.2, linetype = "dashed") +
      ggplot2::annotate("text", x = 0.75, y = 0.15,
                        label = paste0(compare_label, " AUC = ", auc_comp),
                        size = 4, color = cols[n_thr+2])
  }
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "threshold_roc.pdf"),
                    plot = p, width = 7, height = 6, dpi = 300)
  }
  return(p)
}

# ── 11. NRI / IDI analysis ─────────────────────────────────────────────────
#' Optimized NRI/IDI Plot (Academic Style)
#'
#' Generates a horizontal bar plot for NRI and IDI with 95% CI.
#' Style matches academic publications (e.g., Figure 5 in Lambert Leong, Figure 3 in Chen et al. 2022).
#'
#' @param thresh_result1 First threshold result.
#' @param thresh_result2 Second threshold result.
#' @param label1,label2 Model labels.
#' @param cutoffs Numeric vector of risk cutoffs.
#' @param save_plot Logical. Save plot?
#' @param save_dir Output directory.
#' @return Invisible list with NRI/IDI results and a plot.
#' @export
CalculateNRI <- function(thresh_result1,
                         thresh_result2,
                         label1 = "Model 1",
                         label2 = "Model 2",
                         cutoffs = c(0.5),
                         save_plot = FALSE,
                         save_dir  = NULL) {
  .check_clinical_pkgs()
  probs1 <- thresh_result1$probabilities
  probs2 <- thresh_result2$probabilities
  true   <- thresh_result1$true
  binary <- as.numeric(true == thresh_result1$positive)
  
  # NRI
  nri_obj <- tryCatch(
    nricens::nribin(event = binary, p.std = probs2, p.new = probs1,
                    cut = cutoffs, updown = "category"),
    error = function(e) NULL
  )
  
  # IDI
  idii <- tryCatch(
    nricens::improveProb(x1 = probs2, x2 = probs1, y = binary),
    error = function(e) NULL
  )
  
  # Extract estimates safely
  nri_plus <- tryCatch(nri_obj$nri["NRI+", "Estimate"], error = function(e) NA_real_)
  nri_se   <- tryCatch(nri_obj$nri["NRI+", "Std.Error"], error = function(e) NA_real_)
  idi_val  <- if (!is.null(idii)) idii$idi else NA_real_
  idi_se   <- if (!is.null(idii)) idii$se else NA_real_
  
  # Build data frame (in percentage)
  df <- data.frame(
    Metric   = c("NRI+", "IDI"),
    Estimate = c(nri_plus * 100, idi_val * 100),
    Lower    = c((nri_plus - 1.96 * nri_se) * 100,
                 (idi_val  - 1.96 * idi_se) * 100),
    Upper    = c((nri_plus + 1.96 * nri_se) * 100,
                 (idi_val  + 1.96 * idi_se) * 100)
  )
  df <- df[!is.na(df$Estimate), ]
  
  if (nrow(df) == 0) {
    cat("NRI/IDI could not be estimated.\n")
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "Not estimable", size = 6) +
      ggplot2::theme_void()
  } else {
    # Create the plot with academic styling
    p <- ggplot2::ggplot(df, ggplot2::aes(x = Estimate, y = Metric)) +
      # Reference line at zero
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey70", linewidth = 0.8) +
      # Error bars
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = Lower, xmax = Upper),
        height = 0.15,
        linewidth = 1.0,
        color = "#2166ac"
      ) +
      # Points
      ggplot2::geom_point(size = 3, color = "#2166ac", shape = 16) +
      # Labels
      ggplot2::labs(
        title = paste("NRI/IDI:", label1, "vs", label2),
        subtitle = "Error bars represent 95% CI",
        x = "Value (%)",
        y = NULL
      ) +
      # Add reference region shading (optional, for clinical significance)
      ggplot2::annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf,
                        fill = "grey90", alpha = 0.3) +
      # Theme
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40", size = 10),
        axis.title.x  = ggplot2::element_text(face = "bold"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.background   = ggplot2::element_rect(fill = "white", color = NA),
        plot.background    = ggplot2::element_rect(fill = "white", color = NA)
      )
  }
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "nri_idi.pdf"),
                    plot = p, width = 6, height = 3, dpi = 300)
  }
  invisible(list(nri = nri_obj, idi = idii, plot = p))
}
#' Apply a Threshold to Predicted Probabilities
#'
#' Uses a specified threshold from \code{CalculateThresholds} (or a custom value)
#' to convert probabilities into binary predictions, then computes a
#' confusion matrix and common performance metrics.
#'
#' @param thresh_result Output of \code{CalculateThresholds}.
#' @param which_threshold Name of the threshold to use (e.g., "Youden", "PPV_Target").
#'   Ignored if \code{custom_threshold} is provided.
#' @param custom_threshold Numeric threshold value to override \code{which_threshold}.
#' @param newdata Optional data frame (if NULL, uses the true labels stored in thresh_result).
#' @param positive_class Positive class name (optional, auto-detected).
#' @return Invisibly returns a list with \code{predictions}, \code{conf_matrix}, \code{metrics}.
#' @export
ApplyThreshold <- function(thresh_result,
                           which_threshold = "Youden",
                           custom_threshold = NULL,
                           newdata = NULL,
                           positive_class = NULL) {
  prob <- thresh_result$probabilities
  true <- thresh_result$true
  positive <- if (!is.null(positive_class)) positive_class else thresh_result$positive
  negative <- setdiff(levels(true), positive)
  
  thr <- if (!is.null(custom_threshold)) custom_threshold else thresh_result$thresholds[which_threshold]
  
  pred_class <- factor(ifelse(prob > thr, positive, negative),
                       levels = c(negative, positive))
  # If newdata is provided, extract true labels from it (using same outcome column as model_obj)
  if (!is.null(newdata)) {
    # Assumes the group column is stored in thresh_result? Actually we need the outcome column.
    # Better: use thresh_result$true directly, but if newdata has different labels, we handle it.
    # For simplicity, we require that the outcome variable name be passed, but here we assume
    # the user wants to use the same true labels as the threshold calculation.
    # We'll just use thresh_result$true.
    warning("newdata argument not fully implemented; using original true labels.")
  }
  
  cm <- caret::confusionMatrix(pred_class, true, positive = positive)
  metrics <- data.frame(
    Threshold = thr,
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"],
    Accuracy = cm$overall["Accuracy"],
    F1 = cm$byClass["F1"]
  )
  cat(sprintf("--- Applied Threshold = %.4f ---\n", thr))
  print(cm$table)
  print(metrics)
  invisible(list(predictions = pred_class, conf_matrix = cm, metrics = metrics))
}

#' Compare Classification Performance at Specified Thresholds (Bar Plot)
#'
#' Applies user‑chosen thresholds (or Youden) to two sets of predicted
#' probabilities, computes confusion matrices and common metrics, and displays
#' a side‑by‑side bar plot.
#'
#' @param thresh_result1 First threshold result (from \code{CalculateThresholds}).
#' @param thresh_result2 Second threshold result.
#' @param thr1 Threshold for model 1. If \code{NULL}, uses \code{"Youden"} from result1.
#' @param thr2 Threshold for model 2. If \code{NULL}, uses \code{"Youden"} from result2.
#' @param label1,label2 Model labels.
#' @param palette_name Wes Anderson palette for the two models (default \code{"Darjeeling1"}).
#' @param save_plot Save plot?
#' @param save_dir Output directory.
#' @return Invisible list with metrics and plots.
#' @export
CompareClassification <- function(thresh_result1,
                                  thresh_result2,
                                  thr1 = NULL,
                                  thr2 = NULL,
                                  label1 = "Model 1",
                                  label2 = "Model 2",
                                  palette_name = "Darjeeling1",
                                  save_plot = FALSE,
                                  save_dir = NULL) {
  .check_clinical_pkgs()
  # Resolve thresholds
  if (is.null(thr1)) thr1 <- thresh_result1$thresholds["Youden"]
  if (is.null(thr2)) thr2 <- thresh_result2$thresholds["Youden"]
  
  cat(sprintf("Threshold %s: %.4f\n", label1, thr1))
  cat(sprintf("Threshold %s: %.4f\n", label2, thr2))
  
  # Apply thresholds
  pred1 <- factor(ifelse(thresh_result1$probabilities > thr1, 
                         thresh_result1$positive, thresh_result1$negative),
                  levels = c(thresh_result1$negative, thresh_result1$positive))
  pred2 <- factor(ifelse(thresh_result2$probabilities > thr2,
                         thresh_result2$positive, thresh_result2$negative),
                  levels = c(thresh_result2$negative, thresh_result2$positive))
  true <- thresh_result1$true   # same labels for both
  
  cm1 <- caret::confusionMatrix(pred1, true, positive = thresh_result1$positive)
  cm2 <- caret::confusionMatrix(pred2, true, positive = thresh_result2$positive)
  
  # Extract metrics
  extract_metrics <- function(cm, model_name) {
    data.frame(Model = model_name,
               Sensitivity = cm$byClass["Sensitivity"],
               Specificity = cm$byClass["Specificity"],
               PPV = cm$byClass["Pos Pred Value"],
               NPV = cm$byClass["Neg Pred Value"],
               Accuracy = cm$overall["Accuracy"],
               F1 = cm$byClass["F1"],
               row.names = NULL)
  }
  metrics_df <- rbind(extract_metrics(cm1, label1),
                      extract_metrics(cm2, label2))
  
  # Long format for bar plot
  long_df <- tidyr::pivot_longer(metrics_df, -Model, names_to = "Metric", values_to = "Value")
  
  cols <- wesanderson::wes_palette(palette_name, 2, type = "discrete")
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = Metric, y = Value, fill = Model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", Value)),
      position = ggplot2::position_dodge(width = 0.8),
      vjust = -0.5, size = 3.5, colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::labs(title = "Classification Performance at Chosen Thresholds",
                  x = NULL, y = "Value") +
    ggprism::theme_prism(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(save_dir, "compare_classification.pdf"),
                    plot = p, width = 7, height = 5, dpi = 300)
  }
  
  invisible(list(metrics = metrics_df, cm1 = cm1, cm2 = cm2, plot = p))
}
# ── 12. Standalone step functions ──────────────────────────────────────────
#' Clinical Correlation Analysis (Standalone)
#' @inheritParams ClinicalAnalysis
#' @export
ClinicalCorrelation <- function(model_obj, clinical_data = NULL, ...) {
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  PlotClinicalCorrelation(model_obj, clinical_data, ...)
}

#' Subgroup Analysis (Standalone)
#' @inheritParams ClinicalAnalysis
#' @export
ClinicalSubgroup <- function(model_obj, clinical_data = NULL,
                             subgroup_vars = NULL, ...) {
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  for (v in subgroup_vars)
    PlotSubgroupForest(model_obj, clinical_data, v, ...)
}

#' Confounder Adjustment (Standalone)
#' @inheritParams ClinicalAnalysis
#' @export
ClinicalConfounder <- function(model_obj, clinical_data = NULL, ...) {
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  PlotConfounderForest(model_obj, clinical_data, ...)
}

#' Threshold Analysis (Standalone) – Dual Mode
#'
#' Can be used in two ways:
#' 1. Provide \code{model_obj} (and optional \code{model_name}, \code{newdata}),
#'    thresholds are calculated internally and then visualized.
#' 2. Provide a pre‑computed \code{thresh} object (from
#'    \code{CalculateThresholds} or \code{CalculateThresholdsFromProbs}),
#'    all plots are generated directly.  \code{model_obj} can be \code{NULL}.
#'
#' @param model_obj  A Train_Model or caret model.  If \code{thresh} is supplied,
#'   this is ignored.
#' @param thresh     Optional pre‑computed threshold result.  When provided,
#'   \code{model_obj}, \code{newdata}, \code{model_name} and all \code{...}
#'   arguments are not used.
#' @param compare_model  Optional second \code{thresh_result} for ROC comparison and NRI.
#' @param compare_label  Label for the comparison model.
#' @param save_plot  Save plots?
#' @param save_dir   Output directory.
#' @param ...        Further arguments passed to \code{CalculateThresholds} when
#'   \code{thresh} is not supplied (e.g., \code{target_ppv}, \code{target_npv}).
#' @return Invisible list with threshold results.
#' @export
ClinicalThreshold <- function(model_obj = NULL,
                              thresh = NULL,
                              newdata = NULL,
                              model_name = NULL,
                              compare_model = NULL,
                              compare_label = "Comparator",
                              save_plot = FALSE,
                              save_dir = NULL,
                              ...) {
  .check_clinical_pkgs()
  
  # --- Case 2: pre‑computed thresholds ---
  if (!is.null(thresh)) {
    if (!is.list(thresh) || is.null(thresh$thresholds))
      stop("'thresh' must be a list returned by CalculateThresholds or CalculateThresholdsFromProbs.")
    cat("Using pre‑computed threshold object...\n")
    # No clinical data needed for plots
  } else {
    # --- Case 1: compute from model object ---
    if (is.null(model_obj))
      stop("Either 'model_obj' or 'thresh' must be provided.")
    if (is.null(newdata) && inherits(model_obj, "Train_Model"))
      newdata <- model_obj@filtered.set$testing
    thresh <- CalculateThresholds(model_obj, newdata, model_name, ...)
  }
  
  # ---- Visualizations (all use the same 'thresh' object) ----
  PlotThresholdAccuracy(thresh, save_plot = save_plot, save_dir = save_dir)
  PlotThresholdDensity(thresh, save_plot = save_plot, save_dir = save_dir)
  for (thr_name in names(thresh$thresholds)) {
    PlotThresholdWaterfall(thresh, thr_name, save_plot = save_plot, save_dir = save_dir)
    PlotThresholdConfusion(thresh, thr_name, save_plot = save_plot, save_dir = save_dir)
  }
  PlotThresholdROC(thresh, compare_model = compare_model,
                   compare_label = compare_label, save_plot = save_plot, save_dir = save_dir)
  if (!is.null(compare_model)) {
    CalculateNRI(thresh, compare_model,
                 label1 = "Model", label2 = compare_label,
                 save_plot = save_plot, save_dir = save_dir)
  }
  
  invisible(thresh)
}
#' Compare Two Models via Threshold ROC and NRI/IDI (Standalone)
#'
#' @param model_obj1 First model (Train_Model or caret train).
#' @param model_obj2 Second model (Train_Model or caret train).
#' @param newdata Data frame for prediction (default: testing set).
#' @param model_name1,model_name2 Optional model names within Train_Model objects.
#' @param label1,label2 Labels for the two models.
#' @param save_plot Logical. Save plots?
#' @param save_dir Output directory.
#' @param ... Additional arguments passed to `CalculateThresholds` for the first model
#'   (e.g., target_ppv, target_npv).
#' @return Invisible list with threshold results and comparison plots.
#' @export
CompareModelThresholds <- function(model_obj1,
                                   model_obj2,
                                   newdata = NULL,
                                   model_name1 = NULL,
                                   model_name2 = NULL,
                                   label1 = "Model 1",
                                   label2 = "Model 2",
                                   save_plot = FALSE,
                                   save_dir = NULL,
                                   ...) {
  .check_clinical_pkgs()
  # 分别计算两个模型的阈值（默认仅计算约登指数和最大准确率）
  thresh1 <- CalculateThresholds(model_obj1, newdata, model_name1, ...)
  thresh2 <- CalculateThresholds(model_obj2, newdata, model_name2)
  
  # 绘制对比ROC曲线（标注两个模型的阈值点）
  PlotThresholdROC(thresh1, compare_model = thresh2,
                   compare_label = label2,
                   save_plot = save_plot, save_dir = save_dir)
  
  # 计算NRI/IDI
  CalculateNRI(thresh1, thresh2,
               label1 = label1, label2 = label2,
               save_plot = save_plot, save_dir = save_dir)
  
  invisible(list(thresh1 = thresh1, thresh2 = thresh2))
}
# ── 13. Master Clinical Analysis Pipeline ──────────────────────────────────
#' Complete Clinical Analysis Pipeline
#'
#' Runs correlation heatmap, subgroup forests, confounder adjustment,
#' threshold analysis (accuracy curve, density, waterfall, confusion matrix,
#' ROC comparison), and optionally NRI/IDI.
#'
#' @param model_obj      A Train_Model or caret model.
#' @param clinical_data  Data frame of clinical variables. If NULL, uses
#'   \code{@process.info$clinical_data} (set by \code{AttachClinicalData}).
#' @param subgroup_vars  Character vector of categorical variables for subgroup
#'   analysis. If NULL, skipped.
#' @param outcome_var    Binary outcome column name.
#' @param newdata        Data for prediction (default: \code{@filtered.set$testing}).
#' @param model_name     Which model to use (NULL = best, "ensemble", model
#'   name, or caret object).
#' @param compare_model  Optional second \code{thresh_result} for ROC comparison and NRI.
#' @param compare_label  Label for the comparison model.
#' @param save_plots     Save all plots?
#' @param save_dir       Output directory.
#' @param ...            Additional arguments passed to \code{CalculateThresholds}
#'   (e.g., \code{target_se}, \code{target_ppv}).
#' @return Invisible list containing threshold results.
#' @export
ClinicalAnalysis <- function(model_obj,
                             clinical_data  = NULL,
                             subgroup_vars  = NULL,
                             outcome_var    = "group",
                             newdata        = NULL,
                             model_name     = NULL,
                             compare_model  = NULL,
                             compare_label  = "Comparator",
                             save_plots     = TRUE,
                             save_dir       = "./ClinicalResults/",
                             ...) {
  .check_clinical_pkgs()
  if (is.null(clinical_data))
    clinical_data <- model_obj@process.info$clinical_data
  if (is.null(newdata) && inherits(model_obj, "Train_Model"))
    newdata <- model_obj@filtered.set$testing
  if (is.null(save_dir) && save_plots) save_dir <- "./ClinicalResults/"
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  cat("--- Clinical Analysis Pipeline ---\n")
  
  # 1. Correlation matrix
  cat("[1/8] Correlation matrix...\n")
  PlotClinicalCorrelation(model_obj, clinical_data, newdata, model_name,
                          save_plot = save_plots, save_dir = save_dir)
  
  # 2. Subgroup forests
  if (!is.null(subgroup_vars)) {
    cat("[2/8] Subgroup analysis...\n")
    for (v in subgroup_vars) {
      PlotSubgroupForest(model_obj, clinical_data, v, newdata, model_name,
                         save_plot = save_plots, save_dir = save_dir)
    }
  } else {
    cat("[2/8] Subgroup analysis skipped.\n")
  }
  
  # 3. Confounder adjustment
  cat("[3/8] Confounder adjustment...\n")
  PlotConfounderForest(model_obj, clinical_data, outcome_var, newdata, model_name,
                       save_plot = save_plots, save_dir = save_dir)
  
  # 4. Threshold calculation
  cat("[4/8] Calculating thresholds...\n")
  thresh <- CalculateThresholds(model_obj, newdata, model_name, ...)
  print(thresh$thresholds)
  
  # 5. Accuracy / PPV / NPV curve
  cat("[5/8] Threshold accuracy curve...\n")
  PlotThresholdAccuracy(thresh, save_plot = save_plots, save_dir = save_dir)
  
  # 6. Density zones
  cat("[6/8] Decision density plot...\n")
  PlotThresholdDensity(thresh, save_plot = save_plots, save_dir = save_dir)
  
  # 7. Waterfall + confusion matrices
  cat("[7/8] Waterfall and confusion plots...\n")
  for (thr_name in names(thresh$thresholds)) {
    PlotThresholdWaterfall(thresh, thr_name, save_plot = save_plots, save_dir = save_dir)
    PlotThresholdConfusion(thresh, thr_name, save_plot = save_plots, save_dir = save_dir)
  }
  
  # 8. ROC comparison + NRI
  cat("[8/8] ROC comparison and NRI...\n")
  PlotThresholdROC(thresh, compare_model = compare_model, compare_label = compare_label,
                   save_plot = save_plots, save_dir = save_dir)
  if (!is.null(compare_model)) {
    CalculateNRI(thresh, compare_model,
                 label1 = "Novel Model", label2 = compare_label,
                 save_plot = save_plots, save_dir = save_dir)
  }
  
  cat("--- Clinical analysis complete. Plots saved to", save_dir, "---\n")
  invisible(list(thresh = thresh))
}