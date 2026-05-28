## ============================================================
##  viz_functions.R
##  Publication-quality visualization functions
##  Organized by module:
##    § 1. Stat       — distribution, correlation, PCA, DEG
##    § 2. Train_Model — ROC, confusion matrix, feature importance, calibration
##    § 3. Subtyping  — dim-reduction, cluster heatmap, silhouette, alluvial
##    § 4. PrognosiX  — KM, forest, time-ROC, RCS, nomogram, calibration, DCA, risk
##
##  Every function:
##    · accepts an S4 object OR raw data frames
##    · returns a ggplot2 / grob object (further customisable)
##    · saves to PDF / PNG / SVG when save_plot = TRUE
##    · uses theme_prism + wesanderson / RColorBrewer palettes by default
## ============================================================

## ── shared helpers ────────────────────────────────────────────────────────────

.pub_theme <- function(base_size = 13) {
  ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            size  = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.title    = ggplot2::element_text(face  = "bold"),
      legend.title  = ggplot2::element_text(face  = "bold"),
      strip.text    = ggplot2::element_text(face  = "bold")
    )
}

.get_palette <- function(palette_name, n) {
  tryCatch(
    as.character(wesanderson::wes_palette(n = n, name = palette_name,
                                          type = if (n > 5) "continuous" else "discrete")),
    error = function(e)
      RColorBrewer::brewer.pal(max(3L, n), palette_name)[seq_len(n)]
  )
}

.save_plot <- function(p, dir, filename, width, height, format = "pdf") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, paste0(tools::file_path_sans_ext(filename), ".", format))
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
  cat("Plot saved to:", path, "\n")
  invisible(path)
}

#' Get Output Directory for Visualization Functions
#' 
#' Internal helper to get the appropriate output directory based on global config
#' @keywords internal
.get_viz_output_dir <- function(type = c("Stat", "Model", "Subtyping", "PrognosiX")) {
  type <- match.arg(type)
  root <- get_output_root()
  if (is.null(root)) {
    viz_dir <- here::here("Figures", type)
  } else {
    viz_dir <- file.path(root, "Figures", type)
  }
  # Create directory on-demand
  if (!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(viz_dir)
}

## ═════════════════════════════════════════════════════════════════════════════
##  § 1  STAT  ── distribution · correlation · PCA · DEG
## ═════════════════════════════════════════════════════════════════════════════

#' Grouped box-violin plot with statistical comparisons
#'
#' @param object   A \code{Stat} object, or a plain data frame.
#' @param features Character vector of feature names to plot. When \code{NULL}
#'   all numeric columns are used (max 12).
#' @param group_col Grouping column. Ignored when \code{object} is a
#'   \code{Stat} and the slot is already set.
#' @param test      Statistical test passed to \code{ggpubr::stat_compare_means}:
#'   \code{"wilcox.test"} (default) or \code{"t.test"}.
#' @param palette_name Wesanderson / RColorBrewer palette.
#' @param ncol      Number of facet columns.
#' @param base_size Base font size.
#' @param save_plot Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir  Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height  Plot dimensions (inches).
#' @param format    \code{"pdf"}, \code{"png"}, or \code{"svg"}.
#' @returns A \code{ggplot} object.
#' @export
PlotGroupedDistribution <- function(object,
                                    features     = NULL,
                                    group_col    = "group",
                                    test         = "wilcox.test",
                                    palette_name = "Royal1",
                                    ncol         = 3,
                                    base_size    = 13,
                                    save_plot    = FALSE,
                                    save_dir     = NULL,
                                    width        = 10,
                                    height       = 6,
                                    format       = "pdf") {
  cat("Generating grouped distribution plot...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  if (inherits(object, "Stat")) {
    df        <- object@clean.data
    group_col <- object@group_col
  } else {
    df <- as.data.frame(object)
  }
  
  num_cols <- names(df)[sapply(df, is.numeric)]
  num_cols <- setdiff(num_cols, group_col)
  if (is.null(features)) features <- head(num_cols, 12)
  features  <- intersect(features, num_cols)
  if (length(features) == 0) stop("No valid numeric features found.")
  
  cols <- .get_palette(palette_name, length(unique(df[[group_col]])))
  
  long_df <- tidyr::pivot_longer(df[, c(group_col, features), drop = FALSE],
                                 cols = -dplyr::all_of(group_col),
                                 names_to  = "Feature",
                                 values_to = "Value")
  long_df[[group_col]] <- as.factor(long_df[[group_col]])
  
  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = .data[[group_col]],
                                    y = Value,
                                    fill = .data[[group_col]])) +
    ggplot2::geom_violin(alpha = 0.5, colour = NA) +
    ggplot2::geom_boxplot(width = 0.2, outlier.shape = 21,
                          outlier.size = 1.2, colour = "grey30") +
    ggpubr::stat_compare_means(method = test,
                               label  = "p.signif",
                               comparisons = utils::combn(
                                 levels(long_df[[group_col]]), 2, simplify = FALSE)) +
    ggplot2::facet_wrap(~ Feature, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = group_col, y = "Value", fill = group_col) +
    .pub_theme(base_size)
  
  if (save_plot) .save_plot(p, save_dir, "grouped_distribution", width, height, format)
  return(p)
}


#' Publication‑ready correlation heatmap with advanced colour palettes (enhanced)
#'
#' @param object A Stat object or data frame.
#' @param features Features to include (NULL = all numeric, max 40).
#' @param method "pearson" (default) or "spearman".
#' @param p_thresh P‑value threshold for masking non‑significant pairs.
#' @param palette_name RColorBrewer diverging palette name (used only when
#'   \code{color_scheme = "brew"}). Default "RdYlBu".
#' @param color_scheme Colour scheme: \code{"viridis"} (colour‑blind friendly,
#'   default), \code{"brew"} (traditional RColorBrewer), \code{"gradient"}
#'   (custom blue‑white‑red gradient).
#' @param cluster_rows Logical. Hierarchical clustering of rows/cols.
#' @param base_size Base font size.
#' @param save_plot Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height Plot dimensions (inches).
#' @param format File format ("pdf", "png", "svg").
#'
#' @returns A ggplot object.
#' @export
PlotCorrelationHeatmap <- function(object,
                                   features      = NULL,
                                   method        = "pearson",
                                   p_thresh      = 0.05,
                                   palette_name  = "RdYlBu",
                                   color_scheme  = c("viridis", "brew", "gradient"),
                                   cluster_rows  = TRUE,
                                   base_size     = 11,
                                   save_plot     = FALSE,
                                   save_dir      = NULL,
                                   width         = 8,
                                   height        = 7,
                                   format        = "pdf") {
  cat("Generating correlation heatmap...\n")
  color_scheme <- match.arg(color_scheme)
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  df <- if (inherits(object, "Stat")) object@clean.data else as.data.frame(object)
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (is.null(features)) features <- head(num_cols, 40)
  features <- intersect(features, num_cols)
  if (length(features) == 0) stop("No numeric features found.")
  
  mat <- df[, features, drop = FALSE]
  mat <- mat[, apply(mat, 2, stats::var, na.rm = TRUE) > 0]
  
  cor_out <- Hmisc::rcorr(as.matrix(mat), type = method)
  r_mat   <- cor_out$r
  p_mat   <- cor_out$P
  
  if (cluster_rows) {
    ord <- hclust(as.dist(1 - r_mat))$order
    r_mat <- r_mat[ord, ord]
    p_mat <- p_mat[ord, ord]
  }
  
  r_mat[upper.tri(r_mat)] <- NA
  p_mat[upper.tri(p_mat)] <- NA
  
  long_r <- reshape2::melt(r_mat, na.rm = TRUE, value.name = "r")
  long_p <- reshape2::melt(p_mat, na.rm = TRUE, value.name = "p")
  long_df <- merge(long_r, long_p, by = c("Var1", "Var2"))
  long_df$sig <- ifelse(long_df$p < p_thresh, "", "×")
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(Var1, Var2, fill = r))
  
  if (color_scheme == "viridis") {
    p <- p + ggplot2::scale_fill_viridis_c(option = "C", limit = c(-1, 1),
                                           name = paste0(tools::toTitleCase(method), "\n r"))
  } else if (color_scheme == "brew") {
    p <- p + ggplot2::scale_fill_distiller(palette = palette_name, limit = c(-1, 1),
                                           name = paste0(tools::toTitleCase(method), "\n r"))
  } else if (color_scheme == "gradient") {
    p <- p + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                           high = "#B2182B", midpoint = 0,
                                           limit = c(-1, 1),
                                           name = paste0(tools::toTitleCase(method), "\n r"))
  }
  
  p <- p +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sig), size = 3, colour = "grey30", na.rm = TRUE) +
    ggplot2::labs(title = "Correlation Heatmap", x = NULL, y = NULL) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, colour = "grey20"),
      axis.text.y   = ggplot2::element_text(colour = "grey20"),
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      panel.grid    = ggplot2::element_blank(),
      legend.title  = ggplot2::element_text(face = "bold")
    )
  
  if (save_plot) .save_plot(p, save_dir, "correlation_heatmap", width, height, format)
  return(p)
}

#' PCA scatter plot coloured by metadata (zero‑variance safe, colour‑named, enhanced)
#'
#' @param object     A \code{Stat} object or numeric data frame.
#' @param color_by   Column in \code{info.data} (or the data frame) used for
#'   point colour. Defaults to \code{group_col}.
#' @param shape_by   Optional second metadata column for point shape.
#' @param pcs        Integer vector length 2 selecting PCs. Default \code{c(1,2)}.
#' @param ellipse    Draw group ellipses (95 \% CI). Default \code{TRUE}.
#' @param label_points Logical. Label sample names. Default \code{FALSE}.
#' @param palette_name Palette.
#' @param base_size  Base font size.
#' @param save_plot  Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir   Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height Inches.
#' @param format     File format.
#' @returns A \code{ggplot} object.
#' @export
PlotPCA <- function(object,
                    color_by      = NULL,
                    shape_by      = NULL,
                    pcs           = c(1, 2),
                    ellipse       = TRUE,
                    label_points  = FALSE,
                    palette_name  = "Darjeeling1",
                    base_size     = 13,
                    save_plot     = FALSE,
                    save_dir      = NULL,
                    width         = 7,
                    height        = 6,
                    format        = "pdf") {
  cat("Generating PCA plot...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  if (inherits(object, "Stat")) {
    mat      <- object@clean.data[, sapply(object@clean.data, is.numeric)]
    info     <- object@info.data
    color_by <- if (is.null(color_by)) object@group_col else color_by
  } else {
    mat  <- as.data.frame(object)[, sapply(as.data.frame(object), is.numeric)]
    info <- data.frame(row.names = rownames(mat))
  }
  
  keep_cols <- apply(mat, 2, var, na.rm = TRUE) > 0
  if (any(!keep_cols)) {
    removed <- colnames(mat)[!keep_cols]
    message("Removed ", length(removed), " zero‑variance column(s): ",
            paste(removed, collapse = ", "))
    mat <- mat[, keep_cols, drop = FALSE]
  }
  if (ncol(mat) < 2) stop("After removing constant columns, fewer than 2 variables remain.")
  
  pca_res  <- stats::prcomp(mat, scale. = TRUE, center = TRUE)
  var_exp  <- round(summary(pca_res)$importance[2, ] * 100, 1)
  scores   <- as.data.frame(pca_res$x[, pcs])
  colnames(scores) <- paste0("PC", pcs)
  
  if (nrow(info) == nrow(scores) && color_by %in% colnames(info)) {
    scores[[color_by]] <- as.factor(info[rownames(scores), color_by])
  } else {
    scores[[color_by]] <- factor("all")
  }
  
  if (!is.null(shape_by) && shape_by %in% colnames(info)) {
    scores[[shape_by]] <- as.factor(info[rownames(scores), shape_by])
  }
  
  n_grp <- length(levels(scores[[color_by]]))
  cols  <- .get_palette(palette_name, n_grp)
  names(cols) <- levels(scores[[color_by]])
  
  xlab <- paste0("PC", pcs[1], " (", var_exp[pcs[1]], "%)")
  ylab <- paste0("PC", pcs[2], " (", var_exp[pcs[2]], "%)")
  
  if (!is.null(shape_by) && shape_by %in% colnames(scores)) {
    aes_map <- ggplot2::aes(.data[[paste0("PC", pcs[1])]],
                            .data[[paste0("PC", pcs[2])]],
                            colour = .data[[color_by]],
                            shape  = .data[[shape_by]])
  } else {
    aes_map <- ggplot2::aes(.data[[paste0("PC", pcs[1])]],
                            .data[[paste0("PC", pcs[2])]],
                            colour = .data[[color_by]])
  }
  
  p <- ggplot2::ggplot(scores, aes_map) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(title = "PCA", x = xlab, y = ylab, colour = color_by) +
    .pub_theme(base_size)
  
  if (ellipse) p <- p + ggplot2::stat_ellipse(level = 0.95, linewidth = 0.7, linetype = 2)
  if (label_points) p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = rownames(scores)),
                                                      size = 2.5, max.overlaps = 20)
  if (save_plot) .save_plot(p, save_dir, "PCA", width, height, format)
  return(p)
}
#' Feature selection plot: AUC vs. –log10(p‑value) (enhanced)
#'
#' Calculates per‑feature ROC AUC and plots it against –log10(p‑value).
#' Points that satisfy both \code{auc_thresh} and \code{p_thresh} are highlighted.
#'
#' @param deg_df       Data frame with columns \code{feature}, \code{logFC}, \code{p.adjust} / \code{p_value}.
#' @param mat_test     Expression matrix containing the grouping column.
#' @param group_col    Name of the grouping column in \code{mat_test}.
#' @param auc_thresh   AUC threshold for selection. Default 0.55.
#' @param p_thresh     P‑value threshold. Default 0.05.
#' @param selected_fill Fill colour for selected features.
#' @param removed_fill  Fill colour for removed features.
#' @param arrow        Logical. Draw arrows from labels to points.
#' @param base_size    Base font size.
#' @param save_plot    Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir     Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height Plot dimensions (inches).
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotAUCPval <- function(deg_df, mat_test,
                        group_col      = "group",
                        auc_thresh     = 0.55,
                        p_thresh       = 0.05,
                        selected_fill  = "#e84118",
                        removed_fill   = "#7f8fa6",
                        arrow          = TRUE,
                        base_size      = 13,
                        save_plot      = FALSE,
                        save_dir       = NULL,
                        width          = 7,
                        height         = 6,
                        format         = "pdf") {
  cat("Generating AUC-P value selection plot...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  df <- as.data.frame(deg_df)
  if (!"feature" %in% colnames(df)) df$feature <- rownames(df)
  p_col <- if ("p.adjust" %in% colnames(df)) "p.adjust" else "p_value"
  df$neg_log10p <- -log10(df[[p_col]] + 1e-300)
  
  mat <- as.data.frame(mat_test)
  if (!group_col %in% colnames(mat)) stop("Group column not found in mat_test.")
  mat[[group_col]] <- as.factor(mat[[group_col]])
  lev <- levels(mat[[group_col]])
  
  auc_vec <- sapply(df$feature, function(feat) {
    if (!feat %in% colnames(mat)) return(NA_real_)
    tem_dat <- mat[, c(feat, group_col)]
    colnames(tem_dat)[1] <- "feature"
    colnames(tem_dat)[2] <- "group"
    tem_dat$feature <- as.numeric(tem_dat$feature)
    roc_obj <- tryCatch(
      pROC::roc(group ~ feature, data = tem_dat, levels = lev, direction = "auto", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) return(NA_real_)
    round(as.numeric(pROC::auc(roc_obj)), 3)
  })
  
  df$auc <- auc_vec
  if (any(is.na(df$auc))) cat("Features with missing AUC:", sum(is.na(df$auc)), "\n")
  
  df$group <- ifelse(!is.na(df$auc) & df$auc >= auc_thresh & df[[p_col]] < p_thresh,
                     "selected", "removed")
  selected_df <- df[df$group == "selected", , drop = FALSE]
  n_sel <- nrow(selected_df)
  cat("Selected features:", n_sel, "\n")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = auc, y = neg_log10p)) +
    ggplot2::geom_vline(xintercept = auc_thresh, linetype = 2,
                        colour = "black", linewidth = 0.8, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = 2,
                        colour = "black", linewidth = 0.8, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(fill = group, size = neg_log10p),
                        colour = "black", stroke = 0.3, shape = 21, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(selected = selected_fill, removed = removed_fill),
                               name = NULL) +
    ggplot2::scale_size(range = c(1, 10), guide = "none") +
    ggrepel::geom_text_repel(
      data = selected_df,
      ggplot2::aes(label = feature),
      size = 3.5, colour = "black", fontface = "bold.italic",
      min.segment.length = 0,
      arrow = if (arrow) ggplot2::arrow(length = ggplot2::unit(0.01, "npc")) else NULL,
      show.legend = FALSE, max.overlaps = 200,
      segment.curvature = -0.05, segment.square = FALSE, segment.inflect = TRUE,
      force = 0.1, nudge_x = 0.02, direction = "y", hjust = 0
    ) +
    ggplot2::annotate("text",
                      x = max(df$auc, na.rm = TRUE) * 0.95,
                      y = max(df$neg_log10p, na.rm = TRUE) * 0.95,
                      label = paste0("Selected: ", n_sel),
                      colour = selected_fill, fontface = "bold", size = 4, hjust = 1) +
    .pub_theme(base_size) +
    ggplot2::labs(title = "Feature Selection by AUC & Significance",
                  x = "ROC AUC", y = expression(-log[10]~"(P‑value)")) +
    ggplot2::theme(legend.position = "right")
  
  if (save_plot) .save_plot(p, save_dir, "PlotAUCPval", width, height, format)
  return(p)
}

#' Boxplot of top differential features with significance stars (enhanced)
#'
#' @param deg_results DEG data frame with columns 'id', 'change', 'logFC'.
#' @param expr_data   Expression matrix (data frame with features + group column).
#' @param group_col   Column name for the group variable.
#' @param top_n       Number of top features to show.
#' @param palette_fill Colour vector of length 2 for groups.
#' @param base_size   Base font size.
#' @param save_plot   Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir    Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height Plot dimensions.
#' @param format      Output format.
#' @return A ggplot object.
#' @export
PlotDegBoxplot <- function(deg_results,
                           expr_data,
                           group_col   = "group",
                           top_n       = 5,
                           palette_fill = c("#edf8b1", "#2c7fb8"),
                           base_size   = 12,
                           save_plot   = FALSE,
                           save_dir    = NULL,
                           width       = 7,
                           height      = 5,
                           format      = "pdf") {
  cat("Generating DEG boxplot...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  sig_df <- deg_results[deg_results$change != "Stable", ]
  sig_df <- sig_df[order(-abs(sig_df$logFC)), ]
  n_up   <- sum(sig_df$change == "Up")
  n_down <- sum(sig_df$change == "Down")
  top_up   <- head(sig_df[sig_df$change == "Up", ], min(ceiling(top_n/2), n_up))
  top_down <- head(sig_df[sig_df$change == "Down", ], min(floor(top_n/2), n_down))
  selected <- rbind(top_up, top_down)
  features <- selected$id
  
  box_df <- expr_data[, c(group_col, features), drop = FALSE]
  box_df[[group_col]] <- as.factor(box_df[[group_col]])
  long_df <- tidyr::pivot_longer(box_df, cols = -dplyr::all_of(group_col),
                                 names_to = "feature", values_to = "value")
  
  p_vals <- long_df %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(p_value = wilcox.test(value ~ .data[[group_col]])$p.value,
                     .groups = "drop") %>%
    dplyr::mutate(significance = dplyr::case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE           ~ ""
    ))
  
  long_df$log_val <- log10(long_df$value + 1)
  y_pos <- long_df %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(y_position = max(log_val, na.rm = TRUE) + 0.2, .groups = "drop")
  p_vals <- dplyr::left_join(p_vals, y_pos, by = "feature")
  long_df <- dplyr::left_join(long_df, p_vals, by = "feature")
  
  long_df$feature <- factor(long_df$feature,
                            levels = p_vals$feature[order(p_vals$y_position, decreasing = TRUE)])
  
  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = feature, y = log_val, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.size = 0.5, outlier.colour = "grey50") +
    ggplot2::geom_text(
      data = dplyr::distinct(long_df, feature, .keep_all = TRUE),
      ggplot2::aes(x = feature, y = y_position, label = significance),
      size = 5, colour = "black", inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = palette_fill) +
    ggplot2::labs(title = paste0("Top ", top_n, " DEGs"),
                  x = NULL, y = expression(log[10]~(Expression + 1)), fill = group_col) +
    .pub_theme(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  print(p)
  if (save_plot) .save_plot(p, save_dir, "deg_boxplot", width, height, format)
  return(p)
}
#' Feature heatmap (top DEGs or selected markers) - enhanced with progress message
#'
#' @param object       A \code{Stat} object or numeric matrix / data frame.
#' @param features     Features to display. Defaults to top 50 by variance.
#' @param group_col    Annotation column for sample grouping bar.
#' @param scale_rows   Z-score scale across samples. Default \code{TRUE}.
#' @param cluster_rows Cluster rows (features). Default \code{TRUE}.
#' @param cluster_cols Cluster columns (samples). Default \code{TRUE}.
#' @param palette_name RColorBrewer diverging palette. Default \code{"RdYlBu"}.
#' @param ann_palette  Named colour vector for the annotation bar. Auto if \code{NULL}.
#' @param show_rownames Show row (feature) names. Default \code{TRUE}.
#' @param show_colnames Show column (sample) names. Default \code{FALSE}.
#' @param base_fontsize Font size in heatmap.
#' @param save_plot    Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir     Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{pheatmap} object (invisibly).
#' @export
PlotFeatureHeatmap <- function(object,
                               features      = NULL,
                               group_col     = NULL,
                               scale_rows    = TRUE,
                               cluster_rows  = TRUE,
                               cluster_cols  = TRUE,
                               palette_name  = "RdYlBu",
                               ann_palette   = NULL,
                               show_rownames = TRUE,
                               show_colnames = FALSE,
                               base_fontsize = 10,
                               save_plot     = FALSE,
                               save_dir      = NULL,
                               width         = 9,
                               height        = 8,
                               format        = "pdf") {
  cat("Generating feature heatmap...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  if (inherits(object, "Stat")) {
    mat      <- t(object@clean.data[, sapply(object@clean.data, is.numeric)])
    info     <- object@info.data
    group_col <- if (is.null(group_col)) object@group_col else group_col
  } else {
    mat   <- t(as.data.frame(object)[, sapply(as.data.frame(object), is.numeric)])
    info  <- data.frame()
  }
  
  if (is.null(features)) {
    var_ord  <- order(apply(mat, 1, var, na.rm = TRUE), decreasing = TRUE)
    features <- rownames(mat)[head(var_ord, 50)]
  }
  features <- intersect(features, rownames(mat))
  mat <- mat[features, , drop = FALSE]
  
  if (scale_rows) {
    mat <- t(scale(t(mat)))
    mat[is.nan(mat)] <- 0
  }
  
  ann_col <- NULL
  ann_colors <- list()
  if (!is.null(group_col) && nrow(info) > 0 && group_col %in% colnames(info)) {
    ann_col <- data.frame(Group = as.factor(info[[group_col]]),
                          row.names = rownames(info))
    grp_lvls <- levels(ann_col$Group)
    grp_cols <- .get_palette("Darjeeling1", length(grp_lvls))
    names(grp_cols) <- grp_lvls
    ann_colors[["Group"]] <- if (!is.null(ann_palette)) ann_palette else grp_cols
  }
  
  color_ramp <- rev(RColorBrewer::brewer.pal(11, palette_name))
  
  ph <- pheatmap::pheatmap(
    mat,
    color            = colorRampPalette(color_ramp)(100),
    cluster_rows     = cluster_rows,
    cluster_cols     = cluster_cols,
    annotation_col   = ann_col,
    annotation_colors= ann_colors,
    show_rownames    = show_rownames,
    show_colnames    = show_colnames,
    fontsize         = base_fontsize,
    fontsize_row     = base_fontsize - 1,
    border_color     = NA,
    silent           = TRUE
  )
  
  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, paste0("feature_heatmap.", format))
    if (format == "pdf") {
      pdf(path, width = width, height = height)
      grid::grid.draw(ph$gtable)
      dev.off()
    } else {
      png(path, width = width * 100, height = height * 100, res = 100)
      grid::grid.draw(ph$gtable)
      dev.off()
    }
    cat("Heatmap saved to:", path, "\n")
  }
  return(ph)
}


## ═════════════════════════════════════════════════════════════════════════════
##  § 2  TRAIN_MODEL  ── ROC · confusion matrix · feature importance · calibration
## ═════════════════════════════════════════════════════════════════════════════

#' Multi-model ROC comparison
#'
#' Plots ROC curves for all models stored in \code{train.models}, with AUC
#' values in the legend. A combined AUC summary table is returned invisibly.
#'
#' @param object       A \code{Train_Model} object.
#' @param test_data    Optional: data frame for evaluation. When \code{NULL}
#'   the \code{split.data$test} slot is used.
#' @param palette_name Colour palette.
#' @param show_ci      Draw 95 \% bootstrap CI bands. Default \code{FALSE}.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotMultiROC <- function(object,
                         test_data    = NULL,
                         palette_name = "Darjeeling1",
                         show_ci      = FALSE,
                         base_size    = 13,
                         save_plot    = FALSE,
                         save_dir = NULL,  # default: auto-detect from config
                         width        = 7,
                         height       = 6,
                         format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Model")

  if (!inherits(object, "Train_Model"))
    stop("'object' must be a Train_Model.")
  if (length(object@train.models) == 0)
    stop("No trained models found in object@train.models.")

  td <- if (!is.null(test_data)) test_data else object@split.data$test
  if (is.null(td)) stop("Provide 'test_data' or populate split.data$test.")

  group_col <- as.character(object@group_col)
  truth     <- factor(td[[group_col]])
  pos_level <- levels(truth)[2]

  roc_list <- lapply(names(object@train.models), function(nm) {
    model <- object@train.models[[nm]]
    probs <- tryCatch(
      stats::predict(model, newdata = td, type = "prob")[, pos_level],
      error = function(e)
        as.numeric(stats::predict(model, newdata = td, type = "raw") == pos_level)
    )
    roc_obj <- pROC::roc(truth, probs, levels = levels(truth), quiet = TRUE)
    auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
    coords  <- pROC::coords(roc_obj, "all", ret = c("specificity", "sensitivity"))
    data.frame(Model       = paste0(nm, " (AUC=", auc_val, ")"),
               Specificity = coords$specificity,
               Sensitivity = coords$sensitivity)
  })

  roc_df <- dplyr::bind_rows(roc_list)
  cols   <- .get_palette(palette_name, length(object@train.models))

  p <- ggplot2::ggplot(roc_df,
    ggplot2::aes(1 - Specificity, Sensitivity, colour = Model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "ROC Curves — Model Comparison",
                  x     = "1 – Specificity (FPR)",
                  y     = "Sensitivity (TPR)",
                  colour = NULL) +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = c(0.72, 0.22))

  if (save_plot)
    .save_plot(p, save_dir, "multi_ROC", width, height, format)
  return(p)
}


#' Confusion matrix heatmap
#'
#' Renders a tile-based confusion matrix with cell counts and row-wise
#' percentages.
#'
#' @param object     A \code{Train_Model} object.
#' @param model_name Name of the model in \code{train.models}. When \code{NULL}
#'   the best model is used.
#' @param test_data  Optional evaluation data frame.
#' @param palette    Two-colour gradient: \code{c(low, high)}.
#' @param base_size  Base font size.
#' @param save_plot  Logical.
#' @param save_dir   Output directory.
#' @param width,height Inches.
#' @param format     File format.
#' @returns A \code{ggplot} object.
#' @export
PlotConfusionMatrix <- function(object,

                                model_name = NULL,
                                test_data  = NULL,
                                palette    = c("#EEF3FA", "#1F5FA6"),
                                base_size  = 13,
                                save_plot  = FALSE,
                                save_dir = NULL,  # default: auto-detect from config
                                width      = 5,
                                height     = 4.5,
                                format     = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Model")

  if (!inherits(object, "Train_Model"))
    stop("'object' must be a Train_Model.")

  nm    <- if (is.null(model_name)) names(object@train.models)[1] else model_name
  model <- object@train.models[[nm]]
  td    <- if (!is.null(test_data)) test_data else object@split.data$test
  if (is.null(td)) stop("Provide test_data or populate split.data$test.")

  group_col <- as.character(object@group_col)
  truth     <- factor(td[[group_col]])
  pred      <- stats::predict(model, newdata = td, type = "raw")

  cm      <- table(Predicted = factor(pred, levels(truth)),
                   Actual    = truth)
  cm_df   <- as.data.frame(cm)
  cm_df   <- cm_df %>%
    dplyr::group_by(Actual) %>%
    dplyr::mutate(Pct = round(Freq / sum(Freq) * 100, 1)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(cm_df,
    ggplot2::aes(Actual, Predicted, fill = Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Freq, "\n(", Pct, "%)")),
                       size = 4, fontface = "bold") +
    ggplot2::scale_fill_gradient(low = palette[1], high = palette[2]) +
    ggplot2::labs(title = paste("Confusion Matrix —", nm),
                  x = "Actual", y = "Predicted", fill = "Count") +
    .pub_theme(base_size)

  if (save_plot)
    .save_plot(p, save_dir, paste0("confusion_", nm), width, height, format)
  return(p)
}


#' Feature importance plot
#'
#' Bar chart of variable importance scores across trained models. Each model
#' produces one facet panel.
#'
#' @param object      A \code{Train_Model} object.
#' @param top_n       Show top N features per model. Default \code{20}.
#' @param palette_name Colour palette.
#' @param base_size   Base font size.
#' @param save_plot   Logical.
#' @param save_dir    Output directory.
#' @param width,height Inches.
#' @param format      File format.
#' @returns A \code{ggplot} object.
#' @export
PlotFeatureImportance <- function(object,

                                  top_n        = 20,
                                  palette_name = "Zissou1",
                                  base_size    = 12,
                                  save_plot    = FALSE,
                                  save_dir = NULL,  # default: auto-detect from config
                                  width        = 10,
                                  height       = 6,
                                  format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Model")

  if (!inherits(object, "Train_Model"))
    stop("'object' must be a Train_Model.")
  if (length(object@train.models) == 0)
    stop("No trained models found.")

  imp_list <- lapply(names(object@train.models), function(nm) {
    model <- object@train.models[[nm]]
    imp   <- tryCatch(
      caret::varImp(model, scale = TRUE)$importance,
      error = function(e) NULL
    )
    if (is.null(imp)) return(NULL)
    imp_df <- data.frame(Feature    = rownames(imp),
                         Importance = imp[, 1],
                         Model      = nm,
                         row.names  = NULL)
    head(imp_df[order(-imp_df$Importance), ], top_n)
  })

  imp_df <- dplyr::bind_rows(imp_list[!sapply(imp_list, is.null)])
  if (nrow(imp_df) == 0) stop("Could not extract variable importance.")

  imp_df <- imp_df %>%
    dplyr::group_by(Model) %>%
    dplyr::mutate(Feature = factor(Feature, levels = rev(unique(Feature)))) %>%
    dplyr::ungroup()

  cols <- .get_palette(palette_name, length(unique(imp_df$Model)))

  p <- ggplot2::ggplot(imp_df,
    ggplot2::aes(Importance, Feature, fill = Model)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::facet_wrap(~ Model, scales = "free") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = paste0("Top ", top_n, " Feature Importance"),
                  x = "Importance (scaled)", y = NULL) +
    .pub_theme(base_size)

  if (save_plot)
    .save_plot(p, save_dir, "feature_importance", width, height, format)
  return(p)
}


#' Calibration curve
#'
#' Plots predicted probability vs. observed event rate (Platt-style), with
#' a perfectly calibrated reference line and Hosmer–Lemeshow test annotation.
#'
#' @param object      A \code{Train_Model} object.
#' @param model_name  Model to use.
#' @param test_data   Test data frame.
#' @param n_bins      Number of probability bins. Default \code{10}.
#' @param palette_name Palette.
#' @param base_size   Base font size.
#' @param save_plot   Logical.
#' @param save_dir    Output directory.
#' @param width,height Inches.
#' @param format      File format.
#' @returns A \code{ggplot} object.
#' @export
PlotCalibration <- function(object,

                            model_name   = NULL,
                            test_data    = NULL,
                            n_bins       = 10,
                            palette_name = "Royal1",
                            base_size    = 13,
                            save_plot    = FALSE,
                            save_dir = NULL,  # default: auto-detect from config
                            width        = 6,
                            height       = 5.5,
                            format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Model")


  if (!inherits(object, "Train_Model"))
    stop("'object' must be a Train_Model.")

  nm    <- if (is.null(model_name)) names(object@train.models)[1] else model_name
  model <- object@train.models[[nm]]
  td    <- if (!is.null(test_data)) test_data else object@split.data$test

  group_col <- as.character(object@group_col)
  truth     <- as.integer(td[[group_col]]) - 1L
  pos_level <- levels(factor(td[[group_col]]))[2]
  probs     <- stats::predict(model, newdata = td, type = "prob")[, pos_level]

  calib_df <- data.frame(truth = truth, prob = probs)
  calib_df$bin <- cut(calib_df$prob, breaks = seq(0, 1, length.out = n_bins + 1),
                      include.lowest = TRUE)
  calib_sum <- calib_df %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(mean_pred = mean(prob),
                     obs_rate  = mean(truth),
                     n         = dplyr::n(), .groups = "drop")

  col <- .get_palette(palette_name, 2)[1]

  p <- ggplot2::ggplot(calib_sum,
    ggplot2::aes(mean_pred, obs_rate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey40", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = col, alpha = 0.85) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, colour = col,
                         linewidth = 0.9, fill = paste0(col, "40")) +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "n") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title    = paste("Calibration Curve —", nm),
                  subtitle = paste0("n_bins = ", n_bins),
                  x = "Mean Predicted Probability",
                  y = "Observed Event Rate") +
    .pub_theme(base_size)

  if (save_plot)
    .save_plot(p, save_dir, paste0("calibration_", nm), width, height, format)
  return(p)
}


## ═════════════════════════════════════════════════════════════════════════════
##  § 3  SUBTYPING  ── dim-reduction · cluster heatmap · silhouette · alluvial
## ═════════════════════════════════════════════════════════════════════════════

#' t-SNE / UMAP scatter coloured by cluster and features
#'
#' Draws 2-D embedding points coloured by cluster label plus optional
#' metadata overlay. Returns a \code{patchwork} composite when
#' \code{overlay_features} is supplied.
#'
#' @param object           A \code{Subtyping} object.
#' @param reduction        \code{"tsne"} or \code{"umap"}.
#' @param color_by         \code{"cluster"} or a column name in
#'   \code{info.data}.
#' @param overlay_features Character vector of features to colour additional
#'   panels.
#' @param point_size       Point size. Default \code{1.8}.
#' @param palette_name     Palette.
#' @param base_size        Base font size.
#' @param save_plot        Logical.
#' @param save_dir         Output directory.
#' @param width,height     Inches.
#' @param format           File format.
#' @returns A \code{ggplot} or \code{patchwork} object.
#' @export
PlotDimReduction <- function(object,
                             reduction        = c("tsne", "umap"),
                             color_by         = "cluster",
                             overlay_features = NULL,
                             point_size       = 1.8,
                             palette_name     = "Darjeeling1",
                             base_size        = 13,
                             save_plot        = FALSE,
                             save_dir = NULL,  # default: auto-detect from config
                             width            = 6,
                             height           = 5,
                             format           = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Subtyping")

  if (!inherits(object, "Subtyping")) stop("'object' must be a Subtyping object.")
  reduction <- match.arg(reduction)

  vis <- object@visualization.results
  key <- paste0(reduction, ".df")
  if (is.null(vis) || (!key %in% names(vis) && !reduction %in% names(vis)))
    stop("Run Sub_tsne_analyse() or Sub_umap_analyse() first.")
  
  if (key %in% names(vis)) {
    coords_df <- vis[[key]]
  } else {
    coords_df <- vis[[reduction]]$coords
    if (is.null(coords_df))
      coords_df <- vis[[reduction]][, 1:2]
  }
  
  if(ncol(coords_df) >= 2) {
    colnames(coords_df)[1:2] <- c("Dim1", "Dim2")
  }

  cdata <- object@clustered.data
  if (!is.null(cdata) && nrow(cdata) > 0) {
    if ("cluster" %in% colnames(cdata)) {
      coords_df$cluster <- as.factor(cdata[rownames(coords_df), "cluster"])
    } else if ("group" %in% colnames(cdata)) {
      coords_df$cluster <- as.factor(cdata[rownames(coords_df), "group"])
    }
  }
  info <- object@info.data
  if (nrow(info) > 0)
    coords_df <- cbind(coords_df, info[rownames(coords_df), , drop = FALSE])

  ax <- toupper(reduction)
  .one_panel <- function(clr_col, is_factor = TRUE, title_extra = "") {
    if (is_factor) {
      n_lv <- length(unique(coords_df[[clr_col]]))
      cols  <- .get_palette(palette_name, n_lv)
      p <- ggplot2::ggplot(coords_df,
        ggplot2::aes(Dim1, Dim2, colour = as.factor(.data[[clr_col]]))) +
        ggplot2::scale_colour_manual(values = cols,
                                     name   = clr_col)
    } else {
      p <- ggplot2::ggplot(coords_df,
        ggplot2::aes(Dim1, Dim2,
                     colour = as.numeric(.data[[clr_col]]))) +
        ggplot2::scale_colour_viridis_c(name = clr_col, option = "D")
    }
    p + ggplot2::geom_point(size = point_size, alpha = 0.8) +
      ggplot2::labs(title = paste0(ax, title_extra),
                    x = paste(ax, "1"), y = paste(ax, "2")) +
      .pub_theme(base_size)
  }

  main_p <- .one_panel(color_by,
                        is_factor = !color_by %in% names(coords_df)[
                          sapply(coords_df, is.numeric)])

  if (!is.null(overlay_features) && length(overlay_features) > 0) {
    feat_panels <- lapply(overlay_features, function(ft) {
      if (!ft %in% colnames(coords_df)) {
        cd <- object@clean.data
        if (ft %in% colnames(cd))
          coords_df[[ft]] <<- cd[rownames(coords_df), ft]
        else return(NULL)
      }
      .one_panel(ft, is_factor = FALSE,
                  title_extra = paste0(" – ", ft))
    })
    feat_panels <- feat_panels[!sapply(feat_panels, is.null)]
    all_panels  <- c(list(main_p), feat_panels)
    out_p       <- patchwork::wrap_plots(all_panels)
  } else {
    out_p <- main_p
  }

  if (save_plot)
    .save_plot(out_p, save_dir, paste0(reduction, "_", color_by),
               width * max(1, 1 + length(overlay_features)), height, format)
  return(out_p)
}


#' Cluster-level heatmap (top marker features per cluster)
#'
#' Averages feature expression per cluster and renders a heatmap, optionally
#' with a side bar showing per-cluster sample counts.
#'
#' @param object       A \code{Subtyping} object with \code{clustered.data}.
#' @param top_n        Top N features per cluster ranked by mean difference.
#'   Default \code{10}.
#' @param scale_rows   Z-score rows. Default \code{TRUE}.
#' @param palette_name RColorBrewer palette.
#' @param base_fontsize Font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{pheatmap} object (invisibly).
#' @export
PlotClusterHeatmap <- function(object,

                               top_n         = 10,
                               scale_rows    = TRUE,
                               palette_name  = "RdYlBu",
                               base_fontsize = 10,
                               save_plot     = FALSE,
                               save_dir = NULL,  # default: auto-detect from config
                               width         = 9,
                               height        = 8,
                               format        = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Subtyping")


  if (!inherits(object, "Subtyping"))
    stop("'object' must be a Subtyping object.")
  cd <- object@clustered.data
  if (is.null(cd) || nrow(cd) == 0 || !"cluster" %in% colnames(cd))
    stop("Run clustering first (clustered.data empty or missing 'cluster' column).")

  num_feats <- setdiff(names(cd)[sapply(cd, is.numeric)], "cluster")
  cd_num    <- cd[, num_feats, drop = FALSE]
  cl_vec    <- as.factor(cd[["cluster"]])

  # select top_n per cluster by mean expression rank
  cl_means  <- aggregate(cd_num, by = list(cluster = cl_vec), FUN = mean)
  row.names(cl_means) <- cl_means$cluster
  cl_means  <- cl_means[, -1, drop = FALSE]

  top_feats <- unique(unlist(lapply(rownames(cl_means), function(cl) {
    ord <- order(-cl_means[cl, ], decreasing = FALSE)
    colnames(cl_means)[head(ord, top_n)]
  })))

  mat <- t(cd_num[, top_feats, drop = FALSE])
  if (scale_rows) {
    mat <- t(scale(t(mat)))
    mat[is.nan(mat)] <- 0
  }

  ann_col <- data.frame(Cluster = cl_vec, row.names = rownames(cd))
  n_cl    <- nlevels(cl_vec)
  cl_cols <- .get_palette("Darjeeling1", n_cl)
  names(cl_cols) <- levels(cl_vec)

  color_ramp <- rev(RColorBrewer::brewer.pal(11, palette_name))

  ph <- pheatmap::pheatmap(
    mat,
    color          = colorRampPalette(color_ramp)(100),
    annotation_col = ann_col,
    annotation_colors = list(Cluster = cl_cols),
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_colnames  = FALSE,
    fontsize       = base_fontsize,
    border_color   = NA,
    silent         = TRUE
  )

  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, paste0("cluster_heatmap.", format))
    # pheatmap object needs special handling - use grid.draw
    if (format == "pdf") {
      pdf(path, width = width, height = height)
      grid::grid.draw(ph$gtable)
      dev.off()
    } else {
      png(path, width = width * 100, height = height * 100, res = 100)
      grid::grid.draw(ph$gtable)
      dev.off()
    }
    cat("Saved:", path, "\n")
  }
  return(ph)
}


#' Silhouette plot
#'
#' Draws a sorted silhouette width bar chart, with per-cluster average lines
#' and an overall average annotation.
#'
#' @param object       A \code{Subtyping} object with \code{clustered.data}.
#' @param dist_method  Distance for silhouette calculation. Default
#'   \code{"euclidean"}.
#' @param palette_name Palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotSilhouette <- function(object,

                           dist_method  = "euclidean",
                           palette_name = "Darjeeling1",
                           base_size    = 13,
                           save_plot    = FALSE,
                           save_dir = NULL,  # default: auto-detect from config
                           width        = 7,
                           height       = 5,
                           format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Subtyping")


  if (!inherits(object, "Subtyping"))
    stop("'object' must be a Subtyping object.")
  cd <- object@clustered.data
  if (is.null(cd) || !"cluster" %in% colnames(cd))
    stop("Run clustering first.")

  num_mat <- as.matrix(cd[, setdiff(names(cd)[sapply(cd, is.numeric)], "cluster")])
  cl_int  <- as.integer(cd[["cluster"]])
  d_mat   <- stats::dist(num_mat, method = dist_method)
  sil     <- cluster::silhouette(cl_int, d_mat)
  sil_df  <- as.data.frame(sil[, ])
  sil_df  <- sil_df[order(sil_df$cluster, sil_df$sil_width), ]
  sil_df$order   <- seq_len(nrow(sil_df))
  sil_df$cluster <- as.factor(sil_df$cluster)
  avg_sil <- round(mean(sil_df$sil_width), 3)

  cols <- .get_palette(palette_name, nlevels(sil_df$cluster))

  p <- ggplot2::ggplot(sil_df,
    ggplot2::aes(order, sil_width, fill = cluster)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::geom_hline(yintercept = avg_sil,
                        linetype = "dashed", colour = "grey20") +
    ggplot2::annotate("text", x = nrow(sil_df) * 0.02,
                      y = avg_sil + 0.03,
                      label = paste0("Avg = ", avg_sil),
                      hjust = 0, size = 3.5) +
    ggplot2::scale_fill_manual(values = cols, name = "Cluster") +
    ggplot2::labs(title = "Silhouette Plot",
                  x = "Samples (sorted by cluster)", y = "Silhouette Width") +
    .pub_theme(base_size)

  if (save_plot)
    .save_plot(p, save_dir, "silhouette", width, height, format)
  return(p)
}


#' Alluvial / Sankey plot — cluster method comparison
#'
#' Compares cluster assignments from two methods (or time-points) as an
#' alluvial diagram, showing how samples flow between categories.
#'
#' @param object    A \code{Subtyping} object, or a plain data frame with
#'   two cluster columns.
#' @param col_from  Column name for the left axis.
#' @param col_to    Column name for the right axis.
#' @param palette_name Palette for \code{col_from} strata.
#' @param base_size Base font size.
#' @param save_plot Logical.
#' @param save_dir  Output directory.
#' @param width,height Inches.
#' @param format    File format.
#' @returns A \code{ggplot} object.
#' @export
PlotAlluvial <- function(object,

                         col_from     = "cluster",
                         col_to       = "cluster2",
                         palette_name = "Darjeeling1",
                         base_size    = 13,
                         save_plot    = FALSE,
                         save_dir = NULL,  # default: auto-detect from config
                         width        = 7,
                         height       = 6,
                         format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Subtyping")

  df <- if (inherits(object, "Subtyping")) {
    cd <- object@clustered.data
    if (col_to %in% colnames(object@info.data))
      cbind(cd, object@info.data[rownames(cd), col_to, drop = FALSE])
    else cd
  } else {
    as.data.frame(object)
  }

  if (!all(c(col_from, col_to) %in% colnames(df)))
    stop("Both '", col_from, "' and '", col_to, "' must be columns in the data.")

  alluvial_df <- df %>%
    dplyr::count(.data[[col_from]], .data[[col_to]], name = "freq") %>%
    dplyr::mutate(across(c(col_from, col_to), as.factor))

  cols <- .get_palette(palette_name, nlevels(alluvial_df[[col_from]]))

  p <- ggplot2::ggplot(alluvial_df,
    ggplot2::aes(axis1 = .data[[col_from]],
                 axis2 = .data[[col_to]],
                 y     = freq)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data[[col_from]]),
                               width = 1/12, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 1/12, fill = "grey80", colour = "white") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = 3.5, fontface = "bold") +
    ggplot2::scale_x_discrete(limits = c(col_from, col_to), expand = c(0.1, 0)) +
    ggplot2::scale_fill_manual(values = cols, name = col_from) +
    ggplot2::labs(title = "Cluster Assignment Alluvial",
                  y = "Sample Count") +
    .pub_theme(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", size = 12))

  if (save_plot)
    .save_plot(p, save_dir, "alluvial", width, height, format)
  return(p)
}


## ═════════════════════════════════════════════════════════════════════════════
##  § 4  PROGNOSIS  ── KM · forest · time-ROC · RCS · nomogram · calibration · DCA · risk
## ═════════════════════════════════════════════════════════════════════════════

#' Kaplan-Meier survival curve
#'
#' Draws KM curves for up to 4 groups with risk table, median survival lines,
#' log-rank p-value, and optional pairwise p-value annotations.
#'
#' @param object       A \code{PrognosiX} object or a survival data frame with
#'   \code{time} and \code{status} columns.
#' @param group_col    Column to stratify on. When \code{NULL} the entire
#'   cohort is shown as one curve.
#' @param time_unit    Label for the x-axis: \code{"days"}, \code{"months"},
#'   \code{"years"}.
#' @param conf_int     Draw confidence bands. Default \code{FALSE}.
#' @param pairwise_p   Add pairwise p-value table. Default \code{FALSE}.
#' @param palette_name Palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggsurvplot} object (list).
#' @export
PlotKaplanMeier <- function(object,

                            group_col    = NULL,
                            time_unit    = "months",
                            conf_int     = FALSE,
                            pairwise_p   = FALSE,
                            palette_name = "Darjeeling1",
                            base_size    = 13,
                            save_plot    = FALSE,
                            save_dir = NULL,  # default: auto-detect from config
                            width        = 7,
                            height       = 8,
                            format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")

  if (inherits(object, "PrognosiX")) {
    surv_df <- object@survival.data
  } else {
    surv_df <- as.data.frame(object)
  }

  form_str <- if (is.null(group_col)) "Surv(time, status) ~ 1"
              else paste0("Surv(time, status) ~ ", group_col)
  fit  <- survival::survfit(as.formula(form_str), data = surv_df)

  n_grp <- if (is.null(group_col)) 1L
            else length(unique(surv_df[[group_col]]))
  cols  <- .get_palette(palette_name, n_grp)

  km <- survminer::ggsurvplot(
    fit,
    data             = surv_df,
    palette          = cols,
    conf.int         = conf_int,
    pval             = !is.null(group_col),
    pval.method      = !is.null(group_col),
    risk.table       = TRUE,
    risk.table.col   = "strata",
    surv.median.line = "hv",
    xlab             = paste0("Time (", time_unit, ")"),
    ylab             = "Overall Survival Probability",
    legend.title     = if (!is.null(group_col)) group_col else "",
    ggtheme          = .pub_theme(base_size),
    tables.theme     = ggplot2::theme_minimal(base_size = base_size - 2) +
                       ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold"))
  )

  if (pairwise_p && !is.null(group_col)) {
    pw <- survminer::pairwise_survdiff(as.formula(form_str), data = surv_df)
    km$plot <- km$plot +
      ggplot2::labs(caption = paste(
        "Pairwise log-rank p:\n",
        paste(capture.output(print(pw$p.value, digits = 3)),
              collapse = "\n")))
  }

  if (save_plot) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir,
                      paste0("KM_", if (!is.null(group_col)) group_col else "overall",
                             ".", format))
    survminer::ggsave(path, km, width = width, height = height, dpi = 300)
    cat("KM plot saved:", path, "\n")
  }
  return(km)
}


#' Forest plot — univariate / multivariate Cox HR
#'
#' Renders a publication-style forest plot with HR, 95 \% CI, and p-value
#' columns alongside the graphical HR display.
#'
#' @param object       A \code{PrognosiX} object or a data frame with columns
#'   \code{Variable}, \code{HR}, \code{CI_lower}, \code{CI_upper},
#'   \code{P_value}, \code{HR_95CI}.
#' @param analysis     \code{"univariate"} (default) or \code{"multivariate"}.
#'   Selects the result stored in \code{univariate.analysis}.
#' @param hr_limit     x-axis range for HR. Default \code{c(0.1, 10)}.
#' @param log_scale    Use log10 x-axis. Default \code{TRUE}.
#' @param palette_name Palette (2 colours: sig / NS).
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotForestPlot <- function(object,

                           analysis     = "univariate",
                           hr_limit     = c(0.1, 10),
                           log_scale    = TRUE,
                           palette_name = "Royal1",
                           base_size    = 12,
                           save_plot    = FALSE,
                           save_dir = NULL,  # default: auto-detect from config
                           width        = 10,
                           height       = 7,
                           format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")

  if (inherits(object, "PrognosiX")) {
    hr_df <- object@univariate.analysis$hr_results
    if (is.null(hr_df) && analysis == "multivariate")
      hr_df <- object@best.model$hr_results
  } else {
    hr_df <- as.data.frame(object)
  }

  req_cols <- c("Variable", "HR", "CI_lower", "CI_upper", "P_value")
  if (!all(req_cols %in% colnames(hr_df)))
    stop("hr_df must contain: ", paste(req_cols, collapse = ", "))

  hr_df <- hr_df %>%
    dplyr::filter(is.finite(HR) & is.finite(CI_lower) & is.finite(CI_upper)) %>%
    dplyr::filter(HR >= hr_limit[1] & HR <= hr_limit[2]) %>%
    dplyr::mutate(
      Sig      = P_value < 0.05,
      P_label  = ifelse(P_value < 0.001, "<0.001",
                        as.character(round(P_value, 3))),
      HR_label = if ("HR_95CI" %in% colnames(.)) HR_95CI
                 else paste0(round(HR, 2), " (",
                             round(CI_lower, 2), "–",
                             round(CI_upper, 2), ")"),
      Variable = factor(Variable, levels = rev(Variable))
    )

  cols <- .get_palette(palette_name, 2)

  p <- ggplot2::ggplot(hr_df,
    ggplot2::aes(y = Variable, x = HR,
                 xmin = CI_lower, xmax = CI_upper,
                 colour = Sig)) +
    ggplot2::geom_vline(xintercept = 1,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_errorbarh(height = 0.25, linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "grey60", `TRUE` = cols[1]),
      labels = c(`FALSE` = "p ≥ 0.05", `TRUE` = "p < 0.05"),
      name   = NULL) +
    ggplot2::geom_text(ggplot2::aes(x = hr_limit[2] * 1.1, label = HR_label),
                       hjust = 0, size = 3, colour = "black") +
    ggplot2::geom_text(ggplot2::aes(x = hr_limit[2] * 1.5, label = P_label),
                       hjust = 0, size = 3, colour = "black") +
    {if (log_scale) ggplot2::scale_x_log10(limits = c(hr_limit[1],
                                                        hr_limit[2] * 2))
     else ggplot2::scale_x_continuous(limits = c(hr_limit[1],
                                                   hr_limit[2] * 2))} +
    ggplot2::labs(title = paste(tools::toTitleCase(analysis), "Cox Forest Plot"),
                  x = "Hazard Ratio (95% CI)", y = NULL) +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = "top")

  if (save_plot)
    .save_plot(p, save_dir, paste0("forest_", analysis), width, height, format)
  return(p)
}


#' Time-dependent ROC curves
#'
#' Draws AUC-annotated ROC curves for multiple time-points using
#' \pkg{timeROC}.
#'
#' @param object       A \code{PrognosiX} object or a survival data frame.
#' @param marker_col   Column to use as the marker / risk score.
#' @param time_points  Numeric vector of time-points. When \code{NULL} the
#'   25th, 50th, 75th percentiles of \code{time} are used.
#' @param time_unit    Unit label. Default \code{"months"}.
#' @param palette_name Palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotTimeROC <- function(object,

                        marker_col   = "risk_score",
                        time_points  = NULL,
                        time_unit    = "months",
                        palette_name = "Darjeeling1",
                        base_size    = 13,
                        save_plot    = FALSE,
                        save_dir = NULL,  # default: auto-detect from config
                        width        = 6,
                        height       = 5.5,
                        format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")


  surv_df <- if (inherits(object, "PrognosiX")) object@survival.data
              else as.data.frame(object)
  if (!marker_col %in% colnames(surv_df))
    stop("Column '", marker_col, "' not found.")

  time   <- surv_df[["time"]]
  status <- as.numeric(surv_df[["status"]])
  marker <- surv_df[[marker_col]]

  tps <- if (is.null(time_points)) {
    round(stats::quantile(time, c(0.25, 0.5, 0.75)), 1)
  } else {
    time_points
  }

  roc_res <- timeROC::timeROC(T = time, delta = status,
                               marker = marker, cause = 1,
                               weighting = "marginal", times = tps)

  auc_vals <- round(roc_res$AUC, 3)
  cols     <- .get_palette(palette_name, length(tps))

  roc_list <- lapply(seq_along(tps), function(i) {
    fp <- roc_res$FP[, i]
    tp <- roc_res$TP[, i]
    data.frame(FPR       = fp,
               TPR       = tp,
               TimePoint = paste0("t=", tps[i], " ", time_unit,
                                  " (AUC=", auc_vals[i], ")"))
  })
  roc_df <- dplyr::bind_rows(roc_list)

  p <- ggplot2::ggplot(roc_df,
    ggplot2::aes(FPR, TPR, colour = TimePoint)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::coord_equal() +
    ggplot2::labs(title    = "Time-Dependent ROC",
                  subtitle = paste("Marker:", marker_col),
                  x = "False Positive Rate",
                  y = "True Positive Rate") +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = c(0.68, 0.22))

  if (save_plot)
    .save_plot(p, save_dir, "time_ROC", width, height, format)
  return(p)
}


#' Restricted Cubic Spline (RCS) dose-response plot
#'
#' Fits an RCS model (linear, logistic, or Cox) and plots the non-linear
#' association with 95 \% CI ribbon.
#'
#' @param object       A \code{PrognosiX} object or data frame.
#' @param x_col        Continuous predictor column.
#' @param y_col        Outcome column (event status for Cox).
#' @param method       \code{"linear"}, \code{"logistic"}, or \code{"cox"}.
#' @param knots        Number of RCS knots. Default \code{4}.
#' @param adjust_vars  Additional adjustment covariates.
#' @param ref_point    Reference value for OR/HR = 1. Default: median.
#' @param palette_name Single colour palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotRCS <- function(object,

                    x_col        = NULL,
                    y_col        = "status",
                    method       = c("cox", "logistic", "linear"),
                    knots        = 4,
                    adjust_vars  = NULL,
                    ref_point    = NULL,
                    palette_name = "Royal1",
                    base_size    = 13,
                    save_plot    = FALSE,
                    save_dir = NULL,  # default: auto-detect from config
                    width        = 6,
                    height       = 5,
                    format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")

  method  <- match.arg(method)
  surv_df <- if (inherits(object, "PrognosiX")) object@survival.data
              else as.data.frame(object)
  if (is.null(x_col)) stop("Specify 'x_col'.")

  ref <- if (is.null(ref_point)) stats::median(surv_df[[x_col]], na.rm = TRUE)
          else ref_point

  dd <- rms::datadist(surv_df)
  options(datadist = "dd")

  fit <- if (method == "cox") {
    rms::cph(as.formula(paste0("Surv(time, status) ~ rms::rcs(", x_col, ",", knots, ")",
                               if (!is.null(adjust_vars)) paste0("+", paste(adjust_vars, collapse = "+")))),
              data = surv_df, x = TRUE, y = TRUE)
  } else if (method == "logistic") {
    rms::lrm(as.formula(paste0(y_col, " ~ rms::rcs(", x_col, ",", knots, ")",
                               if (!is.null(adjust_vars)) paste0("+", paste(adjust_vars, collapse = "+")))),
              data = surv_df, x = TRUE, y = TRUE)
  } else {
    rms::ols(as.formula(paste0(y_col, " ~ rms::rcs(", x_col, ",", knots, ")",
                               if (!is.null(adjust_vars)) paste0("+", paste(adjust_vars, collapse = "+")))),
              data = surv_df, x = TRUE, y = TRUE)
  }

  pred <- rms::Predict(fit, !!rlang::sym(x_col), ref.zero = TRUE, fun = exp)
  pred_df <- as.data.frame(pred)

  col <- .get_palette(palette_name, 1)
  y_label <- switch(method, cox = "Hazard Ratio", logistic = "Odds Ratio", "β")

  hist_df <- data.frame(x = surv_df[[x_col]])

  p <- ggplot2::ggplot(pred_df,
    ggplot2::aes(x = .data[[x_col]], y = yhat)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = paste0(col, "30"), colour = NA) +
    ggplot2::geom_line(colour = col, linewidth = 1.2) +
    ggplot2::geom_rug(data = hist_df, ggplot2::aes(x = x, y = NULL),
                      sides = "b", alpha = 0.3, colour = "grey30") +
    ggplot2::labs(title    = paste("RCS —", x_col, "vs", y_label),
                  subtitle = paste0(knots, " knots, ref = ", round(ref, 2)),
                  x = x_col, y = y_label) +
    .pub_theme(base_size)

  if (save_plot)
    .save_plot(p, save_dir, paste0("RCS_", x_col), width, height, format)
  return(p)
}


#' Risk score distribution plot
#'
#' Draws a dot plot of per-sample risk scores sorted by value, coloured by
#' high / low risk group, with a KM curve and time-ROC in a composite layout.
#'
#' @param object       A \code{PrognosiX} object or a data frame with
#'   \code{risk_score}, \code{time}, \code{status}.
#' @param score_col    Column name for the risk score. Default
#'   \code{"risk_score"}.
#' @param cutoff_method Method to determine high / low boundary:
#'   \code{"median"} (default), \code{"mean"}, or a numeric value.
#' @param palette_name Palette (2 colours: high / low).
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{patchwork} composite.
#' @export
PlotRiskScore <- function(object,

                          score_col      = "risk_score",
                          cutoff_method  = "median",
                          palette_name   = "Royal1",
                          base_size      = 12,
                          save_plot      = FALSE,
                          save_dir = NULL,  # default: auto-detect from config
                          width          = 10,
                          height         = 4,
                          format         = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")


  surv_df <- if (inherits(object, "PrognosiX")) object@survival.data
              else as.data.frame(object)
  if (!score_col %in% colnames(surv_df))
    stop("Column '", score_col, "' not found.")

  cutoff <- if (is.numeric(cutoff_method)) cutoff_method
             else if (cutoff_method == "mean") mean(surv_df[[score_col]], na.rm = TRUE)
             else stats::median(surv_df[[score_col]], na.rm = TRUE)

  surv_df <- surv_df[order(surv_df[[score_col]]), ]
  surv_df$Risk  <- ifelse(surv_df[[score_col]] >= cutoff, "High", "Low")
  surv_df$Index <- seq_len(nrow(surv_df))

  cols <- .get_palette(palette_name, 2)
  names(cols) <- c("Low", "High")

  # panel A – dot plot
  p_score <- ggplot2::ggplot(surv_df,
    ggplot2::aes(Index, .data[[score_col]], colour = Risk)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = which(surv_df$Risk == "High")[1] - 0.5,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(title = "Risk Score Distribution",
                  x = "Patients (ranked)", y = "Risk Score",
                  colour = "Risk Group") +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = "right")

  # panel B – survival status tiles
  p_status <- ggplot2::ggplot(surv_df,
    ggplot2::aes(Index, 1,
                 colour = as.factor(as.numeric(surv_df[["status"]])))) +
    ggplot2::geom_point(shape = 3, size = 0.6) +
    ggplot2::scale_colour_manual(values = c("0" = "grey70", "1" = "black"),
                                  labels = c("0" = "Censored", "1" = "Event"),
                                  name   = NULL) +
    ggplot2::labs(x = NULL, y = "Event") +
    ggplot2::scale_y_continuous(breaks = NULL) +
    .pub_theme(base_size - 2) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  out_p <- patchwork::wrap_plots(p_score, p_status, ncol = 1,
                                  heights = c(4, 1))

  if (save_plot)
    .save_plot(out_p, save_dir, "risk_score", width, height, format)
  return(out_p)
}


#' Decision Curve Analysis (DCA) plot
#'
#' Plots net benefit vs. threshold probability for one or more models
#' alongside the "treat all" and "treat none" reference lines.
#'
#' @param object       A \code{PrognosiX} object, or a data frame with outcome
#'   and predictor columns.
#' @param predictors   Character vector of predictor / risk score column names.
#' @param outcome_col  Binary outcome column. Default \code{"status"}.
#' @param thresholds   Sequence of decision thresholds. Default
#'   \code{seq(0, 0.5, by = 0.01)}.
#' @param palette_name Palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotDCA <- function(object,

                    predictors   = "risk_score",
                    outcome_col  = "status",
                    thresholds   = seq(0, 0.5, by = 0.01),
                    palette_name = "Darjeeling1",
                    base_size    = 13,
                    save_plot    = FALSE,
                    save_dir = NULL,  # default: auto-detect from config
                    width        = 7,
                    height       = 5.5,
                    format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")

  df <- if (inherits(object, "PrognosiX")) object@survival.data
         else as.data.frame(object)
  df[[outcome_col]] <- as.numeric(df[[outcome_col]])
  n <- nrow(df)
  prev <- mean(df[[outcome_col]], na.rm = TRUE)

  dca_list <- lapply(thresholds, function(pt) {
    treat_all_nb <- prev - (1 - prev) * pt / (1 - pt)
    rows <- lapply(predictors, function(pr) {
      if (!pr %in% colnames(df)) return(NULL)
      test_pos   <- df[[pr]] >= pt
      true_pos   <- sum(test_pos & df[[outcome_col]] == 1, na.rm = TRUE)
      false_pos  <- sum(test_pos & df[[outcome_col]] == 0, na.rm = TRUE)
      net_benefit <- (true_pos / n) - (false_pos / n) * (pt / (1 - pt))
      data.frame(Threshold = pt, Model = pr, NetBenefit = net_benefit)
    })
    ref_row <- data.frame(Threshold  = pt,
                          Model      = "Treat All",
                          NetBenefit = treat_all_nb)
    dplyr::bind_rows(c(rows, list(ref_row)))
  })

  dca_df <- dplyr::bind_rows(dca_list)
  dca_df <- dplyr::bind_rows(
    dca_df,
    data.frame(Threshold  = thresholds,
               Model      = "Treat None",
               NetBenefit = 0)
  )

  model_levels <- c(predictors, "Treat All", "Treat None")
  dca_df$Model <- factor(dca_df$Model, levels = model_levels)

  n_cols <- length(model_levels)
  cols   <- .get_palette(palette_name, max(3, length(predictors)))
  cols   <- c(cols[seq_along(predictors)], "grey40", "grey70")
  names(cols) <- model_levels

  ltys <- c(rep("solid", length(predictors)), "dashed", "dotted")
  names(ltys) <- model_levels

  p <- ggplot2::ggplot(dca_df,
    ggplot2::aes(Threshold, NetBenefit, colour = Model, linetype = Model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = cols)  +
    ggplot2::scale_linetype_manual(values = ltys) +
    ggplot2::coord_cartesian(ylim = c(-0.05, max(dca_df$NetBenefit, na.rm = TRUE) + 0.05)) +
    ggplot2::labs(title    = "Decision Curve Analysis",
                  x = "Threshold Probability",
                  y = "Net Benefit",
                  colour   = NULL,
                  linetype = NULL) +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = c(0.75, 0.75))

  if (save_plot)
    .save_plot(p, save_dir, "DCA", width, height, format)
  return(p)
}
