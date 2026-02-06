#' Survival Analysis for Single Variable
#'
#' Performs survival analysis and generates Kaplan-Meier plots for a specific variable.
#'
#' @param data Data frame.
#' @param time_col Time column name.
#' @param status_col Status column name.
#' @param var_col Variable to stratify by.
#' @param palette_name Palette name.
#' @param save_plots Logical.
#' @param save_dir Save directory.
#' @param plot_width Plot width.
#' @param plot_height Plot height.
#' @param base_size Base font size.
#'
#' @return A list containing the KM object and the combined plot.
#' @export
survival_analysis_var_plot <- function(data,
                                       time_col = "time",
                                       status_col = "status",
                                       var_col=NULL,
                                       palette_name = "Dark2",
                                       save_plots = TRUE,
                                       save_dir = here::here('PrognosiX', "univariate_analysis"),
                                       plot_width = 5,
                                       plot_height = 5,
                                       base_size = 14) {
  cat("Starting survival analysis...\n")

  if (!(var_col %in% names(data))) {
    stop(paste("Variable", var_col, "does not exist in the data frame"))
  }
  if (!(time_col %in% names(data))) {
    stop(paste("Variable", time_col, "does not exist in the data frame"))
  }
  if (!(status_col %in% names(data))) {
    stop(paste("Variable", status_col, "does not exist in the data frame"))
  }

  if (!is.numeric(data[[time_col]])) {
    data[[time_col]] <- suppressWarnings(as.numeric(data[[time_col]]))
  }
  status_vec <- data[[status_col]]
  if (is.factor(status_vec)) status_vec <- as.character(status_vec)
  if (is.logical(status_vec)) status_vec <- as.integer(status_vec)
  if (is.character(status_vec)) {
    status_vec_trim <- trimws(status_vec)
    if (all(na.omit(status_vec_trim) %in% c("0", "1"))) {
      status_vec <- as.numeric(status_vec_trim)
    } else {
      u <- sort(unique(na.omit(status_vec_trim)))
      if (length(u) == 2) {
        status_vec <- ifelse(status_vec_trim == u[1], 0, 1)
      } else {
        stop("Status column is not binary.")
      }
    }
  }
  if (is.numeric(status_vec)) {
    u <- sort(unique(na.omit(status_vec)))
    if (length(u) == 2 && all(u %in% c(1, 2))) {
      status_vec <- ifelse(status_vec == 1, 0, 1)
    } else if (length(u) == 2 && !all(u %in% c(0, 1))) {
       # Try to map min to 0, max to 1
       status_vec <- ifelse(status_vec == min(u), 0, 1)
    }
  }
  data[[status_col]] <- status_vec

  data <- data %>% dplyr::filter(!is.na(data[[var_col]]) & !is.na(data[[status_col]]) & !is.na(data[[time_col]]))
  if (nrow(data) == 0) {
    stop("Filtered data has no valid rows.")
  }

  # Ensure passed data is a pure data frame to avoid random class issues
  data <- as.data.frame(data)

  cat("Data validation successful. Number of valid rows: ", nrow(data),"\n")

  cat("Fitting survival model...\n")
  # Use non-evaluation formula construction for robustness
  f <- as.formula(paste("survival::Surv(", time_col, ", ", status_col, ") ~ ", var_col))
  fit <- survival::survfit(f, data = data)
  # Inject formula into call to avoid symbol lookup issues in ggsurvplot
  fit$call$formula <- f
  
  legend.title <- paste("Risk Group", var_col)
  pal <- wes_palette(palette_name, n = max(2, length(unique(data[[var_col]]))), type = "discrete")
  
  # Debug prints
  cat("DEBUG: Calling ggsurvplot\n")
  print(class(fit))
  print(fit$call)
  
  km_var <- survminer::ggsurvplot(
    fit,
    data = data,
    conf.int = FALSE,
    pval = TRUE,
    pval.method = TRUE,
    title = sprintf("Kaplan-Meier Survival Curve by %s", var_col),
    surv.median.line = "hv",
    risk.table = TRUE,
    xlab = "Follow up time (days)",
    legend = c(0.8, 0.75),
    legend.title = legend.title,
    break.x.by = 100,
    palette = pal
  )

  cat("Survival curve plotted.\n")

  surv_plot <- km_var$plot +
    ggplot2::theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = c(0.8, 0.8)
    )

  risk_table <- km_var$table +
    ggplot2::theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )

  combined_plot <- gridExtra::grid.arrange(surv_plot, risk_table, ncol = 1, heights = c(1.8, 1))

  cat("Displaying combined plot...\n")
  print(combined_plot)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat(paste("Created directory:", save_dir,"\n"))
    }

    ggsave(file.path(save_dir, paste0("survival_time_distribution_by_", var_col,".pdf")),
           plot = combined_plot, width = plot_width, height = plot_height,device = "pdf")

    ggsave(file.path(save_dir, paste0("survival_curve_by_", var_col,".pdf")),
           plot = surv_plot, width = plot_width, height = plot_height,device = "pdf")
    ggsave(file.path(save_dir, paste0("risk_table_by_", var_col,".pdf")),
           plot = risk_table, width = plot_width, height = plot_height,device = "pdf")
    
    cat("Plot saved to:", save_dir, "\n")
  }

  cat("Survival analysis function execution completed.\n")
  cat("plot saved successfully at:", save_dir, "\n")
  return(list(km_var = km_var, combined_plot = combined_plot))
}


#' Plot Kaplan-Meier Curve for PrognosiX Object
#'
#' Wrapper to plot KM curve for a specific variable using PrognosiX object or data frame.
#'
#' @param object PrognosiX object or data frame.
#' @param time_col Time column name.
#' @param status_col Status column name.
#' @param var_col Variable column name.
#' @param palette_name Palette name.
#' @param save_plots Logical.
#' @param save_dir Save directory.
#' @param plot_width Plot width.
#' @param plot_height Plot height.
#' @param base_size Base font size.
#' @param use_subgroup_data Logical.
#'
#' @return Updated PrognosiX object or plot object.
#' @export
plot_var_kaplan_meier <- function(object,
                                  time_col = "time",
                                  status_col = "status",
                                  var_col = NULL,
                                  palette_name = "Dark2",
                                  save_plots = TRUE,
                                  save_dir = here::here('PrognosiX', "univariate_analysis"),
                                  plot_width = 5,
                                  plot_height = 5,
                                  base_size = 14,
                                  use_subgroup_data = FALSE) {


  if (inherits(object, 'PrognosiX')) {
    if (use_subgroup_data) {
      data <- methods::slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
    } else {
      data <- methods::slot(object, "survival.data")
      cat("Using original survival data...\n")
    }
    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")

    if (is.null(data) || nrow(data) == 0) {
      stop("The survival.data in the PrognosiX object is empty.")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  if (is.null(var_col) || !(var_col %in% colnames(data))) {
    stop("The specified var_col is either NULL or not found in the dataset.")
  }

  cat("Starting Kaplan-Meier plotting...\n")

  combined_plot <- survival_analysis_var_plot(data,
                                              time_col = time_col,
                                              status_col = status_col,
                                              var_col = var_col,
                                              palette_name = palette_name,
                                              plot_width = plot_width,
                                              plot_height = plot_height,
                                              save_plots = save_plots,
                                              save_dir = save_dir)
  result_col<- paste0("km_","results_",var_col )
  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object...\n")

    object@univariate.analysis[["km_results"]][[result_col]]<-combined_plot
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }

  cat("Kaplan-Meier plot completed.\n")
  return(combined_plot)
}

#' Plot Survival Time Distribution
#'
#' Plots the distribution of survival time for a given variable.
#'
#' @param data Data frame.
#' @param time_col Time column.
#' @param var_col Variable column.
#' @param palette_name Palette name.
#' @param binwidth Bin width for histogram.
#' @param save_plots Logical.
#' @param save_dir Save directory.
#' @param plot_width Plot width.
#' @param plot_height Plot height.
#' @param base_size Base font size.
#'
#' @return The ggplot object.
#' @export
plot_survival_time_distribution <- function(data,
                                            time_col,
                                            var_col,
                                            palette_name = "AsteroidCity1",
                                            binwidth = 10,
                                            save_plots = TRUE,
                                            save_dir = here::here('PrognosiX', "univariate_analysis"),
                                            plot_width = 5,
                                            plot_height = 5,
                                            base_size = 14) {
  data[[var_col]] <- as.factor(data[[var_col]])

  mean_df <- data %>%
    dplyr::group_by(!!rlang::sym(var_col)) %>%
    dplyr::summarise(mean_time = mean(!!rlang::sym(time_col), na.rm = TRUE),
              se_time = sd(!!rlang::sym(time_col), na.rm = TRUE) / sqrt(dplyr::n())) %>%
    dplyr::ungroup()

  label_df <- data %>%
    dplyr::group_by(!!rlang::sym(var_col)) %>%
    dplyr::summarise(label = paste0(round(mean(!!rlang::sym(time_col), na.rm = TRUE), 1), "±",
                             round(sd(!!rlang::sym(time_col), na.rm = TRUE) / sqrt(dplyr::n()), 1), " days")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(y = seq(12, 12 - 5 * (dplyr::n() - 1), length.out = dplyr::n()))

  colors <- wes_palette(palette_name, n = max(4, length(unique(data[[var_col]]))), type = "discrete")
  p <- ggplot2::ggplot(data, aes(x = !!rlang::sym(time_col), fill = !!rlang::sym(var_col))) +
    geom_histogram(binwidth = binwidth, color = "black", alpha = 0.7) +
    geom_vline(data = mean_df, aes(xintercept = mean_time), color = colors[4], linetype = "solid", size = 1) +
    geom_vline(data = mean_df, aes(xintercept = mean_time - se_time), color = colors[2], linetype = "dashed", size = 0.6) +
    geom_vline(data = mean_df, aes(xintercept = mean_time + se_time), color = colors[2], linetype = "dashed", size = 0.6) +
    scale_fill_manual(values = wes_palette(palette_name, n = max(2, length(unique(data[[var_col]]))), type = "discrete")) +
    labs(title = paste("Distribution of Survival Time by", var_col),
         x = "Survival Time (days)",
         y = "Frequency") +
    facet_grid(paste0(var_col, " ~ .")) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none",
      strip.background = element_rect(fill = "white"),
      panel.grid = element_line(linetype = "dotted", color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.y.left = element_line(color = "black", size = 0.5),
      axis.line.y.right = element_line(color = "black", size = 0.5))

  cat("Plotting survival time distribution...\n")
  print(p)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
    }
    ggsave(file.path(save_dir, paste0("survival_distribution_by_", var_col,".pdf")),
           plot = p, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", file.path(save_dir, paste0("survival_distribution_by_", var_col,".pdf"), "\n"))
    }

  return(p)
}



#' Wrapper for Plotting Survival Time Distribution
#'
#' Wrapper function to plot survival time distribution for PrognosiX object or data frame.
#'
#' @param object PrognosiX object or data frame.
#' @param time_col Time column.
#' @param status_col Status column.
#' @param var_col Variable column.
#' @param binwidth Bin width.
#' @param palette_name Palette name.
#' @param save_plots Logical.
#' @param save_dir Save directory.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param use_subgroup_data Logical.
#'
#' @return Updated PrognosiX object or plot object.
#' @export
plot_var_survival_time <- function(object,
                                   time_col = "time",
                                   status_col = "status",
                                   var_col=NULL,
                                   binwidth = 10,
                                   palette_name = "AsteroidCity1",
                                   save_plots = TRUE,
                                   save_dir = here::here('PrognosiX', "univariate_analysis"),
                                   plot_width = 5,
                                   plot_height = 5,
                                   base_size = 14,
                                   use_subgroup_data = FALSE) {

  if (inherits(object, 'PrognosiX')) {
    if (use_subgroup_data) {
      data <- methods::slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
    } else {
      data <- methods::slot(object, "survival.data")
      cat("Using original survival data...\n")
    }
    status_col <- methods::slot(object, "status_col")
    time_col <- methods::slot(object, "time_col")

    if (is.null(data) || nrow(data) == 0) {
      stop("The survival.data in the PrognosiX object is empty.")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }
  data <- convert_to_numeric(data)

  p <- plot_survival_time_distribution(data=data,
                                       time_col=time_col,
                                       var_col=var_col,
                                       palette_name=palette_name,
                                       binwidth=binwidth,
                                       save_plots=save_plots,
                                       save_dir=save_dir,
                                       plot_width=plot_width,
                                       plot_height=plot_height,
                                       base_size=base_size)
  result_col<- paste0("distribution_time_",var_col)
  if (inherits(object, 'PrognosiX')) {
    cat("Updating 'PrognosiX' object...\n")

    object@univariate.analysis[["distribution_time"]][[result_col]]<-  p
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }

  cat("Plotting function execution completed.\n")
  return(p)
}
