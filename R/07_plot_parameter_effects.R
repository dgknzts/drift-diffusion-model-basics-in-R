# R/utils/plot_parameter_effects.R

#' Create Parameter Effect Plots
#'
#' This function creates visualizations comparing DDM outputs
#' across different parameter values or conditions.
#'
#' @param data_list List of data frames from DDM simulations
#' @param param_name Character string naming the parameter being varied
#' @param param_values Numeric vector of parameter values
#' @param plot_type Character string: "comprehensive", "summary", or "diagnostic"
#' @param rt_xlim Numeric vector of length 2 for RT axis limits
#' @param color_palette Character string naming the color palette to use
#'
#' @return A ggplot object or patchwork composition
#' @export
plot_parameter_effects <- function(data_list,
                                   param_name,
                                   param_values,
                                   plot_type = "comprehensive",
                                   rt_xlim = c(0, 2),
                                   color_palette = "viridis") {

  # Combine data with parameter labels
  if (!is.list(data_list)) {
    stop("data_list must be a list of data frames")
  }

  # Add parameter value to each dataset
  for (i in seq_along(data_list)) {
    data_list[[i]]$param_value <- param_values[i]
    data_list[[i]]$param_label <- paste0(param_name, " = ", param_values[i])
  }

  combined_data <- dplyr::bind_rows(data_list)

  # Filter valid data
  valid_data <- combined_data %>%
    dplyr::filter(!is.na(choice), !is.na(rt), rt > 0)

  if (nrow(valid_data) == 0) {
    warning("No valid data after filtering")
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "No Valid Data"))
  }

  # Choose plot type
  if (plot_type == "comprehensive") {
    return(create_comprehensive_plot(valid_data, param_name, rt_xlim, color_palette))
  } else if (plot_type == "summary") {
    return(create_summary_plot(valid_data, param_name, color_palette))
  } else if (plot_type == "diagnostic") {
    return(create_diagnostic_plot(valid_data, param_name, color_palette))
  } else {
    stop("plot_type must be 'comprehensive', 'summary', or 'diagnostic'")
  }
}

#' Create Comprehensive Parameter Effect Plot
#' @keywords internal
create_comprehensive_plot <- function(data, param_name, rt_xlim, color_palette) {

  # 1. RT distributions
  p1 <- ggplot2::ggplot(data,
                        ggplot2::aes(x = rt, fill = param_label)) +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::facet_wrap(~ factor(choice, labels = c("Lower", "Upper"))) +
    ggplot2::coord_cartesian(xlim = rt_xlim) +
    ggplot2::scale_fill_viridis_d(option = color_palette) +
    ggplot2::labs(
      title = "RT Distributions",
      x = "Reaction Time (s)",
      y = "Density",
      fill = param_name
    ) +
    ggplot2::theme_minimal()

  # 2. Choice proportions
  choice_summary <- data %>%
    dplyr::group_by(param_value, param_label) %>%
    dplyr::summarise(
      prop_upper = mean(choice == 1),
      n_trials = dplyr::n(),
      .groups = 'drop'
    )

  p2 <- ggplot2::ggplot(choice_summary,
                        ggplot2::aes(x = param_value, y = prop_upper)) +
    ggplot2::geom_line(size = 2, color = "darkblue") +
    ggplot2::geom_point(size = 4, color = "darkblue") +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(
      title = "Choice Probability",
      x = param_name,
      y = "P(Upper Boundary)"
    ) +
    ggplot2::theme_minimal()

  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(p1 + p2 +
             patchwork::plot_annotation(
               title = paste("Effects of", param_name, "on DDM Behavior"),
               theme = ggplot2::theme(
                 plot.title = ggplot2::element_text(size = 16, face = "bold")
               )
             ))
  } else {
    return(p1)  # Return first plot if patchwork not available
  }
}

#' Create Summary Parameter Effect Plot
#' @keywords internal
create_summary_plot <- function(data, param_name, color_palette) {

  # Calculate summary statistics
  summary_stats <- data %>%
    dplyr::group_by(param_value, param_label) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      prop_upper = mean(choice == 1),
      mean_rt_all = mean(rt),
      mean_rt_upper = mean(rt[choice == 1]),
      mean_rt_lower = mean(rt[choice == 0]),
      rt_variability = sd(rt),
      rt_skew = (mean(rt) - median(rt)) / sd(rt),
      .groups = 'drop'
    )

  # Create multi-metric plot
  metrics_long <- summary_stats %>%
    dplyr::select(param_value, prop_upper, mean_rt_all,
                  rt_variability, rt_skew) %>%
    tidyr::pivot_longer(
      cols = -param_value,
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric_label = factor(
        metric,
        levels = c("prop_upper", "mean_rt_all", "rt_variability", "rt_skew"),
        labels = c("P(Upper)", "Mean RT", "RT Variability", "RT Skew")
      )
    )

  ggplot2::ggplot(metrics_long,
                  ggplot2::aes(x = param_value, y = value)) +
    ggplot2::geom_line(size = 2, color = "darkblue") +
    ggplot2::geom_point(size = 4, color = "darkblue") +
    ggplot2::facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = paste("Summary Effects of", param_name),
      x = param_name,
      y = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    )
}

#' Create Diagnostic Parameter Plot
#' @keywords internal
create_diagnostic_plot <- function(data, param_name, color_palette) {

  # Calculate diagnostic metrics
  diagnostics <- data %>%
    dplyr::group_by(param_value, param_label) %>%
    dplyr::summarise(
      n_trials = dplyr::n(),
      n_upper = sum(choice == 1),
      n_lower = sum(choice == 0),
      min_rt = min(rt),
      max_rt = max(rt),
      prop_fast = mean(rt < quantile(rt, 0.1)),
      prop_slow = mean(rt > quantile(rt, 0.9)),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      balance_ratio = n_upper / (n_upper + n_lower),
      rt_range = max_rt - min_rt
    )

  # Create diagnostic visualization
  p1 <- ggplot2::ggplot(diagnostics,
                        ggplot2::aes(x = param_value)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = min_rt, ymax = max_rt),
      alpha = 0.3, fill = "gray50"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = min_rt),
      color = "darkred", size = 2, linetype = "dashed"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = max_rt),
      color = "darkgreen", size = 2, linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "RT Range",
      x = param_name,
      y = "RT (s)"
    ) +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(diagnostics,
                        ggplot2::aes(x = param_value, y = balance_ratio)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_line(size = 2, color = "purple") +
    ggplot2::geom_point(size = 4, color = "purple") +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(
      title = "Response Balance",
      x = param_name,
      y = "Proportion Upper"
    ) +
    ggplot2::theme_minimal()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(p1 + p2 +
             patchwork::plot_annotation(
               title = paste("Diagnostic Plots for", param_name),
               theme = ggplot2::theme(
                 plot.title = ggplot2::element_text(size = 16, face = "bold")
               )
             ))
  } else {
    return(p1)
  }
}

#' Create Interactive Parameter Explorer
#'
#' This function creates an interactive visualization for exploring
#' parameter effects using plotly (if available).
#'
#' @param data DDM simulation data
#' @param param_name Name of the parameter being explored
#' @param param_values Vector of parameter values
#'
#' @return A plotly object if available, otherwise a static ggplot
#' @export
create_interactive_explorer <- function(data, param_name, param_values) {

  # Check if plotly is available
  if (!requireNamespace("plotly", quietly = TRUE)) {
    message("plotly not available. Creating static plot instead.")
    return(plot_parameter_effects(
      list(data), param_name, param_values, "summary"
    ))
  }

  # Create interactive plot with plotly
  # [Implementation would go here if needed]

  message("Interactive explorer not yet implemented")
  return(NULL)
}
