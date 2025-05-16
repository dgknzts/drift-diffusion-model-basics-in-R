# R/utils/plot_parameter_comparison.R

#' Plot a Comparison of DDM Simulation Outputs for Parameter Exploration
#'
#' This function takes two sets of DDM simulation data (typically a baseline
#' and a version with one parameter modified) and generates a combined plot
#' showing RT density distributions and choice proportions for comparison.
#'
#' @param data_sim1 Data frame. Simulation output for the first (e.g., baseline) condition.
#'   Must contain 'rt', 'choice', and a 'model_label' column (e.g., "Baseline v=0.15").
#' @param data_sim2 Data frame. Simulation output for the second (e.g., modified) condition.
#'   Must contain 'rt', 'choice', and a 'model_label' column (e.g., "Modified v=0.25").
#' @param param_varied Character. A string describing the parameter that was varied
#'   and its values, for the plot title (e.g., "Effect of Mean Drift Rate (mean_v)").
#' @param rt_xlim Numeric vector of length 2, optional. X-axis limits for the RT density plots.
#'   If NULL (default), limits are determined automatically from the data.
#' @param density_alpha Numeric. Alpha transparency for density plots. Default 0.6.
#' @param hist_binwidth Numeric. Binwidth for optional histograms under densities. Default 0.05.
#'   Set to NULL if no histogram desired.
#'
#' @return A `patchwork` ggplot object.
#'
#' @details
#' The function generates:
#' 1. Overlaid RT density plots for choice 0 and choice 1, comparing `data_sim1` and `data_sim2`.
#' 2. Bar plots (or text) showing choice proportions for each simulation.
#' These are arranged using the `patchwork` package. It assumes `choice` is coded as 0 and 1.
#'
#' @export
#' @examples
#' # This is a conceptual example, actual data needs to be simulated first.
#' # Ensure your DDM simulation functions (e.g., simulate_diffusion_experiment_variable)
#' # are defined and sourced if running this example standalone.
#'
#' if (interactive() &&
#'     requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE) &&
#'     exists("simulate_diffusion_experiment_variable")) { # Check for one of your sim functions
#'
#'   # library(ggplot2) # Loaded by function
#'   # library(dplyr)   # Loaded by function
#'   # library(patchwork)# Loaded by function
#'
#'   # Simulate baseline data
#'   set.seed(123)
#'   baseline_data <- simulate_diffusion_experiment_variable(
#'     n_trials = 500, mean_v = 0.15, a = 1, mean_z = 0.5, s = 0.35, mean_ter = 0.2,
#'     sv = 0, sz = 0, st0 = 0, dt = 0.01 # No across-trial variability for baseline
#'   )
#'   baseline_data$model_label <- "Baseline (v=0.15, sv=0)"
#'
#'   # Simulate data with modified sv
#'   set.seed(123) # Use same seed for comparison of sv effect
#'   modified_data_sv <- simulate_diffusion_experiment_variable(
#'     n_trials = 500, mean_v = 0.15, a = 1, mean_z = 0.5, s = 0.35, mean_ter = 0.2,
#'     sv = 0.2, sz = 0, st0 = 0, dt = 0.01 # sv is now 0.2
#'   )
#'   modified_data_sv$model_label <- "sv = 0.2"
#'
#'   # Create the comparison plot
#'   # p_sv_comp <- plot_ddm_parameter_comparison(
#'   #   data_sim1 = baseline_data,
#'   #   data_sim2 = modified_data_sv,
#'   #   param_varied = "Effect of Drift Rate Variability (sv)",
#'   #   rt_xlim = c(0, 3) # Optional x-limit
#'   # )
#'   # print(p_sv_comp)
#' }
plot_ddm_parameter_comparison <- function(data_sim1,
                                          data_sim2,
                                          param_varied = "Parameter Comparison",
                                          rt_xlim = NULL,
                                          density_alpha = 0.6,
                                          hist_binwidth = 0.05) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Packages 'ggplot2', 'dplyr', and 'patchwork' are needed. Please install them.", call. = FALSE)
  }
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  # --- 1. Combine Data and Prepare ---
  if (!"model_label" %in% names(data_sim1) || !"model_label" %in% names(data_sim2)) {
    stop("Both data_sim1 and data_sim2 must have a 'model_label' column.")
  }
  
  combined_data <- bind_rows(data_sim1, data_sim2) %>%
    filter(!is.na(rt) & !is.na(choice)) %>% # Ensure valid RTs and choices for densities
    mutate(
      choice_label = factor(choice, levels = c(0, 1), labels = c("Choice 0", "Choice 1")),
      model_label = factor(model_label) # Ensure model_label is a factor for consistent ordering
    )
  
  if (nrow(combined_data) == 0) {
    warning("No valid data after filtering NAs from RT and choice.")
    return(ggplot() + labs(title = "No data to plot") + theme_void())
  }
  
  # Determine common x-axis limits if not provided
  current_rt_xlim <- rt_xlim
  if (is.null(current_rt_xlim)) {
    if (nrow(combined_data) > 0 && any(is.finite(combined_data$rt))) {
      q99 <- quantile(combined_data$rt, 0.99, na.rm = TRUE)
      current_rt_xlim <- c(0, q99 * 1.05)
    } else {
      current_rt_xlim <- c(0, 1) # Fallback
    }
  }
  if(any(is.na(current_rt_xlim)) || any(!is.finite(current_rt_xlim)) || current_rt_xlim[2] <= current_rt_xlim[1]){
    current_rt_xlim <- c(0,1) # Robust fallback
  }
  
  
  # --- 2. Create RT Density Plots ---
  # Ensure enough data for density estimation within each group
  data_for_density <- combined_data %>% group_by(model_label, choice_label) %>% filter(n() > 1) %>% ungroup()
  
  p_rt_densities <- ggplot() + theme_void() # Placeholder
  if(nrow(data_for_density) > 0) {
    p_rt_densities <- ggplot(data_for_density, aes(x = rt, fill = model_label, color = model_label)) +
      geom_density(alpha = density_alpha, linewidth = 0.7)
    
    if (!is.null(hist_binwidth) && is.numeric(hist_binwidth) && hist_binwidth > 0) {
      p_rt_densities <- p_rt_densities +
        geom_histogram(aes(y = ..density..), # Plot density on y-axis to match geom_density
                       binwidth = hist_binwidth,
                       alpha = density_alpha / 2, # Make histogram more transparent
                       position = "identity") # Important for overlaid histograms
    }
    
    p_rt_densities <- p_rt_densities +
      facet_wrap(~choice_label, ncol = 1) + # Separate rows for choice 0 and 1
      coord_cartesian(xlim = current_rt_xlim) +
      labs(x = "Reaction Time (s)", y = "Density", fill = "Simulation:", color = "Simulation:") +
      theme_bw(base_size = 10) +
      theme(legend.position = "top",
            strip.text = element_text(face = "bold", size=10),
            plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
  }
  
  
  # --- 3. Calculate and Plot Choice Proportions ---
  choice_props <- combined_data %>% # Use combined_data before RT filtering for proportions
    group_by(model_label, choice_label) %>%
    summarise(N = n(), .groups = "drop") %>%
    group_by(model_label) %>%
    mutate(Proportion = N / sum(N)) %>%
    ungroup()
  
  p_choice_props <- ggplot() + theme_void() # Placeholder
  if(nrow(choice_props) > 0){
    p_choice_props <- ggplot(choice_props, aes(x = choice_label, y = Proportion, fill = model_label)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth=0.3) +
      geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
                position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
      facet_wrap(~model_label, ncol=1) + # Separate plot for each model for clarity
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult=c(0, 0.05))) +
      labs(x = "Choice Outcome", y = "Proportion") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", # Legend is already with density plot
            strip.text = element_text(face = "bold", size=10),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
  }
  
  # --- 4. Combine Plots with Patchwork ---
  final_plot <- ggplot() + labs(title = "Not enough data for one or more plots") + theme_void()
  
  if (is.ggplot(p_rt_densities) && !"theme_void" %in% class(p_rt_densities$theme) && 
      is.ggplot(p_choice_props) && !"theme_void" %in% class(p_choice_props$theme)) {
    
    final_plot <- p_rt_densities + p_choice_props +
      plot_layout(widths = c(2, 1), guides = "collect") + # Density plot wider
      plot_annotation(
        title = param_varied,
        theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                      legend.position = "top")
      ) & theme(legend.position = 'top') # Apply legend position to collected legends
    
  } else if (is.ggplot(p_rt_densities) && !"theme_void" %in% class(p_rt_densities$theme)) {
    final_plot <- p_rt_densities + 
      plot_annotation(title = param_varied, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  } else if (is.ggplot(p_choice_props) && !"theme_void" %in% class(p_choice_props$theme)) {
    final_plot <- p_choice_props +
      plot_annotation(title = param_varied, theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  }
  
  
  return(final_plot)
}