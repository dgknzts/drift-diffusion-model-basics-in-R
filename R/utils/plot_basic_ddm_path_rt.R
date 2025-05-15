# R/utils/plot_basic_ddm_path_rt.R 

#' Plot DDM Evidence Paths with Flanking RT Histograms
#'
#' Creates a multi-panel plot visualizing Diffusion Decision Model (DDM) dynamics.
#' The central panel displays simulated evidence accumulation paths. Flanking this
#' panel are histograms of decision times: one above for trials reaching the
#' upper boundary, and one below (flipped) for trials reaching the lower boundary.
#' Requires the 'ggplot2', 'dplyr', and 'patchwork' packages.
#'
#' @param trials_data_list A list. Each element must be the complete output
#'   from `simulate_diffusion_trial_with_path()` or a compatible function
#'   that returns `path_data`, `choice`, and `decision_time`.
#' @param v_drift Numeric. The mean drift rate used in the simulations (for the average drift line).
#' @param z_start Numeric. The mean starting point `z` used in the simulations.
#' @param a_threshold Numeric. The upper decision threshold `a`.
#' @param max_time_to_plot Numeric, optional. Sets the maximum limit for the
#'   x-axis (decision time in seconds). If `NULL` (default), the function
#'   auto-calculates a limit based on observed data.
#' @param hist_binwidth Numeric, optional. Binwidth for RT histograms (seconds). Default is `0.05`.
#' @param main_plot_title Character, optional. Overall title for the combined plot.
#'
#' @return A `patchwork` ggplot object.
#' @export
plot_ddm_paths_with_histograms <- function(trials_data_list,
                                           v_drift,
                                           z_start,
                                           a_threshold,
                                           max_time_to_plot = NULL,
                                           hist_binwidth = 0.05,
                                           main_plot_title = "DDM Simulation Output") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Packages 'ggplot2', 'dplyr', and 'patchwork' are needed. Please install them.", call. = FALSE)
  }
  
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  # --- Determine FINAL x-axis limit for all plots (final_plot_xlim) ---
  # This must be calculated first.
  final_plot_xlim <- NULL
  if (!is.null(max_time_to_plot) && is.finite(max_time_to_plot) && max_time_to_plot > 0) {
    final_plot_xlim <- max_time_to_plot
  } else {
    temp_all_path_times <- unlist(lapply(trials_data_list, function(x) if(!is.null(x$path_data) && nrow(x$path_data)>0) x$path_data$time_s else numeric(0)))
    temp_all_decision_times <- sapply(trials_data_list, function(x) if(!is.null(x$decision_time)) x$decision_time else NA_real_) # Handle NULL decision_time
    
    all_finite_times_for_auto_limit <- c(
      temp_all_path_times[is.finite(temp_all_path_times)],
      temp_all_decision_times[is.finite(temp_all_decision_times)]
    )
    if (length(all_finite_times_for_auto_limit) > 0 && any(is.finite(all_finite_times_for_auto_limit))) {
      final_plot_xlim <- max(all_finite_times_for_auto_limit, na.rm = TRUE) * 1.05
    } else {
      final_plot_xlim <- 1.0
    }
  }
  if (is.na(final_plot_xlim) || !is.finite(final_plot_xlim) || final_plot_xlim <= 0) {
    final_plot_xlim <- 1.0
  }
  
  # --- 1. Prepare Path Data: Process each trial individually for outcome and truncation ---
  processed_paths_list <- list() # Use a new list
  
  for (i in seq_along(trials_data_list)) {
    trial_output <- trials_data_list[[i]]
    
    if (is.null(trial_output$path_data) || nrow(trial_output$path_data) == 0) {
      next # Skip if no path data
    }
    
    path_df_current_trial <- as.data.frame(trial_output$path_data) # Ensure it's a data frame
    
    # Determine the display outcome for THIS trial based on final_plot_xlim
    display_outcome_char <- NA_character_
    original_choice <- trial_output$choice
    original_dt_val <- trial_output$decision_time # Can be NA
    
    if (is.na(original_choice)) { # Was an original timeout
      display_outcome_char <- "Timeout"
    } else if (!is.na(original_dt_val) && original_dt_val > final_plot_xlim) {
      # Made a choice, but its true decision time is beyond the current display x-limit
      display_outcome_char <- "Timeout" # Treat as timeout for display color
    } else {
      # Made a choice within display limit (or original_dt_val is NA but original_choice is not - less likely)
      display_outcome_char <- as.character(original_choice) # "0" or "1"
    }
    
    # Truncate path data for plotting up to final_plot_xlim
    path_df_truncated <- path_df_current_trial %>%
      filter(time_s <= final_plot_xlim)
    
    if (nrow(path_df_truncated) > 0) {
      path_df_truncated$trial_id <- paste("Trial", i)
      path_df_truncated$display_outcome_factor <- factor(display_outcome_char, levels = c("0", "1", "Timeout"))
      processed_paths_list[[length(processed_paths_list) + 1]] <- path_df_truncated
    }
  }
  
  if (length(processed_paths_list) == 0) {
    warning("No valid path data remains after filtering by max_time_to_plot.")
    return(ggplot() + labs(title = "No path data to plot") + theme_void())
  }
  plot_df_paths <- bind_rows(processed_paths_list)
  
  
  # --- 2. Prepare RT Data for HISTOGRAMS ---
  rt_data_for_hist <- data.frame(
    decision_time = sapply(trials_data_list, function(x) x$decision_time),
    choice = sapply(trials_data_list, function(x) x$choice)
  ) %>%
    filter(!is.na(choice), !is.na(decision_time), is.finite(decision_time),
           decision_time >= 0, decision_time <= final_plot_xlim)
  
  rt_upper <- rt_data_for_hist %>% filter(choice == 1)
  rt_lower <- rt_data_for_hist %>% filter(choice == 0)
  
  # --- 3b. Determine common y-axis (count) limit for HISTOGRAMS ---
  # (This logic remains the same)
  max_count <- 0
  hist_breaks <- seq(0, final_plot_xlim + hist_binwidth, by = hist_binwidth)
  if (nrow(rt_upper) > 1) {
    hist_upper_counts <- hist(rt_upper$decision_time, breaks = hist_breaks, plot = FALSE)$counts
    if (length(hist_upper_counts) > 0) max_count <- max(max_count, hist_upper_counts, na.rm = TRUE)
  }
  if (nrow(rt_lower) > 1) {
    hist_lower_counts <- hist(rt_lower$decision_time, breaks = hist_breaks, plot = FALSE)$counts
    if (length(hist_lower_counts) > 0) max_count <- max(max_count, hist_lower_counts, na.rm = TRUE)
  }
  y_hist_limit <- max_count * 1.1
  if (y_hist_limit == 0 || !is.finite(y_hist_limit)) y_hist_limit <- 1
  
  # --- 4. Create Main Path Plot (p_main) ---
  min_evidence_val <- min(0, min(plot_df_paths$evidence, na.rm = TRUE) - 0.05 * a_threshold, na.rm = TRUE)
  max_evidence_val <- max(a_threshold, max(plot_df_paths$evidence, na.rm = TRUE) + 0.05 * a_threshold, na.rm = TRUE)
  if (is.infinite(min_evidence_val) || is.na(min_evidence_val) || length(min_evidence_val)==0) min_evidence_val <- -0.1 * a_threshold # Added length check
  if (is.infinite(max_evidence_val) || is.na(max_evidence_val) || length(max_evidence_val)==0) max_evidence_val <- a_threshold + 0.1 * a_threshold # Added length check
  
  
  p_main <- ggplot(plot_df_paths, aes(x = time_s, y = evidence, group = trial_id, color = display_outcome_factor)) + # CRITICAL: use display_outcome_factor
    geom_line(alpha = 0.25, linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_hline(yintercept = a_threshold, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_hline(yintercept = z_start, linetype = "dotted", color = "grey15", linewidth = 0.8) +
    annotate("segment", x = 0, y = z_start, xend = final_plot_xlim, yend = z_start + v_drift * final_plot_xlim,
             color = "black", linetype = "solid", linewidth = 1, alpha = 0.7) +
    scale_color_manual(name = "Displayed Path Outcome",
                       values = c("0" = "salmon", "1" = "steelblue", "Timeout" = "grey50"),
                       labels = c("0" = "Path to Lower (within limit)",
                                  "1" = "Path to Upper (within limit)",
                                  "Timeout" = "Path Timed Out (at display limit)"),
                       drop = FALSE) + # drop = FALSE ensures all levels appear in legend
    labs(x = NULL, y = "Accumulated Evidence") +
    coord_cartesian(xlim = c(0, final_plot_xlim), ylim = c(min_evidence_val, max_evidence_val), expand = FALSE, clip = "off") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", axis.title.y = element_text(size = 14),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.box = "horizontal", plot.margin = margin(t = 5, r = 15, b = 0, l = 15))
  p_main <- p_main +
    annotate("text", x = final_plot_xlim, y = a_threshold, label = paste(" a =", round(a_threshold, 2)), hjust = 1.1, vjust = -0.5, size = 3.5, color = "black") +
    annotate("text", x = final_plot_xlim, y = 0, label = " 0", hjust = 1.1, vjust = 1.5, size = 3.5, color = "black")
  
  # --- 5. Create Upper RT Histogram (p_upper_hist) ---
  # ... (no changes here) ...
  p_upper_hist <- ggplot() + theme_void() + theme(plot.margin = margin(t = 5, r = 15, b = 0, l = 15))
  if (nrow(rt_upper) > 1) {
    p_upper_hist <- ggplot(rt_upper, aes(x = decision_time)) +
      geom_histogram(binwidth = hist_binwidth, fill = "steelblue", color = "black", alpha = 0.7, boundary = 0) +
      scale_x_continuous(limits = c(0, final_plot_xlim), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, y_hist_limit), expand = expansion(mult = c(0, 0.05))) +
      labs(title = NULL, x = NULL, y = "Count") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
            axis.title.y = element_text(size = 10), plot.margin = margin(t = 5, r = 15, b = 0, l = 15))
  }
  
  # --- 6. Create Lower RT Histogram (p_lower_hist) ---
  # ... (no changes here) ...
  p_lower_hist <- ggplot() + theme_void() + theme(plot.margin = margin(t = 0, r = 15, b = 5, l = 15))
  if (nrow(rt_lower) > 1) {
    p_lower_hist <- ggplot(rt_lower, aes(x = decision_time)) +
      geom_histogram(binwidth = hist_binwidth, fill = "salmon", color = "black", alpha = 0.7, boundary = 0) +
      scale_x_continuous(limits = c(0, final_plot_xlim), expand = c(0, 0)) +
      scale_y_reverse(limits = c(y_hist_limit, 0), expand = expansion(mult = c(0.05, 0))) +
      labs(title = NULL, x = "Decision Time (s)", y = "Count") +
      theme_minimal(base_size = 10) +
      theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 10),
            plot.margin = margin(t = 0, r = 15, b = 5, l = 15))
  }
  
  # --- 7. Combine Plots using Patchwork ---
  # ... (patchwork combination logic remains the same) ...
  plot_list <- list()
  if(is.ggplot(p_upper_hist) && !"theme_void" %in% class(p_upper_hist$theme)) plot_list <- c(plot_list, list(p_upper_hist))
  if(is.ggplot(p_main)) plot_list <- c(plot_list, list(p_main)) 
  if(is.ggplot(p_lower_hist) && !"theme_void" %in% class(p_lower_hist$theme)) plot_list <- c(plot_list, list(p_lower_hist))
  
  height_ratios <- c()
  if(is.ggplot(p_upper_hist) && !"theme_void" %in% class(p_upper_hist$theme)) height_ratios <- c(height_ratios, 1)
  if(is.ggplot(p_main)) height_ratios <- c(height_ratios, 3)
  if(is.ggplot(p_lower_hist) && !"theme_void" %in% class(p_lower_hist$theme)) height_ratios <- c(height_ratios, 1)
  
  combined_final_plot <- ggplot() + labs(title = "Not enough data or plots to combine") + theme_void() 
  if(length(plot_list) > 0 && length(height_ratios) == length(plot_list) && length(plot_list) >=1) { 
    for(i in seq_along(plot_list)) {
      if(is.ggplot(plot_list[[i]])) plot_list[[i]] <- plot_list[[i]] + labs(title=NULL)
    }
    
    combined_final_plot <- Reduce(`/`, plot_list) +
      plot_layout(heights = height_ratios, guides = 'collect') &
      theme(legend.position = 'bottom')
    combined_final_plot <- combined_final_plot + 
      plot_annotation(title = main_plot_title,
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  } else if (length(plot_list) == 1 && is.ggplot(plot_list[[1]])) { 
    combined_final_plot <- plot_list[[1]] + 
      labs(title = main_plot_title) + 
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  }
  
  return(combined_final_plot)
}