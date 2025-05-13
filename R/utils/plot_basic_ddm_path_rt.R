# R/utils/plot_basic_ddm_path_rt.R (or your chosen file name)

# ... (Roxygen comments as before, ensure they still match the simplified x-axis logic if you update them) ...

plot_ddm_paths_with_histograms <- function(trials_data_list,
                                           v_drift,
                                           z_start,
                                           a_threshold,
                                           max_time_to_plot = NULL, # User can specify, or it's auto
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
  
  # --- 1. Prepare Path Data ---
  all_paths_df_list <- vector("list", length(trials_data_list))
  for (i in seq_along(trials_data_list)) {
    trial_output <- trials_data_list[[i]]
    if (!is.null(trial_output$path_data) && nrow(trial_output$path_data) > 0) {
      path_df <- trial_output$path_data
      path_df$trial_id <- paste("Trial", i)
      path_df$final_choice <- trial_output$choice
      all_paths_df_list[[i]] <- path_df
    }
  }
  plot_df_paths <- bind_rows(all_paths_df_list)
  
  if (nrow(plot_df_paths) == 0) {
    warning("No valid path data to plot.")
    return(ggplot() + labs(title = "No path data") + theme_void())
  }
  
  # --- 2. Prepare RT Data for Histograms (using decision_time) ---
  rt_data_for_hist <- data.frame(
    decision_time = sapply(trials_data_list, function(x) if(!is.na(x$choice) && !is.na(x$decision_time)) x$decision_time else NA_real_),
    choice = sapply(trials_data_list, function(x) x$choice)
  ) %>%
    filter(!is.na(choice), !is.na(decision_time))
  
  rt_upper <- rt_data_for_hist %>% filter(choice == 1)
  rt_lower <- rt_data_for_hist %>% filter(choice == 0)
  
  # --- 3. Determine common x-axis limit for all plots ---
  # If max_time_to_plot is not provided by user, calculate from data
  current_max_time_to_plot <- max_time_to_plot # Use user value if provided
  if (is.null(current_max_time_to_plot)) {
    all_finite_times <- c(
      if(nrow(plot_df_paths) > 0) plot_df_paths$time_s[is.finite(plot_df_paths$time_s)] else numeric(0),
      if(nrow(rt_data_for_hist) > 0) rt_data_for_hist$decision_time[is.finite(rt_data_for_hist$decision_time)] else numeric(0)
    )
    if (length(all_finite_times) > 0) {
      current_max_time_to_plot <- max(all_finite_times, na.rm = TRUE) * 1.05 # Small buffer
    } else {
      current_max_time_to_plot <- 1.0 # Fallback if no valid times at all
    }
  }
  # Ensure it's a single finite positive number
  if (is.na(current_max_time_to_plot) || !is.finite(current_max_time_to_plot) || current_max_time_to_plot <= 0) {
    current_max_time_to_plot <- 1.0
  }
  
  
  # --- 3b. Determine common y-axis (count) limit for HISTOGRAMS ---
  max_count <- 0
  # Define breaks for consistent binning when calculating max_count
  hist_breaks <- seq(0, current_max_time_to_plot + hist_binwidth, by = hist_binwidth)
  
  if (nrow(rt_upper) > 1) { # Need >1 for hist() to give counts typically
    hist_upper_counts <- hist(rt_upper$decision_time, breaks = hist_breaks, plot = FALSE)$counts
    if (length(hist_upper_counts) > 0) max_count <- max(max_count, hist_upper_counts, na.rm = TRUE)
  }
  if (nrow(rt_lower) > 1) {
    hist_lower_counts <- hist(rt_lower$decision_time, breaks = hist_breaks, plot = FALSE)$counts
    if (length(hist_lower_counts) > 0) max_count <- max(max_count, hist_lower_counts, na.rm = TRUE)
  }
  y_hist_limit <- max_count * 1.1
  if (y_hist_limit == 0 || !is.finite(y_hist_limit)) y_hist_limit <- 1 # Prevent limit of 0 or NA/Inf
  
  
  # --- 4. Create Main Path Plot (p_main) ---
  min_evidence_val <- min(0, min(plot_df_paths$evidence, na.rm = TRUE) - 0.05 * a_threshold, na.rm = TRUE)
  max_evidence_val <- max(a_threshold, max(plot_df_paths$evidence, na.rm = TRUE) + 0.05 * a_threshold, na.rm = TRUE)
  if (is.infinite(min_evidence_val) || is.na(min_evidence_val)) min_evidence_val <- -0.1 * a_threshold
  if (is.infinite(max_evidence_val) || is.na(max_evidence_val)) max_evidence_val <- a_threshold + 0.1 * a_threshold
  
  # Define custom breaks for the y-axis
  y_axis_breaks <- sort(unique(c(0, z_start, a_threshold))) 
  # We can add more breaks if the range is large, or let ggplot2 add some too.
  # For more control, you could define them more explicitly, e.g.,
  # y_axis_breaks <- pretty(c(min_evidence_val, max_evidence_val), n=5) # ggplot2's default-like breaks
  # y_axis_breaks <- sort(unique(c(y_axis_breaks, 0, z_start, a_threshold))) # Combine with key points
  
  p_main <- ggplot(plot_df_paths, aes(x = time_s, y = evidence, group = trial_id, color = factor(final_choice, levels=c(0,1,NA)))) +
    geom_line(alpha = 0.25, linewidth = 0.4) +
    # Horizontal lines are still useful for emphasis even with breaks
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_hline(yintercept = a_threshold, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_hline(yintercept = z_start, linetype = "dotted", color = "grey15", linewidth = 0.8) +
    scale_color_manual(values = c("0" = "salmon", "1" = "steelblue", "NA" = "grey70"),
                       labels = c("0" = "Path to Lower", "1" = "Path to Upper", "NA" = "Path Timeout"),
                       name = "Path Outcome", drop = FALSE,
                       guide = guide_legend(override.aes = list(alpha = 1, linewidth=1.5))) +
    # ADDED/MODIFIED scale_y_continuous for custom breaks
    scale_y_continuous(breaks = y_axis_breaks, # Use our defined breaks
                       labels = scales::number_format(accuracy = 0.01)) + # Format labels if needed
    labs(x = NULL, y = "Accumulated Evidence") +
    coord_cartesian(xlim = c(0, current_max_time_to_plot), ylim = c(min_evidence_val, max_evidence_val), expand=FALSE, clip="off") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size=16, face="bold"), # Will be overridden
          axis.title.y = element_text(size=14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey80"), # Ensure major y grid lines are visible
          panel.grid.minor.y = element_blank(), # Optionally remove minor y grid lines
          legend.box = "horizontal",
          plot.margin = margin(t=5, r=15, b=0, l=5))
  
  # Add boundary labels to the main plot - on the RIGHT side
  p_main <- p_main +
    annotate("text",
             x = current_max_time_to_plot,
             y = a_threshold,
             label = paste("a =", round(a_threshold,2)),
             hjust = 1.1,
             vjust = -0.5,
             size=3.5, color="black") +
    annotate("text",
             x = current_max_time_to_plot,
             y = 0,
             label = "0",
             hjust = 1.1,
             vjust = 1.5,
             size=3.5, color="black")
  
  # --- 5. Create Upper RT Histogram (p_upper_hist) ---
  p_upper_hist <- ggplot() + theme_void() + theme(plot.margin = margin(t=5, r=15, b=0, l=15))
  if (nrow(rt_upper) > 1) {
    p_upper_hist <- ggplot(rt_upper, aes(x = decision_time)) +
      geom_histogram(binwidth = hist_binwidth, fill = "steelblue", color = "black", alpha = 0.7,
                     breaks = hist_breaks) + # Use common breaks
      scale_x_continuous(limits = c(0, current_max_time_to_plot), expand = c(0,0)) +
      scale_y_continuous(limits = c(0, y_hist_limit), expand = expansion(mult = c(0, 0.05))) +
      labs(title = NULL, x = NULL, y = "Count") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.title.y = element_text(size=10),
            plot.margin = margin(t=5, r=15, b=0, l=15))
  }
  
  # --- 6. Create Lower RT Histogram (p_lower_hist) ---
  p_lower_hist <- ggplot() + theme_void() + theme(plot.margin = margin(t=0, r=15, b=5, l=15))
  if (nrow(rt_lower) > 1) {
    p_lower_hist <- ggplot(rt_lower, aes(x = decision_time)) +
      geom_histogram(binwidth = hist_binwidth, fill = "salmon", color = "black", alpha = 0.7,
                     breaks = hist_breaks) + # Use common breaks
      scale_x_continuous(limits = c(0, current_max_time_to_plot), expand = c(0,0)) +
      scale_y_reverse(limits = c(y_hist_limit, 0), expand = expansion(mult = c(0.05, 0))) +
      labs(title = NULL, x = "Decision Time (s)", y = "Count") +
      theme_minimal(base_size = 10) +
      theme(axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=10),
            plot.margin = margin(t=0, r=15, b=5, l=15))
  }
  
  # --- 7. Combine Plots using Patchwork ---
  plot_list <- list()
  if(is.ggplot(p_upper_hist) && !"theme_void" %in% class(p_upper_hist$theme)) plot_list <- c(plot_list, list(p_upper_hist))
  if(is.ggplot(p_main)) plot_list <- c(plot_list, list(p_main)) # Main plot should always exist if data
  if(is.ggplot(p_lower_hist) && !"theme_void" %in% class(p_lower_hist$theme)) plot_list <- c(plot_list, list(p_lower_hist))
  
  height_ratios <- c()
  if(is.ggplot(p_upper_hist) && !"theme_void" %in% class(p_upper_hist$theme)) height_ratios <- c(height_ratios, 1)
  if(is.ggplot(p_main)) height_ratios <- c(height_ratios, 3)
  if(is.ggplot(p_lower_hist) && !"theme_void" %in% class(p_lower_hist$theme)) height_ratios <- c(height_ratios, 1)
  
  if(length(plot_list) > 0) { # Should be at least 1 (p_main)
    # Remove titles from individual plots before combining, as plot_annotation will add one
    for(i in seq_along(plot_list)) {
      if(is.ggplot(plot_list[[i]])) plot_list[[i]] <- plot_list[[i]] + labs(title=NULL)
    }
    
    final_plot <- Reduce(`/`, plot_list) +
      plot_layout(heights = height_ratios, guides = 'collect') &
      theme(legend.position = 'bottom')
    final_plot <- final_plot + plot_annotation(title = main_plot_title, 
                                               theme = theme(plot.title = element_text(hjust = 0.5, size=16, face="bold")))
  } else {
    final_plot <- ggplot() + labs(title = "Not enough data for full plot") + theme_void()
  }
  
  return(final_plot)
}