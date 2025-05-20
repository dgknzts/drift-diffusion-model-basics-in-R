# R/utils/plot_qpp.R

#' Plot Quantile Probability Plot with RT Density Ridges and Filled Quantile Regions
#'
#' For each condition (faceted) and response type (Correct/Error, positioned by
#' their probability), this function draws the full RT distribution as a density
#' ridge. The regions between specified quantiles within each ridge are filled
#' with different colors. Text labels indicate "Correct" and "Error" ridges.
#'
#' @param data_list Named list of data.frames (one per condition). Each
#'   data.frame must have `rt` and `choice` columns.
#' @param correct_choice_value Value in the `choice` column that denotes a
#'   correct/upper-boundary response. Default `1`.
#' @param error_choice_value Value in the `choice` column that denotes an
#'   error/lower-boundary response. Default `0`.
#' @param plot_title Title for the plot.
#' @param rt_lim Optional numeric vector [min, max] for RT (y-axis after `coord_flip()`).
#' @param condition_numbers Logical; if `TRUE` conditions are labelled by a
#'   simple integer index (1, 2, …). If `FALSE`, the names in `data_list` are
#'   used verbatim. Default `TRUE`.
#' @param num_quantiles_fill Integer. Number of quantile regions to fill (e.g.,
#'   `4` for quartiles, resulting in 4 filled regions bounded by 0, Q1, Median, Q3, Max).
#'   Passed to `ggridges::stat_density_ridges(quantiles = ...)`. Default `4`.
#' @param ridge_scale Passed to `ggridges::geom_density_ridges(scale = ...)`. Default `0.9`.
#' @param alpha_quantile_fill Alpha transparency applied to quantile region fills. Default `0.7`.
#' @param show_median_line Logical. If TRUE, a line is drawn for the median. Default `TRUE`.
#' @param median_line_color Color for the median line. Default `"black"`.
#' @param median_line_size Size for the median line. Default `0.5`.
#' @param text_label_size Numeric. Font size for "Correct"/"Error" labels. Default `3.5`.
#' @param text_label_offset_factor Numeric. Factor to control offset of "Correct"/"Error"
#'   labels from their respective probability lines (on the probability axis after `coord_flip()`).
#'   Default `0.03`.
#'
#' @return A `ggplot` object showing the ridge-style quantile-probability plot.
#' @export
#'
#' @examples
#' if (interactive() &&
#'     requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("dplyr",  quietly = TRUE) &&
#'     requireNamespace("ggridges", quietly = TRUE) &&
#'     requireNamespace("scales", quietly = TRUE) &&
#'     exists("simulate_diffusion_experiment_variable")) { # Check for your sim function
#'
#'   # Simulate some data for example
#'   set.seed(123)
#'   data_condA <- simulate_diffusion_experiment_variable(
#'     n_trials=300, mean_v=0.1, a=1, mean_z=0.5, s=0.1, mean_ter=0.2, dt=0.001
#'   )
#'   set.seed(456)
#'   data_condB <- simulate_diffusion_experiment_variable(
#'     n_trials=300, mean_v=0.2, a=1, mean_z=0.5, s=0.1, mean_ter=0.2, dt=0.001
#'   )
#'   example_data_list <- list("Condition A (v=0.1)" = data_condA,
#'                             "Condition B (v=0.2)" = data_condB)
#'
#'   # p <- plot_qpp(
#'   #   data_list = example_data_list,
#'   #   plot_title = "Example QPP with Filled Quantiles",
#'   #   num_quantiles_fill = 4, # Quartiles
#'   #   rt_lim = c(0.2, 1.5)
#'   # )
#'   # print(p)
#' }
plot_qpp <- function(
    data_list,
    correct_choice_value = 1,
    error_choice_value = 0,
    plot_title = "RT Density Ridges with Filled Quantile Regions",
    rt_lim = NULL,
    condition_numbers = TRUE,
    num_quantiles_fill = 4,
    ridge_scale = 0.9,
    alpha_quantile_fill = 0.7,
    show_median_line = TRUE,
    median_line_color = "black",
    median_line_size = 0.5,
    text_label_size = 3.5,
    text_label_offset_factor = 0.03,
    vertical_label_offset = 0.01
) {
  # ---- Package checks ----
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr",  quietly = TRUE) ||
      !requireNamespace("ggridges", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) {
    stop("Packages 'ggplot2', 'dplyr', 'ggridges', and 'scales' are required.", call. = FALSE)
  }
  library(ggplot2)
  library(dplyr)
  library(ggridges)
  library(scales)

  if (!is.list(data_list) || is.null(names(data_list))) {
    stop("'data_list' must be a *named* list of data frames.")
  }

  # ---- Build long format data for plotting ----
  plot_data_list  <- list()
  cond_index  <- 0L

  for (cond_name in names(data_list)) {
    cond_index <- cond_index + 1L
    dat <- data_list[[cond_name]]

    if (!all(c("rt", "choice") %in% names(dat))) {
      warning(sprintf("Cond '%s' missing 'rt'/'choice' – skipped.", cond_name)); next
    }
    dat_filtered <- dplyr::filter(dat, !is.na(rt) & !is.na(choice) & is.finite(rt))
    if (nrow(dat_filtered) < 2) {
      warning(sprintf("Cond '%s' has < 2 valid trials – skipped.", cond_name)); next
    }
    n_tot_valid_condition <- nrow(dat_filtered)
    label_this_cond <- if (condition_numbers) sprintf("Cond %d", cond_index) else cond_name

    correct_df <- dplyr::filter(dat_filtered, choice == correct_choice_value)
    if(nrow(correct_df) >= 2) {
      p_correct  <- nrow(correct_df) / n_tot_valid_condition
      plot_data_list[[length(plot_data_list) + 1L]] <- correct_df %>%
        mutate(probability_group = p_correct, response_type = "Correct",
               condition_label_facet = label_this_cond,
               ridge_group = paste(label_this_cond, "Correct"))
    }
    error_df <- dplyr::filter(dat_filtered, choice == error_choice_value)
    if(nrow(error_df) >= 2) {
      p_error    <- nrow(error_df) / n_tot_valid_condition
      plot_data_list[[length(plot_data_list) + 1L]] <- error_df %>%
        mutate(probability_group = p_error, response_type = "Error",
               condition_label_facet = label_this_cond,
               ridge_group = paste(label_this_cond, "Error"))
    }
  }

  if (length(plot_data_list) == 0) {
    warning("No suitable data to plot."); return(ggplot() + theme_void() + labs(title="No data"))
  }
  plot_df  <- dplyr::bind_rows(plot_data_list) %>%
    mutate(response_type = factor(response_type, levels=c("Error", "Correct")))

  # ---- Prepare data for text labels ----
  text_label_positions <- plot_df %>%
    group_by(condition_label_facet, response_type, probability_group) %>%
    summarise(rt_for_label_pos = if(n() > 1) max(rt, na.rm = TRUE) else if(n()==1) rt else NA_real_,
              .groups = "drop") %>%
    filter(!is.na(rt_for_label_pos)) %>%
    mutate(
      prob_axis_label_pos = ifelse(response_type == "Correct",
                                   probability_group - text_label_offset_factor,
                                   probability_group - text_label_offset_factor),
      prob_axis_label_pos = case_when(
        prob_axis_label_pos > 1.0 ~ 1.0 - (text_label_offset_factor / 2),
        prob_axis_label_pos < 0.0 ~ 0.0 + (text_label_offset_factor / 2),
        TRUE ~ prob_axis_label_pos
      )
    )

  # ---- Build plot ----
  plt <- ggplot(plot_df,
                aes(x = rt, y = probability_group, group = ridge_group)) + # Main aes for x,y,group
    ggridges::stat_density_ridges( # Layer for filled quantile regions
      aes(fill = factor(stat(quantile))), # Fill is specific to this layer
      geom = "density_ridges_gradient",
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.01, height = 0.001),
      point_shape = '-', point_size = 4, point_alpha = 1,
      calc_ecdf = TRUE,
      quantiles = num_quantiles_fill,
      alpha = alpha_quantile_fill,
      scale = ridge_scale,
      rel_min_height = 0.01,
      color = "black", # Outline of the main ridge shape
      linewidth = 0.3
    )

  if(show_median_line && num_quantiles_fill >=2 ){
    plt <- plt + ggridges::stat_density_ridges( # Layer for median line
      aes(x = rt, y = probability_group, group = ridge_group), # Ensure aes is specified or inherited correctly
      geom = "density_ridges",
      calc_ecdf = FALSE,
      quantiles = 2, # Median
      quantile_lines = TRUE,
      alpha = 0, # Transparent fill for this layer
      scale = ridge_scale,
      color = median_line_color, # Color of the median line
      linewidth = median_line_size,
      linetype = "solid"
    )
  }

  if(nrow(text_label_positions) > 0) { # Layer for text labels
    plt <- plt +
      geom_text(data = text_label_positions,
                aes(x = rt_for_label_pos + vertical_label_offset, y = prob_axis_label_pos, label = response_type),
                inherit.aes = FALSE, # CRITICAL: Do not inherit fill from main aes
                fontface = "bold", size = text_label_size,
                hjust = 0.5,
                color = "black")
  }

  plt <- plt +
    facet_wrap(~condition_label_facet, ncol = 1, scales = "fixed") +
    coord_flip() +
    scale_x_continuous(limits = rt_lim, name = "Reaction Time (s)") +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.2),
                       labels = scales::percent_format(accuracy = 1), name = "Response Probability") +
    scale_fill_viridis_d(name = paste0(num_quantiles_fill, "-Quantile Regions"), guide = guide_legend(reverse = TRUE)) +
    labs(title = plot_title) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face="bold"),
          legend.position = "top")

  return(plt)
}
