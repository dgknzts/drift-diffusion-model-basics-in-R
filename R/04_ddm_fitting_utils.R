# R/04_ddm_fitting_utils.R

#' Calculate Summary Statistics for DDM Fitting
#'
#' Calculates choice proportions and RT quantiles for correct and error responses,
#' suitable for use in an objective function for DDM parameter estimation.
#'
#' @param data A data frame with 'choice' and 'rt' columns.
#' @param quantiles_to_calc Numeric vector of probabilities for RT quantiles.
#'   Default: `c(0.1, 0.3, 0.5, 0.7, 0.9)`.
#' @param correct_choice_value Value for correct/upper response. Default: `1`.
#' @param error_choice_value Value for error/lower response. Default: `0`.
#' @param min_n_for_quantiles Minimum number of observations per response type
#'   to calculate quantiles (otherwise returns NAs for quantiles). Default: `10`.
#'
#' @return A named numeric vector of summary statistics. Names will be like
#'   "p_correct", "rt_correct_q10", "rt_error_q50", etc.
#' @export
calculate_ddm_summary_stats <- function(data,
                                        quantiles_to_calc = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                        correct_choice_value = 1,
                                        error_choice_value = 0,
                                        min_n_for_quantiles = 10) {

  if (!all(c("rt", "choice") %in% names(data))) {
    stop("Data must contain 'rt' and 'choice' columns.")
  }
  data_filtered <- dplyr::filter(data, !is.na(rt) & !is.na(choice) & is.finite(rt))
  if (nrow(data_filtered) == 0) {
    warning("No valid trials in data for summary stats.")
    # Return a structure of NAs
    out_names <- c("p_correct", "n_correct", "n_error")
    for (q in quantiles_to_calc) {
      out_names <- c(out_names, paste0("rt_correct_q", q * 100), paste0("rt_error_q", q * 100))
    }
    stats <- rep(NA_real_, length(out_names))
    names(stats) <- out_names
    return(stats)
  }

  n_total_valid <- nrow(data_filtered)

  # Correct responses
  correct_trials <- data_filtered %>% dplyr::filter(choice == correct_choice_value)
  n_correct <- nrow(correct_trials)
  p_correct <- n_correct / n_total_valid
  rt_correct_quantiles <- rep(NA_real_, length(quantiles_to_calc))
  if (n_correct >= min_n_for_quantiles) {
    rt_correct_quantiles <- quantile(correct_trials$rt, probs = quantiles_to_calc, na.rm = TRUE, type = 7)
  }
  names(rt_correct_quantiles) <- paste0("rt_correct_q", quantiles_to_calc * 100)

  # Error responses
  error_trials <- data_filtered %>% dplyr::filter(choice == error_choice_value)
  n_error <- nrow(error_trials)
  # p_error <- n_error / n_total_valid # Not strictly needed if we have p_correct
  rt_error_quantiles <- rep(NA_real_, length(quantiles_to_calc))
  if (n_error >= min_n_for_quantiles) {
    rt_error_quantiles <- quantile(error_trials$rt, probs = quantiles_to_calc, na.rm = TRUE, type = 7)
  }
  names(rt_error_quantiles) <- paste0("rt_error_q", quantiles_to_calc * 100)

  # Combine all stats
  all_stats <- c(
    p_correct = p_correct,
    n_correct = n_correct,
    n_error = n_error,
    rt_correct_quantiles,
    rt_error_quantiles
  )
  return(all_stats)
}

#' DDM Objective Function for Parameter Estimation
#'
#' Simulates DDM data with a given set of parameters, calculates summary
#' statistics, and returns a discrepancy score (e.g., sum of squared errors)
#' against a target set of summary statistics.
#'
#' @param params_to_test A named numeric vector of DDM parameters to simulate with.
#'   Names must match arguments of `simulate_diffusion_experiment_variable`
#'   (e.g., `c(mean_v = 0.1, a = 1.0, ...)`).
#' @param target_stats A named numeric vector of target summary statistics,
#'   as produced by `calculate_ddm_summary_stats`.
#' @param n_sim_per_eval Integer. Number of trials to simulate for each
#'   evaluation of the objective function. Default `2000`.
#' @param fixed_params A named list of DDM parameters that are held constant
#'   during optimization (e.g., `s`, `dt`).
#' @param quantiles_for_obj_func Numeric vector of quantiles used in `target_stats`
#'   and for calculating stats from simulated data.
#' @param weight_p_correct Numeric. Weight for the choice proportion error. Default `1`.
#'   Can be increased if fitting accuracy is more important.
#' @param verbose Logical. If TRUE, print iteration details. Default `FALSE`.
#'
#' @return A single numeric value representing the discrepancy (error).
#' @export
ddm_objective_function <- function(params_to_test, # Vector of parameters being optimized
                                   target_stats,
                                   param_names_optim, # Names of parameters in params_to_test
                                   n_sim_per_eval = 2000,
                                   fixed_params = list(s = 0.1, dt = 0.001), # Default fixed
                                   quantiles_for_obj_func = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                   weight_p_correct = 1.0, # Can increase if accuracy fit is poor
                                   verbose = FALSE,
                                   constrain_z_to_a_div_2 = FALSE) {

  current_iter_params <- fixed_params
  for(i in seq_along(param_names_optim)){
    current_iter_params[[param_names_optim[i]]] <- params_to_test[i]
  }

  # ---- NEW: Handle mean_z = a/2 constraint ----
  if (constrain_z_to_a_div_2) {
    if (!"a" %in% names(current_iter_params)) {
      stop("If 'constrain_z_to_a_div_2' is TRUE, 'a' must be either a fixed or optimized parameter.")
    }
    current_iter_params[["mean_z"]] <- current_iter_params[["a"]] / 2
    # Remove mean_z if it was accidentally passed in params_to_test or fixed_params
    # when this constraint is active, to avoid conflict.
    # However, param_names_optim should NOT include mean_z if this constraint is active.
  }

  # Ensure all necessary mean parameters are present for simulate_diffusion_experiment_variable
  # If some variability params are not being optimized, they default to 0 in the simulator
  # This assumes simulate_diffusion_experiment_variable handles missing sv, sz, st0 by defaulting them
  required_means <- c("mean_v", "a", "mean_z", "mean_ter")
  for(p_name in required_means){
    if(!p_name %in% names(current_iter_params)){
      # This case should ideally not happen if fixed_params and param_names_optim cover things
      # Or if the simulator has defaults for these. For safety, could stop or use a global default.
      stop(paste("Required parameter", p_name, "not found for simulation."))
    }
  }

  # Simulate data with current candidate parameters
  sim_args <- c(list(n_trials = n_sim_per_eval), current_iter_params)

  # Catch errors from simulation (e.g., bad parameter combinations)
  current_sim_data <- tryCatch({
    do.call(simulate_diffusion_experiment_variable, sim_args)
  }, error = function(e) {
    if(verbose) cat("Error in simulation with params:", paste(names(params_to_test), round(params_to_test,3), collapse=", "), "\n", conditionMessage(e), "\n")
    return(NULL) # Return NULL if simulation fails
  })

  if (is.null(current_sim_data) || nrow(current_sim_data) == 0) {
    return(1e6) # Return a large error if simulation failed or produced no data
  }

  # Calculate summary stats for the simulated data
  current_sim_stats <- calculate_ddm_summary_stats(
    current_sim_data,
    quantiles_to_calc = quantiles_for_obj_func
  )

  # Calculate discrepancy (e.g., Sum of Squared Errors - SSE)
  # Ensure stats are aligned and handle NAs
  # Only compare stats that are present in both and not NA
  common_names <- intersect(names(target_stats), names(current_sim_stats))

  target_common <- target_stats[common_names]
  sim_common <- current_sim_stats[common_names]

  valid_comparison_idx <- !is.na(target_common) & !is.na(sim_common)

  if(sum(valid_comparison_idx) == 0) {
    if(verbose) cat("No valid stats to compare for params:", paste(names(params_to_test), round(params_to_test,3), collapse=", "), "\n")
    return(1e6) # Large error if no overlap or all NAs
  }

  target_valid <- target_common[valid_comparison_idx]
  sim_valid <- sim_common[valid_comparison_idx]

  # Weight choice proportion error more if desired
  # This is a simple SSE. Chi-square is often better as it accounts for variance of stats.
  error_terms <- (target_valid - sim_valid)^2

  # Apply weight to choice proportion if it's part of the comparison
  if ("p_correct" %in% names(target_valid)) {
    p_correct_idx <- which(names(target_valid) == "p_correct")
    error_terms[p_correct_idx] <- error_terms[p_correct_idx] * weight_p_correct
    # Could also weight n_correct/n_error if using those
  }

  total_error <- sum(error_terms)

  # Penalty for too few error or correct responses if they are expected
  # This helps guide the optimizer if it lands in a region with no errors/corrects
  # when the target data has them.
  min_expected_n <- 5 # Arbitrary small number
  if(target_stats["n_correct"] > min_expected_n && current_sim_stats["n_correct"] < min_expected_n) total_error <- total_error + 100
  if(target_stats["n_error"] > min_expected_n && current_sim_stats["n_error"] < min_expected_n) total_error <- total_error + 100


  if (verbose) {
    cat("Iter: Params = [", paste(sprintf("%.3f", params_to_test), collapse = ", "),
        "], Error = ", sprintf("%.4f", total_error), "\n")
  }

  # optim minimizes, so return positive error
  return(total_error)
}
