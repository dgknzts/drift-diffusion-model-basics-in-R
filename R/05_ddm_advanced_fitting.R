# R/05_ddm_advanced_fitting.R

# Load required libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required but not available")
}

#' Calculate Binned RT Proportions for DDM Fitting
#'
#' Calculates the proportion of responses falling into predefined RT bins
#' for correct and error choices.
#'
#' @param data A data frame with 'choice' and 'rt' columns.
#' @param rt_bins A numeric vector defining the boundaries of RT bins (e.g., `seq(0, 3, by = 0.1)`).
#' @param correct_choice_value Value for correct/upper response. Default: `1`.
#' @param error_choice_value Value for error/lower response. Default: `0`.
#' @param min_trials_per_bin_for_calc If a bin has fewer than this many trials in the target data,
#'   it might be excluded from likelihood calculation to avoid instability, or given less weight.
#'   For simplicity now, we will use all bins but this is a consideration.
#'
#' @return A data frame with columns: response_type ("Correct", "Error"),
#'   rt_bin_midpoint, observed_n (count in bin), observed_prop (proportion in bin).
#' @export
calculate_binned_rt_proportions <- function(data,
                                            rt_bins,
                                            correct_choice_value = 1,
                                            error_choice_value = 0) {
  if (!all(c("rt", "choice") %in% names(data))) {
    stop("Data must contain 'rt' and 'choice' columns.")
  }
  data_filtered <- dplyr::filter(data, !is.na(rt) & !is.na(choice) & is.finite(rt) & rt >= 0)
  if (nrow(data_filtered) == 0) {
    return(data.frame(response_type = character(), rt_bin_label = character(),
                      rt_bin_midpoint = numeric(), observed_n = integer(),
                      observed_prop = numeric(), stringsAsFactors = FALSE))
  }

  # Ensure rt_bins start at 0 if min(rt) can be 0
  if (rt_bins[1] > 0 && any(data_filtered$rt < rt_bins[1])) {
    rt_bins <- c(0, rt_bins)
  }
  # Ensure rt_bins cover the max RT
  max_rt_observed <- max(data_filtered$rt, na.rm = TRUE)
  if (max_rt_observed > tail(rt_bins, 1)) {
    rt_bins <- c(rt_bins, max_rt_observed + (rt_bins[2]-rt_bins[1])) # Add one more bin
  }
  rt_bins <- unique(sort(rt_bins)) # Ensure sorted and unique

  binned_data_list <- list()

  # Correct responses
  correct_trials <- data_filtered %>% dplyr::filter(choice == correct_choice_value)
  if (nrow(correct_trials) > 0) {
    correct_cuts <- cut(correct_trials$rt, breaks = rt_bins, include.lowest = TRUE, right = FALSE)
    correct_counts <- as.data.frame(table(correct_cuts), stringsAsFactors = FALSE)
    colnames(correct_counts) <- c("rt_bin_label", "observed_n")
    correct_counts$response_type <- "Correct"
    correct_counts$observed_prop <- correct_counts$observed_n / nrow(correct_trials) # Prop within corrects
    binned_data_list[["correct"]] <- correct_counts
  }

  # Error responses
  error_trials <- data_filtered %>% dplyr::filter(choice == error_choice_value)
  if (nrow(error_trials) > 0) {
    error_cuts <- cut(error_trials$rt, breaks = rt_bins, include.lowest = TRUE, right = FALSE)
    error_counts <- as.data.frame(table(error_cuts), stringsAsFactors = FALSE)
    colnames(error_counts) <- c("rt_bin_label", "observed_n")
    error_counts$response_type <- "Error"
    error_counts$observed_prop <- error_counts$observed_n / nrow(error_trials) # Prop within errors
    binned_data_list[["error"]] <- error_counts
  }

  if(length(binned_data_list) == 0) {
    return(data.frame(response_type = character(), rt_bin_label = character(),
                      rt_bin_midpoint = numeric(), observed_n = integer(),
                      observed_prop = numeric(), stringsAsFactors = FALSE))
  }

  binned_df <- dplyr::bind_rows(binned_data_list)

  # Add midpoints for easier plotting/matching (optional)
  # More robust parsing of bin labels to avoid NA warnings
  binned_df$rt_bin_midpoint <- sapply(binned_df$rt_bin_label, function(label) {
    # Try to extract numbers from labels like "[0,0.15)" or "[0.15,0.3)"
    label_clean <- as.character(label)
    
    # More robust regex to extract decimal numbers
    numbers <- suppressWarnings(as.numeric(unlist(regmatches(label_clean, gregexpr("[0-9]*\\.?[0-9]+", label_clean)))))
    
    # Remove any NAs and take first two valid numbers
    numbers <- numbers[!is.na(numbers)]
    
    # If we found at least 2 numbers, return their midpoint
    if(length(numbers) >= 2) {
      return((numbers[1] + numbers[2]) / 2)
    } else if(length(numbers) == 1) {
      # If only one number found, use it as approximation
      return(numbers[1])
    } else {
      # If no numbers found, return 0 (will be handled below)
      return(0)
    }
  })
  
  return(binned_df)
}


#' DDM Objective Function using Binned Log-Likelihood
#'
#' Calculates the negative log-likelihood of observing target binned RT proportions
#' given DDM parameters. Aims to be MINIMIZED by optim().
#'
#' @param params_to_test A named numeric vector of DDM parameters to simulate with.
#' @param target_binned_props A data frame of target binned RT proportions,
#'   as produced by `calculate_binned_rt_proportions`. Must include columns:
#'   `response_type`, `rt_bin_label`, `observed_n` (from target data).
#' @param param_names_optim Names of parameters in `params_to_test`.
#' @param n_sim_per_eval Number of DDM trials to simulate per evaluation.
#' @param fixed_params List of DDM parameters held constant.
#' @param rt_bins Numeric vector defining RT bin boundaries (must match those used for target_binned_props).
#' @param constrain_z_to_a_div_2 Logical, if TRUE, mean_z is set to a/2.
#' @param small_constant_for_ll A small constant added to predicted counts to avoid log(0). Default 1e-6.
#' @param verbose Logical.
#' @param debug Logical.
#'
#' @return Negative log-likelihood value (to be minimized).
#' @export
ddm_binned_likelihood_objective <- function(params_to_test,
                                            target_binned_props, # This now includes target Ns
                                            param_names_optim,
                                            n_sim_per_eval = 2000,
                                            fixed_params = list(s = 0.1, dt = 0.001),
                                            rt_bins,
                                            constrain_z_to_a_div_2 = FALSE,
                                            small_constant_for_ll = 1e-6, # To avoid log(0)
                                            verbose = FALSE,
                                            debug = FALSE) {

  if(debug) cat("DEBUG: Starting objective function with params:", paste(round(params_to_test, 4), collapse=", "), "\n")

  current_iter_params <- fixed_params
  for(i in seq_along(param_names_optim)){
    current_iter_params[[param_names_optim[i]]] <- params_to_test[i]
  }
  
  if(debug) cat("DEBUG: Current params after assignment:", paste(names(current_iter_params), "=", round(unlist(current_iter_params), 4), collapse=", "), "\n")
  
  if (constrain_z_to_a_div_2) {
    if (!"a" %in% names(current_iter_params)) {
      if(debug) cat("DEBUG: ERROR - 'a' not found for z constraint\n")
      stop("'a' must be defined if constraining z.")
    }
    current_iter_params[["mean_z"]] <- current_iter_params[["a"]] / 2
    if(debug) cat("DEBUG: Set mean_z =", current_iter_params[["mean_z"]], "\n")
  }
  
  # Validate required parameters for simulation
  required_params <- c("mean_v", "a", "mean_z", "s", "dt", "mean_ter")
  missing_params <- required_params[!required_params %in% names(current_iter_params)]
  if(length(missing_params) > 0) {
    if(debug) cat("DEBUG: Missing required parameters:", paste(missing_params, collapse=", "), "\n")
    return(1e7)
  }
  
  # Set default values for optional variability parameters if not present
  if(!"sv" %in% names(current_iter_params)) current_iter_params[["sv"]] <- 0
  if(!"sz" %in% names(current_iter_params)) current_iter_params[["sz"]] <- 0
  if(!"st0" %in% names(current_iter_params)) current_iter_params[["st0"]] <- 0

  # Separate simulation parameters from binning parameters
  sim_param_names <- c("n_trials", "mean_v", "a", "mean_z", "s", "dt", "mean_ter", "sv", "sz", "st0", "max_decision_time")
  sim_args <- current_iter_params[names(current_iter_params) %in% sim_param_names]
  sim_args[["n_trials"]] <- n_sim_per_eval  # Add n_trials
  
  if(debug) cat("DEBUG: Simulation args:", paste(names(sim_args), "=", round(unlist(sim_args), 4), collapse=", "), "\n")
  
  current_sim_data <- tryCatch({
    do.call(simulate_diffusion_experiment_variable, sim_args)
  }, error = function(e) { 
    if(debug) cat("DEBUG: Simulation failed with error:", conditionMessage(e), "\n")
    return(NULL) 
  })

  if (is.null(current_sim_data) || nrow(current_sim_data) == 0) {
    if(debug) cat("DEBUG: No simulation data returned, returning large error\n")
    return(1e7) # Large error
  }
  
  if(debug) cat("DEBUG: Simulation successful, got", nrow(current_sim_data), "trials\n")

  # Calculate binned proportions for the *current simulated data*
  sim_binned_props <- calculate_binned_rt_proportions(
    current_sim_data,
    rt_bins = rt_bins, # Use same bins as target
    correct_choice_value = fixed_params$correct_choice_value %||% 1, # Get from fixed or default
    error_choice_value = fixed_params$error_choice_value %||% 0
  )

  # Merge with target counts to ensure alignment and get target N for each bin
  # Target_binned_props should have: response_type, rt_bin_label, observed_n (from target)
  # Sim_binned_props has: response_type, rt_bin_label, observed_prop (from sim)

  # Calculate overall P(Correct) and P(Error) from simulation
  sim_n_total <- nrow(dplyr::filter(current_sim_data, !is.na(choice) & !is.na(rt) & is.finite(rt)))
  sim_p_correct <- sum(current_sim_data$choice == (fixed_params$correct_choice_value %||% 1), na.rm=TRUE) / sim_n_total
  sim_p_error   <- sum(current_sim_data$choice == (fixed_params$error_choice_value %||% 0), na.rm=TRUE) / sim_n_total

  # If sim_p_correct or sim_p_error is 0, and target has counts there, this is bad.
  # Add small constant to avoid p=0 if target has counts.
  if(sim_p_correct == 0 && any(target_binned_props$response_type == "Correct" & target_binned_props$observed_n >0)) sim_p_correct = small_constant_for_ll
  if(sim_p_error == 0 && any(target_binned_props$response_type == "Error" & target_binned_props$observed_n >0)) sim_p_error = small_constant_for_ll


  # Join target Ns with simulated proportions
  comparison_df <- target_binned_props %>%
    select(response_type, rt_bin_label, target_n = observed_n) %>% # observed_n from target
    left_join(sim_binned_props %>% select(response_type, rt_bin_label, sim_prop_within_resp = observed_prop),
              by = c("response_type", "rt_bin_label")) %>%
    mutate(sim_prop_within_resp = ifelse(is.na(sim_prop_within_resp), 0, sim_prop_within_resp), # Bins with 0 sim count
           # Calculate overall probability for this bin
           sim_prob_overall = ifelse(response_type == "Correct",
                                     sim_prop_within_resp * sim_p_correct,
                                     sim_prop_within_resp * sim_p_error)
    )

  # Add small constant to predicted probabilities to avoid log(0)
  comparison_df$sim_prob_overall <- pmax(comparison_df$sim_prob_overall, small_constant_for_ll)

  # Ensure probabilities are properly normalized (should sum to 1)
  total_prob <- sum(comparison_df$sim_prob_overall, na.rm = TRUE)
  if(abs(total_prob - 1.0) > 0.1) { # If far from 1, normalize
    comparison_df$sim_prob_overall <- comparison_df$sim_prob_overall / total_prob
    # Re-apply minimum to avoid log(0) after normalization
    comparison_df$sim_prob_overall <- pmax(comparison_df$sim_prob_overall, small_constant_for_ll)
  }

  # Calculate log-likelihood (multinomial)
  # LL = sum(N_observed_in_bin * log(P_predicted_for_bin))
  # We want to MINIMIZE negative LL
  neg_log_likelihood <- -sum(comparison_df$target_n * log(comparison_df$sim_prob_overall), na.rm = TRUE)

  # Penalty for not producing choices of a type present in target
  target_has_corrects <- any(target_binned_props$response_type == "Correct" & target_binned_props$observed_n > 0)
  target_has_errors   <- any(target_binned_props$response_type == "Error" & target_binned_props$observed_n > 0)

  if(target_has_corrects && sim_p_correct < small_constant_for_ll*10) neg_log_likelihood <- neg_log_likelihood + 1000
  if(target_has_errors && sim_p_error < small_constant_for_ll*10) neg_log_likelihood <- neg_log_likelihood + 1000


  if (verbose) {
    cat("Iter: Params = [", paste(sprintf("%.3f", params_to_test), collapse = ", "),
        "], NegLL = ", sprintf("%.4f", neg_log_likelihood), "\n")
  }
  if(!is.finite(neg_log_likelihood)) return(1e7) # Return large error if not finite

  return(neg_log_likelihood)
}

# Helper for %||% (from rlang, but define here if not using rlang)
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Run DDM Optimization Multiple Times with Random Starts
#'
#' @param n_starts Integer. Number of random starting points for optimization.
#' @param target_stats_or_binned_props Target statistics (either from `calculate_ddm_summary_stats`
#'   or `calculate_binned_rt_proportions` depending on objective_fn_name).
#' @param param_names_optim Character vector of parameter names to optimize.
#' @param initial_guesses_means Named numeric vector for mean of initial guess distribution.
#' @param initial_guesses_sds Named numeric vector for SD of initial guess distribution.
#' @param lower_bounds Named numeric vector of lower bounds for parameters.
#' @param upper_bounds Named numeric vector of upper bounds for parameters.
#' @param objective_fn_name Character. Name of the objective function to use
#'   (e.g., "ddm_objective_function" or "ddm_binned_likelihood_objective").
#' @param n_sim_per_eval Integer. Trials per objective function evaluation.
#' @param fixed_params List of fixed DDM parameters.
#' @param optim_maxit Integer. Max iterations for `optim()`.
#' @param ... Additional arguments to pass to the objective function (e.g., `rt_bins`, `constrain_z_to_a_div_2`).
#'
#' @return A list containing the best `optim` result and a data frame of all runs.
#' @export
run_ddm_optimization_multi_start <- function(n_starts = 5,
                                             target_stats_or_binned_props,
                                             param_names_optim,
                                             initial_guesses_means, # e.g., c(mean_v=0.1, a=1.0)
                                             initial_guesses_sds,   # e.g., c(mean_v=0.05, a=0.2)
                                             lower_bounds,
                                             upper_bounds,
                                             objective_fn_name = "ddm_binned_likelihood_objective",
                                             n_sim_per_eval = 1000,
                                             fixed_params = list(s = 0.1, dt = 0.001),
                                             optim_maxit = 100,
                                             ...) { # Pass ... to objective function

  all_optim_results <- list()
  best_value <- Inf
  best_result_obj <- NULL

  objective_function <- get(objective_fn_name) # Get the actual function

  # Ensure names align for means, sds, bounds, and optim params
  if (!all(param_names_optim %in% names(initial_guesses_means)) ||
      !all(param_names_optim %in% names(initial_guesses_sds)) ||
      !all(param_names_optim %in% names(lower_bounds)) ||
      !all(param_names_optim %in% names(upper_bounds))) {
    stop("Parameter names in param_names_optim must match names in initial_guesses_means/sds and bounds.")
  }

  # Reorder if necessary, though optim takes named 'par'
  initial_guesses_means <- initial_guesses_means[param_names_optim]
  initial_guesses_sds <- initial_guesses_sds[param_names_optim]
  lower_bounds_ordered <- lower_bounds[param_names_optim]
  upper_bounds_ordered <- upper_bounds[param_names_optim]


  cat(paste("Starting optimization with", n_starts, "random starts...\n"))

  for (i_start in 1:n_starts) {
    cat(paste("--- Optimization Run:", i_start, "/", n_starts, "---\n"))
    # Generate random starting parameters based on means and sds, respecting bounds
    current_initial_guesses <- numeric(length(param_names_optim))
    names(current_initial_guesses) <- param_names_optim
    for (p_name in param_names_optim) {
      guess <- rnorm(1, mean = initial_guesses_means[[p_name]], sd = initial_guesses_sds[[p_name]])
      # Clip to bounds
      current_initial_guesses[[p_name]] <- max(lower_bounds_ordered[[p_name]],
                                               min(guess, upper_bounds_ordered[[p_name]]))
    }
    cat("Random Initial Guesses for this run:\n")
    print(round(current_initial_guesses,3))

    optim_run <- tryCatch({
      optim(
        par = current_initial_guesses,
        fn = objective_function,
        # Args for objective_function:
        target_binned_props = target_stats_or_binned_props, # Renamed argument
        param_names_optim = param_names_optim,
        n_sim_per_eval = n_sim_per_eval,
        fixed_params = fixed_params,
        ..., # Pass through other args like rt_bins, constrain_z
        # optim settings:
        method = "L-BFGS-B",
        lower = lower_bounds_ordered,
        upper = upper_bounds_ordered,
        control = list(maxit = optim_maxit, trace = 0, # trace=0 for less verbose during multi-start
                       parscale = abs(current_initial_guesses) + 0.05, factr=1e4) # Reduced factr for better precision
      )
    }, error = function(e) {
      cat("Error during optim run ", i_start, ": ", conditionMessage(e), "\n")
      return(list(value = Inf, par=current_initial_guesses, convergence = -1, message=conditionMessage(e))) # Return a dummy error list
    })

    all_optim_results[[i_start]] <- list(
      initial_params = current_initial_guesses,
      estimated_params = optim_run$par,
      obj_value = optim_run$value,
      convergence = optim_run$convergence,
      message = optim_run$message
    )

    if (is.finite(optim_run$value) && optim_run$value < best_value) {
      best_value <- optim_run$value
      best_result_obj <- optim_run
      cat(paste("New best objective value found:", round(best_value, 4), "\n"))
    }
    cat(paste("Run", i_start, "finished. Obj value:", round(optim_run$value, 4), "Convergence:", optim_run$convergence, "\n"))
  }

  # Prepare a summary data frame of all runs
  summary_df_all_runs <- do.call(rbind, lapply(seq_along(all_optim_results), function(i) {
    res <- all_optim_results[[i]]
    df <- data.frame(run = i, obj_value = res$obj_value, convergence = res$convergence)
    # Add initial and estimated params
    for (p_name in param_names_optim) {
      df[[paste0("initial_", p_name)]] <- res$initial_params[p_name]
      df[[paste0("estimated_", p_name)]] <- res$estimated_params[p_name]
    }
    df
  }))

  cat("\n--- Overall Best Result ---\n")
  if(!is.null(best_result_obj)){
    cat("Best Objective Value:", best_result_obj$value, "\n")
    cat("Best Parameters:\n")
    print(round(best_result_obj$par,4))
  } else {
    cat("No successful optimization run found.\n")
  }


  return(list(best_optim_result = best_result_obj, all_runs_summary = summary_df_all_runs))
}


#' Diagnostic Function for DDM Fitting Issues
#'
#' Tests the objective function with true parameters to help diagnose fitting problems.
#'
#' @param true_params Named list of true DDM parameters
#' @param target_binned_props Target binned proportions
#' @param objective_fn_name Name of objective function to test
#' @param n_sim_per_eval Number of simulation trials
#' @param fixed_params Fixed parameters
#' @param rt_bins RT bins for likelihood calculation
#' @param constrain_z_to_a_div_2 Whether to constrain mean_z = a/2
#' @param param_names_optim Names of parameters being optimized
#'
#' @return List with diagnostic information
#' @export
diagnose_ddm_fitting <- function(true_params,
                                 target_binned_props,
                                 objective_fn_name = "ddm_binned_likelihood_objective",
                                 n_sim_per_eval = 2000,
                                 fixed_params = list(dt = 0.001),
                                 rt_bins,
                                 constrain_z_to_a_div_2 = TRUE,
                                 param_names_optim) {
  
  cat("=== DDM Fitting Diagnostics ===\n")
  
  # Extract true values for optimized parameters
  true_values <- sapply(param_names_optim, function(p) true_params[[p]])
  
  cat("True parameter values:\n")
  print(round(true_values, 4))
  
  # Test objective function with true parameters
  objective_function <- get(objective_fn_name)
  
  true_obj_value <- objective_function(
    params_to_test = true_values,
    target_binned_props = target_binned_props,
    param_names_optim = param_names_optim,
    n_sim_per_eval = n_sim_per_eval,
    fixed_params = fixed_params,
    rt_bins = rt_bins,
    constrain_z_to_a_div_2 = constrain_z_to_a_div_2,
    verbose = FALSE,
    debug = TRUE  # Enable debugging for diagnostics
  )
  
  cat("\nObjective function value with TRUE parameters:", round(true_obj_value, 4), "\n")
  
  # Test with slightly perturbed parameters
  perturbed_values <- true_values * (1 + rnorm(length(true_values), 0, 0.1))
  perturbed_obj_value <- objective_function(
    params_to_test = perturbed_values,
    target_binned_props = target_binned_props,
    param_names_optim = param_names_optim,
    n_sim_per_eval = n_sim_per_eval,
    fixed_params = fixed_params,
    rt_bins = rt_bins,
    constrain_z_to_a_div_2 = constrain_z_to_a_div_2,
    verbose = FALSE,
    debug = TRUE  # Enable debugging for diagnostics
  )
  
  cat("Objective function value with PERTURBED parameters:", round(perturbed_obj_value, 4), "\n")
  cat("Difference (should be positive if objective function working):", round(perturbed_obj_value - true_obj_value, 4), "\n")
  
  # Check target data properties
  cat("\nTarget data summary:\n")
  cat("Number of bins:", nrow(target_binned_props), "\n")
  cat("Response types:", unique(target_binned_props$response_type), "\n")
  cat("Total target counts:", sum(target_binned_props$observed_n), "\n")
  
  return(list(
    true_obj_value = true_obj_value,
    perturbed_obj_value = perturbed_obj_value,
    difference = perturbed_obj_value - true_obj_value,
    target_summary = target_binned_props
  ))
}

# Add improved fitting function at the end

#' Improved DDM Objective Function with Independent mean_z
#'
#' This version treats mean_z as an independent parameter rather than constraining it to a/2.
#' This allows for better parameter recovery and avoids circular dependencies.
#'
#' @param params_to_test A named numeric vector of DDM parameters to simulate with.
#' @param target_binned_props A data frame of target binned RT proportions.
#' @param param_names_optim Names of parameters in `params_to_test`.
#' @param n_sim_per_eval Number of DDM trials to simulate per evaluation.
#' @param fixed_params List of DDM parameters held constant.
#' @param rt_bins Numeric vector defining RT bin boundaries.
#' @param small_constant_for_ll A small constant added to predicted counts to avoid log(0).
#' @param verbose Logical.
#' @param debug Logical.
#'
#' @return Negative log-likelihood value (to be minimized).
#' @export
ddm_improved_likelihood_objective <- function(params_to_test,
                                             target_binned_props,
                                             param_names_optim,
                                             n_sim_per_eval = 2000,
                                             fixed_params = list(s = 0.1, dt = 0.001),
                                             rt_bins,
                                             small_constant_for_ll = 1e-6,
                                             verbose = FALSE,
                                             debug = FALSE) {

  if(debug) cat("DEBUG: Starting improved objective function with params:", paste(round(params_to_test, 4), collapse=", "), "\n")

  current_iter_params <- fixed_params
  for(i in seq_along(param_names_optim)){
    current_iter_params[[param_names_optim[i]]] <- params_to_test[i]
  }
  
  if(debug) cat("DEBUG: Current params after assignment:", paste(names(current_iter_params), "=", round(unlist(current_iter_params), 4), collapse=", "), "\n")
  
  # Validate that mean_z is within bounds (0, a)
  if("mean_z" %in% names(current_iter_params) && "a" %in% names(current_iter_params)) {
    if(current_iter_params[["mean_z"]] <= 0 || current_iter_params[["mean_z"]] >= current_iter_params[["a"]]) {
      if(debug) cat("DEBUG: mean_z outside valid bounds relative to a\n")
      return(1e7)
    }
  }
  
  # Validate sz parameter doesn't exceed sensible bounds
  if("sz" %in% names(current_iter_params) && "a" %in% names(current_iter_params)) {
    if(current_iter_params[["sz"]] >= current_iter_params[["a"]] * 0.8) {
      if(debug) cat("DEBUG: sz too large relative to a\n")
      return(1e7)
    }
  }
  
  # Validate required parameters for simulation
  required_params <- c("mean_v", "a", "mean_z", "s", "dt", "mean_ter")
  missing_params <- required_params[!required_params %in% names(current_iter_params)]
  if(length(missing_params) > 0) {
    if(debug) cat("DEBUG: Missing required parameters:", paste(missing_params, collapse=", "), "\n")
    return(1e7)
  }
  
  # Set default values for optional variability parameters if not present
  if(!"sv" %in% names(current_iter_params)) current_iter_params[["sv"]] <- 0
  if(!"sz" %in% names(current_iter_params)) current_iter_params[["sz"]] <- 0
  if(!"st0" %in% names(current_iter_params)) current_iter_params[["st0"]] <- 0

  # Separate simulation parameters from binning parameters
  sim_param_names <- c("n_trials", "mean_v", "a", "mean_z", "s", "dt", "mean_ter", "sv", "sz", "st0", "max_decision_time")
  sim_args <- current_iter_params[names(current_iter_params) %in% sim_param_names]
  sim_args[["n_trials"]] <- n_sim_per_eval
  
  if(debug) cat("DEBUG: Simulation args:", paste(names(sim_args), "=", round(unlist(sim_args), 4), collapse=", "), "\n")
  
  current_sim_data <- tryCatch({
    do.call(simulate_diffusion_experiment_variable, sim_args)
  }, error = function(e) { 
    if(debug) cat("DEBUG: Simulation failed with error:", conditionMessage(e), "\n")
    return(NULL) 
  })

  if (is.null(current_sim_data) || nrow(current_sim_data) == 0) {
    if(debug) cat("DEBUG: No simulation data returned, returning large error\n")
    return(1e7)
  }
  
  if(debug) cat("DEBUG: Simulation successful, got", nrow(current_sim_data), "trials\n")

  # Calculate binned proportions for the current simulated data
  sim_binned_props <- calculate_binned_rt_proportions(
    current_sim_data,
    rt_bins = rt_bins,
    correct_choice_value = fixed_params$correct_choice_value %||% 1,
    error_choice_value = fixed_params$error_choice_value %||% 0
  )

  # Calculate overall P(Correct) and P(Error) from simulation
  sim_n_total <- nrow(dplyr::filter(current_sim_data, !is.na(choice) & !is.na(rt) & is.finite(rt)))
  sim_p_correct <- sum(current_sim_data$choice == (fixed_params$correct_choice_value %||% 1), na.rm=TRUE) / sim_n_total
  sim_p_error   <- sum(current_sim_data$choice == (fixed_params$error_choice_value %||% 0), na.rm=TRUE) / sim_n_total

  # Add small constant to avoid p=0 if target has counts
  if(sim_p_correct == 0 && any(target_binned_props$response_type == "Correct" & target_binned_props$observed_n >0)) sim_p_correct = small_constant_for_ll
  if(sim_p_error == 0 && any(target_binned_props$response_type == "Error" & target_binned_props$observed_n >0)) sim_p_error = small_constant_for_ll

  # Join target Ns with simulated proportions
  comparison_df <- target_binned_props %>%
    dplyr::select(response_type, rt_bin_label, target_n = observed_n) %>%
    dplyr::left_join(sim_binned_props %>% dplyr::select(response_type, rt_bin_label, sim_prop_within_resp = observed_prop),
              by = c("response_type", "rt_bin_label")) %>%
    dplyr::mutate(sim_prop_within_resp = ifelse(is.na(sim_prop_within_resp), 0, sim_prop_within_resp),
           sim_prob_overall = ifelse(response_type == "Correct",
                                     sim_prop_within_resp * sim_p_correct,
                                     sim_prop_within_resp * sim_p_error)
    )

  # Add small constant to predicted probabilities to avoid log(0)
  comparison_df$sim_prob_overall <- pmax(comparison_df$sim_prob_overall, small_constant_for_ll)

  # Ensure probabilities are properly normalized
  total_prob <- sum(comparison_df$sim_prob_overall, na.rm = TRUE)
  if(abs(total_prob - 1.0) > 0.1) {
    comparison_df$sim_prob_overall <- comparison_df$sim_prob_overall / total_prob
    comparison_df$sim_prob_overall <- pmax(comparison_df$sim_prob_overall, small_constant_for_ll)
  }

  # Calculate log-likelihood (multinomial)
  neg_log_likelihood <- -sum(comparison_df$target_n * log(comparison_df$sim_prob_overall), na.rm = TRUE)

  # Penalty for not producing choices of a type present in target
  target_has_corrects <- any(target_binned_props$response_type == "Correct" & target_binned_props$observed_n > 0)
  target_has_errors   <- any(target_binned_props$response_type == "Error" & target_binned_props$observed_n > 0)

  if(target_has_corrects && sim_p_correct < small_constant_for_ll*10) neg_log_likelihood <- neg_log_likelihood + 1000
  if(target_has_errors && sim_p_error < small_constant_for_ll*10) neg_log_likelihood <- neg_log_likelihood + 1000

  if (verbose) {
    cat("Iter: Params = [", paste(sprintf("%.3f", params_to_test), collapse = ", "),
        "], NegLL = ", sprintf("%.4f", neg_log_likelihood), "\n")
  }
  if(!is.finite(neg_log_likelihood)) return(1e7)

  return(neg_log_likelihood)
}


#' Create Improved RT Bins
#'
#' Creates finer RT bins for better likelihood estimation
#'
#' @param target_data Data frame with rt column
#' @param bin_width Bin width in seconds (default: 0.05)
#' @param min_rt_quantile Lower quantile for bin range (default: 0.01)
#' @param max_rt_quantile Upper quantile for bin range (default: 0.99)
#'
#' @return Numeric vector of bin boundaries
#' @export
create_improved_rt_bins <- function(target_data, bin_width = 0.05, min_rt_quantile = 0.01, max_rt_quantile = 0.99) {
  valid_rts <- target_data$rt[!is.na(target_data$rt) & is.finite(target_data$rt) & target_data$rt > 0]
  
  if(length(valid_rts) == 0) {
    return(seq(0, 2, by = bin_width))
  }
  
  min_rt <- max(0, quantile(valid_rts, min_rt_quantile, na.rm = TRUE))
  max_rt <- quantile(valid_rts, max_rt_quantile, na.rm = TRUE)
  
  # Extend slightly beyond observed range
  min_rt <- max(0, min_rt - bin_width)
  max_rt <- max_rt + bin_width
  
  return(seq(min_rt, max_rt, by = bin_width))
}

#' Recommended DDM Parameter Recovery Function
#'
#' This is the recommended approach for DDM parameter fitting based on empirical testing.
#' It uses the improved likelihood objective with independent mean_z optimization,
#' fine RT bins, and optimized settings for better parameter recovery.
#'
#' @param target_data Data frame with 'choice' and 'rt' columns from your experiment
#' @param param_names_to_fit Character vector of parameter names to optimize. 
#'   Default includes mean_z independently: c("mean_v", "a", "mean_z", "s", "mean_ter", "sv", "sz", "st0")
#' @param n_starts Number of random optimization starts (default: 6)
#' @param n_sim_per_eval Number of simulations per objective function evaluation (default: 2500)
#' @param bin_width RT bin width in seconds (default: 0.05)
#' @param optim_maxit Maximum optimization iterations (default: 200)
#' @param dt Time step for simulation (default: 0.001)
#' @param verbose Whether to show optimization progress (default: TRUE)
#' @param custom_bounds Optional list with 'lower' and 'upper' named vectors for custom bounds
#' @param custom_initial_means Optional named vector for custom starting values
#'
#' @return List with optimization results and parameter recovery summary
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit DDM to your data
#' results <- fit_ddm_recommended(your_data, verbose = TRUE)
#' 
#' # View results
#' print(results$recovery_summary)
#' print(results$best_parameters)
#' }
fit_ddm_recommended <- function(target_data,
                               param_names_to_fit = c("mean_v", "a", "mean_z", "s", "mean_ter", "sv", "sz", "st0"),
                               n_starts = 6,
                               n_sim_per_eval = 2500,
                               bin_width = 0.05,
                               optim_maxit = 200,
                               dt = 0.001,
                               verbose = TRUE,
                               custom_bounds = NULL,
                               custom_initial_means = NULL) {
  
  if (verbose) cat("=== DDM Parameter Fitting (Recommended Approach) ===\n")
  
  # Validate data
  if (!all(c("choice", "rt") %in% names(target_data))) {
    stop("target_data must contain 'choice' and 'rt' columns")
  }
  
  # Create optimized RT bins
  rt_bins <- create_improved_rt_bins(target_data, bin_width = bin_width)
  target_binned <- calculate_binned_rt_proportions(target_data, rt_bins = rt_bins)
  
  if (verbose) {
    cat("Data: ", nrow(target_data), "trials\n")
    cat("RT bins: ", length(rt_bins)-1, "bins (", bin_width, "s resolution)\n")
    cat("Parameters to fit: ", paste(param_names_to_fit, collapse = ", "), "\n")
  }
  
  # Set up default bounds and starting values
  if (is.null(custom_bounds)) {
    lower_bounds <- c(
      mean_v = -0.5, a = 0.3, mean_z = 0.05, s = 0.05, mean_ter = 0.05,
      sv = 0.0, sz = 0.0, st0 = 0.0
    )[param_names_to_fit]
    
    upper_bounds <- c(
      mean_v = 0.6, a = 1.5, mean_z = 0.8, s = 0.3, mean_ter = 0.25,
      sv = 0.4, sz = 0.1, st0 = 0.08
    )[param_names_to_fit]
  } else {
    lower_bounds <- custom_bounds$lower[param_names_to_fit]
    upper_bounds <- custom_bounds$upper[param_names_to_fit]
  }
  
  if (is.null(custom_initial_means)) {
    initial_means <- c(
      mean_v = 0.2, a = 0.8, mean_z = 0.4, s = 0.1, mean_ter = 0.1,
      sv = 0.1, sz = 0.03, st0 = 0.02
    )[param_names_to_fit]
  } else {
    initial_means <- custom_initial_means[param_names_to_fit]
  }
  
  # Set conservative SDs for initial guess generation
  initial_sds <- pmax(abs(initial_means) * 0.15, 0.01)
  names(initial_sds) <- param_names_to_fit
  
  if (verbose) cat("Starting optimization with", n_starts, "random starts...\n")
  
  # Run optimization
  results <- run_ddm_optimization_multi_start(
    n_starts = n_starts,
    target_stats_or_binned_props = target_binned,
    param_names_optim = param_names_to_fit,
    initial_guesses_means = initial_means,
    initial_guesses_sds = initial_sds,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    objective_fn_name = "ddm_improved_likelihood_objective",
    n_sim_per_eval = n_sim_per_eval,
    fixed_params = list(dt = dt, correct_choice_value = 1, error_choice_value = 0),
    optim_maxit = optim_maxit,
    rt_bins = rt_bins,
    verbose = verbose
  )
  
  # Create summary
  if (!is.null(results$best_optim_result)) {
    best_params <- results$best_optim_result$par
    names(best_params) <- param_names_to_fit
    
    convergence_status <- switch(as.character(results$best_optim_result$convergence),
                                "0" = "Successful convergence",
                                "1" = "Maximum iterations reached", 
                                "10" = "Degeneracy in optimization",
                                "51" = "Warning from L-BFGS-B",
                                "52" = "Error from L-BFGS-B",
                                paste("Convergence code:", results$best_optim_result$convergence))
    
    summary_info <- list(
      best_parameters = best_params,
      convergence_status = convergence_status,
      final_objective_value = results$best_optim_result$value,
      optimization_successful = results$best_optim_result$convergence %in% c(0, 1),
      all_runs_summary = results$all_runs_summary,
      fitting_approach = "Independent mean_z (recommended)",
      rt_bins_used = length(rt_bins) - 1,
      simulations_per_evaluation = n_sim_per_eval
    )
    
    if (verbose) {
      cat("\n=== FITTING RESULTS ===\n")
      cat("Status:", convergence_status, "\n")
      cat("Final objective value:", round(summary_info$final_objective_value, 2), "\n")
      cat("Best parameters:\n")
      print(round(best_params, 4))
    }
    
  } else {
    summary_info <- list(
      best_parameters = NULL,
      convergence_status = "All optimization runs failed",
      optimization_successful = FALSE
    )
    if (verbose) cat("ERROR: All optimization runs failed\n")
  }
  
  return(summary_info)
}
