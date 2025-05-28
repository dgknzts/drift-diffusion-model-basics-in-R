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

# Helper for %||% (from rlang, but define here if not using rlang)
`%||%` <- function(a, b) if (is.null(a)) b else a

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
                                            target_binned_props,
                                            param_names_optim,
                                            n_sim_per_eval = 2000,
                                            fixed_params = list(s = 0.1, dt = 0.001),
                                            rt_bins,
                                            constrain_z_to_a_div_2 = FALSE,
                                            small_constant_for_ll = 1e-6,
                                            verbose = FALSE,
                                            debug = FALSE) {

  # Quick parameter setup
  current_iter_params <- c(fixed_params, setNames(as.list(params_to_test), param_names_optim))
  
  if (constrain_z_to_a_div_2) {
    current_iter_params[["mean_z"]] <- current_iter_params[["a"]] / 2
  }
  
  # Validate required parameters
  required_params <- c("mean_v", "a", "mean_z", "s", "dt", "mean_ter")
  if (!all(required_params %in% names(current_iter_params))) {
    return(1e7)
  }
  
  # Set default values for optional parameters
  current_iter_params$sv <- current_iter_params$sv %||% 0
  current_iter_params$sz <- current_iter_params$sz %||% 0
  current_iter_params$st0 <- current_iter_params$st0 %||% 0

  # Simulate data
  sim_data <- tryCatch({
    do.call(simulate_diffusion_experiment_variable, 
            c(list(n_trials = n_sim_per_eval), 
              current_iter_params[names(current_iter_params) %in% 
                c("mean_v", "a", "mean_z", "s", "dt", "mean_ter", "sv", "sz", "st0", "max_decision_time")]))
  }, error = function(e) NULL)

  if (is.null(sim_data) || nrow(sim_data) == 0) return(1e7)

  # Calculate binned proportions efficiently
  sim_binned_props <- calculate_binned_rt_proportions(
    sim_data,
    rt_bins = rt_bins,
    correct_choice_value = fixed_params$correct_choice_value %||% 1,
    error_choice_value = fixed_params$error_choice_value %||% 0
  )

  # Calculate overall probabilities efficiently
  sim_n_total <- sum(!is.na(sim_data$choice) & !is.na(sim_data$rt) & is.finite(sim_data$rt))
  sim_p_correct <- sum(sim_data$choice == (fixed_params$correct_choice_value %||% 1), na.rm=TRUE) / sim_n_total
  sim_p_error <- sum(sim_data$choice == (fixed_params$error_choice_value %||% 0), na.rm=TRUE) / sim_n_total

  # Handle zero probabilities
  if(sim_p_correct == 0 && any(target_binned_props$response_type == "Correct" & target_binned_props$observed_n > 0)) {
    sim_p_correct <- small_constant_for_ll
  }
  if(sim_p_error == 0 && any(target_binned_props$response_type == "Error" & target_binned_props$observed_n > 0)) {
    sim_p_error <- small_constant_for_ll
  }

  # Efficient comparison calculation
  comparison_df <- merge(
    target_binned_props[, c("response_type", "rt_bin_label", "observed_n")],
    sim_binned_props[, c("response_type", "rt_bin_label", "observed_prop")],
    by = c("response_type", "rt_bin_label"),
    all.x = TRUE
  )
  
  # Calculate probabilities
  comparison_df$sim_prop_within_resp <- ifelse(is.na(comparison_df$observed_prop), 0, comparison_df$observed_prop)
  comparison_df$sim_prob_overall <- ifelse(
    comparison_df$response_type == "Correct",
    comparison_df$sim_prop_within_resp * sim_p_correct,
    comparison_df$sim_prop_within_resp * sim_p_error
  )

  # Ensure minimum probability and normalize
  comparison_df$sim_prob_overall <- pmax(comparison_df$sim_prob_overall, small_constant_for_ll)
  total_prob <- sum(comparison_df$sim_prob_overall, na.rm = TRUE)
  if(abs(total_prob - 1.0) > 0.1) {
    comparison_df$sim_prob_overall <- comparison_df$sim_prob_overall / total_prob
    comparison_df$sim_prob_overall <- pmax(comparison_df$sim_prob_overall, small_constant_for_ll)
  }

  # Calculate negative log-likelihood
  neg_log_likelihood <- -sum(comparison_df$observed_n * log(comparison_df$sim_prob_overall), na.rm = TRUE)

  # Add penalties for missing response types
  if(any(target_binned_props$response_type == "Correct" & target_binned_props$observed_n > 0) && 
     sim_p_correct < small_constant_for_ll*10) {
    neg_log_likelihood <- neg_log_likelihood + 1000
  }
  if(any(target_binned_props$response_type == "Error" & target_binned_props$observed_n > 0) && 
     sim_p_error < small_constant_for_ll*10) {
    neg_log_likelihood <- neg_log_likelihood + 1000
  }

  if (verbose) {
    cat("Iter: Params = [", paste(sprintf("%.3f", params_to_test), collapse = ", "),
        "], NegLL = ", sprintf("%.4f", neg_log_likelihood), "\n")
  }

  if(!is.finite(neg_log_likelihood)) return(1e7)
  return(neg_log_likelihood)
}

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
