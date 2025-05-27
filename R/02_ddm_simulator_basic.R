# R/02_ddm_simulator_basic.R

#' Validate DDM Parameters
#'
#' Internal function to validate DDM parameters and provide informative error messages.
#'
#' @param v Numeric. Drift rate.
#' @param a Numeric. Threshold separation.
#' @param z Numeric. Starting point.
#' @param s Numeric. Noise standard deviation.
#' @param dt Numeric. Time step.
#' @param ter Numeric. Non-decision time.
#' @param max_decision_time Numeric. Maximum decision time.
#' @param function_name Character. Name of calling function for error messages.
#'
#' @return NULL (invisible). Throws informative errors if validation fails.
#' @keywords internal
validate_ddm_parameters <- function(v, a, z, s = 0.1, dt = 0.001, 
                                   ter = 0.1, max_decision_time = 5.0,
                                   function_name = "DDM function") {
  
  # Check for required parameters
  if (missing(v)) stop(paste(function_name, ": 'v' (drift rate) is required."))
  if (missing(a)) stop(paste(function_name, ": 'a' (threshold) is required."))
  if (missing(z)) stop(paste(function_name, ": 'z' (starting point) is required."))
  
  # Validate parameter types and ranges
  if (!is.numeric(v) || length(v) != 1) {
    stop(paste(function_name, ": 'v' must be a single numeric value."))
  }
  
  if (!is.numeric(a) || length(a) != 1 || a <= 0) {
    stop(paste(function_name, ": 'a' must be a single positive numeric value."))
  }
  
  if (!is.numeric(z) || length(z) != 1 || z <= 0 || z >= a) {
    stop(paste(function_name, ": 'z' must be between 0 and 'a' (exclusive). Current: z =", z, ", a =", a))
  }
  
  if (!is.numeric(s) || length(s) != 1 || s <= 0) {
    stop(paste(function_name, ": 's' must be a single positive numeric value."))
  }
  
  if (!is.numeric(dt) || length(dt) != 1 || dt <= 0) {
    stop(paste(function_name, ": 'dt' must be a single positive numeric value."))
  }
  
  if (!is.numeric(ter) || length(ter) != 1 || ter < 0) {
    stop(paste(function_name, ": 'ter' must be a single non-negative numeric value."))
  }
  
  if (!is.numeric(max_decision_time) || length(max_decision_time) != 1 || max_decision_time <= 0) {
    stop(paste(function_name, ": 'max_decision_time' must be a single positive numeric value."))
  }
  
  # Warnings for potentially problematic parameter combinations
  if (dt > 0.1) {
    warning(paste(function_name, ": Large time step (dt =", dt, ") may affect simulation accuracy. Consider dt <= 0.01."))
  }
  
  if (abs(v) > 5) {
    warning(paste(function_name, ": Very large drift rate (|v| =", abs(v), ") may lead to very fast decisions."))
  }
  
  if (z/a < 0.1 || z/a > 0.9) {
    warning(paste(function_name, ": Extreme starting point bias (z/a =", round(z/a, 3), "). Consider 0.1 < z/a < 0.9."))
  }
  
  invisible(NULL)
}

#' Simulate a single trial of the Diffusion Decision Model (DDM).
#'
#' This function simulates the accumulation of evidence over time until one of
#' two decision boundaries (0 or 'a') is reached. The simulation uses an
#' Euler-Maruyama approximation for the stochastic differential equation.
#'
#' @param v Numeric. Drift rate. Represents the average rate of evidence
#'   accumulation towards the upper boundary 'a'. Negative values indicate
#'   drift towards the lower boundary 0.
#' @param a Numeric. Threshold separation. The upper decision boundary.
#'   The lower boundary is fixed at 0. Must be positive.
#' @param z Numeric. Starting point of evidence accumulation. Must be between
#'   0 and 'a' (exclusive, i.e., 0 < z < a). Represents initial bias.
#'   A common unbiased starting point is a/2.
#' @param s Numeric. Standard deviation of the within-trial normally distributed
#'   increments (noise). Often fixed to 0.1 in many DDM applications to make
#'   other parameters scale appropriately. Default is 0.1.
#' @param dt Numeric. Time step for the simulation (in seconds). Smaller values
#'   increase accuracy but also computation time. Default is 0.001.
#' @param ter Numeric. Non-decision time (in seconds). Represents time for
#'   perceptual encoding and response execution, added to the decision time.
#'   Default is 0.1.
#' @param max_decision_time Numeric. Maximum allowed time for the decision process
#'   (in seconds) before the trial is considered a non-decision. This prevents
#'   infinite loops for very slow or non-terminating processes. Default is 5.0.
#'
#' @return A list containing:
#'   \item{choice}{Integer: 1 if the upper boundary 'a' was reached,
#'     0 if the lower boundary 0 was reached, or NA if `max_decision_time`
#'     was exceeded.}
#'   \item{rt}{Numeric: The total reaction time (decision_time + ter) in seconds,
#'     or NA if `max_decision_time` was exceeded.}
#'   \item{decision_time}{Numeric: The time taken for the evidence accumulation
#'     process to reach a boundary, in seconds (or NA).}
#' @export
#'
#' @examples
#' set.seed(123)
#' # Simulate a trial with positive drift, starting in the middle
#' trial_result <- simulate_diffusion_trial(v = 0.2, a = 1.0, z = 0.5, ter = 0.15)
#' print(trial_result)
#' #> $choice
#' #> [1] 1
#' #> $rt
#' #> [1] 0.611
#' #> $decision_time
#' #> [1] 0.461

#' # Simulate a trial with negative drift
#' trial_result_neg_drift <- simulate_diffusion_trial(v = -0.1, a = 1.2, z = 0.6, ter = 0.2)
#' print(trial_result_neg_drift)
#' #> $choice
#' #> [1] 0
#' #> $rt
#' #> [1] 1.213
#' #> $decision_time
#' #> [1] 1.013
simulate_diffusion_trial <- function(v,
                                     a,
                                     z,
                                     s = 0.1,
                                     dt = 0.001,
                                     ter = 0.1,
                                     max_decision_time = 5.0) {
  
  # Parameter validation
  validate_ddm_parameters(v, a, z, s, dt, ter, max_decision_time, "simulate_diffusion_trial")
  
  evidence <- z
  current_decision_time <- 0
  max_steps <- max_decision_time / dt # Convert max_decision_time to number of steps
  
  # Simulate the evidence accumulation process
  time_steps_taken <- 0
  while (evidence > 0 && evidence < a && time_steps_taken < max_steps) {
    # Generate random increment based on Wiener process approximation
    # Increment = deterministic_part + random_part
    # deterministic_part = v * dt
    # random_part = s * sqrt(dt) * N(0,1) where N(0,1) is a standard normal deviate
    increment <- rnorm(n = 1, mean = v * dt, sd = s * sqrt(dt))
    
    evidence <- evidence + increment
    time_steps_taken <- time_steps_taken + 1
  }
  
  decision_time <- time_steps_taken * dt
  
  # Determine outcome
  if (time_steps_taken >= max_steps) {
    # Max decision time reached without hitting a boundary
    choice <- NA
    rt <- NA
    decision_time_output <- NA # Explicitly set for return
  } else {
    # A boundary was hit
    if (evidence >= a) {
      choice <- 1 # Upper boundary 'a'
    } else { # evidence <= 0
      choice <- 0 # Lower boundary 0
    }
    rt <- decision_time + ter
    decision_time_output <- decision_time
  }
  
  return(list(choice = choice, rt = rt, decision_time = decision_time_output))
}


#' Simulate multiple trials of the Diffusion Decision Model (DDM).
#'
#' This function calls `simulate_diffusion_trial()` multiple times and
#' compiles the results into a data frame.
#'
#' @param n_trials Integer. The number of trials to simulate. Default is 100.
#' @param v Numeric. Drift rate.
#' @param a Numeric. Threshold separation.
#' @param z Numeric. Starting point.
#' @param s Numeric. Noise standard deviation. Default is 0.1.
#' @param dt Numeric. Time step. Default is 0.001.
#' @param ter Numeric. Non-decision time. Default is 0.1.
#' @param max_decision_time Numeric. Maximum decision time. Default is 5.0.
#' @param verbose Logical. If TRUE, prints progress information. Default is FALSE.
#'
#' @return A data frame with `n_trials` rows and columns:
#'   \item{trial}{Integer: The trial number.}
#'   \item{choice}{Integer: The choice made on that trial (1 for upper, 0 for lower, or NA).}
#'   \item{rt}{Numeric: The total reaction time (in seconds) for that trial (or NA).}
#'   \item{decision_time}{Numeric: The decision time component (in seconds) for that trial (or NA).}
#' @export
#'
#' @examples
#' set.seed(456)
#' ddm_experiment_data <- simulate_diffusion_experiment(
#'   n_trials = 5,
#'   v = 0.15,
#'   a = 0.8,
#'   z = 0.4,
#'   ter = 0.2
#' )
#' print(ddm_experiment_data)
#' #>   trial choice        rt decision_time
#' #> 1     1      1 0.5090001     0.3090001
#' #> 2     2      1 0.4250001     0.2250001
#' #> 3     3      0 0.3830001     0.1830001
#' #> 4     4      1 0.7940001     0.5940001
#' #> 5     5      0 0.3940001     0.1940001
simulate_diffusion_experiment <- function(n_trials = 100, 
                                          v,
                                          a, 
                                          z, 
                                          s = 0.1,
                                          dt = 0.001,
                                          ter = 0.1,
                                          max_decision_time = 5.0,
                                          verbose = FALSE) {
  
  # Validate parameters
  validate_ddm_parameters(v, a, z, s, dt, ter, max_decision_time, "simulate_diffusion_experiment")
  
  if (!is.numeric(n_trials) || length(n_trials) != 1 || n_trials <= 0 || n_trials != round(n_trials)) {
    stop("'n_trials' must be a single positive integer.")
  }
  
  if (verbose) {
    cat("Simulating", n_trials, "DDM trials with parameters:\n")
    cat("  v =", v, ", a =", a, ", z =", z, ", s =", s, ", ter =", ter, ", dt =", dt, "\n")
    cat("  Progress: ")
  }
  
  # Create an empty list to store results
  results_list <- vector("list", n_trials)
  
  # Run the simulation n_trials times and store results
  for (i in 1:n_trials) {
    results_list[[i]] <- simulate_diffusion_trial(
      v = v,
      a = a,
      z = z,
      s = s,
      dt = dt,
      ter = ter,
      max_decision_time = max_decision_time
    )
    
    # Progress indicator
    if (verbose && i %% max(1, round(n_trials/10)) == 0) {
      cat(round(100 * i / n_trials), "% ")
    }
  }
  
  if (verbose) cat("\nCompleted!\n")
  
  # Convert the list of lists into a data frame
  df_results <- data.frame(
    trial = 1:n_trials,
    choice = sapply(results_list, function(trial_output) trial_output$choice),
    rt = sapply(results_list, function(trial_output) trial_output$rt),
    decision_time = sapply(results_list, function(trial_output) trial_output$decision_time)
  )
  
  return(df_results)
}

#' Simulate a single DDM trial and store the evidence accumulation path.
#'
#' This version of the DDM trial simulator records the evidence level at each
#' time step, allowing for visualization of the accumulation path.
#'
#' @param v Numeric. Drift rate.
#' @param a Numeric. Threshold separation.
#' @param z Numeric. Starting point of evidence accumulation.
#' @param s Numeric. Standard deviation of within-trial noise. Default is 0.1.
#' @param dt Numeric. Time step for the simulation. Default is 0.001.
#' @param ter Numeric. Non-decision time. Default is 0.1.
#' @param max_decision_time Numeric. Maximum allowed decision time. Default is 5.0.
#'
#' @return A list containing:
#'   \item{choice}{Integer: 1 for upper, 0 for lower, NA for timeout.}
#'   \item{rt}{Numeric: Total reaction time (or NA).}
#'   \item{decision_time}{Numeric: Decision time component (or NA).}
#'   \item{path_data}{Data frame: Columns for 'time_s' (time in seconds) and
#'     'evidence'. Stores the trajectory.}
#'   \item{a}{Numeric: The threshold value (for convenience in plotting).}
#' @export
#'
#' @examples
#' set.seed(1)
#' path_trial <- simulate_diffusion_trial_with_path(v = 0.2, a = 1, z = 0.5)
#' # head(path_trial$path_data)
#' # plot(path_trial$path_data$time_s, path_trial$path_data$evidence, type='l',
#' #      xlab="Time (s)", ylab="Evidence", main="Example DDM Path")
#' # abline(h=c(0, path_trial$a), col="red", lty=2)
simulate_diffusion_trial_with_path <- function(v,
                                               a,
                                               z,
                                               s = 0.1,
                                               dt = 0.001,
                                               ter = 0.1,
                                               max_decision_time = 5.0) {
  # Parameter validation
  validate_ddm_parameters(v, a, z, s, dt, ter, max_decision_time, "simulate_diffusion_trial_with_path")
  
  evidence <- z
  current_decision_time <- 0
  max_steps <- max_decision_time / dt
  
  # Pre-allocate for path data for efficiency (estimate max size)
  # Add a little buffer to max_steps for the case where it exactly hits max_steps
  estimated_path_length <- as.integer(max_steps + 10)
  path_time_s <- numeric(estimated_path_length)
  path_evidence <- numeric(estimated_path_length)
  
  time_steps_taken <- 0
  path_idx <- 1 # Index for storing path data
  
  # Store initial state
  path_time_s[path_idx] <- current_decision_time
  path_evidence[path_idx] <- evidence
  path_idx <- path_idx + 1
  
  while (evidence > 0 && evidence < a && time_steps_taken < max_steps) {
    increment <- rnorm(n = 1, mean = v * dt, sd = s * sqrt(dt))
    evidence <- evidence + increment
    time_steps_taken <- time_steps_taken + 1
    current_decision_time <- time_steps_taken * dt
    
    # Store path data
    if (path_idx <= estimated_path_length) {
      path_time_s[path_idx] <- current_decision_time
      path_evidence[path_idx] <- evidence
    } else {
      # This case should ideally not be hit if estimated_path_length is sufficient
      # Fallback: dynamically grow (less efficient)
      path_time_s <- c(path_time_s, current_decision_time)
      path_evidence <- c(path_evidence, evidence)
    }
    path_idx <- path_idx + 1
  }
  
  decision_time <- current_decision_time # or time_steps_taken * dt
  
  # Trim unused parts of pre-allocated vectors
  actual_path_length <- path_idx - 1
  path_df <- data.frame(
    time_s = path_time_s[1:actual_path_length],
    evidence = path_evidence[1:actual_path_length]
  )
  
  if (time_steps_taken >= max_steps) {
    choice <- NA
    rt <- NA
    decision_time_output <- NA
  } else {
    if (evidence >= a) {
      choice <- 1
    } else { # evidence <= 0
      choice <- 0
    }
    rt <- decision_time + ter
    decision_time_output <- decision_time
  }
  
  return(list(choice = choice,
              rt = rt,
              decision_time = decision_time_output,
              path_data = path_df,
              a = a)) # Return 'a' to make plotting thresholds easier
}

#' Calculate comprehensive summary statistics for DDM experiment data
#'
#' Computes accuracy, mean/median RTs, choice proportions, and other useful
#' statistics for DDM simulation results.
#'
#' @param ddm_data Data frame. Output from `simulate_diffusion_experiment()`.
#' @param correct_response Integer. Which choice (0 or 1) is considered "correct" 
#'   for accuracy calculations. If NULL, accuracy is not calculated. Default is 1.
#'
#' @return A list containing summary statistics:
#'   \item{n_trials}{Total number of trials}
#'   \item{n_valid}{Number of valid trials (non-timeout)}
#'   \item{n_timeout}{Number of timeout trials}
#'   \item{timeout_rate}{Proportion of timeout trials}
#'   \item{choice_proportions}{Table of choice proportions}
#'   \item{accuracy}{Overall accuracy (if correct_response specified)}
#'   \item{rt_overall}{RT statistics for all valid trials}
#'   \item{rt_by_choice}{RT statistics split by choice}
#'   \item{decision_time_overall}{Decision time statistics for all valid trials}
#'   \item{decision_time_by_choice}{Decision time statistics split by choice}
#' @export
#'
#' @examples
#' set.seed(123)
#' data <- simulate_diffusion_experiment(n_trials = 1000, v = 0.2, a = 1.0, z = 0.5)
#' summary_stats <- summarize_ddm_data(data, correct_response = 1)
#' print(summary_stats$rt_overall)
summarize_ddm_data <- function(ddm_data, correct_response = 1) {
  
  if (!is.data.frame(ddm_data)) {
    stop("'ddm_data' must be a data frame from simulate_diffusion_experiment()")
  }
  
  required_cols <- c("trial", "choice", "rt", "decision_time")
  if (!all(required_cols %in% names(ddm_data))) {
    stop("'ddm_data' must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  n_trials <- nrow(ddm_data)
  n_valid <- sum(!is.na(ddm_data$choice))
  n_timeout <- sum(is.na(ddm_data$choice))
  timeout_rate <- n_timeout / n_trials
  
  # Choice proportions
  valid_data <- ddm_data[!is.na(ddm_data$choice), ]
  choice_proportions <- if (nrow(valid_data) > 0) {
    prop.table(table(valid_data$choice))
  } else {
    table(numeric(0))
  }
  
  # Accuracy
  accuracy <- if (!is.null(correct_response) && nrow(valid_data) > 0) {
    mean(valid_data$choice == correct_response)
  } else {
    NA
  }
  
  # RT statistics helper function
  calc_rt_stats <- function(rt_values) {
    if (length(rt_values) == 0 || all(is.na(rt_values))) {
      return(list(n = 0, mean = NA, median = NA, sd = NA, min = NA, max = NA,
                  q25 = NA, q75 = NA, skewness = NA))
    }
    rt_clean <- rt_values[!is.na(rt_values)]
    if (length(rt_clean) == 0) {
      return(list(n = 0, mean = NA, median = NA, sd = NA, min = NA, max = NA,
                  q25 = NA, q75 = NA, skewness = NA))
    }
    
    list(
      n = length(rt_clean),
      mean = mean(rt_clean),
      median = median(rt_clean),
      sd = sd(rt_clean),
      min = min(rt_clean),
      max = max(rt_clean),
      q25 = quantile(rt_clean, 0.25),
      q75 = quantile(rt_clean, 0.75),
      skewness = (mean(rt_clean) - median(rt_clean)) / sd(rt_clean)
    )
  }
  
  # Overall RT statistics
  rt_overall <- calc_rt_stats(valid_data$rt)
  decision_time_overall <- calc_rt_stats(valid_data$decision_time)
  
  # RT statistics by choice
  rt_by_choice <- list()
  decision_time_by_choice <- list()
  
  for (choice_val in unique(valid_data$choice)) {
    if (!is.na(choice_val)) {
      choice_data <- valid_data[valid_data$choice == choice_val, ]
      rt_by_choice[[paste0("choice_", choice_val)]] <- calc_rt_stats(choice_data$rt)
      decision_time_by_choice[[paste0("choice_", choice_val)]] <- calc_rt_stats(choice_data$decision_time)
    }
  }
  
  return(list(
    n_trials = n_trials,
    n_valid = n_valid,
    n_timeout = n_timeout,
    timeout_rate = timeout_rate,
    choice_proportions = choice_proportions,
    accuracy = accuracy,
    rt_overall = rt_overall,
    rt_by_choice = rt_by_choice,
    decision_time_overall = decision_time_overall,
    decision_time_by_choice = decision_time_by_choice
  ))
}

#' Create a parameter set for DDM simulations
#'
#' Convenience function to create a named list of DDM parameters with validation.
#' Useful for systematic parameter exploration and documentation.
#'
#' @param v Numeric. Drift rate.
#' @param a Numeric. Threshold separation.
#' @param z Numeric. Starting point.
#' @param s Numeric. Noise standard deviation. Default is 0.1.
#' @param dt Numeric. Time step. Default is 0.001.
#' @param ter Numeric. Non-decision time. Default is 0.1.
#' @param max_decision_time Numeric. Maximum decision time. Default is 5.0.
#' @param name Character. Optional name for the parameter set.
#'
#' @return A list containing the validated parameters with class 'ddm_params'.
#' @export
#'
#' @examples
#' # Create parameter sets for different conditions
#' easy_params <- create_ddm_params(v = 0.3, a = 1.0, z = 0.5, name = "Easy Condition")
#' hard_params <- create_ddm_params(v = 0.1, a = 1.0, z = 0.5, name = "Hard Condition")
#' 
#' print(easy_params)
create_ddm_params <- function(v, a, z, s = 0.1, dt = 0.001, ter = 0.1, 
                              max_decision_time = 5.0, name = NULL) {
  
  # Validate parameters
  validate_ddm_parameters(v, a, z, s, dt, ter, max_decision_time, "create_ddm_params")
  
  params <- list(
    v = v,
    a = a, 
    z = z,
    s = s,
    dt = dt,
    ter = ter,
    max_decision_time = max_decision_time,
    name = name
  )
  
  class(params) <- "ddm_params"
  return(params)
}

#' Print method for DDM parameters
#' @param x A ddm_params object
#' @param ... Additional arguments (unused)
#' @export
print.ddm_params <- function(x, ...) {
  cat("DDM Parameters")
  if (!is.null(x$name)) {
    cat(": ", x$name)
  }
  cat("\n")
  cat("  Drift rate (v):", x$v, "\n")
  cat("  Threshold (a):", x$a, "\n")
  cat("  Starting point (z):", x$z, sprintf(" (z/a = %.3f)", x$z/x$a), "\n")
  cat("  Noise (s):", x$s, "\n")
  cat("  Non-decision time (ter):", x$ter, "\n")
  cat("  Time step (dt):", x$dt, "\n")
  cat("  Max decision time:", x$max_decision_time, "\n")
  invisible(x)
}

#' Generate a systematic parameter grid for DDM exploration
#'
#' Creates a data frame with all combinations of specified parameter values,
#' useful for systematic exploration of DDM parameter space.
#'
#' @param v_values Numeric vector. Drift rate values to explore.
#' @param a_values Numeric vector. Threshold values to explore.
#' @param z_values Numeric vector. Starting point values to explore. Can also
#'   be specified as proportions of threshold (e.g., 0.5 for z = a/2).
#' @param s_values Numeric vector. Noise values. Default is 0.1.
#' @param ter_values Numeric vector. Non-decision time values. Default is 0.1.
#' @param z_as_proportion Logical. If TRUE, z_values are treated as proportions
#'   of threshold (z = z_value * a). Default is FALSE.
#'
#' @return A data frame where each row represents a parameter combination.
#' @export
#'
#' @examples
#' # Create a parameter grid for exploration
#' param_grid <- create_parameter_grid(
#'   v_values = c(-0.2, 0, 0.2),
#'   a_values = c(0.8, 1.0, 1.2),
#'   z_values = c(0.3, 0.5, 0.7),
#'   z_as_proportion = TRUE
#' )
#' print(head(param_grid))
create_parameter_grid <- function(v_values, a_values, z_values,
                                  s_values = 0.1, ter_values = 0.1,
                                  z_as_proportion = FALSE) {
  
  if (z_as_proportion) {
    # Create grid with proportional z values
    base_grid <- expand.grid(
      v = v_values,
      a = a_values,
      z_prop = z_values,
      s = s_values,
      ter = ter_values,
      stringsAsFactors = FALSE
    )
    base_grid$z <- base_grid$z_prop * base_grid$a
    base_grid$z_prop <- NULL
  } else {
    # Create grid with absolute z values
    base_grid <- expand.grid(
      v = v_values,
      a = a_values,
      z = z_values,
      s = s_values,
      ter = ter_values,
      stringsAsFactors = FALSE
    )
  }
  
  # Validate that all z values are between 0 and a
  valid_rows <- base_grid$z > 0 & base_grid$z < base_grid$a
  if (!all(valid_rows)) {
    n_invalid <- sum(!valid_rows)
    warning(paste("Removed", n_invalid, "parameter combinations where z was not between 0 and a."))
    base_grid <- base_grid[valid_rows, ]
  }
  
  # Add parameter set identifiers
  base_grid$param_set <- paste0("set_", seq_len(nrow(base_grid)))
  
  return(base_grid)
}

# --- Example Usage (you can run this interactively in RStudio) ---
# if (interactive()) {
#   # Test a single trial
#   set.seed(789) # for reproducibility
#   single_ddm_trial <- simulate_diffusion_trial(
#     v = 0.25,    # Drift rate
#     a = 1.2,     # Threshold separation
#     z = 0.6,     # Starting point (a/2 for unbiased)
#     s = 0.1,     # Noise (default)
#     dt = 0.001,  # Time step (default)
#     ter = 0.1,   # Non-decision time
#     max_decision_time = 5.0 # Max decision time (default)
#   )
#   cat("Single DDM Trial Result:\n")
#   print(single_ddm_trial)
#   
#   # Test a small experiment
#   set.seed(101112)
#   ddm_experiment_results <- simulate_diffusion_experiment(
#     n_trials = 1000,
#     v = 0.2,   # Moderate positive drift
#     a = 1.0,
#     z = 0.5,   # Start in the middle
#     s = 0.1,
#     ter = 0.1
#   )
#   cat("\nDDM Experiment Results (first few trials):\n")
#   print(head(ddm_experiment_results))
#   
#   # Example summary statistics
#   summary_stats <- summarize_ddm_data(ddm_experiment_results, correct_response = 1)
#   cat("\nSummary Statistics:\n")
#   print(summary_stats$rt_overall)
#   
#   # Example parameter creation
#   params <- create_ddm_params(v = 0.2, a = 1.0, z = 0.5, name = "Example")
#   print(params)
#   
#   # Example for plotting RT distributions (requires ggplot2 and dplyr)
#   if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
#     library(ggplot2)
#     library(dplyr)
#     
#     set.seed(131415)
#     large_ddm_data <- simulate_diffusion_experiment(
#       n_trials = 2000,
#       v = 0.15,      # Moderate drift
#       a = 1.2,       # Reasonable threshold
#       z = 0.6,       # Unbiased start (a/2)
#       s = 0.1,       # Standard noise
#       ter = 0.2,     # Reasonable non-decision time
#       dt = 0.001,    # Standard time step
#       max_decision_time = 4 # Default
#     )
#     
#     # Filter out NA RTs for plotting
#     plot_data <- large_ddm_data %>% filter(!is.na(rt))
#     
#     if (nrow(plot_data) > 0) {
#       p <- ggplot(plot_data, aes(x = rt, fill = factor(choice))) +
#         geom_histogram(binwidth = 0.05, alpha = 0.7, position = "identity") +
#         facet_wrap(~factor(choice, labels = c("Lower Boundary (0)", "Upper Boundary (1)"))) +
#         labs(title = "DDM RT Distributions",
#              x = "Reaction Time (s)",
#              y = "Frequency",
#              fill = "Choice") +
#         theme_minimal()
#       print(p)
#       
#       # Choice proportions
#       cat("\nChoice Proportions (1=Upper, 0=Lower):\n")
#       print(prop.table(table(plot_data$choice, useNA = "ifany")))
#       
#     } else {
#       cat("\nNo valid trials to plot (all may have hit max_decision_time).\n")
#     }
#   } else {
#     cat("\nInstall ggplot2 and dplyr packages to see example plots (run: install.packages(c(\"ggplot2\", \"dplyr\")) ).\n")
#   }
# }


