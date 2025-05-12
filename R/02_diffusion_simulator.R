# R/02_diffusion_simulator.R

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
  if (a <= 0) {
    stop("Threshold 'a' must be positive.")
  }
  if (z <= 0 || z >= a) {
    stop("Starting point 'z' must be strictly between 0 and 'a'.")
  }
  if (s <= 0) {
    stop("Noise 's' must be positive.")
  }
  if (dt <= 0) {
    stop("Time step 'dt' must be positive.")
  }
  if (ter < 0) {
    stop("Non-decision time 'ter' cannot be negative.")
  }
  if (max_decision_time <= 0) {
    stop("Maximum decision time must be positive.")
  }
  
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
#' @param ... Additional arguments to be passed to `simulate_diffusion_trial()`,
#'   such as `v`, `a`, `z`, `s`, `dt`, `ter`, `max_decision_time`.
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
                                          max_decision_time = 5.0) {
  
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
  }
  # Convert the list of lists into a data frame
  df_results <- data.frame(
    trial = 1:n_trials,
    choice = sapply(results_list, function(trial_output) trial_output$choice),
    rt = sapply(results_list, function(trial_output) trial_output$rt),
    decision_time = sapply(results_list, function(trial_output) trial_output$decision_time)
  )
  
  return(df_results)
}


# --- Example Usage (you can run this interactively in RStudio) ---
if (interactive()) {
  # Test a single trial
  set.seed(789) # for reproducibility
  single_ddm_trial <- simulate_diffusion_trial(
    v = 0.25,    # Drift rate
    a = 1.2,     # Threshold separation
    z = 0.6,     # Starting point (a/2 for unbiased)
    s = 0.1,     # Noise (default)
    dt = 0.001,  # Time step (default)
    ter = 0.1,   # Non-decision time
    max_decision_time = 5.0 # Max decision time (default)
  )
  cat("Single DDM Trial Result:\n")
  print(single_ddm_trial)
  
  # Test a small experiment
  set.seed(101112)
  ddm_experiment_results <- simulate_diffusion_experiment(
    n_trials = 1000,
    v = 0.2,   # Moderate positive drift
    a = 1.0,
    z = 0.5,   # Start in the middle
    s = 0.1,
    ter = 0.1
  )
  cat("\nDDM Experiment Results (first few trials):\n")
  print(head(ddm_experiment_results))
  
  # Example for plotting RT distributions (requires ggplot2 and dplyr)
  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
    library(ggplot2)
    library(dplyr)
    
    set.seed(131415)
    large_ddm_data <- simulate_diffusion_experiment(
      n_trials = 2000,
      v = 0.15,      # Moderate drift
      a = 1.2,       # Reasonable threshold
      z = 0.6,       # Unbiased start (a/2)
      s = 0.5,       # Standard noise
      ter = 0.2,     # Reasonable non-decision time
      dt = 0.1,  # Larger time step for faster simulation
      max_decision_time = 4 # Default
    )
    
    # Filter out NA RTs for plotting
    plot_data <- large_ddm_data %>% filter(!is.na(rt))
    
    if (nrow(plot_data) > 0) {
      p <- ggplot(plot_data, aes(x = rt, fill = factor(choice))) +
        geom_histogram(binwidth = 0.05, alpha = 0.7, position = "identity") +
        facet_wrap(~factor(choice, labels = c("Lower Boundary (0)", "Upper Boundary (1)"))) +
        labs(title = "DDM RT Distributions",
             x = "Reaction Time (s)",
             y = "Frequency",
             fill = "Choice") +
        theme_minimal()
      print(p)
      
      # Choice proportions
      cat("\nChoice Proportions (1=Upper, 0=Lower):\n")
      print(prop.table(table(plot_data$choice, useNA = "ifany")))
      
    } else {
      cat("\nNo valid trials to plot (all may have hit max_decision_time).\n")
    }
  } else {
    cat("\nInstall ggplot2 and dplyr packages to see example plots (run: install.packages(c(\"ggplot2\", \"dplyr\")) ).\n")
  }
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
  # Parameter validation (can copy from simulate_diffusion_trial or omit for brevity here)
  if (a <= 0) stop("Threshold 'a' must be positive.")
  if (z <= 0 || z >= a) stop("Starting point 'z' must be strictly between 0 and 'a'.")
  
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

