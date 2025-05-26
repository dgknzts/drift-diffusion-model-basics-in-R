# R/01_random_walk_simulator.R

#' Simulate a single trial of a simple random walk model.
#'
#' In this model, evidence accumulates in discrete steps until it crosses
#' one of two (symmetric) decision thresholds. This is a discrete-time
#' approximation of the continuous-time diffusion decision model.
#'
#' The mathematical model can be expressed as:
#' E(t+1) = E(t) + drift + N(0, sd_step^2)
#' where E(t) is evidence at time t, and the process continues until
#' |E(t)| >= threshold.
#'
#' @param start_point Numeric. The starting value of the evidence accumulator.
#'   Default is 0.
#' @param threshold Numeric. The positive value of the decision threshold.
#'   The decision is made when `abs(evidence) >= threshold`. Default is 10.
#' @param drift Numeric. The average amount added to the evidence at each step.
#'   Default is 0.1.
#' @param sd_step Numeric. The standard deviation of the amount added to the
#'   evidence at each step (reflecting noise). Must be positive. Default is 1.
#' @param max_steps Integer. The maximum number of steps allowed for a trial
#'   before it's considered a non-decision (to prevent infinite loops).
#'   Default is 10000.
#' @param return_path Logical. If TRUE, returns the full evidence path.
#'   Default is FALSE.
#'
#' @return A list containing:
#'   \item{choice}{Integer: 1 if the upper threshold was crossed, -1 if the
#'     lower threshold was crossed, or NA if `max_steps` was reached before
#'     a threshold.}
#'   \item{rt}{Integer: The number of steps taken to reach a threshold (reaction time),
#'     or NA if `max_steps` was reached.}
#'   \item{evidence_path}{Numeric vector: The full evidence path (only if return_path = TRUE).}
#' @export
#'
#' @examples
#' set.seed(123)
#' trial_result <- simulate_random_walk_trial(threshold = 5, drift = 0.2, sd_step = 1.5)
#' print(trial_result)
#' #> $choice
#' #> [1] 1
#' #> $rt
#' #> [1] 16
#' 
#' # Get full evidence path
#' trial_with_path <- simulate_random_walk_trial(threshold = 5, drift = 0.2, 
#'                                              return_path = TRUE)
#' plot(trial_with_path$evidence_path, type = "l", 
#'      main = "Evidence Accumulation Path")
#' abline(h = c(-5, 5), lty = 2, col = "red")

simulate_random_walk_trial <- function(start_point = 0,
                                       threshold = 10,
                                       drift = 0.1,
                                       sd_step = 1,
                                       max_steps = 10000,
                                       return_path = FALSE) {
  
  # Input validation
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0) {
    stop("Threshold must be a single positive numeric value.")
  }
  
  if (!is.numeric(sd_step) || length(sd_step) != 1 || sd_step <= 0) {
    stop("sd_step must be a single positive numeric value.")
  }
  
  if (!is.numeric(start_point) || length(start_point) != 1) {
    stop("start_point must be a single numeric value.")
  }
  
  if (!is.numeric(drift) || length(drift) != 1) {
    stop("drift must be a single numeric value.")
  }
  
  if (!is.numeric(max_steps) || length(max_steps) != 1 || max_steps <= 0 || max_steps != round(max_steps)) {
    stop("max_steps must be a single positive integer.")
  }
  
  if (!is.logical(return_path) || length(return_path) != 1) {
    stop("return_path must be a single logical value.")
  }
  
  # Check if start_point already exceeds threshold
  if (abs(start_point) >= threshold) {
    choice <- ifelse(start_point >= threshold, 1, -1)
    rt <- 0
    evidence_path <- if (return_path) start_point else NULL
    
    result <- list(choice = choice, rt = rt)
    if (return_path) result$evidence_path <- evidence_path
    return(result)
  }
  
  evidence <- start_point
  time_steps <- 0
  
  # Initialize evidence path if requested
  if (return_path) {
    evidence_path <- numeric(max_steps + 1)
    evidence_path[1] <- start_point
  }
  
  while (abs(evidence) < threshold && time_steps < max_steps) {
    # Generate a random step from a normal distribution
    step_increment <- rnorm(n = 1, mean = drift, sd = sd_step)
    
    # Update evidence
    evidence <- evidence + step_increment
    
    # Increment time
    time_steps <- time_steps + 1
    
    # Store evidence path if requested
    if (return_path) {
      evidence_path[time_steps + 1] <- evidence
    }
  }
  
  # Determine the outcome
  if (time_steps == max_steps) {
    # Max steps reached without a decision
    choice <- NA
    rt <- NA
  } else {
    # A threshold was crossed
    rt <- time_steps
    if (evidence >= threshold) {
      choice <- 1 # Upper threshold
    } else {
      choice <- -1 # Lower threshold (since abs(evidence) >= threshold was met)
    }
  }
  
  # Prepare return value
  result <- list(choice = choice, rt = rt)
  
  if (return_path) {
    if (is.na(rt)) {
      result$evidence_path <- evidence_path[1:(max_steps + 1)]
    } else {
      result$evidence_path <- evidence_path[1:(rt + 1)]
    }
  }
  
  return(result)
}


#' Simulate multiple trials of a simple random walk model.
#'
#' This function calls `simulate_random_walk_trial()` multiple times and
#' compiles the results into a data frame. For large simulations, progress
#' is displayed.
#'
#' @param n_trials Integer. The number of trials to simulate. Default is 100.
#' @param show_progress Logical. Whether to show progress for large simulations
#'   (n_trials > 1000). Default is TRUE.
#' @param ... Additional arguments to be passed to `simulate_random_walk_trial()`,
#'   such as `start_point`, `threshold`, `drift`, `sd_step`, `max_steps`.
#'   Note: `return_path` is not supported in this function for memory efficiency.
#'
#' @return A data frame with `n_trials` rows and columns:
#'   \item{trial}{Integer: The trial number.}
#'   \item{choice}{Integer: The choice made on that trial (1, -1, or NA).}
#'   \item{rt}{Integer: The reaction time (number of steps) for that trial (or NA).}
#' @export
#'
#' @examples
#' set.seed(456)
#' experiment_data <- simulate_random_walk_experiment(n_trials = 5,
#'                                                   threshold = 8,
#'                                                   drift = 0.05,
#'                                                   sd_step = 1.2)
#' print(experiment_data)
#' #>   trial choice rt
#' #> 1     1      1 45
#' #> 2     2      1 52
#' #> 3     3     -1 24
#' #> 4     4      1 70
#' #> 5     5     -1 51

simulate_random_walk_experiment <- function(n_trials = 100, show_progress = TRUE, ...) {
  
  # Input validation
  if (!is.numeric(n_trials) || length(n_trials) != 1 || n_trials <= 0 || n_trials != round(n_trials)) {
    stop("n_trials must be a single positive integer.")
  }
  
  if (!is.logical(show_progress) || length(show_progress) != 1) {
    stop("show_progress must be a single logical value.")
  }
  
  # Capture additional arguments and ensure return_path is FALSE for efficiency
  trial_args <- list(...)
  trial_args$return_path <- FALSE
  
  # Show progress for large simulations
  if (show_progress && n_trials > 1000) {
    message(sprintf("Simulating %d trials...", n_trials))
    
    # Use replicate with progress indication
    progress_interval <- max(1, floor(n_trials / 20))
    results_list <- vector("list", n_trials)
    
    for (i in 1:n_trials) {
      if (i %% progress_interval == 0 || i == n_trials) {
        message(sprintf("Progress: %d/%d trials completed (%.1f%%)", 
                       i, n_trials, 100 * i / n_trials))
      }
      results_list[[i]] <- do.call(simulate_random_walk_trial, trial_args)
    }
  } else {
    # Use replicate for smaller simulations
    results_list <- replicate(n_trials, 
                             do.call(simulate_random_walk_trial, trial_args), 
                             simplify = FALSE)
  }
  
  # Convert the list of lists into a data frame more efficiently
  df_results <- data.frame(
    trial = 1:n_trials,
    choice = vapply(results_list, function(x) x$choice, numeric(1)),
    rt = vapply(results_list, function(x) x$rt, numeric(1))
  )
  
  return(df_results)
}


#' Plot the evidence path of a single random walk trial
#'
#' @param evidence_path Numeric vector. The evidence path from a trial.
#' @param threshold Numeric. The decision threshold for plotting boundaries.
#' @param main Character. Title for the plot.
#' @param ... Additional arguments passed to plot().
#'
#' @return No return value, creates a plot.
#' @export
#'
#' @examples
#' set.seed(123)
#' trial <- simulate_random_walk_trial(threshold = 5, return_path = TRUE)
#' plot_evidence_path(trial$evidence_path, threshold = 5)

plot_evidence_path <- function(evidence_path, threshold, main = "Evidence Accumulation Path", ...) {
  
  if (!is.numeric(evidence_path) || !is.numeric(threshold)) {
    stop("evidence_path and threshold must be numeric.")
  }
  
  plot(0:(length(evidence_path) - 1), evidence_path, 
       type = "l", lwd = 2,
       xlab = "Time Steps", ylab = "Evidence",
       main = main, ...)
  
  # Add threshold lines
  abline(h = c(-threshold, threshold), lty = 2, col = "red", lwd = 2)
  abline(h = 0, lty = 3, col = "gray50")
  
  # Add threshold labels
  text(length(evidence_path) * 0.1, threshold + threshold * 0.1, 
       "Upper Threshold", col = "red", cex = 0.8)
  text(length(evidence_path) * 0.1, -threshold - threshold * 0.1, 
       "Lower Threshold", col = "red", cex = 0.8)
}

# --- Example Usage (you can run this interactively in RStudio) ---
if (interactive()) {
  # Test a single trial with path
  set.seed(101) # for reproducibility
  single_trial <- simulate_random_walk_trial(drift = 0.2, threshold = 5, return_path = TRUE)
  cat("Single Trial Result:\n")
  print(list(choice = single_trial$choice, rt = single_trial$rt))
  
  # Plot the evidence path
  plot_evidence_path(single_trial$evidence_path, threshold = 5)
  
  # Test a small experiment
  set.seed(202)
  experiment_results <- simulate_random_walk_experiment(n_trials = 10,
                                                        drift = -0.05,
                                                        threshold = 7,
                                                        sd_step = 1.5,
                                                        show_progress = FALSE)
  cat("\nExperiment Results (first few trials):\n")
  print(head(experiment_results))
  
  # RT distribution for a larger experiment
  set.seed(303)
  larger_experiment <- simulate_random_walk_experiment(n_trials = 1000,
                                                       drift = 0.1,
                                                       threshold = 10,
                                                       sd_step = 1,
                                                       show_progress = FALSE)
  hist(larger_experiment$rt[!is.na(larger_experiment$rt)],
       breaks = 50,
       main = "RT Distribution for Random Walk",
       xlab = "RT (number of steps)")

  # Calculate choice proportions
  cat("\nChoice Proportions (1=Upper, -1=Lower):\n")
  print(prop.table(table(larger_experiment$choice, useNA = "ifany")))
}