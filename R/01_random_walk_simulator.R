# R/01_random_walk_simulator.R

#' Simulate a single trial of a simple random walk model.
#'
#' In this model, evidence accumulates in discrete steps until it crosses
#' one of two (symmetric) decision thresholds.
#'
#' @param start_point Numeric. The starting value of the evidence accumulator.
#'   Default is 0.
#' @param threshold Numeric. The positive value of the decision threshold.
#'   The decision is made when `abs(evidence) >= threshold`. Default is 10.
#' @param drift Numeric. The average amount added to the evidence at each step.
#'   Default is 0.1.
#' @param sd_step Numeric. The standard deviation of the amount added to the
#'   evidence at each step (reflecting noise). Default is 1.
#' @param max_steps Integer. The maximum number of steps allowed for a trial
#'   before it's considered a non-decision (to prevent infinite loops).
#'   Default is 10000.
#'
#' @return A list containing:
#'   \item{choice}{Integer: 1 if the upper threshold was crossed, -1 if the
#'     lower threshold was crossed, or NA if `max_steps` was reached before
#'     a threshold.}
#'   \item{rt}{Integer: The number of steps taken to reach a threshold (reaction time),
#'     or NA if `max_steps` was reached.}
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



simulate_random_walk_trial <- function(start_point = 0,
                                       threshold = 10,
                                       drift = 0.1,
                                       sd_step = 1,
                                       max_steps = 10000) {
  
  if (threshold <= 0) {
    stop("Threshold must be a positive value.")
  }
  
  evidence <- start_point
  time_steps <- 0
  
  while (abs(evidence) < threshold && time_steps < max_steps) {
    # Generate a random step from a normal distribution
    step_increment <- rnorm(n = 1, mean = drift, sd = sd_step)
    
    # Update evidence
    evidence <- evidence + step_increment
    
    # Increment time
    time_steps <- time_steps + 1
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
  
  return(list(choice = choice, rt = rt))
}


#' Simulate multiple trials of a simple random walk model.
#'
#' This function calls `simulate_random_walk_trial()` multiple times and
#' compiles the results into a data frame.
#'
#' @param n_trials Integer. The number of trials to simulate. Default is 100.
#' @param ... Additional arguments to be passed to `simulate_random_walk_trial()`,
#'   such as `start_point`, `threshold`, `drift`, `sd_step`, `max_steps`.
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


simulate_random_walk_experiment <- function(n_trials = 100, ...) {
  # Use replicate to run the single trial simulation n_trials times
  # simplify = FALSE ensures that we get a list of lists
  results_list <- replicate(n_trials, simulate_random_walk_trial(...), simplify = FALSE)
  
  # Convert the list of lists into a data frame
  # We can extract 'choice' and 'rt' from each element of results_list
  df_results <- data.frame(
    trial = 1:n_trials,
    choice = sapply(results_list, function(trial_output) trial_output$choice),
    rt = sapply(results_list, function(trial_output) trial_output$rt)
  )
  
  return(df_results)
}

# --- Example Usage (you can run this interactively in RStudio) ---
if (interactive()) {
  # Test a single trial
  set.seed(101) # for reproducibility
  single_trial <- simulate_random_walk_trial(drift = 0.2, threshold = 5)
  cat("Single Trial Result:\n")
  print(single_trial)
  
  # Test a small experiment
  set.seed(202)
  experiment_results <- simulate_random_walk_experiment(n_trials = 10,
                                                        drift = -0.05,
                                                        threshold = 7,
                                                        sd_step = 1.5)
  cat("\nExperiment Results (first few trials):\n")
  print(head(experiment_results))
  
  # RT distribution for a larger experiment
  set.seed(303)
  larger_experiment <- simulate_random_walk_experiment(n_trials = 1000,
                                                       drift = 0.1,
                                                       threshold = 10,
                                                       sd_step = 1)
  hist(larger_experiment$rt[!is.na(larger_experiment$rt)],
       breaks = 50,
       main = "RT Distribution for Random Walk",
       xlab = "RT (number of steps)")

  # Calculate choice proportions
  cat("\nChoice Proportions (1=Upper, -1=Lower):\n")
  print(prop.table(table(larger_experiment$choice, useNA = "ifany")))
}