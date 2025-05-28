#' Simulate a Single DDM Trial with Across-Trial Parameter Variability
#'
#' This function simulates one trial of the Diffusion Decision Model (DDM),
#' incorporating trial-to-trial variability in key parameters: drift rate (`v`),
#' starting point (`z`), and non-decision time (`ter`). This allows for more
#' realistic modeling of behavioral data compared to a DDM with fixed parameters.
#' On each call, values for `v`, `z`, and `ter` for the current trial are sampled
#' from specified distributions before the evidence accumulation process begins.
#'
#' @details
#' The evidence accumulation process itself follows the standard DDM:
#' evidence starts at `z_trial` and accumulates with a mean rate of `v_trial`
#' and within-trial noise `s` until it reaches either the lower boundary (0)
#' or the upper boundary (`a`). The time step `dt` controls the granularity
#' of the simulation.
#'
#' Parameter Variability:
#' \itemize{
#'   \item **Drift Rate (`v_trial`):** Sampled from a Normal distribution:
#'     `v_trial ~ Normal(mean_v, sv)`. `sv` is the standard deviation of this
#'     across-trial distribution.
#'   \item **Starting Point (`z_trial`):** Sampled from a Uniform distribution:
#'     `z_trial ~ Uniform(mean_z - sz/2, mean_z + sz/2)`. `sz` is the total
#'     range of this uniform distribution. The sampled `z_trial` is subsequently
#'     clipped to ensure it remains strictly between 0 and `a` to prevent
#'     simulation errors.
#'   \item **Non-Decision Time (`ter_trial`):** Sampled from a Uniform distribution:
#'     `ter_trial ~ Uniform(mean_ter - st0/2, mean_ter + st0/2)`. `st0` is the
#'     total range. The sampled `ter_trial` is clipped to be non-negative.
#' }
#' If the variability parameter (e.g., `sv`) is set to 0, the corresponding
#' model parameter (e.g., `v_trial`) will simply take its mean value (`mean_v`).
#' The threshold `a` and within-trial noise `s` are currently assumed to be
#' constant across trials in this implementation.
#'
#' @param mean_v Numeric. The mean of the distribution from which the trial-specific
#'   drift rate (`v_trial`) is sampled.
#' @param a Numeric. The threshold separation (upper boundary). The lower boundary
#'   is fixed at 0. This parameter is assumed to be constant across trials here.
#'   Must be positive.
#' @param mean_z Numeric. The mean of the distribution from which the trial-specific
#'   starting point (`z_trial`) is sampled. Must be strictly between 0 and `a`.
#' @param s Numeric, optional. The standard deviation of the within-trial Gaussian
#'   noise in the evidence accumulation process. This is assumed to be constant
#'   across trials. Default is `0.1`.
#' @param dt Numeric, optional. The time step (in seconds) used for the discrete
#'   approximation of the continuous diffusion process. Smaller values lead to
#'   more accurate simulations but increase computation time. Default is `0.001`.
#' @param mean_ter Numeric, optional. The mean of the distribution from which the
#'   trial-specific non-decision time (`ter_trial`) is sampled. Default is `0.1`.
#' @param sv Numeric, optional. The standard deviation of the across-trial
#'   Normal distribution for drift rate (`v`). If `sv = 0` (default), drift rate
#'   does not vary across trials (i.e., `v_trial = mean_v`).
#' @param sz Numeric, optional. The total range of the across-trial Uniform
#'   distribution for the starting point (`z`). If `sz = 0` (default), the
#'   starting point does not vary (i.e., `z_trial = mean_z`).
#' @param st0 Numeric, optional. The total range of the across-trial Uniform
#'   distribution for non-decision time (`ter`). If `st0 = 0` (default),
#'   non-decision time does not vary (i.e., `ter_trial = mean_ter`).
#' @param max_decision_time Numeric, optional. The maximum duration (in seconds)
#'   allowed for the evidence accumulation process. If a boundary is not reached
#'   within this time, the trial is considered a timeout (choice=NA, rt=NA).
#'   Default is `5.0`.
#'
#' @return A list containing:
#'   \item{choice}{Integer: `1` if the upper boundary `a` was reached,
#'     `0` if the lower boundary 0 was reached, or `NA` if `max_decision_time`
#'     was exceeded.}
#'   \item{rt}{Numeric: The total reaction time (decision_time + `ter_trial`)
#'     in seconds, or `NA` if a timeout occurred.}
#'   \item{decision_time}{Numeric: The time taken (in seconds) for the evidence
#'     accumulation process to reach a boundary, or `NA` if a timeout occurred.}
#'   \item{v_trial}{Numeric: The specific drift rate sampled and used for this trial.}
#'   \item{z_trial}{Numeric: The specific starting point sampled and used for this trial.}
#'   \item{ter_trial}{Numeric: The specific non-decision time sampled and used for this trial.}
#' @export
#'
#' @seealso \code{\link{simulate_diffusion_trial}} for a DDM simulation without across-trial variability,
#'   \code{\link{simulate_diffusion_experiment_variable}} for simulating multiple trials with variability.
#'
#' @examples
#' set.seed(123)
#' # Single trial with variability in v
#' trial_var_v <- simulate_diffusion_trial_variable(
#'   mean_v = 0.2, a = 1.0, mean_z = 0.5, mean_ter = 0.15,
#'   sv = 0.1, sz = 0, st0 = 0
#' )
#' print(trial_var_v)
#'
#' # Single trial with variability in v, z, and ter
#' trial_all_var <- simulate_diffusion_trial_variable(
#'   mean_v = 0.1, a = 1.2, mean_z = 0.6, mean_ter = 0.2,
#'   sv = 0.05, sz = 0.1, st0 = 0.04
#' )
#' print(trial_all_var)
simulate_diffusion_trial_variable <- function(mean_v,
                                              a,
                                              mean_z,
                                              s = 0.1,
                                              dt = 0.001,
                                              mean_ter = 0.1,
                                              sv = 0,
                                              sz = 0,
                                              st0 = 0,
                                              max_decision_time = 5.0) {
  # --- 1. Sample trial-specific parameters ---
  v_trial <- rnorm(1, mean = mean_v, sd = sv)
  
  if (sz > 0) {
    z_trial <- runif(1, min = mean_z - sz / 2, max = mean_z + sz / 2)
  } else {
    z_trial <- mean_z
  }
  epsilon_z <- 1e-6 
  z_trial_clipped <- max(epsilon_z, min(z_trial, a - epsilon_z))
  if (z_trial_clipped <=0 || z_trial_clipped >=a || is.na(z_trial_clipped)) { 
    z_trial_clipped <- mean_z 
    if (z_trial_clipped <=0 || z_trial_clipped >=a || is.na(z_trial_clipped)) {
      z_trial_clipped = a/2 
    }
  }
  z_trial <- z_trial_clipped
  
  
  if (st0 > 0) {
    ter_trial <- runif(1, min = mean_ter - st0 / 2, max = mean_ter + st0 / 2)
  } else {
    ter_trial <- mean_ter
  }
  ter_trial <- max(0, ter_trial)
  
  
  # --- 2. Run the core diffusion process ---
  evidence <- z_trial 
  current_decision_time <- 0
  max_steps <- max_decision_time / dt
  time_steps_taken <- 0
  
  while (evidence > 0 && evidence < a && time_steps_taken < max_steps) {
    increment <- rnorm(n = 1, mean = v_trial * dt, sd = s * sqrt(dt)) 
    evidence <- evidence + increment
    time_steps_taken <- time_steps_taken + 1
  }
  
  decision_time <- time_steps_taken * dt
  
  choice_val <- NA
  rt_val <- NA
  decision_time_out <- NA
  
  if (time_steps_taken >= max_steps) {
    choice_val <- NA
    rt_val <- NA
    decision_time_out <- NA
  } else {
    if (evidence >= a) {
      choice_val <- 1
    } else { 
      choice_val <- 0
    }
    rt_val <- decision_time + ter_trial 
    decision_time_out <- decision_time
  }
  
  return(list(
    choice = choice_val,
    rt = rt_val,
    decision_time = decision_time_out,
    v_trial = v_trial, 
    z_trial = z_trial,
    ter_trial = ter_trial
  ))
}


#' Simulate Multiple DDM Trials with Across-Trial Parameter Variability.
#'
#' This function serves as a wrapper around `simulate_diffusion_trial_variable`
#' to conduct a full DDM experiment by simulating a specified number of trials.
#' It collects the outcomes (choice, RT, decision time) as well as the
#' trial-specific parameters (`v_trial`, `z_trial`, `ter_trial`) that were
#' sampled and used for each individual trial.
#'
#' @param n_trials Integer. The total number of DDM trials to simulate for the
#'   experiment. Default is `100`.
#' @param mean_v Numeric. The mean of the distribution for `v_trial`.
#' @param a Numeric. The threshold separation (constant across trials).
#' @param mean_z Numeric. The mean of the distribution for `z_trial`.
#' @param s Numeric, optional. Within-trial noise SD. Default `0.1`.
#' @param dt Numeric, optional. Simulation time step. Default `0.001`.
#' @param mean_ter Numeric, optional. The mean of the distribution for `ter_trial`. Default `0.1`.
#' @param sv Numeric, optional. SD for across-trial `v` variability. Default `0`.
#' @param sz Numeric, optional. Range for across-trial `z` variability. Default `0`.
#' @param st0 Numeric, optional. Range for across-trial `ter` variability. Default `0`.
#' @param max_decision_time Numeric, optional. Max decision time. Default `5.0`.
#'
#' @return A data frame where each row represents a simulated trial. Columns include:
#'   trial, choice, rt, decision_time, v_trial, z_trial, ter_trial.
#' @export
simulate_diffusion_experiment_variable <- function(n_trials = 100,
                                                   mean_v,
                                                   a,
                                                   mean_z,
                                                   s = 0.1,
                                                   dt = 0.001,
                                                   mean_ter = 0.1,
                                                   sv = 0,
                                                   sz = 0,
                                                   st0 = 0,
                                                   max_decision_time = 5.0) {
  
  # Pre-allocate vectors for trial-specific parameters
  v_trials <- rnorm(n_trials, mean = mean_v, sd = sv)
  
  if (sz > 0) {
    z_trials <- runif(n_trials, min = mean_z - sz/2, max = mean_z + sz/2)
    # Clip z values to valid range
    z_trials <- pmax(1e-6, pmin(z_trials, a - 1e-6))
  } else {
    z_trials <- rep(mean_z, n_trials)
  }
  
  if (st0 > 0) {
    ter_trials <- runif(n_trials, min = mean_ter - st0/2, max = mean_ter + st0/2)
    ter_trials <- pmax(0, ter_trials)
  } else {
    ter_trials <- rep(mean_ter, n_trials)
  }
  
  # Pre-allocate results vectors
  choices <- integer(n_trials)
  rts <- numeric(n_trials)
  decision_times <- numeric(n_trials)
  
  # Vectorized simulation
  max_steps <- max_decision_time / dt
  for (i in 1:n_trials) {
    evidence <- z_trials[i]
    time_steps_taken <- 0
    
    while (evidence > 0 && evidence < a && time_steps_taken < max_steps) {
      increment <- rnorm(1, mean = v_trials[i] * dt, sd = s * sqrt(dt))
      evidence <- evidence + increment
      time_steps_taken <- time_steps_taken + 1
    }
    
    decision_time <- time_steps_taken * dt
    
    if (time_steps_taken >= max_steps) {
      choices[i] <- NA
      rts[i] <- NA
      decision_times[i] <- NA
    } else {
      choices[i] <- if (evidence >= a) 1L else 0L
      rts[i] <- decision_time + ter_trials[i]
      decision_times[i] <- decision_time
    }
  }
  
  # Create data frame directly from vectors
  data.frame(
    trial = 1:n_trials,
    choice = choices,
    rt = rts,
    decision_time = decision_times,
    v_trial = v_trials,
    z_trial = z_trials,
    ter_trial = ter_trials
  )
}


#' Simulate a single DDM trial with across-trial parameter variability AND store path.
#'
#' Combines features of simulate_diffusion_trial_variable and
#' simulate_diffusion_trial_with_path.
#' (Add more detailed Roxygen comments here similar to simulate_diffusion_trial_variable)
#'
#' @param mean_v Numeric. Mean drift rate.
#' @param a Numeric. Threshold.
#' @param mean_z Numeric. Mean starting point.
#' @param s Numeric. Within-trial noise. Default 0.1.
#' @param dt Numeric. Time step. Default 0.001.
#' @param mean_ter Numeric. Mean non-decision time. Default 0.1.
#' @param sv Numeric. SD of v variability. Default 0.
#' @param sz Numeric. Range of z variability. Default 0.
#' @param st0 Numeric. Range of ter variability. Default 0.
#' @param max_decision_time Numeric. Max decision time. Default 5.0.
#'
#' @return A list containing choice, rt, decision_time, v_trial, z_trial, ter_trial,
#'   and path_data (data frame with time_s, evidence).
#' @export
simulate_diffusion_trial_variable_with_path <- function(mean_v,
                                                        a,
                                                        mean_z,
                                                        s = 0.1,
                                                        dt = 0.001,
                                                        mean_ter = 0.1,
                                                        sv = 0,
                                                        sz = 0,
                                                        st0 = 0,
                                                        max_decision_time = 5.0) {
  
  # --- 1. Sample trial-specific parameters ---
  v_trial <- rnorm(1, mean = mean_v, sd = sv)
  
  if (sz > 0) {
    z_samp <- runif(1, min = mean_z - sz / 2, max = mean_z + sz / 2)
  } else {
    z_samp <- mean_z
  }
  epsilon_z <- 1e-6
  z_trial <- max(epsilon_z, min(z_samp, a - epsilon_z))
  if (z_trial <= 0 || z_trial >= a || is.na(z_trial)) {
    z_trial <- mean_z
    if (z_trial <= 0 || z_trial >= a || is.na(z_trial)) {
      z_trial <- a / 2
    }
  }
  
  if (st0 > 0) {
    ter_samp <- runif(1, min = mean_ter - st0 / 2, max = mean_ter + st0 / 2)
  } else {
    ter_samp <- mean_ter
  }
  ter_trial <- max(0, ter_samp)
  
  # --- 2. Run the core diffusion process using trial-specific parameters & store path ---
  evidence <- z_trial
  current_decision_time <- 0
  max_steps <- max_decision_time / dt
  
  estimated_path_length <- as.integer(max_steps + 10)
  path_time_s <- numeric(estimated_path_length)
  path_evidence <- numeric(estimated_path_length)
  
  time_steps_taken <- 0
  path_idx <- 1
  
  path_time_s[path_idx] <- current_decision_time
  path_evidence[path_idx] <- evidence
  path_idx <- path_idx + 1
  
  while (evidence > 0 && evidence < a && time_steps_taken < max_steps) {
    increment <- rnorm(n = 1, mean = v_trial * dt, sd = s * sqrt(dt))
    evidence <- evidence + increment
    time_steps_taken <- time_steps_taken + 1
    current_decision_time <- time_steps_taken * dt
    
    if (path_idx <= estimated_path_length) {
      path_time_s[path_idx] <- current_decision_time
      path_evidence[path_idx] <- evidence
    } else {
      path_time_s <- c(path_time_s, current_decision_time) # Fallback grow
      path_evidence <- c(path_evidence, evidence)
    }
    path_idx <- path_idx + 1
  }
  
  decision_time_val <- current_decision_time
  actual_path_length <- path_idx - 1
  path_df_out <- data.frame(
    time_s = path_time_s[1:actual_path_length],
    evidence = path_evidence[1:actual_path_length]
  )
  
  choice_val <- NA
  rt_val <- NA
  
  if (time_steps_taken >= max_steps) {
    # Path data still recorded up to this point
  } else {
    if (evidence >= a) {
      choice_val <- 1
    } else { # evidence <= 0
      choice_val <- 0
    }
    rt_val <- decision_time_val + ter_trial
  }
  
  return(list(
    choice = choice_val,
    rt = rt_val,
    decision_time = if(is.na(choice_val)) NA else decision_time_val, # decision_time is NA if timeout
    v_trial = v_trial,
    z_trial = z_samp, # Return original sampled z before clipping for inspection
    ter_trial = ter_samp, # Return original sampled ter
    path_data = path_df_out 
  ))
}