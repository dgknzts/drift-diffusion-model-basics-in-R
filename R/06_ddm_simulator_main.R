#' Simulate a Single DDM Trial with Full Parameter Set
#'
#' This function simulates one trial of the Diffusion Decision Model (DDM),
#' incorporating trial-to-trial variability in key parameters: drift rate (`v`),
#' starting point (`z`), and non-decision time (`ter`). It also includes the
#' response execution bias parameter (`d`) which allows for different non-decision
#' times between upper and lower threshold responses.
#'
#' @details
#' The evidence accumulation process follows the standard DDM:
#' evidence starts at `z_trial` and accumulates with a mean rate of `v_trial`
#' and within-trial noise `s` until it reaches either the lower boundary (0)
#' or the upper boundary (`a`). The time step `dt` controls the granularity
#' of the simulation.
#'
#' Parameter Variability:
#' \itemize{
#'   \item **Drift Rate (`v_trial`):** Sampled from a Normal distribution:
#'     `v_trial ~ Normal(mean_v, sv)`.
#'   \item **Starting Point (`z_trial`):** Sampled from a Uniform distribution:
#'     `z_trial ~ Uniform(mean_z - sz/2, mean_z + sz/2)`.
#'   \item **Non-Decision Time (`ter_trial`):** Base value sampled from a Uniform 
#'     distribution: `ter_base ~ Uniform(mean_ter - st0/2, mean_ter + st0/2)`.
#'     Then adjusted by `d`: `ter_upper = ter_base - 0.5*d`, `ter_lower = ter_base + 0.5*d`.
#' }
#'
#' @param mean_v Numeric. The mean of the distribution from which the trial-specific
#'   drift rate is sampled.
#' @param a Numeric. The threshold separation (upper boundary). Must be positive.
#' @param mean_z Numeric. The mean of the distribution from which the trial-specific
#'   starting point is sampled. Must be strictly between 0 and `a`.
#' @param s Numeric, optional. The standard deviation of the within-trial Gaussian
#'   noise. Default is `0.1`.
#' @param dt Numeric, optional. The time step (in seconds) for the discrete
#'   approximation. Default is `0.001`.
#' @param mean_ter Numeric, optional. The mean of the distribution from which the
#'   base non-decision time is sampled. Default is `0.1`.
#' @param d Numeric, optional. The difference in non-decision time between responses
#'   at the upper vs lower threshold. Default is `0`.
#' @param sv Numeric, optional. The standard deviation of the across-trial Normal
#'   distribution for drift rate. Default is `0`.
#' @param sz Numeric, optional. The total range of the across-trial Uniform
#'   distribution for starting point. Default is `0`.
#' @param st0 Numeric, optional. The total range of the across-trial Uniform
#'   distribution for non-decision time. Default is `0`.
#' @param max_decision_time Numeric, optional. The maximum duration (in seconds)
#'   allowed for the evidence accumulation process. Default is `5.0`.
#'
#' @return A list containing:
#'   \item{choice}{Integer: `1` if upper boundary reached, `0` if lower, `NA` if timeout.}
#'   \item{rt}{Numeric: The total reaction time (decision_time + ter_trial) in seconds.}
#'   \item{decision_time}{Numeric: The time for evidence accumulation in seconds.}
#'   \item{v_trial}{Numeric: The drift rate used for this trial.}
#'   \item{z_trial}{Numeric: The starting point used for this trial.}
#'   \item{ter_trial}{Numeric: The actual non-decision time used for this trial.}
#' @export
#'
#' @examples
#' set.seed(123)
#' # Single trial with d parameter (upper responses 25ms faster)
#' trial <- simulate_diffusion_trial(
#'   mean_v = 0.2, a = 1.0, mean_z = 0.5, mean_ter = 0.15,
#'   d = 0.05, sv = 0.1, sz = 0, st0 = 0
#' )
#' print(trial)
simulate_diffusion_trial <- function(mean_v,
                                     a,
                                     mean_z,
                                     s = 0.1,
                                     dt = 0.001,
                                     mean_ter = 0.1,
                                     d = 0,
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
  if (z_trial_clipped <= 0 || z_trial_clipped >= a || is.na(z_trial_clipped)) { 
    z_trial_clipped <- mean_z 
    if (z_trial_clipped <= 0 || z_trial_clipped >= a || is.na(z_trial_clipped)) {
      z_trial_clipped <- a/2 
    }
  }
  z_trial <- z_trial_clipped
  
  # Sample base non-decision time
  if (st0 > 0) {
    ter_base <- runif(1, min = mean_ter - st0 / 2, max = mean_ter + st0 / 2)
  } else {
    ter_base <- mean_ter
  }
  ter_base <- max(0, ter_base)
  
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
  ter_trial <- NA
  
  if (time_steps_taken >= max_steps) {
    choice_val <- NA
    rt_val <- NA
    decision_time_out <- NA
  } else {
    # Apply d parameter to determine response-specific non-decision time
    if (evidence >= a) {
      choice_val <- 1
      ter_trial <- ter_base - 0.5 * d  # Upper threshold
    } else { 
      choice_val <- 0
      ter_trial <- ter_base + 0.5 * d  # Lower threshold
    }
    ter_trial <- max(0, ter_trial)  # Ensure non-negative
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


#' Simulate Multiple DDM Trials with Full Parameter Set
#'
#' This function serves as a wrapper around `simulate_diffusion_trial`
#' to conduct a full DDM experiment by simulating a specified number of trials.
#' It collects the outcomes (choice, RT, decision time) as well as the
#' trial-specific parameters that were sampled and used for each individual trial.
#'
#' @param n_trials Integer. The total number of DDM trials to simulate. Default is `100`.
#' @param mean_v Numeric. The mean drift rate.
#' @param a Numeric. The threshold separation.
#' @param mean_z Numeric. The mean starting point.
#' @param s Numeric, optional. Within-trial noise SD. Default `0.1`.
#' @param dt Numeric, optional. Simulation time step. Default `0.001`.
#' @param mean_ter Numeric, optional. The mean non-decision time. Default `0.1`.
#' @param d Numeric, optional. Difference in non-decision time. Default `0`.
#' @param sv Numeric, optional. SD for across-trial drift variability. Default `0`.
#' @param sz Numeric, optional. Range for across-trial starting point variability. Default `0`.
#' @param st0 Numeric, optional. Range for across-trial non-decision time variability. Default `0`.
#' @param max_decision_time Numeric, optional. Max decision time. Default `5.0`.
#'
#' @return A data frame where each row represents a simulated trial. Columns include:
#'   trial, choice, rt, decision_time, v_trial, z_trial, ter_trial.
#' @export
#'
#' @examples
#' set.seed(456)
#' # Simulate experiment with d parameter
#' data <- simulate_diffusion_experiment(
#'   n_trials = 200,
#'   mean_v = 0.15, a = 1.0, mean_z = 0.5, mean_ter = 0.2,
#'   d = 0.02, sv = 0.1, sz = 0.05, st0 = 0.03
#' )
#' # Check mean RTs by response
#' aggregate(rt ~ choice, data = data, mean)
simulate_diffusion_experiment <- function(n_trials = 100,
                                          mean_v,
                                          a,
                                          mean_z,
                                          s = 0.1,
                                          dt = 0.001,
                                          mean_ter = 0.1,
                                          d = 0,
                                          sv = 0,
                                          sz = 0,
                                          st0 = 0,
                                          max_decision_time = 5.0) {
  
  # Pre-allocate results list
  results <- vector("list", n_trials)
  
  # Run trials
  for (i in 1:n_trials) {
    trial_result <- simulate_diffusion_trial(
      mean_v = mean_v,
      a = a,
      mean_z = mean_z,
      s = s,
      dt = dt,
      mean_ter = mean_ter,
      d = d,
      sv = sv,
      sz = sz,
      st0 = st0,
      max_decision_time = max_decision_time
    )
    
    results[[i]] <- data.frame(
      trial = i,
      choice = trial_result$choice,
      rt = trial_result$rt,
      decision_time = trial_result$decision_time,
      v_trial = trial_result$v_trial,
      z_trial = trial_result$z_trial,
      ter_trial = trial_result$ter_trial
    )
  }
  
  # Combine all trials into single data frame
  do.call(rbind, results)
}


#' Simulate a Single DDM Trial with Path Recording
#'
#' This function extends simulate_diffusion_trial by also recording the
#' complete evidence accumulation path. Useful for visualization and
#' understanding the dynamics of the decision process.
#'
#' @inheritParams simulate_diffusion_trial
#'
#' @return A list containing all outputs from simulate_diffusion_trial plus:
#'   \item{path_data}{Data frame with columns time_s and evidence showing the accumulation path.}
#' @export
#'
#' @examples
#' set.seed(789)
#' # Single trial with path recording
#' trial_with_path <- simulate_diffusion_trial_with_path(
#'   mean_v = 0.2, a = 1.0, mean_z = 0.5, mean_ter = 0.15,
#'   d = 0.05, sv = 0.1
#' )
#' plot(trial_with_path$path_data$time_s, trial_with_path$path_data$evidence, 
#'      type = "l", xlab = "Time (s)", ylab = "Evidence")
#' abline(h = c(0, 1.0), lty = 2)
simulate_diffusion_trial_with_path <- function(mean_v,
                                               a,
                                               mean_z,
                                               s = 0.1,
                                               dt = 0.001,
                                               mean_ter = 0.1,
                                               d = 0,
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
    ter_base <- runif(1, min = mean_ter - st0 / 2, max = mean_ter + st0 / 2)
  } else {
    ter_base <- mean_ter
  }
  ter_base <- max(0, ter_base)
  
  # --- 2. Run the core diffusion process with path recording ---
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
      path_time_s <- c(path_time_s, current_decision_time)
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
  ter_trial <- NA
  
  if (time_steps_taken >= max_steps) {
    # Timeout
  } else {
    # Apply d parameter
    if (evidence >= a) {
      choice_val <- 1
      ter_trial <- ter_base - 0.5 * d
    } else {
      choice_val <- 0
      ter_trial <- ter_base + 0.5 * d
    }
    ter_trial <- max(0, ter_trial)
    rt_val <- decision_time_val + ter_trial
  }
  
  return(list(
    choice = choice_val,
    rt = rt_val,
    decision_time = if(is.na(choice_val)) NA else decision_time_val,
    v_trial = v_trial,
    z_trial = z_samp,
    ter_trial = ter_trial,
    path_data = path_df_out
  ))
}