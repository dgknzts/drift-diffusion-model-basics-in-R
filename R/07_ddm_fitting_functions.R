# DDM Fitting Functions: Maximum Likelihood, Chi-Square, and Kolmogorov-Smirnov Methods
# This script provides functions for fitting DDM parameters using three different
# optimization criteria, as described in Voss, Voss & Lerche (2015)

#' Calculate DDM Density for Maximum Likelihood
#'
#' Calculates the probability density for a given RT and response using the
#' Navarro & Fuss (2009) method for numerical stability.
#'
#' @param rt Numeric. Response time (must be > ter).
#' @param response Integer. Response (0 = lower, 1 = upper threshold).
#' @param v Numeric. Drift rate.
#' @param a Numeric. Threshold separation.
#' @param z Numeric. Starting point.
#' @param s Numeric. Diffusion coefficient (usually 0.1 or 1).
#' @param ter Numeric. Non-decision time.
#' @param d Numeric. Difference in non-decision time.
#' @param sv Numeric. Inter-trial variability in drift.
#' @param sz Numeric. Inter-trial variability in starting point.
#' @param st0 Numeric. Inter-trial variability in non-decision time.
#'
#' @return Numeric. Probability density.
#' @export
ddm_density <- function(rt, response, v, a, z, s = 0.1, ter = 0.1, d = 0,
                        sv = 0, sz = 0, st0 = 0) {
  
  # Adjust ter based on response and d parameter
  if (response == 1) {
    ter_adj <- ter - 0.5 * d
  } else {
    ter_adj <- ter + 0.5 * d
  }
  
  # Decision time
  t <- rt - ter_adj
  
  # Return very small density if RT is impossible
  if (t <= 0) return(1e-10)
  
  # For lower threshold, flip parameters
  if (response == 0) {
    v <- -v
    z <- a - z
  }
  
  # Basic density without variability
  if (sv == 0 && sz == 0 && st0 == 0) {
    # Use Navarro & Fuss (2009) switching method
    if (t < 0.5) {
      # Small-time expansion
      k_max <- max(2, ceiling(sqrt(-2 * log(1e-10) * s^2 * t) / a))
      dens <- 0
      for (k in -k_max:k_max) {
        w <- (2 * k * a + z)
        dens <- dens + w * exp(-w^2 / (2 * s^2 * t))
      }
      dens <- dens * (s / sqrt(2 * pi * t^3)) * exp((v * z - v^2 * s^2 * t / 2) / s^2)
    } else {
      # Large-time expansion
      k_max <- max(1, ceiling(sqrt(-2 * log(1e-10) * a^2 / (pi^2 * s^2 * t))))
      dens <- 0
      for (k in 1:k_max) {
        dens <- dens + k * sin(k * pi * z / a) * exp(-k^2 * pi^2 * s^2 * t / (2 * a^2))
      }
      dens <- dens * (pi * s^2 / a^2) * exp((v * z - v^2 * s^2 * t / 2) / s^2)
    }
    
    return(max(dens, 1e-10))
  }
  
  # With variability - numerical integration
  # This is a simplified implementation
  n_points <- 15
  dens_sum <- 0
  
  if (sv > 0) {
    # Integrate over drift variability
    v_points <- seq(v - 3*sv, v + 3*sv, length.out = n_points)
    v_weights <- dnorm(v_points, v, sv)
    v_weights <- v_weights / sum(v_weights)
    
    for (i in 1:n_points) {
      dens_sum <- dens_sum + v_weights[i] * 
        ddm_density(rt, response, v_points[i], a, z, s, ter, d, 0, sz, st0)
    }
    return(dens_sum)
  }
  
  if (sz > 0) {
    # Integrate over starting point variability
    z_points <- seq(max(1e-6, z - sz/2), min(a - 1e-6, z + sz/2), length.out = n_points)
    z_weights <- rep(1/n_points, n_points)
    
    for (i in 1:n_points) {
      dens_sum <- dens_sum + z_weights[i] * 
        ddm_density(rt, response, v, a, z_points[i], s, ter, d, 0, 0, st0)
    }
    return(dens_sum)
  }
  
  if (st0 > 0) {
    # Integrate over non-decision time variability
    ter_points <- seq(max(0, ter - st0/2), ter + st0/2, length.out = n_points)
    ter_weights <- rep(1/n_points, n_points)
    
    for (i in 1:n_points) {
      dens_sum <- dens_sum + ter_weights[i] * 
        ddm_density(rt, response, v, a, z, s, ter_points[i], d, 0, 0, 0)
    }
    return(dens_sum)
  }
}


#' Maximum Likelihood Objective Function
#'
#' Calculates negative log-likelihood for DDM parameters given data.
#'
#' @param params Numeric vector of parameters being optimized.
#' @param data Data frame with 'rt' and 'choice' columns.
#' @param param_names Character vector of parameter names.
#' @param fixed_params Named list of fixed parameters.
#' @param constrain_z_to_a_div_2 Logical. Constrain z = a/2.
#'
#' @return Numeric. Negative log-likelihood.
#' @export
ddm_objective_ml <- function(params, data, param_names, fixed_params = list(),
                             constrain_z_to_a_div_2 = FALSE) {
  
  # Combine parameters
  all_params <- fixed_params
  for (i in seq_along(param_names)) {
    all_params[[param_names[i]]] <- params[i]
  }
  
  # Handle constraints
  if (constrain_z_to_a_div_2 && "a" %in% names(all_params)) {
    all_params[["z"]] <- all_params[["a"]] / 2
  }
  
  # Ensure we have mean_v instead of v, etc.
  param_mapping <- list(
    "mean_v" = "v", "mean_z" = "z", "mean_ter" = "ter"
  )
  
  for (long_name in names(param_mapping)) {
    short_name <- param_mapping[[long_name]]
    if (long_name %in% names(all_params) && !(short_name %in% names(all_params))) {
      all_params[[short_name]] <- all_params[[long_name]]
    }
  }
  
  # Set defaults
  defaults <- list(v = 0, a = 1, z = 0.5, s = 0.1, ter = 0.1, 
                   d = 0, sv = 0, sz = 0, st0 = 0)
  for (p in names(defaults)) {
    if (!(p %in% names(all_params))) {
      all_params[[p]] <- defaults[[p]]
    }
  }
  
  # Calculate log-likelihood
  ll <- 0
  valid_trials <- 0
  
  for (i in 1:nrow(data)) {
    if (!is.na(data$rt[i]) && !is.na(data$choice[i])) {
      dens <- ddm_density(
        rt = data$rt[i],
        response = data$choice[i],
        v = all_params$v,
        a = all_params$a,
        z = all_params$z,
        s = all_params$s,
        ter = all_params$ter,
        d = all_params$d,
        sv = all_params$sv,
        sz = all_params$sz,
        st0 = all_params$st0
      )
      
      ll <- ll + log(max(dens, 1e-10))
      valid_trials <- valid_trials + 1
    }
  }
  
  # Return negative for minimization
  if (valid_trials < 10) return(1e6)
  return(-ll)
}


#' Chi-Square Objective Function
#'
#' Calculates chi-square statistic comparing observed and predicted RT quantiles.
#' Uses the method from Ratcliff & Tuerlinckx (2002).
#'
#' @inheritParams ddm_objective_ml
#' @param n_sim Number of trials to simulate for predictions.
#' @param quantiles Quantiles to use for binning (default: 0.1, 0.3, 0.5, 0.7, 0.9).
#'
#' @return Numeric. Chi-square statistic.
#' @export
ddm_objective_cs <- function(params, data, param_names, fixed_params = list(),
                             constrain_z_to_a_div_2 = FALSE, n_sim = 2000,
                             quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  
  # Combine parameters (same as ML)
  all_params <- fixed_params
  for (i in seq_along(param_names)) {
    all_params[[param_names[i]]] <- params[i]
  }
  
  if (constrain_z_to_a_div_2 && "a" %in% names(all_params)) {
    all_params[["mean_z"]] <- all_params[["a"]] / 2
  }
  
  # Ensure all needed parameters exist with proper names
  param_mapping <- list(
    "v" = "mean_v", "z" = "mean_z", "ter" = "mean_ter"
  )
  
  for (short_name in names(param_mapping)) {
    long_name <- param_mapping[[short_name]]
    if (short_name %in% names(all_params) && !(long_name %in% names(all_params))) {
      all_params[[long_name]] <- all_params[[short_name]]
    }
  }
  
  # Simulate data
  sim_args <- c(list(n_trials = n_sim), all_params)
  sim_data <- tryCatch({
    do.call(simulate_diffusion_experiment, sim_args)
  }, error = function(e) NULL)
  
  if (is.null(sim_data)) return(1e6)
  
  # Calculate observed bins
  obs_data <- data[!is.na(data$rt) & !is.na(data$choice), ]
  obs_0 <- obs_data$rt[obs_data$choice == 0]
  obs_1 <- obs_data$rt[obs_data$choice == 1]
  
  # Need minimum responses for each type
  if (length(obs_0) < 12 || length(obs_1) < 12) return(1e6)
  
  # Get quantile boundaries
  q_bounds_0 <- quantile(obs_0, probs = quantiles)
  q_bounds_1 <- quantile(obs_1, probs = quantiles)
  
  # Count observations in bins
  obs_counts_0 <- hist(obs_0, breaks = c(-Inf, q_bounds_0, Inf), plot = FALSE)$counts
  obs_counts_1 <- hist(obs_1, breaks = c(-Inf, q_bounds_1, Inf), plot = FALSE)$counts
  
  # Count predictions in bins
  sim_0 <- sim_data$rt[sim_data$choice == 0 & !is.na(sim_data$rt)]
  sim_1 <- sim_data$rt[sim_data$choice == 1 & !is.na(sim_data$rt)]
  
  if (length(sim_0) < 5 || length(sim_1) < 5) return(1e6)
  
  pred_counts_0 <- hist(sim_0, breaks = c(-Inf, q_bounds_0, Inf), plot = FALSE)$counts
  pred_counts_1 <- hist(sim_1, breaks = c(-Inf, q_bounds_1, Inf), plot = FALSE)$counts
  
  # Scale predictions to match number of observations
  pred_counts_0 <- pred_counts_0 * length(obs_0) / length(sim_0)
  pred_counts_1 <- pred_counts_1 * length(obs_1) / length(sim_1)
  
  # Calculate chi-square
  cs <- 0
  for (i in 1:length(obs_counts_0)) {
    if (pred_counts_0[i] > 0) {
      cs <- cs + (obs_counts_0[i] - pred_counts_0[i])^2 / pred_counts_0[i]
    }
  }
  for (i in 1:length(obs_counts_1)) {
    if (pred_counts_1[i] > 0) {
      cs <- cs + (obs_counts_1[i] - pred_counts_1[i])^2 / pred_counts_1[i]
    }
  }
  
  return(cs)
}


#' Kolmogorov-Smirnov Objective Function
#'
#' Calculates KS statistic comparing observed and predicted CDFs.
#' Uses the combined CDF method from Voss & Voss (2007).
#'
#' @inheritParams ddm_objective_cs
#'
#' @return Numeric. Negative log of KS p-value (for minimization).
#' @export
ddm_objective_ks <- function(params, data, param_names, fixed_params = list(),
                             constrain_z_to_a_div_2 = FALSE, n_sim = 2000) {
  
  # Combine parameters (same as CS)
  all_params <- fixed_params
  for (i in seq_along(param_names)) {
    all_params[[param_names[i]]] <- params[i]
  }
  
  if (constrain_z_to_a_div_2 && "a" %in% names(all_params)) {
    all_params[["mean_z"]] <- all_params[["a"]] / 2
  }
  
  # Ensure proper parameter names
  param_mapping <- list(
    "v" = "mean_v", "z" = "mean_z", "ter" = "mean_ter"
  )
  
  for (short_name in names(param_mapping)) {
    long_name <- param_mapping[[short_name]]
    if (short_name %in% names(all_params) && !(long_name %in% names(all_params))) {
      all_params[[long_name]] <- all_params[[short_name]]
    }
  }
  
  # Simulate data
  sim_args <- c(list(n_trials = n_sim), all_params)
  sim_data <- tryCatch({
    do.call(simulate_diffusion_experiment, sim_args)
  }, error = function(e) NULL)
  
  if (is.null(sim_data)) return(1e6)
  
  # Create combined distributions (negative RTs for lower threshold)
  obs_combined <- numeric(nrow(data))
  sim_combined <- numeric(nrow(sim_data))
  
  obs_valid <- !is.na(data$rt) & !is.na(data$choice)
  obs_combined[obs_valid] <- ifelse(data$choice[obs_valid] == 0, 
                                    -data$rt[obs_valid], 
                                    data$rt[obs_valid])
  
  sim_valid <- !is.na(sim_data$rt) & !is.na(sim_data$choice)
  sim_combined[sim_valid] <- ifelse(sim_data$choice[sim_valid] == 0, 
                                    -sim_data$rt[sim_valid], 
                                    sim_data$rt[sim_valid])
  
  # Remove NAs
  obs_combined <- obs_combined[obs_valid]
  sim_combined <- sim_combined[sim_valid]
  
  if (length(obs_combined) < 10 || length(sim_combined) < 10) return(1e6)
  
  # Calculate KS statistic
  ks_result <- ks.test(obs_combined, sim_combined)
  
  # Return negative log p-value for minimization
  return(-log(max(ks_result$p.value, 1e-10)))
}


#' Fit DDM Using Specified Method
#'
#' General fitting function that can use ML, CS, or KS methods.
#'
#' @param data Data frame with 'rt' and 'choice' columns.
#' @param method Character. One of "ML", "CS", or "KS".
#' @param params_to_fit Character vector of parameters to estimate.
#' @param initial_values Named vector of starting values.
#' @param fixed_params Named list of fixed parameters.
#' @param lower_bounds Named vector of lower bounds.
#' @param upper_bounds Named vector of upper bounds.
#' @param constrain_z_to_a_div_2 Logical. Constrain z = a/2.
#' @param n_sim For CS and KS methods, number of simulations.
#' @param control List of control parameters for optim().
#'
#' @return List with optimization results and method-specific information.
#' @export
fit_ddm <- function(data, 
                    method = c("ML", "CS", "KS"),
                    params_to_fit = c("mean_v", "a", "mean_ter"),
                    initial_values = NULL,
                    fixed_params = list(s = 0.1, d = 0, sv = 0, sz = 0, st0 = 0),
                    lower_bounds = NULL,
                    upper_bounds = NULL,
                    constrain_z_to_a_div_2 = TRUE,
                    n_sim = 2000,
                    control = list(maxit = 100, trace = 1)) {
  
  method <- match.arg(method)
  
  # Set defaults for initial values
  if (is.null(initial_values)) {
    defaults <- list(
      mean_v = 0.1, a = 1.0, mean_z = 0.5, mean_ter = 0.3,
      d = 0, sv = 0.1, sz = 0.1, st0 = 0.05
    )
    initial_values <- unlist(defaults[params_to_fit])
  }
  
  # Set defaults for bounds
  if (is.null(lower_bounds)) {
    lower_defaults <- list(
      mean_v = -2, a = 0.1, mean_z = 0.01, mean_ter = 0,
      d = -0.2, sv = 0, sz = 0, st0 = 0
    )
    lower_bounds <- unlist(lower_defaults[params_to_fit])
  }
  
  if (is.null(upper_bounds)) {
    upper_defaults <- list(
      mean_v = 2, a = 3, mean_z = 0.99, mean_ter = 1,
      d = 0.2, sv = 1, sz = 0.5, st0 = 0.5
    )
    upper_bounds <- unlist(upper_defaults[params_to_fit])
  }
  
  # Select objective function
  obj_fn <- switch(method,
                   ML = ddm_objective_ml,
                   CS = function(...) ddm_objective_cs(..., n_sim = n_sim),
                   KS = function(...) ddm_objective_ks(..., n_sim = n_sim)
  )
  
  # Run optimization
  start_time <- Sys.time()
  
  result <- optim(
    par = initial_values,
    fn = obj_fn,
    data = data,
    param_names = params_to_fit,
    fixed_params = fixed_params,
    constrain_z_to_a_div_2 = constrain_z_to_a_div_2,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = control
  )
  
  end_time <- Sys.time()
  
  # Add names to parameters
  names(result$par) <- params_to_fit
  
  # Add method-specific information
  result$method <- method
  result$fitting_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Calculate fit statistics based on method
  n_obs <- nrow(data[!is.na(data$rt) & !is.na(data$choice), ])
  n_params <- length(params_to_fit)
  
  if (method == "ML") {
    result$log_likelihood <- -result$value
    result$AIC <- 2 * result$value + 2 * n_params
    result$BIC <- 2 * result$value + n_params * log(n_obs)
  } else if (method == "CS") {
    result$chi_square <- result$value
    result$df <- 12 - n_params  # 6 bins x 2 responses - parameters
  } else if (method == "KS") {
    result$ks_p_value <- exp(-result$value)
  }
  
  class(result) <- c("ddm_fit", class(result))
  return(result)
}


#' Print Method for DDM Fit Results
#'
#' @param x Object of class ddm_fit.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.ddm_fit <- function(x, ...) {
  cat("\nDDM Fit Results\n")
  cat("===============\n")
  cat("Method:", x$method, "\n")
  cat("Convergence:", ifelse(x$convergence == 0, "Success", "Failed"), "\n")
  cat("Fitting time:", round(x$fitting_time, 2), "seconds\n\n")
  
  cat("Parameter Estimates:\n")
  print(round(x$par, 4))
  
  cat("\nFit Statistics:\n")
  if (x$method == "ML") {
    cat("Log-likelihood:", round(x$log_likelihood, 2), "\n")
    cat("AIC:", round(x$AIC, 2), "\n")
    cat("BIC:", round(x$BIC, 2), "\n")
  } else if (x$method == "CS") {
    cat("Chi-square:", round(x$chi_square, 2), "\n")
    cat("df:", x$df, "\n")
  } else if (x$method == "KS") {
    cat("KS p-value:", round(x$ks_p_value, 4), "\n")
  }
}


#' Compare Multiple Fitting Methods
#'
#' Fits the same data using all three methods for comparison.
#'
#' @inheritParams fit_ddm
#'
#' @return List containing results from all three methods.
#' @export
compare_fitting_methods <- function(data, 
                                    params_to_fit = c("mean_v", "a", "mean_ter"),
                                    initial_values = NULL,
                                    fixed_params = list(s = 0.1, d = 0, sv = 0, sz = 0, st0 = 0),
                                    constrain_z_to_a_div_2 = TRUE,
                                    n_sim = 2000,
                                    control = list(maxit = 100, trace = 0)) {
  
  methods <- c("ML", "CS", "KS")
  results <- list()
  
  cat("Comparing DDM fitting methods...\n")
  
  for (method in methods) {
    cat("\nFitting with", method, "method...\n")
    
    results[[method]] <- fit_ddm(
      data = data,
      method = method,
      params_to_fit = params_to_fit,
      initial_values = initial_values,
      fixed_params = fixed_params,
      constrain_z_to_a_div_2 = constrain_z_to_a_div_2,
      n_sim = n_sim,
      control = control
    )
  }
  
  # Create comparison table
  comparison <- data.frame(
    Method = methods,
    Convergence = sapply(results, function(x) x$convergence == 0),
    Time_seconds = sapply(results, function(x) round(x$fitting_time, 2))
  )
  
  # Add parameter estimates
  for (param in params_to_fit) {
    comparison[[param]] <- sapply(results, function(x) round(x$par[param], 4))
  }
  
  cat("\n\nComparison Summary:\n")
  print(comparison)
  
  return(list(
    results = results,
    comparison = comparison
  ))
}