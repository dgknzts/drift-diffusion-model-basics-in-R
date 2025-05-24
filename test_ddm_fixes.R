# Test script for DDM fitting fixes
# Load required libraries
library(dplyr)

# Source the functions
source("R/03_ddm_simulator_variable.R")
source("R/05_ddm_advanced_fitting.R")

cat("=== Testing DDM Fitting Fixes ===\n")

# Set up test parameters
true_params <- list(
  mean_v   = 0.18,
  a        = 0.5,
  s        = 0.2,
  mean_ter = 0.08,
  sv       = 0.20,
  sz       = 0.05,
  st0      = 0.04
)

# Generate test data
set.seed(42)
target_data <- simulate_diffusion_experiment_variable(
  n_trials = 1000,
  mean_v = true_params$mean_v,
  a = true_params$a,
  mean_z = true_params$a / 2,
  s = true_params$s,
  dt = 0.001,
  mean_ter = true_params$mean_ter,
  sv = true_params$sv,
  sz = true_params$sz,
  st0 = true_params$st0
)

cat("Generated", nrow(target_data), "trials\n")
cat("Proportion correct:", mean(target_data$choice == 1, na.rm = TRUE), "\n")

# Define RT bins
max_rt <- max(target_data$rt, na.rm = TRUE)
rt_bins <- seq(0, ceiling(max_rt * 1.1), by = 0.15)

# Calculate target binned proportions  
target_binned_props <- calculate_binned_rt_proportions(
  target_data, 
  rt_bins = rt_bins
)

cat("Number of bins:", nrow(target_binned_props), "\n")

# Test the diagnostic function
param_names <- c("mean_v", "a", "s", "mean_ter", "sv", "sz", "st0")
fixed_params <- list(dt = 0.001, correct_choice_value = 1, error_choice_value = 0)

cat("\n=== Running Diagnostic with Debug Info ===\n")
diagnostic_results <- diagnose_ddm_fitting(
  true_params = true_params,
  target_binned_props = target_binned_props,
  objective_fn_name = "ddm_binned_likelihood_objective",
  n_sim_per_eval = 500,  # Small for testing
  fixed_params = fixed_params,
  rt_bins = rt_bins,
  constrain_z_to_a_div_2 = TRUE,
  param_names_optim = param_names
)

# Also test a single simulation to see if the basic simulator works
cat("\n=== Testing Basic Simulation ===\n")
test_sim <- tryCatch({
  simulate_diffusion_experiment_variable(
    n_trials = 10,
    mean_v = 0.18,
    a = 0.5,
    mean_z = 0.25,
    s = 0.2,
    dt = 0.001,
    mean_ter = 0.08,
    sv = 0.20,
    sz = 0.05,
    st0 = 0.04
  )
}, error = function(e) {
  cat("Simulation failed:", conditionMessage(e), "\n")
  return(NULL)
})

if(!is.null(test_sim)) {
  cat("Basic simulation works! Generated", nrow(test_sim), "trials\n")
} else {
  cat("Basic simulation failed!\n")
}

cat("\n=== Test Results ===\n")
if(diagnostic_results$difference > 0) {
  cat("SUCCESS: Objective function working correctly!\n")
  cat("True params give better (lower) objective value than perturbed params.\n")
} else {
  cat("WARNING: Objective function may have issues.\n")
  cat("Difference should be positive, got:", diagnostic_results$difference, "\n")
}

cat("\n=== Test Complete ===\n") 