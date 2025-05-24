# Test the DDM objective function specifically
cat("=== Testing DDM Objective Function (UPDATED) ===\n")

# Load required components
library(dplyr)
source("R/03_ddm_simulator_variable.R")
source("R/05_ddm_advanced_fitting.R")  # Reload with fixes

# Generate target data with known parameters
set.seed(123)
true_params <- list(
  mean_v = 0.2,
  a = 0.8,
  s = 0.15,
  mean_ter = 0.1,
  sv = 0.1,
  sz = 0.05,
  st0 = 0.02
)

cat("1. Generating target data...\n")
target_data <- simulate_diffusion_experiment_variable(
  n_trials = 200,  # Smaller for testing
  mean_v = true_params$mean_v,
  a = true_params$a,
  mean_z = true_params$a / 2,  # Constraint: mean_z = a/2
  s = true_params$s,
  dt = 0.001,
  mean_ter = true_params$mean_ter,
  sv = true_params$sv,
  sz = true_params$sz,
  st0 = true_params$st0
)

cat("   Generated", nrow(target_data), "trials\n")
cat("   Proportion correct:", round(mean(target_data$choice == 1, na.rm = TRUE), 3), "\n")

# Calculate target binned proportions
rt_bins <- seq(0, max(target_data$rt, na.rm = TRUE) * 1.2, by = 0.15)
target_binned_props <- calculate_binned_rt_proportions(target_data, rt_bins = rt_bins)

cat("   Number of bins:", nrow(target_binned_props), "\n")

# Test objective function with TRUE parameters
cat("2. Testing objective function with TRUE parameters (FIXED VERSION)...\n")
param_names <- c("mean_v", "a", "s", "mean_ter", "sv", "sz", "st0")
true_values <- sapply(param_names, function(p) true_params[[p]])
fixed_params <- list(dt = 0.001, correct_choice_value = 1, error_choice_value = 0)

obj_value_true <- ddm_binned_likelihood_objective(
  params_to_test = true_values,
  target_binned_props = target_binned_props,
  param_names_optim = param_names,
  n_sim_per_eval = 200,  # Small for testing
  fixed_params = fixed_params,
  rt_bins = rt_bins,
  constrain_z_to_a_div_2 = TRUE,
  debug = TRUE  # Enable debugging
)

cat("   Objective value with TRUE parameters:", obj_value_true, "\n")

# Only test perturbed if the true parameters worked
if(obj_value_true < 1e6) {
  # Test with slightly worse parameters
  cat("3. Testing with perturbed parameters...\n")
  perturbed_values <- true_values * c(1.2, 0.9, 1.1, 1.15, 0.8, 1.3, 1.4)  # Make them worse
  
  obj_value_perturbed <- ddm_binned_likelihood_objective(
    params_to_test = perturbed_values,
    target_binned_props = target_binned_props,
    param_names_optim = param_names,
    n_sim_per_eval = 200,
    fixed_params = fixed_params,
    rt_bins = rt_bins,
    constrain_z_to_a_div_2 = TRUE,
    debug = FALSE  # Less verbose for second test
  )
  
  cat("   Objective value with PERTURBED parameters:", obj_value_perturbed, "\n")
  
  # Check if objective function is working correctly
  cat("4. Results analysis...\n")
  difference <- obj_value_perturbed - obj_value_true
  cat("   Difference (perturbed - true):", round(difference, 4), "\n")
  
  if(difference > 0) {
    cat("   ✅ SUCCESS: Objective function working correctly!\n")
    cat("      True parameters give better (lower) objective value\n")
  } else if(difference < 0) {
    cat("   ⚠️  WARNING: True parameters give worse objective value than perturbed\n")
  } else {
    cat("   ⚠️  WARNING: Both give identical objective values\n")
  }
} else {
  cat("   ❌ PROBLEM: True parameters still returning large error value\n")
  cat("   Need to investigate further...\n")
}

cat("\n=== Objective Function Test Complete ===\n") 