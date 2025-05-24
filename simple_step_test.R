# Simple step-by-step test
cat("=== Simple Step-by-Step Test ===\n")

library(dplyr)
source("R/03_ddm_simulator_variable.R")
source("R/05_ddm_advanced_fitting.R")

cat("Step 1: Test basic simulation directly\n")
basic_sim <- simulate_diffusion_experiment_variable(
  n_trials = 10,
  mean_v = 0.2,
  a = 0.8,
  mean_z = 0.4,
  s = 0.15,
  dt = 0.001,
  mean_ter = 0.1,
  sv = 0.1,
  sz = 0.05,
  st0 = 0.02
)
cat("   Basic simulation result: NROW =", nrow(basic_sim), "\n")

cat("Step 2: Test binned proportions calculation\n")
rt_bins <- seq(0, 2, by = 0.2)
binned_result <- calculate_binned_rt_proportions(basic_sim, rt_bins = rt_bins)
cat("   Binned proportions result: NROW =", nrow(binned_result), "\n")

cat("Step 3: Test objective function with simple parameters\n")
# Use the same parameters as simulation for perfect match
test_params <- c(mean_v = 0.2, a = 0.8, s = 0.15, mean_ter = 0.1, sv = 0.1, sz = 0.05, st0 = 0.02)
fixed_params_simple <- list(dt = 0.001)  # Only essential fixed param

obj_result <- ddm_binned_likelihood_objective(
  params_to_test = test_params,
  target_binned_props = binned_result,
  param_names_optim = names(test_params),
  n_sim_per_eval = 50,  # Very small for testing
  fixed_params = fixed_params_simple,
  rt_bins = rt_bins,
  constrain_z_to_a_div_2 = TRUE,
  debug = TRUE
)

cat("   Objective function result:", obj_result, "\n")

if(obj_result < 1e6) {
  cat("✅ SUCCESS: Objective function is working!\n")
} else {
  cat("❌ PROBLEM: Still getting large error value\n")
}

cat("=== Simple Test Complete ===\n") 