# Minimal test for DDM fitting
cat("=== Minimal DDM Test ===\n")

# Test 1: Load libraries
cat("1. Loading dplyr...\n")
if(require(dplyr, quietly = TRUE)) {
  cat("   dplyr loaded successfully\n")
} else {
  cat("   ERROR: dplyr not available\n")
}

# Test 2: Load simulator
cat("2. Loading simulator...\n")
tryCatch({
  source("R/03_ddm_simulator_variable.R")
  cat("   Simulator loaded successfully\n")
}, error = function(e) {
  cat("   ERROR loading simulator:", conditionMessage(e), "\n")
})

# Test 3: Test basic simulation
cat("3. Testing basic simulation...\n")
test_sim <- tryCatch({
  simulate_diffusion_experiment_variable(
    n_trials = 5,
    mean_v = 0.2,
    a = 0.8,
    mean_z = 0.4,
    s = 0.1,
    dt = 0.001,
    mean_ter = 0.1
  )
}, error = function(e) {
  cat("   ERROR in simulation:", conditionMessage(e), "\n")
  return(NULL)
})

if(!is.null(test_sim) && nrow(test_sim) > 0) {
  cat("   Basic simulation WORKS! Generated", nrow(test_sim), "trials\n")
  cat("   Sample data:\n")
  print(head(test_sim, 2))
} else {
  cat("   Basic simulation FAILED\n")
}

# Test 4: Load advanced functions  
cat("4. Loading advanced fitting functions...\n")
tryCatch({
  source("R/05_ddm_advanced_fitting.R")
  cat("   Advanced functions loaded successfully\n")
}, error = function(e) {
  cat("   ERROR loading advanced functions:", conditionMessage(e), "\n")
})

# Test 5: Test binned proportions calculation
if(!is.null(test_sim) && nrow(test_sim) > 0) {
  cat("5. Testing binned proportions calculation...\n")
  rt_bins <- seq(0, 2, by = 0.2)
  
  binned_props <- tryCatch({
    calculate_binned_rt_proportions(test_sim, rt_bins = rt_bins)
  }, error = function(e) {
    cat("   ERROR in binned proportions:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if(!is.null(binned_props) && nrow(binned_props) > 0) {
    cat("   Binned proportions calculation WORKS!\n")
    cat("   Number of bins:", nrow(binned_props), "\n")
  } else {
    cat("   Binned proportions calculation FAILED\n")
  }
}

cat("\n=== Minimal Test Complete ===\n") 