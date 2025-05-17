# =============================================
# DDM_Basics_R Project Validation Script
# =============================================
# This script checks all components of the DDM_Basics_R package
# including R functions, utilities, and vignettes

# Helper function for colorful output
print_status <- function(message, status) {
  status_color <- switch(status,
                         "OK" = "\033[32m", # green
                         "WARNING" = "\033[33m", # yellow
                         "ERROR" = "\033[31m", # red
                         "INFO" = "\033[36m") # cyan
  
  cat(status_color, sprintf("[%s] ", status), "\033[0m", message, "\n", sep = "")
}

# Start validation
cat("\n============================================\n")
cat("VALIDATING DDM_BASICS_R PROJECT\n")
cat("============================================\n\n")

# =============================================
# 1. Check Required Packages
# =============================================
print_status("Checking required packages...", "INFO")

required_packages <- c("shiny", "ggplot2", "dplyr", "patchwork", 
                      "knitr", "rmarkdown", "devtools")
missing_packages <- c()

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    print_status(sprintf("Package '%s' is installed", pkg), "OK")
  } else {
    print_status(sprintf("Package '%s' is NOT installed", pkg), "ERROR")
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("\nMissing packages. Install them with:\n")
  cat(sprintf('install.packages(c("%s"))\n', paste(missing_packages, collapse = '", "')))
}

# =============================================
# 2. Check R Source Files
# =============================================
print_status("\nChecking R source files...", "INFO")

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
print_status(sprintf("Found %d R files", length(r_files)), "INFO")

# Try to source each R file
for (file in r_files) {
  tryCatch({
    source(file)
    print_status(sprintf("Sourced '%s'", file), "OK")
  }, error = function(e) {
    print_status(sprintf("Error in '%s': %s", file, e$message), "ERROR")
  }, warning = function(w) {
    print_status(sprintf("Warning in '%s': %s", file, w$message), "WARNING")
  })
}

# =============================================
# 3. Basic Functionality Test
# =============================================
print_status("\nTesting basic functionality...", "INFO")

# Test functions if they exist
test_functions <- function() {
  # Test random walk simulator
  if (exists("simulate_random_walk_trial")) {
    tryCatch({
      result <- simulate_random_walk_trial(
        drift_rate = 0.1, 
        threshold = 1, 
        start_point = 0.5, 
        noise_sd = 0.1, 
        dt = 0.001, 
        max_decision_time = 1.0
      )
      print_status("Random walk simulation works", "OK")
    }, error = function(e) {
      print_status(sprintf("Random walk simulation error: %s", e$message), "ERROR")
    })
  } else {
    print_status("Function 'simulate_random_walk_trial' not found", "WARNING")
  }
  
  # Test basic diffusion model
  if (exists("simulate_diffusion_trial")) {
    tryCatch({
      result <- simulate_diffusion_trial(
        v = 0.15, 
        a = 1.0, 
        z = 0.5, 
        s = 0.1, 
        ter = 0.2, 
        dt = 0.001, 
        max_decision_time = 1.0
      )
      print_status("Basic diffusion model works", "OK")
    }, error = function(e) {
      print_status(sprintf("Basic diffusion model error: %s", e$message), "ERROR")
    })
  } else {
    print_status("Function 'simulate_diffusion_trial' not found", "WARNING")
  }
  
  # Test variable diffusion model
  if (exists("simulate_diffusion_trial_variable")) {
    tryCatch({
      result <- simulate_diffusion_trial_variable(
        mean_v = 0.15, 
        a = 1.0, 
        mean_z = 0.5, 
        s = 0.1, 
        mean_ter = 0.2, 
        sv = 0.1, 
        sz = 0.1, 
        st0 = 0.05, 
        dt = 0.001, 
        max_decision_time = 1.0
      )
      print_status("Variable diffusion model works", "OK")
    }, error = function(e) {
      print_status(sprintf("Variable diffusion model error: %s", e$message), "ERROR")
    })
  } else {
    print_status("Function 'simulate_diffusion_trial_variable' not found", "WARNING")
  }
}

# Run function tests in a safer way
tryCatch({
  test_functions()
}, error = function(e) {
  print_status(sprintf("Error in function tests: %s", e$message), "ERROR")
})

# =============================================
# 4. Check Vignettes
# =============================================
print_status("\nChecking vignettes...", "INFO")

vignette_files <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
print_status(sprintf("Found %d vignette files", length(vignette_files)), "INFO")

# Check if we can parse each vignette
if (requireNamespace("knitr", quietly = TRUE)) {
  for (file in vignette_files) {
    tryCatch({
      # Just parse the vignette structure, don't render fully
      knitr::knit(file, output = tempfile(), tangle = TRUE, quiet = TRUE)
      print_status(sprintf("Vignette '%s' parsed successfully", basename(file)), "OK")
    }, error = function(e) {
      print_status(sprintf("Error in vignette '%s': %s", basename(file), e$message), "ERROR")
    }, warning = function(w) {
      print_status(sprintf("Warning in vignette '%s': %s", basename(file), w$message), "WARNING")
    })
  }
} else {
  print_status("Skipping vignette parsing as 'knitr' is not installed", "WARNING")
}

# =============================================
# 5. Check pkgdown configuration (if available)
# =============================================
print_status("\nChecking pkgdown configuration...", "INFO")

if (file.exists("_pkgdown.yml")) {
  tryCatch({
    yaml_content <- readLines("_pkgdown.yml")
    print_status("_pkgdown.yml exists and can be read", "OK")
    
    # Check if yaml is valid if yaml package is available
    if (requireNamespace("yaml", quietly = TRUE)) {
      yaml::yaml.load_file("_pkgdown.yml")
      print_status("_pkgdown.yml contains valid YAML", "OK")
    }
  }, error = function(e) {
    print_status(sprintf("Error in _pkgdown.yml: %s", e$message), "ERROR")
  })
} else {
  print_status("_pkgdown.yml file not found", "WARNING")
}

# =============================================
# 6. Verify project structure
# =============================================
print_status("\nVerifying project structure...", "INFO")

expected_dirs <- c("R", "vignettes", "docs", "data")
for (dir in expected_dirs) {
  if (dir.exists(dir)) {
    print_status(sprintf("Directory '%s/' exists", dir), "OK")
  } else {
    print_status(sprintf("Directory '%s/' is missing", dir), "WARNING")
  }
}

# =============================================
# Summary
# =============================================
cat("\n============================================\n")
cat("VALIDATION COMPLETE\n")
cat("============================================\n")
cat("\nRun this script in R or RStudio with:\n")
cat("  source('check_project.R')\n\n")
cat("If there are any errors, fix them before proceeding.\n")
cat("If there are warnings, review them and address if necessary.\n\n") 