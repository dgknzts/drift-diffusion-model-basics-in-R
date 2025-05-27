# R/utils/check_dependencies.R

#' Check and Install Required Packages for DDM Analysis
#'
#' Checks if required packages are installed and optionally installs missing ones.
#' Useful for ensuring the vignette examples will run properly.
#'
#' @param packages Character vector. Package names to check. Default includes
#'   all packages used in the DDM vignette.
#' @param install_missing Logical. Whether to install missing packages.
#'   Default is FALSE (just report missing packages).
#' @param quiet Logical. Whether to suppress messages. Default is FALSE.
#'
#' @return Logical vector indicating which packages are available.
#' @export
#'
#' @examples
#' # Check if required packages are available
#' check_ddm_dependencies()
#' 
#' # Install missing packages (if desired)
#' check_ddm_dependencies(install_missing = TRUE)
check_ddm_dependencies <- function(packages = c("ggplot2", "dplyr", "tidyr", 
                                                "gridExtra", "patchwork", "knitr"),
                                   install_missing = FALSE,
                                   quiet = FALSE) {
  
  available <- sapply(packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })
  
  missing <- packages[!available]
  
  if (!quiet) {
    if (length(missing) == 0) {
      cat("✓ All required packages are available!\n")
    } else {
      cat("✗ Missing packages:", paste(missing, collapse = ", "), "\n")
      
      if (install_missing) {
        cat("Installing missing packages...\n")
        install.packages(missing)
        
        # Re-check after installation
        available_after <- sapply(packages, function(pkg) {
          requireNamespace(pkg, quietly = TRUE)
        })
        
        still_missing <- packages[!available_after]
        if (length(still_missing) == 0) {
          cat("✓ All packages successfully installed!\n")
        } else {
          cat("✗ Failed to install:", paste(still_missing, collapse = ", "), "\n")
        }
        
        return(available_after)
      } else {
        cat("Run with install_missing = TRUE to install automatically.\n")
        cat("Or install manually with: install.packages(c(", 
            paste0('"', missing, '"', collapse = ", "), "))\n")
      }
    }
  }
  
  return(available)
}

#' Load DDM Analysis Environment
#'
#' Loads all required packages and sources DDM functions for analysis.
#' Convenience function for setting up the complete DDM analysis environment.
#'
#' @param check_packages Logical. Whether to check package availability first.
#'   Default is TRUE.
#' @param quiet Logical. Whether to suppress loading messages. Default is FALSE.
#'
#' @return NULL (invisible). Function is called for side effects.
#' @export
#'
#' @examples
#' # Set up complete DDM analysis environment
#' setup_ddm_environment()
setup_ddm_environment <- function(check_packages = TRUE, quiet = FALSE) {
  
  if (check_packages) {
    deps_ok <- check_ddm_dependencies(quiet = quiet)
    if (!all(deps_ok)) {
      stop("Some required packages are missing. Run check_ddm_dependencies(install_missing = TRUE) first.")
    }
  }
  
  # Load packages
  packages_to_load <- c("ggplot2", "dplyr", "tidyr", "gridExtra", "patchwork")
  
  for (pkg in packages_to_load) {
    if (!quiet) cat("Loading", pkg, "...\n")
    library(pkg, character.only = TRUE)
  }
  
  # Source DDM functions (assuming they're in the package)
  # This would be automatic if this were a proper R package
  
  if (!quiet) {
    cat("✓ DDM analysis environment ready!\n")
    cat("Available functions:\n")
    cat("  - simulate_diffusion_trial()\n")
    cat("  - simulate_diffusion_experiment()\n") 
    cat("  - simulate_diffusion_trial_with_path()\n")
    cat("  - create_ddm_params()\n")
    cat("  - create_parameter_grid()\n")
    cat("  - summarize_ddm_data()\n")
    cat("  - plot_ddm_paths_with_histograms()\n")
    cat("  - plot_rt_comparison()\n")
    cat("  - plot_parameter_effects()\n")
    cat("  - plot_rt_qq()\n")
    cat("  - create_summary_table()\n")
  }
  
  invisible(NULL)
} 