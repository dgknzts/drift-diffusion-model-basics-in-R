# .Rprofile for DDM_Basics_R project

# Load commonly used packages
if (interactive()) {
  suppressMessages({
    library(devtools)
    library(pkgdown)
  })
  
  # Custom helper functions
  build_site <- function() {
    source("scripts/update_site.R")
  }
  
  launch_app <- function() {
    source("R/run_app.R")
    run_ddm_app()
  }
  
  cat("DDM_Basics_R project .Rprofile loaded\n")
  cat("Helper functions available:\n")
  cat("  build_site() - Build the pkgdown site\n")
  cat("  launch_app() - Launch the Shiny app\n")
} 