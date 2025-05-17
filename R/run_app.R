#' Launch the DDM Interactive Demonstration App
#'
#' This function launches a Shiny app that allows users to
#' interactively explore Diffusion Decision Model parameters
#' and their effects on decision times and choice probabilities.
#'
#' @return The app object
#' @export
#'
#' @examples
#' \dontrun{
#' run_ddm_app()
#' }
run_ddm_app <- function() {
  app_dir <- system.file("shiny-app", package = "DDM_Basics_R")
  if (app_dir == "") {
    # For development when not installed as a package
    app_dir <- "inst/shiny-app"
    if (!dir.exists(app_dir)) {
      stop("App directory not found. Try installing the package first.")
    }
  }
  
  shiny::runApp(app_dir, display.mode = "normal")
} 