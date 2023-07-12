#' Start the shiny peak picking app from the R console
#' 
#' The current working directory will become the working directory for the app.
#' The app will display as a browser pop-up.
#' 
#' @export
runPeakPicking <- function() {
  appDir <- system.file("shiny", "non_target_app", package = "ntsworkflow")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `ntsworkflow`.", call. = FALSE)
  }
  globalwd <<- getwd()
  shiny::runApp(appDir, display.mode = "normal")
}