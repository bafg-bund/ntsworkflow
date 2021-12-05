
#' @export
runPeakPicking <- function() {
  appDir <- system.file("shiny", "non_target_app", package = "ntsworkflow")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `ntsworkflow`.", call. = FALSE)
  }
  globalwd <<- getwd()
  shiny::runApp(appDir, display.mode = "normal")
}