# Copyright 2016-2023 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow
# ntsworkflow is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any 
# later version.
# 
# ntsworkflow is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with ntsworkflow. If not, see <https://www.gnu.org/licenses/>.


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