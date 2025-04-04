
if (suppressMessages(suppressWarnings(require(xcms))))
  message("test passed")

detach("package:xcms", unload = T)

if (suppressMessages(suppressWarnings(require(ntsworkflow))))
  message("test passed")

detach("package:ntsworkflow", unload = T)