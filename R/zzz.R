# Package admin functions


#' @useDynLib ntsworkflow
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nPlease cite ntsworkflow if you use it: see citation('ntsworkflow') for details.\n")
}

