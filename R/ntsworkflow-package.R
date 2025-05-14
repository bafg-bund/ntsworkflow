
#' @keywords internal 
"_PACKAGE"

#' ntsworkflow: A package for non-target-screening data processing.
#'
#'
#' @name ntsworkflow
#' @useDynLib ntsworkflow, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @import foreach
#' @import dplyr
NULL


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nPlease cite ntsworkflow if you use it: see citation('ntsworkflow') for details.\n")
}


# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow

