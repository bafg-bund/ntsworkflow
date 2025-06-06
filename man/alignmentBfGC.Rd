% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{alignmentBfGC}
\alias{alignmentBfGC}
\title{Align peaks in separate peak lists by m/z and retention time tolerance}
\usage{
alignmentBfGC(peaklistR, mz_dev, DeltaRT, mz_dev_unit)
}
\arguments{
\item{peaklistR}{list of peak tables}

\item{mz_dev}{m/z tolerance}

\item{DeltaRT}{RT tolerance (s)}

\item{mz_dev_unit}{Units for the m/z tolerance 1 = ppm, 2 = mDa}
}
\value{
A matrix of class Rcpp::integerMatrix giving the .
}
\description{
This function will align the peaks in multiple peak lists by
comparing the m/z and RT.
}
