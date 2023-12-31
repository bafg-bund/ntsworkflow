% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenFormR_modified.R
\name{fixNoLossSpace}
\alias{fixNoLossSpace}
\title{Fix Whitespace Error in GenForm Loss Annotation}
\usage{
fixNoLossSpace(anno_MS_file, replace=TRUE)
}
\arguments{
\item{anno_MS_file}{File name for annotated cleaned MSMS text file output from \code{\link{annotateCleanMSMS}}}

\item{replace}{If \code{TRUE}, replace the input file. If \code{FALSE}, saved under the new name
\code{"*_ed.txt"}.}
}
\value{
File name of the "fixed" file.
}
\description{
This takes a filename with e.g. the output from \code{\link{annotateCleanMSMS}} and
removes a space in the \code{"no loss"} annotation to allow for easier handling.
}
\seealso{
\code{\link{annotateCleanMSMS}}
}
\author{
Emma Schymanski (\code{\link{emma.schymanski@uni.lu}})
}
