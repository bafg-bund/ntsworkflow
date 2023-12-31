% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/componentization.R
\name{alig_componentisation}
\alias{alig_componentisation}
\title{Second stage componentisation using the alignment table}
\usage{
alig_componentisation(
  altable,
  rttols = 3,
  fracComponMatch = 0.5,
  mztol = 0.005,
  pearsonCorr = 0.5,
  pol = "pos",
  numcores = 6
)
}
\arguments{
\item{altable}{alignment table matrix ("grouped" in the app)}

\item{rttols}{rt tolerance in s}

\item{fracComponMatch}{fraction of samples with matching peak shape}

\item{mztol}{mz tolerance}

\item{pearsonCorr}{minimum pearson's r for comparing intensity trends}

\item{pol}{polarity, must be either "pos" or "neg"}

\item{numcores}{number of cores to use for parallel distance matrix computations}
}
\value{
Alignment table with the group column "Gruppe" replaced with the new values
}
\description{
Based on four criteria: rt, peak shape, intensity correlation and common mz differences.
Will replace group column in the alignment table. Distance matrices for all 
four criteria are combined and then DBSCAN clusters the features.
}
