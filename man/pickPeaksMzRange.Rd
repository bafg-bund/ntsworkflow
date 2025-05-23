% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peak-picking.R
\name{pickPeaksMzRange}
\alias{pickPeaksMzRange}
\title{Peak-picking for an m/z range}
\usage{
pickPeaksMzRange(
  daten,
  mz_min,
  mz_max,
  mz_step,
  rt_min,
  rt_max,
  sn,
  int_threshold,
  peak_NoiseScans,
  precursormzTol,
  peakwidth_min,
  peakwidth_max,
  maxPeaksPerSignal
)
}
\arguments{
\item{daten}{Measurement data of class `xcms::xcmsRAW`}

\item{mz_min}{m/z range to peak-picking (lower) (Da)}

\item{mz_max}{m/z range to peak-picking (upper) (Da)}

\item{mz_step}{m/z width for extracted ion chromatogram (Da)}

\item{rt_min}{Retention time range minimum in which to look for peaks (s)}

\item{rt_max}{Retention time range maximum in which to look for peaks (s)}

\item{sn}{Minimum signal-to-noise ratio (apex peak height over noise spread before and after peak)}

\item{int_threshold}{Minimum peak intensity (at peak apex)}

\item{peak_NoiseScans}{Number of scans before and after peak to measure noise}

\item{precursormzTol}{m/z tolerance for linking MS2 fragment spectra by the precursor m/z (ppm)}

\item{peakwidth_min}{Minimum peak width given (s)}

\item{peakwidth_max}{Maximum peak width given (s)}

\item{maxPeaksPerSignal}{Maximum number of sub-peaks within a peak (direction changes) within a peak}
}
\value{
`data.frame` of the peak inventory (peak table)
}
\description{
Function will pick chromatographic peaks in the raw data file by binning and 
analyzing the extracted ion chromatograms.
}
