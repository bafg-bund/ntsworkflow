% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/suspect_screening.R
\name{ms2_search}
\alias{ms2_search}
\title{Suspect search without peak searching, looking only at MS2 spectra}
\usage{
ms2_search(
  data_path,
  db_path,
  rttolm = 1,
  mztolu = 0.5,
  mztolu_fine = 0.005,
  chromatography = "dx.doi.org/10.1016/j.chroma.2015.11.014",
  pol = "pos",
  CE_s = 30:40,
  CES_s = 0:15,
  instr = "LC-ESI-QTOF TripleTOF 5600 SCIEX",
  ceunit = c("V", "eV"),
  comparison = "dot_product",
  threshold = 400,
  rt_res = 1.5,
  rtoffset = 0,
  ndp_m = 2,
  ndp_n = 1,
  mztolu_ms2 = 0.015,
  compounds = NULL
)
}
\arguments{
\item{data_path}{can be a vector of file locations or a list of xcmsRaw objects}

\item{rt_res}{resolution of two peaks in chromatography}

\item{ndp_m}{Peak intensity weighting factor for dot-product}

\item{ndp_n}{m/z weighting factor for dot-product}

\item{mztolu_ms2}{m/z window for dot-product}

\item{compounds}{character vector of compounds which should be processed}
}
\description{
Suspect search without peak searching, looking only at MS2 spectra
}