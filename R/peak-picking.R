


#' Peak-picking for an m/z range
#'
#' @description Function will pick chromatographic peaks in the raw data file by binning and 
#' analyzing the extracted ion chromatograms. 
#' 
#' @param daten Measurement data of class `xcms::xcmsRAW`
#' @param mz_min m/z range to peak-picking (lower) (Da)
#' @param mz_max m/z range to peak-picking (upper) (Da)
#' @param mz_step m/z width for extracted ion chromatogram (Da)
#' @param rt_min Retention time range minimum in which to look for peaks (s)
#' @param rt_max Retention time range maximum in which to look for peaks (s)
#' @param sn Minimum signal-to-noise ratio (apex peak height over noise spread before and after peak)
#' @param int_threshold Minimum peak intensity (at peak apex)
#' @param peak_NoiseScans Number of scans before and after peak to measure noise
#' @param precursormzTol m/z tolerance for linking MS2 fragment spectra by the precursor m/z (ppm)
#' @param peakwidth_min Minimum peak width given (s)
#' @param peakwidth_max Maximum peak width given (s)
#' @param maxPeaksPerSignal Maximum number of sub-peaks within a peak (direction changes) within a peak
#'
#' @returns `data.frame` of the peak inventory (peak table)
#' @export
pickPeaksMzRange <- function(
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
) {
  
  
  cn <- c(
    "mz",
    "RT",
    "Intensity",
    "XIC_Intensity",
    "Scan",
    "LeftendRT",
    "RightendRT",
    "Leftendscan",
    "Rightendscan",
    "NoiseDeviation",
    "Area",
    "FWHM_left",
    "FWHM_right",
    "Baseline",
    "XIC-step",
    "MS2scan",
    "SN",
    "InSourceFragmentOf"
  )
  
  
  RT_Tol_scan <- 3
  mz_Tol_ppm <- 20
  
  
  sctLarger <- daten@scantime > rt_min
  sctSmaller <- daten@scantime < rt_max
  if (all(!sctLarger) || all(!sctLarger))
    return(as.data.frame(matrix(NA, 0, 18, dimnames = list(NULL, cn))))
  rt_min_scan <- min(which(sctLarger))
  rt_max_scan <- max(which(sctSmaller))
  
  
  massrange <- c(mz_min, mz_min+round((mz_max-mz_min)/10))
  pl <- NULL
  
  # Run peak picking algorithm
  plrows <- lapply(
    seq(mz_min, mz_max, by = mz_step*0.5), 
    pickPeaksOneEic,
    rawData = daten, mz_step = mz_step,
    rt_min_scan = rt_min_scan,
    rt_max_scan = rt_max_scan, 
    sn = sn,
    int_threshold = int_threshold,
    NoiseScans = peak_NoiseScans,
    precursormzTol = precursormzTol,
    peakwidth_min = peakwidth_min,
    peakwidth_max = peakwidth_max,
    maxPeaksPerSignal = maxPeaksPerSignal
  )
  
  if (all(vapply(plrows, is.null, logical(1))))
    return(as.data.frame(matrix(NA, 0, 18, dimnames = list(NULL, cn))))
  
  pl <- do.call("rbind", plrows)
  
  # if nothing found, return empty peaklist
  if (is.null(pl)) 
    return(as.data.frame(matrix(NA, 0, 18, dimnames = list(NULL, cn))))
  
  if (nrow(pl) >= 1) {
    
    pl <- pl[order(pl[,1]), , drop = F]
    
    # Delete shoulder peaks which are inside the RT range of another peak:
    for (i in 1:nrow(pl)) {
      doppelte <- which((abs(pl[,1]-pl[i,1]) < pl[i,1]/1000000*mz_Tol_ppm) & 
                          (pl[,2] >= pl[i,5]) & (pl[,2] <= pl[i,6])) 
      doppelte <- doppelte[doppelte != i]
      pl[doppelte, 1] <- 0  
    }
    
    pl <- pl[pl[,1] > 0, , drop = F]
    
  }
  
  # Calculation of S/N
  if (nrow(pl) > 0)
    pl <- cbind(pl,pl[,3]/pl[,10],0)
  
  colnames(pl) <- cn
  
  # Remove rows which are below intensity threshold and not within RT window
  pl <- pl[pl[, "Intensity"] >= int_threshold, , drop = F]
  pl <- pl[pl[, "RT"] >= rt_min, , drop = F]
  pl <- pl[pl[, "RT"] <= rt_max, , drop = F]
  
  as.data.frame(pl)
}


#' Peakpeaking algorithm using C++
#' 
#' This function wraps the c++ peakpicking algorithm. The function picks peaks
#' in a specific mz range (defined by i and mz_step). Used by the function
#' FindPeaks_BfG.
#'
#' @param i Lower m/z for extracted ion chromatogram in Da
#' @param rawData Measurement data of class `xcms::xcmsRAW`
#' @param mz_step m/z width for extracted ion chromatogram in Da
#' @param rt_min_scan Retention time range minimum in which to look for peaks, in scans
#' @param rt_max_scan Retention time range maximum in which to look for peaks, in scans
#' @param sn Minimum signal-to-noise ratio (apex peak height over noise spread before and after peak)
#' @param int_threshold Minimum peak intensity (at peak apex)
#' @param NoiseScans Number of scans before and after peak to measure noise
#' @param peakwidth_min Minimum peak width given in seconds
#' @param peakwidth_max Maximum peak width given in seconds
#' @param precursormzTol m/z tolerance for linking MS2 fragment spectra by the precursor m/z, in ppm
#' @param maxPeaksPerSignal Maximum number of sub-peaks within a peak (direction changes) within a peak
#'
#' @returns 'matrix' of peaks detected in the extracted ion chromatogram, each peak
#' is one row with the columns representing various parameters for the peak.
pickPeaksOneEic <- function(
    i, 
    rawData,
    mz_step, 
    rt_min_scan,
    rt_max_scan, 
    sn,
    int_threshold,
    NoiseScans,
    peakwidth_min,
    peakwidth_max,
    precursormzTol,
    maxPeaksPerSignal
  ) {
  
  maxima <- NULL
  
  XIC <- xcms::rawEIC(rawData, mzrange = c(i, i+mz_step))
  XIC <- XIC$intensity
  stopifnot(
    is.double(XIC), 
    length(XIC) > 0
  )
  
  maxima <- pickPeaksOneEicCpp(
    mz = i, mz_step = mz_step, XIC, scantime = rawData@scantime, 
    min_intensity = int_threshold, sn = sn, noisescans = NoiseScans, 
    peakwidth_min = peakwidth_min, peakwidth_max = peakwidth_max, 
    maxPeaksPerSignal = maxPeaksPerSignal
  )


  if (nrow(maxima) > 0) {
    for (j in 1:nrow(maxima)) {
      mass_spectrum <- xcms::getScan(rawData, maxima[j,5], mzrange = c(i,i+mz_step))
      exactmass <- 0
      if (nrow(mass_spectrum) > 0) {
        exactmass <- mass_spectrum[which.max(mass_spectrum[,2]),1]
        maxima[j,3] <- max(mass_spectrum[,2])
        if (nrow(mass_spectrum) == 1) 
          maxima[j,3] <- maxima[j,3]-maxima[j,14]
      }
      if ((exactmass == 0) | (exactmass < i+mz_step/4) | (exactmass > i+mz_step/4*3)) 
        exactmass <- 0 
      maxima[j,1] <- exactmass
      ms2 <- which(
        (rawData@msnRt > maxima[j,6]) & 
        (rawData@msnRt < maxima[j,7]) & 
        (abs(rawData@msnPrecursorMz-exactmass) <= exactmass/1000000*precursormzTol)
      )
      if (length(ms2) > 0) 
        maxima[j,16] <- ms2[which.min(abs(rawData@msnRt[ms2]-maxima[j,2]))]
    }
  }
  maxima <- maxima[maxima[,1] > 0,,drop = FALSE]
  if (nrow(maxima) > 0) 
    return(maxima)
}

#' Find Pearson's correlation between intensity trends
#' 
#' @description Wrapper for the cpp function correlates_with. Takes a matrix of intensities
#' (each row representing one feature) and a row, returns which rows have
#' correlation coefficients within threshold
#' 
#' @param aligned_intensities Matrix of intensities
#' @param zeile Reference row of matrix to which all other are compared to
#' @param koeffizient Pearson correlation coefficient (minimum)
#' 
#' @returns Boolean. True for those rows of the matrix which have a correlation 
#' coefficient higher than the threshold. 
#' 
#' @export
correlates_with_r <- function(aligned_intensities, zeile, koeffizient){
  return(correlates_with(aligned_intensities, zeile, koeffizient))
}




#' Summarize componentization groups into one column 
#' 
#' @description Will take the componentization information from each sample in the
#' alignment table and attempt to summerize this into one column named "Gruppe".
#' 
#' @param alig Alignment table
#' 
#' @returns Alignment table with the additional column Gruppe
#' @export
summarize_groups <- function(alig) {
  gruppenzaehler <- 1
  
  newGruppe <- numeric(nrow(alig))
  
  intCols <- alig[, grep("^Int_", colnames(alig))]
  grupCols <- alig[, grep("^gruppe_", colnames(alig))]
  
  maxIntCol <- max.col(intCols)
  maxSamp <- stringr::str_match(colnames(intCols)[maxIntCol], "Int_(\\d+)$")[,2]
  # the mapply function fails if more only one sample present
  maxGrup <- mapply(function(x, y) grupCols[x, y], seq_len(nrow(alig)), maxIntCol)
  alig <- as.data.frame(alig)  # so that you can use the much faster [[ function
  for (i in seq_along(newGruppe)) {  #i<-2
    if (newGruppe[i] != 0)
      next
    maxGrupCol <- paste0("gruppe_", maxSamp[i])
    # b = rows with Gruppe == 0 AND same groups found in first samples of a
    b <- which(newGruppe == 0 & alig[[maxGrupCol]] == maxGrup[i])  
    newGruppe[b] <- gruppenzaehler  
    gruppenzaehler <- gruppenzaehler + 1
  }
  
  # Eliminieren von ein-Peak-Gruppen  (Gruppe wird 0)
  for (i in seq_len(max(newGruppe))) {
    thisGruppe <- which(newGruppe == i)
    if (length(thisGruppe) == 1)
      newGruppe[thisGruppe] <- 0
  }
  
  # Schließen der entstandenen Zwischenräume
  V1 <- sort(unique(newGruppe)[-(which(unique(newGruppe) == 0))])
  V2 <- seq_along(V1) 
  for (i in 1:length(V1)) {
    newGruppe[which(newGruppe == V1[i])] <- V2[i]
  }
  alig$Gruppe <- newGruppe
  alig <- as.matrix(alig)
  alig
}


#' Get a recommendation for the m/z-step parameter (binning width)
#' 
#' @description Calculates the optimum mz step based on mass differences. Takes a sample
#' of 100 spectra to make the calculation
#' 
#' @param daten raw data as list of xcmsRAW objects
#' @param probs quantile to use (see quantile function)
#' 
#' @returns Recommended m/z width (Da)
#' 
#' @export
optimumMzStep <- function(daten, probs) {
  get_diff <- function(scan) {
    msspektrum <- xcms::getScan(daten, scan)
    -quantile(-diff(msspektrum[, 1], lag = 1), probs = probs)
  }
  ns <- ifelse(length(daten@scanindex) < 100, length(daten@scanindex), 100)
  diffs <- vapply(sample(seq_along(daten@scanindex), ns), get_diff, numeric(1))
  round(min(diffs), 2)
}

# Copyright 2016-2024 Bundesanstalt für Gewässerkunde
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
