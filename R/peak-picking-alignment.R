


#' Peakpeaking algorithm using c++
#' 
#' This function wraps the c++ peakpicking algorithm. The function picks peaks
#' in a specific mz range (defined by i and mz_step). Used by the function
#' FindPeaks_BfG.
#'
#' @param i lower mz for extracted ion chromatogram
#' @param rawData 
#' @param mz_step mz width for extracted ion chromatogram
#' @param rt_min_scan 
#' @param rt_max_scan 
#' @param sn 
#' @param int_threshold 
#' @param NoiseScans 
#' @param peakwidth_min 
#' @param peakwidth_max 
#' @param precursormzTol 
#' @param maxPeaksPerSignal 
#'
#' @return matrix of peaks detected in the extracted ion chromatogram, each peak
#' is one row with the columns representing various parameters for the peak.
#' @export
peakpicking_BfG_cpp <- function(
    i, rawData, mz_step, rt_min_scan, rt_max_scan, 
    sn, int_threshold, NoiseScans, peakwidth_min, peakwidth_max,
    precursormzTol, maxPeaksPerSignal) {
  
  maxima <- NULL
  
  XIC <- xcms::rawEIC(rawData, mzrange = c(i,i+mz_step))
  XIC <- XIC$intensity
  
  maxima <- peakPickingBfGC(
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
#' Wrapper for the cpp function correlates_with. Takes a matrix of intensities
#' (each row representing one feature) and a row, returns which rows have
#' correlation coefficients within threshold
#' 
#' @param aligned_intensities Matrix of intensities
#' @param zeile Reference row of matrix
#' @param koeffizient Pearson correlation coefficient 
#' 
#' @return Boolean. True for those rows of the matrix which have a correlation 
#' coefficient higher than the threshold. 
#' 
#' @export
correlates_with_r <- function(aligned_intensities, zeile, koeffizient){
  return(correlates_with(aligned_intensities, zeile, koeffizient))
}


#' Alignment of peaks from difference samples
#' 
#' This a wrapper for the cpp function alignmentBfGC. It takes a list of peaklists
#' and produces an alignment table. 
#' 
#' @param peaklists list of peaklists
#' @param mz_dev mz tolerance
#' @param DeltaRT rt tolerance in s
#' @param mz_dev_unit units for mz tolerance (ppm or mDa)
#' 
#' @return matrix of aligned peaks, one aligned feature per row
#' 
#' @export
alignment_BfG_cpp <- function(peaklists, mz_dev, DeltaRT, mz_dev_unit){
  peaklistC <- list()
  for (i in 1:length(peaklists)) {
    peaklistC[[i]] <- as.matrix(peaklists[[i]][,1:3])
  }
  # set mz_dev_unit to integer 1 or 2, 1 for ppm, 2 for mDa
  stopifnot(mz_dev_unit %in% c("ppm", "mDa"))
  mz_dev_unit <- switch(mz_dev_unit, ppm = 1, mDa = 2)
  
  alignedtable <- alignmentBfGC(peaklistC, mz_dev, DeltaRT, mz_dev_unit)
  grouped2 <- matrix(0,nrow=nrow(alignedtable-1),ncol=length(peaklists)*6+4)
  
  Gruppe <- numeric(nrow(grouped2))
  
  spaltennamen <- c("mean_mz","mean_RT","MS2Fit", "Gruppe")
  
  for (i in 1:ncol(alignedtable)) {
    spaltennamen <- c(
      spaltennamen,
      paste0("PeakID_",as.character(i)),
      paste0("mz_",as.character(i)),
      paste0("RT_",as.character(i)),
      paste0("Int_",as.character(i)),
      paste0("ms2scan_",as.character(i)),
      paste0("gruppe_",as.character(i))
    )
    
    for (j in 1:(nrow(alignedtable))) {
      if (alignedtable[j,i] > 0) {
        grouped2[j,i*6-1] <- peaklists[[i]]$peak_id_all[alignedtable[j,i]]
        grouped2[j,i*6] <- peaklists[[i]]$mz[alignedtable[j,i]]
        grouped2[j,i*6+1] <- peaklists[[i]]$RT[alignedtable[j,i]]
        grouped2[j,i*6+2] <- peaklists[[i]]$Intensity[alignedtable[j,i]]
        grouped2[j,i*6+3] <- peaklists[[i]]$MS2scan[alignedtable[j,i]]
        grouped2[j,i*6+4] <- peaklists[[i]]$Gruppe[alignedtable[j,i]]
      }
      
      if (i == (ncol(alignedtable))) {
        grouped2[j,1] <- mean(grouped2[j,which(alignedtable[j,] > 0)*6])
        grouped2[j,2] <- mean(grouped2[j,which(alignedtable[j,] > 0)*6+1])
      }
    }
  }
  
  colnames(grouped2) <- spaltennamen

  grouped2
}

#' Summarize componentization groups into one column 
#' 
#' Will take the componentization information from each sample in the
#' alignment table and attempt to summerize this into one column named "Gruppe".
#' 
#' @param alig Alignment table
#' 
#' @return Matrix. Alignment table with additional column Gruppe
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


#' Align peaks with parallel processing
#' 
#' To align very large datasets it is necessary to split the alignment process
#' into parallel jobs. This function slips the peak lists and runs the
#' alignment_BfG_cpp on each chunk in parallel. The results are combined. Uses
#' forking therefore not available on Windows.
#'
#' @param peaklists list of peaklists
#' @param ppm_dev mz tolerance
#' @param DeltaRT rt tolerance in s
#' @param mz_dev_unit units for mz tolerance (ppm or mDa)
#' @param mDa_split mz gap between chunks
#' @param numSplits Number parallel chunks to create
#'
#' @return matrix of aligned peaks, one aligned feature per row
#' @export
#'
alignment_BfG_cpp_par <- function(peaklists, ppm_dev, DeltaRT, mz_dev_unit, 
                                  mDa_split = 100, numSplits = 16) {
  # split the peaks
  plTemp <- do.call("rbind", peaklists)
  plTemp <- plTemp[order(plTemp$mz), ]
  splitMz <- plTemp$mz[(plTemp$mz - dplyr::lag(plTemp$mz))*1000 > mDa_split][-1]
  # this will create too many splits, assuming we have 16 cores, we want 16 splits in total each
  # a similar number of features to compare
  mzGroups <- findInterval(plTemp$mz, splitMz, all.inside = T)
  numEachGroup <- cumsum(table(mzGroups))
  # around about how many should each group have?
  goal <- round(nrow(plTemp) / numSplits)
  # alright these are the numbers you need to find in the cumulative sums
  findme <- goal * 1:numSplits  
  ind <- vapply(findme, function(x) which.min(abs(numEachGroup - x)), FUN.VALUE = numeric(1))
  ind <- c(1, ind)
  if (ind[numSplits+1] != max(mzGroups))  # the last group needs to be in this index
    ind[numSplits+1] <- max(mzGroups)
  # collect new groups together

  for (i in 1:numSplits) { 
    if (i == 1) {
      mzGroups[mzGroups %in% ind[i]:ind[i+1]] <- i
    } else {
      mzGroups[mzGroups %in% (ind[i]+1):ind[i+1]] <- i
    }
  }
  plnew <- split(plTemp, mzGroups)
  plnew2 <- lapply(plnew, function(x) split(x, x$sample_id))
  
  # there are some missing peak list tables in some of the mass groups
  # You need to find these and fill them with a dummy pl table with mz set to 0
  dummyPltable <- plnew2[[1]][[1]][1,]
  dummyPltable[1,"mz"] <- 0 
  
  for (mzGroup in seq_along(plnew2)) {
    for (samp in plTemp$sample_id) {
      if (is.null(plnew2[[mzGroup]][[as.character(samp)]])) {
        plnew2[[mzGroup]][[as.character(samp)]] <- dummyPltable
      }
    }
  }
  
  # sort the peaklists into the right order
  for (mzGroup in seq_along(plnew2)) {
    plnew2[[mzGroup]] <- plnew2[[mzGroup]][as.character(sort(as.numeric(names(plnew2[[2]]))))]
  }
  # browser()

  # okay, now we are ready to align...in parallel! WAha ha ha haaah!!
  result <- parallel::mclapply(plnew2, function(x) {
    ntsworkflow::alignment_BfG_cpp(x, ppm_dev, DeltaRT, mz_dev_unit)
  })
  result <- do.call("rbind", result)

  # remove the dummy rows
  result <- result[result[, "mean_mz"] != 0, ]
  
  # just need to reorder and then done!
  result[order(result[, 1]), ]
}





#' Peak finding algorithm
#'
#' Function will pick chromatographic peaks in the raw data file by binning and 
#' analyzing the extracted ion chromatograms. The peak finding is done by
#' taking the derivative of the chromatogram and looking for points where it 
#' crosses 0. More information can be found in Dietrich, C., Wick, A., & 
#' Ternes, T. A. (2021). Open source feature detection for non‐target LC‐MS 
#' analytics. Rapid Communications in Mass Spectrometry, e9206. doi:10.1002/rcm.9206  
#' 
#' 
#' @param daten from xcms
#' @param mz_min 
#' @param mz_max 
#' @param mz_step 
#' @param rt_min 
#' @param rt_max 
#' @param sn 
#' @param int_threshold 
#' @param peak_NoiseScans 
#' @param precursormzTol 
#' @param peakwidth_min 
#' @param peakwidth_max 
#' @param maxPeaksPerSignal 
#'
#' @return data.frame of the peak inventory list
#' @export
FindPeaks_BfG <- function(daten,  
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
                          maxPeaksPerSignal) {
  
  RT_Tol_scan <- 3
  mz_Tol_ppm <- 20
  
  rt_min_scan <- min(which(daten@scantime > rt_min))
  rt_max_scan <- max(which(daten@scantime < rt_max))
  if (is.infinite(rt_min_scan)) rt_min_scan <- 1
  if (is.infinite(rt_max_scan)) rt_max_scan <- length(daten@scantime)
 
 
  massrange <- c(mz_min, mz_min+round((mz_max-mz_min)/10))
  pl <- NULL
  
  # Run peak picking algorithm
  plrows <- lapply(seq(mz_min, mz_max, by = mz_step*0.5), peakpicking_BfG_cpp,
               rawData = daten, mz_step = mz_step,
               rt_min_scan = rt_min_scan,
               rt_max_scan = rt_max_scan, 
               sn = sn,
               int_threshold = int_threshold,
               NoiseScans = peak_NoiseScans,
               precursormzTol = precursormzTol,
               peakwidth_min = peakwidth_min,
               peakwidth_max = peakwidth_max,
               maxPeaksPerSignal = maxPeaksPerSignal)
  
  pl <- do.call("rbind", plrows)

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
  
  # if nothing found, return empty peaklist
  if (is.null(pl)) 
    return(matrix(NA, 0, 18, dimnames = list(NULL, cn)))
  
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


#' Add additional peak and sample ID information to peaklists
#' 
#' This needs to be run before alignment
#'
#' @param peaklist 
#' @param datenList 
#'
#' @return New peaklist list in the same form as before but with additional columns
#' sample_id and peak_id_all. If these were already present, they are deleted and re-written
#' @export
new_peaklist <- function(peaklist, datenList) {
  peaklisttemp <- Map(function(pl, sid) transform(pl, sample_id = sid, peak_id_all = NA), 
                                  peaklist, seq_along(datenList))
  new_pl <- do.call("rbind", peaklisttemp)
  new_pl$peak_id_all <- seq_len(nrow(new_pl))
  split(new_pl, new_pl$sample_id)
}





#' Get a recommendation for the mz step parameter (binning width)
#' 
#' Calculates the optimum mz step based on mass differences. Takes a sample
#' of 100 spectra to make the calculation
#' 
#' @param daten raw data as list of xcmsRAW objects
#' @param probs quantile to use (see quantile function)
#' 
#' @return numeric. Width in Da
#' 
#' @export
optimumMzStep <- function(daten, probs) {
  get_diff <- function(scan) {
    msspektrum <- xcms::getScan(daten, scan)
    -quantile(-diff(msspektrum[, 1], lag = 1), probs = probs)
  }
  diffs <- vapply(sample(seq_along(daten@scanindex), 100), get_diff, numeric(1))
  round(min(diffs), 2)
}

