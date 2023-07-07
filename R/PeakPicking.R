


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



#### I AM HERE ####

#' Peakpicking algorithm
#' 
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
  
  #browser()
  # if nothing found, return empty peaklist
  if (is.null(pl)) 
    return(matrix(NA, 0, 18, dimnames = list(NULL, cn)))
  
  if (nrow(pl) >= 1) {
    
  pl <- pl[order(pl[,1]), , drop = F]
  
  #delete peaks smaller than peakwidth_min:
  #peaklist <- peaklist[((peaklist[,6]-peaklist[,5]) > peakwidth_min),,drop=FALSE]
  
  #delete peaks broader than peakwidth_max:
  #we will convert peakwidth_max into a maximum FWHM, in order to not delete tailing peaks
  #for this we assume (as an approximation) that FWHM is 50% of the full (baseline) peak width
  #peaklist <- peaklist[((peaklist[,12]-peaklist[,11]) < peakwidth_max/2),,drop=FALSE]
  
  #delete shoulder peaks which are inside the RT range of another peak:
  #browser()
  for (i in 1:nrow(pl)) {
    doppelte <- which((abs(pl[,1]-pl[i,1]) < pl[i,1]/1000000*mz_Tol_ppm) & 
                        (pl[,2] >= pl[i,5]) & (pl[,2] <= pl[i,6])) 
    doppelte <- doppelte[doppelte != i]
    pl[doppelte, 1] <- 0  
  }
  
  pl <- pl[pl[,1] > 0, , drop = F]
  
  }
  
  # calculation of SN
  if (nrow(pl) > 0)
    pl <- cbind(pl,pl[,3]/pl[,10],0)

  colnames(pl) <- cn
  
  # edit KJ 2018-12-07 ####
  # remove rows which are below intensity threshold and not within RT window
  pl <- pl[pl[, "Intensity"] >= int_threshold, , drop = F]
  pl <- pl[pl[, "RT"] >= rt_min, , drop = F]
  pl <- pl[pl[, "RT"] <= rt_max, , drop = F]
  
  as.data.frame(pl)
}


#' Correct peaklists for parallel computing
#'
#' @param peaklists 
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



alignment_BfG <- function(peaklists, ppm_dev, DeltaRT, ms2ppm, datenList){
  #browser()
  
  peaklist2 <- Map(function(pl, sid) transform(subset(pl, , mz:Gruppe), sample_id = sid), peaklist, seq_along(datenList))
  peaklists <- Reduce(rbind,peaklist2)
  
  ende <- FALSE
  spaltenlaenge <- max(table(peaklists$sample_id)) #nrow der längsten Peaklist
  SuchStart <- 1
  
  ergebnis <- matrix(, nrow = 0, ncol = max(peaklists$sample_id)*6)
  masse <- 0
  
  ii <- 0
  ms2vergleich <- vector()
  
 
  peaklists$Scan <- 0

  
  while (ende == FALSE) {
    ende <- TRUE
    iii <- SuchStart - 1
    jjj <<- iii
    
    ii <- ii + 1
    ergebnis <- rbind(ergebnis,0)
    abbruch <- FALSE
    ms2spectra <- list()
    
    # KJ: gets the first peak from the first peaklist as mz and rt
    while (abbruch == FALSE) {
      iii <- iii + 1
      for (i in 1:max(peaklists$sample_id)) {
        peaklists_i <- peaklists[peaklists$sample_id == i,]
        if (iii <= nrow(peaklists_i)) {
          #if ((peaklists[[i]][iii, "Scan"] == 0) & (abbruch == FALSE)) {
          if ((peaklists_i[iii, "Scan"] == 0) & (abbruch == FALSE)) {
            masse <- peaklists_i[iii, "mz"]
            RT <- peaklists_i[iii, "RT"]
            abbruch <- TRUE
            ende <- FALSE
          }  
        }
        if (iii > spaltenlaenge) 
          abbruch <- TRUE
      }
      SuchStart <- iii
      
    }#end while abbruch == false
    
    # get that mz and rt of the "highest" peak in all files
    kandidaten <- list()
    highest_i <- 1
    for (i in 1:max(peaklists$sample_id)) {
      kandidaten[[i]] <- which((abs(peaklists$mz - masse) < (masse*ppm_dev/1000000)) & 
                                 (abs(peaklists$RT - RT) < DeltaRT) & 
                                 (peaklists$Scan == 0) & (peaklists$sample_id == i))
      #order by intensity:
      kandidaten[[i]] <- kandidaten[[i]][order(peaklists[kandidaten[[i]], "Intensity"], decreasing=TRUE)]
      if (length(kandidaten[[highest_i]]) > 0) {
        if (length(kandidaten[[i]]) > 0) {
          if (peaklists$Intensity[kandidaten[[i]][1]] > peaklists[[highest_i]]$Intensity[kandidaten[[highest_i]][1]]) 
            highest_i <- i
        } 
      } else {
        highest_i <- i
      }
    }  
    
    masse <- peaklists[[highest_i]][kandidaten[[highest_i]][1],1]
    RT <- peaklists[[highest_i]][kandidaten[[highest_i]][1],2]
    
    #browser(expr = !is.na(masse) && isTRUE(all.equal(masse, 237.1016, tolerance = 0.005, scale = 1)))
    # KJ only check files in which this peak was found
    zuPruefende_i <- which(lengths(kandidaten) > 0)
    
    for (i in zuPruefende_i) {
      # if this is not the sample with the highest peak, check candidates again according to new mz and rt
      if (i != highest_i) {
        kandidaten[[i]] <- which((abs(peaklists[[i]]$mz-masse) < (masse*ppm/1000000)) & 
                                   (abs(peaklists[[i]]$RT-RT)<DeltaRT) & (peaklists[[i]]$Scan == 0))
      }  
      
      if (length(kandidaten[[i]]) > 0) {
        passendeZeile <- kandidaten[[i]][which.min(abs(peaklists[[i]]$RT[kandidaten[[i]]] - 
                                                         peaklists[[highest_i]]$RT[kandidaten[[highest_i]][1]]))]
        ergebnis[ii,i*6-5] <- peaklists[[i]]$peak_id_all[passendeZeile]  # peak ID from peaklist
        ergebnis[ii,i*6-4] <- peaklists[[i]]$mz[passendeZeile]
        ergebnis[ii,i*6-3] <- peaklists[[i]]$RT[passendeZeile]
        ergebnis[ii,i*6-2] <- peaklists[[i]]$Intensity[passendeZeile]
        ergebnis[ii,i*6-1] <- peaklists[[i]]$MS2scan[passendeZeile]
        ergebnis[ii,i*6] <- peaklists[[i]]$Gruppe[passendeZeile]
        peaklists[[i]]$Scan[passendeZeile] <- 1
      }
      
      #if ((ergebnis[ii,i*6-1] > 0)) {
      #  RTMS2 <- datenList[[i]]@msnRt[peaklists[[i]][passendeZeile, "MS2scan"]]
      #  if ((RTMS2 > peaklists[[i]][passendeZeile, "LeftendRT"]) & 
      #      (RTMS2 < peaklists[[i]][passendeZeile, "RightendRT"])) {
      #    
      #    ms2spectra[[i]] <- xcms::getMsnScan(datenList[[i]], ergebnis[ii,i*6-1])
      #    if (any(is.na(ms2spectra[[i]]))) 
      #      ms2spectra[[i]] <- NULL
      #  } else {
      #    ms2spectra[[i]] <- NULL
      #    
      #    ergebnis[ii,i*6-1] <- 0
      #  }
      #}  
    }  
    
    #remove empty ms2 spectra:
    #ms2spectra <- Filter(Negate(function(x) is.null(unlist(x))), ms2spectra)
    
    #if there are at least two samples with ms2 spectra, compare them:
    #if (length(ms2spectra) > 1) {
    #  ms2vergleich[ii] <- spektraVergleichen(ms2spectra,ms2ppm)
    #} else {
    #  ms2vergleich[ii] <- 0
    #}
    
    
  }#end while ende == false
  
  
  mean_mz <- vector()
  mean_RT <- vector()
  Gruppe <- numeric(nrow(ergebnis))
  for (i in 1:nrow(ergebnis)) {
    mean_mz[i] <- mean(ergebnis[i,which(ergebnis[i,seq(2,ncol(ergebnis),by=6)] > 0)*6-4])
    mean_RT[i] <- mean(ergebnis[i,which(ergebnis[i,seq(3,ncol(ergebnis),by=6)] > 0)*6-3])
  }
  #ergebnis <- cbind(mean_mz,mean_RT,ms2vergleich,Gruppe,ergebnis)
  ergebnis <- cbind(mean_mz,mean_RT,0,Gruppe,ergebnis)
  
  spaltennamen <- c("mean_mz","mean_RT","MS2Fit", "Gruppe")
  for (i in 1:(ncol((ergebnis)-4)/6)) {
    spaltennamen <- c(spaltennamen,
                      paste0("peak_id_all_",as.character(i)),
                      paste0("mz_",as.character(i)),
                      paste0("RT_",as.character(i)),
                      paste0("Int_",as.character(i)),
                      paste0("ms2scan_",as.character(i)),
                      paste0("gruppe_",as.character(i)))  # componentization information
  }
  colnames(ergebnis) <- spaltennamen
  
  # Vereinheitlichung der Gruppen
  gruppenzaehler <- 1
  
  #for (i in 1:nrow(ergebnis)){
  #  if (ergebnis[i,"Gruppe"]==0){
  #    # max_Int = col with highest Int in row i
  #    max_Int <- grep("Int_", colnames(ergebnis)) [which.max(ergebnis[i,grep("Int_", colnames(ergebnis))])]
  #    max_gruppe <- paste0("gruppe_",stringr::str_match(colnames(ergebnis)[max_Int], "Int_(\\d+)$")[,2])
  #    # b = rows with Gruppe == 0 AND same groups found in first samples of a
  #    b <- which(ergebnis[,"Gruppe"]==0 & ergebnis[,max_gruppe]==ergebnis[i,max_gruppe])
  #    ergebnis[b,"Gruppe"] <- gruppenzaehler
  #    gruppenzaehler <- gruppenzaehler+1}}
  
  # Eliminieren von ein-Peak-Gruppen  (Gruppe wird 0)
  #for (i in 1:max(ergebnis[,"Gruppe"])){
  #  if (length(which(ergebnis[,"Gruppe"]==i))==1) ergebnis[which(ergebnis[,"Gruppe"]==i),"Gruppe"] <-0}
  
  # Schließen der entstandenen Zwischenräume
  #V1 <- sort(unique(ergebnis[,"Gruppe"])[-(which(unique(ergebnis[,"Gruppe"])==0))])
  #V2 <- seq_along(V1) 
  #for (i in 1:length(V1)){
  #  ergebnis[which(ergebnis[,"Gruppe"]==V1[i]),"Gruppe"] <- V2[i]}
  
  
  ergebnis <- ergebnis[!is.na(ergebnis[,1]), ]
  
  neue_gruppen <- matrix(0,nrow = nrow(ergebnis), ncol = ncol((ergebnis)-4)/6)
  #ii <- 1
  gruppennr <- 1
  for (ii in 1:nrow(ergebnis)) {
    if (any(neue_gruppen[ii,] == 0)) {
      for (i in 1:(ncol((ergebnis)-4)/6)) {
        if ((ergebnis[ii,i*6+4] > 0) & (neue_gruppen[ii,i]) == 0) {
          neue_gruppen[which(ergebnis[,i*6+4] == ergebnis[ii,i*6+4]),i] <- gruppennr
        }
      }
      gruppennr <- gruppennr+1    
    }
  }
  
  ungleich <- NULL
  for (i in 1:nrow(neue_gruppen)) {
    if (length(unique(neue_gruppen[i,neue_gruppen[i,] != 0])) > 1) ungleich <- c(ungleich, i)
  }
  
  fehlstellen <- NULL
  for (i in 1:nrow(neue_gruppen)) {
    if (any(neue_gruppen[i,] == 0)) fehlstellen <- c(fehlstellen, i)
  }
  
  return(ergebnis)
}


#' @export
grouping_BfG <- function(peaklists, ppm, DeltaRT, ms2ppm, datenList){
  #browser()
  ende <- FALSE
  spaltenlaenge <- 1 #nrow der l?ngsten Peaklist
  SuchStart <- 1

  ergebnis <- matrix(, nrow = 0, ncol = length(peaklists)*6)
  masse <- 0
  
  ii <- 0
  ms2vergleich <- vector()
  
  for (i in 1:length(peaklists)) {
    if (spaltenlaenge < nrow(peaklists[[i]])) 
      spaltenlaenge <- nrow(peaklists[[i]])
    peaklists[[i]]$Scan <- 0
  }
  
  while (ende == FALSE) {
    ende <- TRUE
    iii <- SuchStart - 1
    jjj <<- iii

    ii <- ii + 1
    ergebnis <- rbind(ergebnis,0)
    abbruch <- FALSE
    ms2spectra <- list()

    # KJ: gets the first peak from the first peaklist as mz and rt
    while (abbruch == FALSE) {
      iii <- iii + 1
      for (i in 1:length(peaklists)) {
        if (iii <= nrow(peaklists[[i]])) {
          if ((peaklists[[i]][iii, "Scan"] == 0) & (abbruch == FALSE)) {
            masse <- peaklists[[i]][iii, "mz"]
            RT <- peaklists[[i]][iii, "RT"]
            abbruch <- TRUE
            ende <- FALSE
          }  
        }
        if (iii > spaltenlaenge) 
          abbruch <- TRUE
      }
      SuchStart <- iii
      
    }#end while abbruch == false
    
    # get that mz and rt of the "highest" peak in all files
    kandidaten <- list()
    highest_i <- 1
    for (i in 1:length(peaklists)) {
      kandidaten[[i]] <- which((abs(peaklists[[i]]$mz - masse) < (masse*ppm/1000000)) & 
                                 (abs(peaklists[[i]]$RT - RT) < DeltaRT) & 
                                 (peaklists[[i]]$Scan == 0))
      #order by intensity:
      kandidaten[[i]] <- kandidaten[[i]][order(peaklists[[i]][kandidaten[[i]], "Intensity"], decreasing=TRUE)]
      if (length(kandidaten[[highest_i]]) > 0) {
        if (length(kandidaten[[i]]) > 0) {
          if (peaklists[[i]]$Intensity[kandidaten[[i]][1]] > peaklists[[highest_i]]$Intensity[kandidaten[[highest_i]][1]]) 
            highest_i <- i
        } 
      } else {
        highest_i <- i
      }
    }  
    
    masse <- peaklists[[highest_i]][kandidaten[[highest_i]][1],1]
    RT <- peaklists[[highest_i]][kandidaten[[highest_i]][1],2]
    
    #browser(expr = !is.na(masse) && isTRUE(all.equal(masse, 237.1016, tolerance = 0.005, scale = 1)))
    # KJ only check files in which this peak was found
    zuPruefende_i <- which(lengths(kandidaten) > 0)
    
    for (i in zuPruefende_i) {
      # if this is not the sample with the highest peak, check candidates again according to new mz and rt
      if (i != highest_i) {
        kandidaten[[i]] <- which((abs(peaklists[[i]]$mz-masse) < (masse*ppm/1000000)) & 
                                   (abs(peaklists[[i]]$RT-RT)<DeltaRT) & (peaklists[[i]]$Scan == 0))
      }  
      
      if (length(kandidaten[[i]]) > 0) {
        passendeZeile <- kandidaten[[i]][which.min(abs(peaklists[[i]]$RT[kandidaten[[i]]] - 
                                                         peaklists[[highest_i]]$RT[kandidaten[[highest_i]][1]]))]
        ergebnis[ii,i*6-5] <- peaklists[[i]]$peak_id_all[passendeZeile]  # peak ID from peaklist
        ergebnis[ii,i*6-4] <- peaklists[[i]]$mz[passendeZeile]
        ergebnis[ii,i*6-3] <- peaklists[[i]]$RT[passendeZeile]
        ergebnis[ii,i*6-2] <- peaklists[[i]]$Intensity[passendeZeile]
        ergebnis[ii,i*6-1] <- peaklists[[i]]$MS2scan[passendeZeile]
        ergebnis[ii,i*6] <- peaklists[[i]]$Gruppe[passendeZeile]
        peaklists[[i]]$Scan[passendeZeile] <- 1
      }
      
      #if ((ergebnis[ii,i*6-1] > 0)) {
      #  RTMS2 <- datenList[[i]]@msnRt[peaklists[[i]][passendeZeile, "MS2scan"]]
      #  if ((RTMS2 > peaklists[[i]][passendeZeile, "LeftendRT"]) & 
      #      (RTMS2 < peaklists[[i]][passendeZeile, "RightendRT"])) {
      #    
      #    ms2spectra[[i]] <- xcms::getMsnScan(datenList[[i]], ergebnis[ii,i*6-1])
      #    if (any(is.na(ms2spectra[[i]]))) 
      #      ms2spectra[[i]] <- NULL
      #  } else {
      #    ms2spectra[[i]] <- NULL
      #    
      #    ergebnis[ii,i*6-1] <- 0
      #  }
      #}  
    }  
    
    #remove empty ms2 spectra:
    #ms2spectra <- Filter(Negate(function(x) is.null(unlist(x))), ms2spectra)
    
    #if there are at least two samples with ms2 spectra, compare them:
    #if (length(ms2spectra) > 1) {
    #  ms2vergleich[ii] <- spektraVergleichen(ms2spectra,ms2ppm)
    #} else {
    #  ms2vergleich[ii] <- 0
    #}
    
    
  }#end while ende == false

  
  mean_mz <- vector()
  mean_RT <- vector()
  Gruppe <- numeric(nrow(ergebnis))
  for (i in 1:nrow(ergebnis)) {
    mean_mz[i] <- mean(ergebnis[i,which(ergebnis[i,seq(2,ncol(ergebnis),by=6)] > 0)*6-4])
    mean_RT[i] <- mean(ergebnis[i,which(ergebnis[i,seq(3,ncol(ergebnis),by=6)] > 0)*6-3])
  }
  #ergebnis <- cbind(mean_mz,mean_RT,ms2vergleich,Gruppe,ergebnis)
  ergebnis <- cbind(mean_mz,mean_RT,0,Gruppe,ergebnis)
 
  spaltennamen <- c("mean_mz","mean_RT","MS2Fit", "Gruppe")
  for (i in 1:(ncol((ergebnis)-4)/6)) {
    spaltennamen <- c(spaltennamen,
                      paste0("peak_id_all_",as.character(i)),
                      paste0("mz_",as.character(i)),
                      paste0("RT_",as.character(i)),
                      paste0("Int_",as.character(i)),
                      paste0("ms2scan_",as.character(i)),
                      paste0("gruppe_",as.character(i)))  # componentization information
  }
  colnames(ergebnis) <- spaltennamen
  
  # Vereinheitlichung der Gruppen
  gruppenzaehler <- 1
  
  #for (i in 1:nrow(ergebnis)){
  #  if (ergebnis[i,"Gruppe"]==0){
  #    # max_Int = col with highest Int in row i
  #    max_Int <- grep("Int_", colnames(ergebnis)) [which.max(ergebnis[i,grep("Int_", colnames(ergebnis))])]
  #    max_gruppe <- paste0("gruppe_",stringr::str_match(colnames(ergebnis)[max_Int], "Int_(\\d+)$")[,2])
  #    # b = rows with Gruppe == 0 AND same groups found in first samples of a
  #    b <- which(ergebnis[,"Gruppe"]==0 & ergebnis[,max_gruppe]==ergebnis[i,max_gruppe])
  #    ergebnis[b,"Gruppe"] <- gruppenzaehler
  #    gruppenzaehler <- gruppenzaehler+1}}
  
  # Eliminieren von ein-Peak-Gruppen  (Gruppe wird 0)
  #for (i in 1:max(ergebnis[,"Gruppe"])){
  #  if (length(which(ergebnis[,"Gruppe"]==i))==1) ergebnis[which(ergebnis[,"Gruppe"]==i),"Gruppe"] <-0}
  
  # Schließen der entstandenen Zwischenräume
  #V1 <- sort(unique(ergebnis[,"Gruppe"])[-(which(unique(ergebnis[,"Gruppe"])==0))])
  #V2 <- seq_along(V1) 
  #for (i in 1:length(V1)){
  #  ergebnis[which(ergebnis[,"Gruppe"]==V1[i]),"Gruppe"] <- V2[i]}
  
  
  ergebnis <- ergebnis[!is.na(ergebnis[,1]), ]
  
  neue_gruppen <- matrix(0,nrow = nrow(ergebnis), ncol = ncol((ergebnis)-4)/6)
  #ii <- 1
  gruppennr <- 1
  for (ii in 1:nrow(ergebnis)) {
    if (any(neue_gruppen[ii,] == 0)) {
  for (i in 1:(ncol((ergebnis)-4)/6)) {
    if ((ergebnis[ii,i*6+4] > 0) & (neue_gruppen[ii,i]) == 0) {
      neue_gruppen[which(ergebnis[,i*6+4] == ergebnis[ii,i*6+4]),i] <- gruppennr
    }
  }
  gruppennr <- gruppennr+1    
  }
  }
  
  ungleich <- NULL
  for (i in 1:nrow(neue_gruppen)) {
    if (length(unique(neue_gruppen[i,neue_gruppen[i,] != 0])) > 1) ungleich <- c(ungleich, i)
  }
  
  fehlstellen <- NULL
  for (i in 1:nrow(neue_gruppen)) {
    if (any(neue_gruppen[i,] == 0)) fehlstellen <- c(fehlstellen, i)
  }
  
  return(ergebnis)
}

#' @export
FindNeutralLoss <- function(selected,loss,ppm=30, RT_Tol=10){
  
  #specificFragmentations        <- c(18.01056,17.0266,43.9898)
  #names(specificFragmentations) <- c("H2O"   ,"NH3"  ,"CO2")
  
  #if (is.numeric(loss)) {
  #  specificFragmentations <- c(specificFragmentations,loss)
  #  names(specificFragmentations) <- c(names(specificFragmentations), "other")
  #  loss <- "other"
  #}
  
  
  ergebnis <- matrix(FALSE,nrow=nrow(peaklist[[selected]]),ncol=length(loss))
  #colnames(ergebnis) <- loss
  colnames(ergebnis) <- names(loss)
  
  #InSourceFragmentation <<- vector()
  #InSourceFragmentation[1:nrow(peaklist[[selected]])] <<- 0
  
  ausgabe <- NULL
  if (length(loss) > 0) {
    for (i in 1:nrow(peaklist[[selected]])) {
      if (peaklist[[selected]][i,15] > 0) {
        mutterion <- datenList[[selected]]@msnPrecursorMz[peaklist[[selected]][i,15]]
        ms2 <- getMsnScan(datenList[[selected]],peaklist[[selected]][i,15])
        ms2[,1] <- ms2[,1]-mutterion
        for (fragment in 1:length(loss)) {
          #ergebnis[i,fragment] <- any(abs(ms2[,1]+specificFragmentations[which(names(specificFragmentations)==loss[fragment])]) < mutterion/1000000*ppm)
          ergebnis[i,fragment] <- any(abs(ms2[,1]+loss[fragment]) < mutterion/1000000*ppm)
        }
      }
      ausgabe <- c(ausgabe,all(ergebnis[i,]))
    }
    
    if (length(loss) == 1) {
      for (i in which(ausgabe)) {
        CheckForInSourceFragmenation <- which((abs(peaklist[[selected]][,1]-(peaklist[[selected]][i,1]-loss)) < peaklist[[selected]][i,1]/1000000*ppm) & (abs(peaklist[[selected]][,2]-peaklist[[selected]][i,2]) < RT_Tol))
        #if (length(CheckForInSourceFragmenation) > 0) InSourceFragmentation[CheckForInSourceFragmenation] <<- i
        if (length(CheckForInSourceFragmenation) > 0) peaklist[[selected]][CheckForInSourceFragmenation,17] <<- i
      }  
    }
    #return(which(ergebnis))
    #return(ergebnis)
    return(which(ausgabe))
  }
  else
  {
    return(NULL)
  }
}

#' @export
FindFragmentIons <- function(selected,fragment_ions,ppm=30, RT_Tol=10){
  ergebnis <- matrix("",nrow=0,ncol=5)
  if (length(fragment_ions) > 0) {
    for (i in 1:nrow(peaklist[[selected]])) {
      if (peaklist[[selected]][i,15] > 0) {
        RT <- datenList[[selected]]@msnRt[peaklist[[selected]][i,15]]
        if ((RT >= peaklist[[selected]][i,5]) & (RT <= peaklist[[selected]][i,6])) {
          mutterion <- datenList[[selected]]@msnPrecursorMz[peaklist[[selected]][i,15]]
          ms2 <- getMsnScan(datenList[[selected]],peaklist[[selected]][i,15])
          ms2 <- ms2[ms2[,2]>max(ms2[,2])*0.1,,drop=FALSE]
          
          for (fragment in 1:nrow(ms2)) {
            fragment_ion_ID <- which(abs(ms2[fragment,1]-fragment_ions[]) < ms2[fragment,1]/1000000*ppm)
            if (length(fragment_ion_ID) == 1) {
              ergebnis <- rbind(ergebnis,c(i, ms2[fragment,1], names(fragment_ion_ID), (ms2[fragment,1]-fragment_ions[fragment_ion_ID])/1000000,fragment_ion_ID))
            }
            CheckForInSourceFragmenation <- which((peaklist[[selected]][,1] < (mutterion-mutterion/1000000*ppm)) & (abs(peaklist[[selected]][,1]-ms2[fragment,1]) < peaklist[[selected]][i,1]/1000000*ppm) & (abs(peaklist[[selected]][,2]-peaklist[[selected]][i,2]) < RT_Tol))
            if (length(CheckForInSourceFragmenation) > 0) peaklist[[selected]][CheckForInSourceFragmenation,17] <<- i
          }  
          
        }
      }
    }
    colnames(ergebnis) <- c("peaklistID","fragment-mz","Fragment-ion","ppm","FragmentIon_ID")
    return(ergebnis)
  }  
  else
  {
    return(NULL)
  }
  
}

#' @export
FindNeutralLoss2 <- function(selected,loss,ppm=30, RT_Tol=10){
  
  ergebnis <- matrix("",nrow=0,ncol=7)
  
  if (length(loss) > 0) {
    for (i in 1:nrow(peaklist[[selected]])) {
      if (peaklist[[selected]][i,15] > 0) {
        RT <- datenList[[selected]]@msnRt[peaklist[[selected]][i,15]]
        if ((RT >= peaklist[[selected]][i,5]) & (RT <= peaklist[[selected]][i,6])) {  
          mutterion <- datenList[[selected]]@msnPrecursorMz[peaklist[[selected]][i,15]]
          ms2 <- getMsnScan(datenList[[selected]],peaklist[[selected]][i,15])
          ms2 <- ms2[ms2[,2]>max(ms2[,2])*0.1,,drop=FALSE]
          
          for (fragment in 1:nrow(ms2)) {
            fragmenttype <- which(abs(ms2[fragment,1]+loss[]-mutterion) < mutterion/1000000*ppm)
            if (length(fragmenttype) == 1) {
              ergebnis <- rbind(ergebnis,c(i,ms2[fragment,1],mutterion-ms2[fragment,1],loss[fragmenttype],names(fragmenttype),fragmenttype,0))
            } 
            if (length(fragmenttype) > 1) {
              for (ii in 1:length(fragmenttype)) {
                ergebnis <- rbind(ergebnis,c(i,ms2[fragment,1],mutterion-ms2[fragment,1],loss[fragmenttype[ii]],names(fragmenttype[ii]),fragmenttype[ii],0))
              } 
            }
            if (length(fragmenttype) == 0) {
              for (fragmenttype in 1:length(loss)) { 
                fragmenttype2 <- which(abs(ms2[fragment,1]+loss[]+loss[fragmenttype]-mutterion) < mutterion/1000000*ppm)
                if (length(fragmenttype2) == 1) {
                  ergebnis <- rbind(ergebnis,c(i,ms2[fragment,1],mutterion-ms2[fragment,1],loss[fragmenttype]+loss[fragmenttype2],paste0(names(loss[fragmenttype]),"/",names(fragmenttype2)),fragmenttype,fragmenttype2))
                }
                if (length(fragmenttype2) > 1) {
                  for (ii in 1:length(fragmenttype2)) {
                    ergebnis <- rbind(ergebnis,c(i,ms2[fragment,1],mutterion-ms2[fragment,1],loss[fragmenttype]+loss[fragmenttype2[ii]],paste0(names(loss[fragmenttype]),"/",names(fragmenttype2[ii])),fragmenttype,fragmenttype2[ii]))
                  } 
                }
              }
            }
          }
        }
      }
    }
    
    colnames(ergebnis) <- c("peaklistID","fragment-mz","delta_mz_real","delta_mz_theory","NeutralLoss","NeutralLossID1","NeutralLossID2")
    return(ergebnis)
  }
  else
  {
    return(NULL)
  }
}

#' @export
hasSpecificPattern <- function(masse, inputspectrum, isotopePattern, mztol = 20, inttol = 0.10) { 
  output <- FALSE #we first assume that masse does not have a Cl or Br pattern
  #isotope <- generateSinglePattern(Cl = Cl, Br = Br)
  #transfer isotope masses from generateSinglePattern (pure ClBr compounds) to mass of "masse"
  for (i in nrow(isotopePattern):1) isotopePattern[i,1] <- isotopePattern[i,1]-isotopePattern[1,1]+masse
  inputspectrum <- inputspectrum[which((inputspectrum[,1] > masse-masse*mztol/1000000) & (inputspectrum[,1] < isotopePattern[nrow(isotopePattern),1]+isotopePattern[nrow(isotopePattern),1]*mztol/1000000)),]
  if (length(inputspectrum) > 2) {  
    #among all masses within mztol around "masse" get the one with highest intensity
    peakmassrow <- which.max(inputspectrum[which((abs(inputspectrum[,1]-masse) < masse*mztol/1000000)),2])
    peakmass <- inputspectrum[peakmassrow,1]
    NrOfIsotopesFound <- 1 #the monoisotopic mass is there for sure, otherwise we wouldn't search for patterns
    intensitaet <- inputspectrum[which.max(inputspectrum[which((abs(inputspectrum[,1]-masse) < masse*mztol/1000000)),2]),2]
    
    iii <- peakmassrow+1 #iii will be the counter to scan inputspectrum for isotopes. we will start our scan at row peakmassrow (the monoisotopic mass)
    
    #neueZeile <- cbind(inputspectrum[ii,1],inputspectrum[ii,2]) #prepare a new row for the output matrix. we will check later if we need it
    for (isotopnr in 2: nrow(isotopePattern)) { #scan all isotopes of the pattern, except the monoisotopic mass
      iii <- iii-1
      while ((inputspectrum[iii,1] < isotopePattern[isotopnr,1]+peakmass*mztol/1000000) & (iii < nrow(inputspectrum))) {
        #check rows of inputspectrum until mass+mztol exceeds isotopic mass
        iii <- iii+1 #go for next mass in inputspectrum
        if ((abs(inputspectrum[iii,1] - isotopePattern[isotopnr,1]) < peakmass*mztol/1000000) & (abs(inputspectrum[iii,2]/intensitaet - isotopePattern[isotopnr,2]/isotopePattern[1,2]) < inttol)) {
          NrOfIsotopesFound <- NrOfIsotopesFound+1
          #if mass and intensity of this row in inputspectrum fits isotopic pattern, another isotope was found
        } #ende if 
      } #ende while
    } # ende for isotopnr  
    
    if (NrOfIsotopesFound == nrow(isotopePattern)) { #only if all isotopes of the pattern were found output is TRUE
      output <- TRUE 
    } 
    
  } #ende if nrow(inputspectrum) > 0
  return(output)
} # ende function

#' @export
IsotopesAndAdducts <- function(peakliste, daten, ppm = 10, RT_Tol = 10){
  C13 <- vector(mode = "numeric", length = nrow(peakliste))
  NaAddukt <- vector(mode = "numeric", length = nrow(peakliste))
  KAddukt <- vector(mode = "numeric", length = nrow(peakliste))
  NH4Addukt <- vector(mode = "numeric", length = nrow(peakliste))
  EthylaminAddukt <- vector(mode = "numeric", length = nrow(peakliste))
  S1 <- vector(mode = "numeric", length = nrow(peakliste))
  S2 <- vector(mode = "numeric", length = nrow(peakliste))
  Cl1 <- vector(mode = "numeric", length = nrow(peakliste))
  Cl2 <- vector(mode = "numeric", length = nrow(peakliste))
  Cl3 <- vector(mode = "numeric", length = nrow(peakliste))
  Cl4 <- vector(mode = "numeric", length = nrow(peakliste))
  Br1 <- vector(mode = "numeric", length = nrow(peakliste))
  Br2 <- vector(mode = "numeric", length = nrow(peakliste))
  Br3 <- vector(mode = "numeric", length = nrow(peakliste))
  Br4 <- vector(mode = "numeric", length = nrow(peakliste))
  inttol <- 0.10 #intensity tolerance for isotope patterns (Cl, Br)
  inttol_S <- 0.02 #intensity tolerance for isotope patterns (Cl, Br)
  
  
  
  #withProgress(message = 'Isotope & Adduct Search', min = 0, max = nrow(peakliste), value = 0, {
    for (i in 1:nrow(peakliste)) {
      #setProgress(value = i) 
      spektrum <- xcms::getScan(daten,peakliste[i,4])
      
      #13C:
      #search for potential 13C entries: (intensity of the proposed 13C peak must be below (monoisotpic mass)/12 as maximum (minimum one H per C => divided by 13; 1.1% 13Cp peak intensity per C => 13/1.1 = 12)) 
      zwischenergebnis <- which((abs(peakliste[,1]-(peakliste[i,1]+1.00335)) < peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (peakliste[,3]/peakliste[i,3]*100 <= peakliste[i,1]/12))
      if (length(zwischenergebnis) > 0) {
        for (j in 1:length(zwischenergebnis)) {
          if (peakliste[i,4] != peakliste[zwischenergebnis[j],4]) { #if the 12C and 13C are found at different scan numbers
            #get the IDs of the 12C and 13C of the spectrum at the first scan number:
            spektrum_monoisotopic <- which(spektrum[,1]==peakliste[i,1])
            spektrum_isotope <- which.min(abs(spektrum[,1]-peakliste[zwischenergebnis[j],1]))
            #get the IDs of the monoisotopic and the isotope of the spectrum at the second scan number:
            spektrum2 <- getScan(daten,peakliste[zwischenergebnis[j],4])
            spektrum2_monoisotopic <- which.min(abs(spektrum2[,1]-peakliste[i,1]))
            spektrum2_isotope <- which(spektrum2[,1]==peakliste[zwischenergebnis[j],1])
            #only keep the 13C if the relative intensities in both scans are within the tolerance:
            if (abs(spektrum[spektrum_isotope,2]/spektrum[spektrum_monoisotopic,2]-spektrum2[spektrum2_isotope,2]/spektrum2[spektrum2_monoisotopic,2]) > inttol) zwischenergebnis[j] <- 0
          } 
        }
        C13[zwischenergebnis[zwischenergebnis > 0]] <- i
      }
      
      
      
      #with 1 sulfur
      isotopePattern <- cbind(c(0.0,1.9958),c(1.0,0.045))
      isotopes <- vector() #we will store the IDs of the isotopes here
      zwischenergebnis <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]*isotopePattern[2,2])) < peakliste[i,3]*inttol_S))
      if (length(zwischenergebnis) > 0) {
        for (j in 1:length(zwischenergebnis)) {
          if (peakliste[i,4] != peakliste[zwischenergebnis[j],4]) { #if the monoisotopic and the isotope are found at different scan numbers
            #get the IDs of the monoisotopic and the isotope of the spectrum at the first scan number:
            spektrum_monoisotopic <- which(spektrum[,1]==peakliste[i,1])
            spektrum_isotope <- which.min(abs(spektrum[,1]-peakliste[zwischenergebnis[j],1]))
            #get the IDs of the monoisotopic and the isotope of the spectrum at the second scan number:
            spektrum2 <- xcms::getScan(daten,peakliste[zwischenergebnis[j],4])
            spektrum2_monoisotopic <- which.min(abs(spektrum2[,1]-peakliste[i,1]))
            spektrum2_isotope <- which(spektrum2[,1]==peakliste[zwischenergebnis[j],1])
            #only keep the isotope if the relative intensities in both scans are within the tolerance:
            if (abs(spektrum[spektrum_isotope,2]/spektrum[spektrum_monoisotopic,2]-spektrum2[spektrum2_isotope,2]/spektrum2[spektrum2_monoisotopic,2]) > inttol_S) zwischenergebnis[j] <- 0
          }
        }
        if (any(zwischenergebnis > 0)) { #if there is one peak in "zwischenergebnis" which fits, remember it. 
          isotopes <- c(isotopes, zwischenergebnis[zwischenergebnis > 0])
        }
        S1[i] <- i
        S1[isotopes] <- i
      }
      
      
      #with 2 sulfur
      isotopePattern <- cbind(c(0.0,1.9958),c(1.0,0.09))
      isotopes <- vector() #we will store the IDs of the isotopes here
      zwischenergebnis <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]*isotopePattern[2,2])) < peakliste[i,3]*inttol_S))
      if (length(zwischenergebnis) > 0) {
        for (j in 1:length(zwischenergebnis)) {
          if (peakliste[i,4] != peakliste[zwischenergebnis[j],4]) { #if the monoisotopic and the isotope are found at different scan numbers
            #get the IDs of the monoisotopic and the isotope of the spectrum at the first scan number:
            spektrum_monoisotopic <- which(spektrum[,1]==peakliste[i,1])
            spektrum_isotope <- which.min(abs(spektrum[,1]-peakliste[zwischenergebnis[j],1]))
            #get the IDs of the monoisotopic and the isotope of the spectrum at the second scan number:
            spektrum2 <- xcms::getScan(daten,peakliste[zwischenergebnis[j],4])
            spektrum2_monoisotopic <- which.min(abs(spektrum2[,1]-peakliste[i,1]))
            spektrum2_isotope <- which(spektrum2[,1]==peakliste[zwischenergebnis[j],1])
            #only keep the isotope if the relative intensities in both scans are within the tolerance:
            if (abs(spektrum[spektrum_isotope,2]/spektrum[spektrum_monoisotopic,2]-spektrum2[spektrum2_isotope,2]/spektrum2[spektrum2_monoisotopic,2]) > inttol_S) zwischenergebnis[j] <- 0
          }
        }
        if (any(zwischenergebnis > 0)) { #if there is one peak in "zwischenergebnis" which fits, remember it. 
          isotopes <- c(isotopes, zwischenergebnis[zwischenergebnis > 0])
        }
        S2[i] <- i
        S2[isotopes] <- i
      }
      
      #Chlorinated (1 time)
      isotopePattern <- cbind(c(0.0,1.99705),c(1.0,0.32))
      #first check if there is an isotopic pattern. Otherwise we might assign peaks due to an only partially fitting isotopic pattern (e.g. an Cl-[M+4] peak without an [M+2])
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        isIsotopicPattern <- TRUE
        isotopes <- vector() #we will store the IDs of the isotopes here
        zwischenergebnis <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]*isotopePattern[2,2])) < peakliste[i,3]*inttol))
        if (length(zwischenergebnis) > 0) {
          for (j in 1:length(zwischenergebnis)) {
            if (peakliste[i,4] != peakliste[zwischenergebnis[j],4]) { #if the monoisotopic and the isotope are found at different scan numbers
              #get the IDs of the monoisotopic and the isotope of the spectrum at the first scan number:
              spektrum_monoisotopic <- which(spektrum[,1]==peakliste[i,1])
              spektrum_isotope <- which.min(abs(spektrum[,1]-peakliste[zwischenergebnis[j],1]))
              #get the IDs of the monoisotopic and the isotope of the spectrum at the second scan number:
              spektrum2 <- xcms::getScan(daten,peakliste[zwischenergebnis[j],4])
              spektrum2_monoisotopic <- which.min(abs(spektrum2[,1]-peakliste[i,1]))
              spektrum2_isotope <- which(spektrum2[,1]==peakliste[zwischenergebnis[j],1])
              #only keep the isotope if the relative intensities in both scans are within the tolerance:
              if (abs(spektrum[spektrum_isotope,2]/spektrum[spektrum_monoisotopic,2]-spektrum2[spektrum2_isotope,2]/spektrum2[spektrum2_monoisotopic,2]) > inttol) zwischenergebnis[j] <- 0
            }
          }
          if (any(zwischenergebnis > 0)) { #if there is one peak in "zwischenergebnis" which fits, remember it. 
            isotopes <- c(isotopes, zwischenergebnis[zwischenergebnis > 0])
          } else { #otherwise this is isotopic pattern
            isIsotopicPattern <- FALSE
          }
        }
        if (isIsotopicPattern) {
          Cl1[i] <- i
          Cl1[isotopes] <- i
        }
      }  
      
      #Chlorinated (2 times)
      isotopePattern <- cbind(c(0.0,1.99705,3.99410),c(1.0,0.64,0.1))
      #first check if there is an isotopic pattern. Otherwise we might assign peaks due to an only partially fitting isotopic pattern (e.g. an Cl-[M+4] peak without an [M+2])
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        isIsotopicPattern <- TRUE
        isotopes <- vector() #we will store the IDs of the isotopes here
        for (isotopeNumber in 2:3) { #search for all isotopes
          zwischenergebnis <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[isotopeNumber,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]*isotopePattern[isotopeNumber,2])) < peakliste[i,3]*inttol))
          if (length(zwischenergebnis) > 0) {
            for (j in 1:length(zwischenergebnis)) {
              if (peakliste[i,4] != peakliste[zwischenergebnis[j],4]) { #if the monoisotopic and the isotope are found at different scan numbers
                #get the IDs of the monoisotopic and the isotope of the spectrum at the first scan number:
                spektrum_monoisotopic <- which(spektrum[,1]==peakliste[i,1])
                spektrum_isotope <- which.min(abs(spektrum[,1]-peakliste[zwischenergebnis[j],1]))
                #get the IDs of the monoisotopic and the isotope of the spectrum at the second scan number:
                spektrum2 <- xcms::getScan(daten,peakliste[zwischenergebnis[j],4])
                spektrum2_monoisotopic <- which.min(abs(spektrum2[,1]-peakliste[i,1]))
                spektrum2_isotope <- which(spektrum2[,1]==peakliste[zwischenergebnis[j],1])
                #only keep the isotope if the relative intensities in both scans are within the tolerance:
                if (abs(spektrum[spektrum_isotope,2]/spektrum[spektrum_monoisotopic,2]-spektrum2[spektrum2_isotope,2]/spektrum2[spektrum2_monoisotopic,2]) > inttol) zwischenergebnis[j] <- 0
              }
            }
            if (any(zwischenergebnis > 0)) { #if there is one peak in "zwischenergebnis" which fits, remember it. 
              isotopes <- c(isotopes, zwischenergebnis[zwischenergebnis > 0])
            } else { #otherwise this is isotopic pattern
              isIsotopicPattern <- FALSE
            }
          }
        }
        if (isIsotopicPattern) {
          Cl2[i] <- i
          Cl2[isotopes] <- i
        }
      } 
      
      #Chlorinated (3 times)
      isotopePattern <- cbind(c(0.0,1.99705,3.99410),c(1.0,0.96,0.31))
      #first check if there is an isotopic pattern. Otherwise we might assign peaks due to an only partially fitting isotopic pattern (e.g. an Cl-[M+4] peak without an [M+2])
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        isIsotopicPattern <- TRUE
        isotopes <- vector() #we will store the IDs of the isotopes here
        for (isotopeNumber in 2:3) { #search for all isotopes
          zwischenergebnis <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[isotopeNumber,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]*isotopePattern[isotopeNumber,2])) < peakliste[i,3]*inttol))
          if (length(zwischenergebnis) > 0) {
            for (j in 1:length(zwischenergebnis)) {
              if (peakliste[i,4] != peakliste[zwischenergebnis[j],4]) { #if the monoisotopic and the isotope are found at different scan numbers
                #get the IDs of the monoisotopic and the isotope of the spectrum at the first scan number:
                spektrum_monoisotopic <- which(spektrum[,1]==peakliste[i,1])
                spektrum_isotope <- which.min(abs(spektrum[,1]-peakliste[zwischenergebnis[j],1]))
                #get the IDs of the monoisotopic and the isotope of the spectrum at the second scan number:
                spektrum2 <- xcms::getScan(daten,peakliste[zwischenergebnis[j],4])
                spektrum2_monoisotopic <- which.min(abs(spektrum2[,1]-peakliste[i,1]))
                spektrum2_isotope <- which(spektrum2[,1]==peakliste[zwischenergebnis[j],1])
                #only keep the isotope if the relative intensities in both scans are within the tolerance:
                if (abs(spektrum[spektrum_isotope,2]/spektrum[spektrum_monoisotopic,2]-spektrum2[spektrum2_isotope,2]/spektrum2[spektrum2_monoisotopic,2]) > inttol) zwischenergebnis[j] <- 0
              }
            }
            if (any(zwischenergebnis > 0)) { #if there is one peak in "zwischenergebnis" which fits, remember it. 
              isotopes <- c(isotopes, zwischenergebnis[zwischenergebnis > 0])
            } else { #otherwise this is isotopic pattern
              isIsotopicPattern <- FALSE
            }
          }
        }
        if (isIsotopicPattern) {
          Cl3[i] <- i
          Cl3[isotopes] <- i
        }
      } 
      
      #Chlorinated (4 times)
      #inputspectrum <- peakliste[abs(peakliste[,2]-peakliste[i,2])<RT_Tol,,drop = FALSE]
      #inputspectrum <- cbind(inputspectrum[,1],inputspectrum[,3])
      isotopePattern <- cbind(c(0.0,1.99705,3.99410,5.99115),c(0.78,1,0.48,0.1))
      #if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = inputspectrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        #enter number (ID) of monoisotopic peak in row of all three isotopes:
        Cl4[i] <- i
        MainPeak <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]/isotopePattern[1,2])) < peakliste[,3]*inttol))
        Cl4[MainPeak] <- i
        MainPeak <- MainPeak[which.min(abs(peakliste[i,4]-peakliste[MainPeak,4]))]
        MainPeak <- MainPeak[which.min(abs(peakliste[i,1]+isotopePattern[2,1]-peakliste[MainPeak,1]))]
        Cl4[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[3,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[3,2])) < peakliste[MainPeak,3]*inttol))] <- i
        Cl4[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[4,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[4,2])) < peakliste[MainPeak,3]*inttol))] <- i
      }
      
      #Brominated (1 time)
      #inputspectrum <- peakliste[abs(peakliste[,2]-peakliste[i,2])<RT_Tol,,drop = FALSE]
      #inputspectrum <- cbind(inputspectrum[,1],inputspectrum[,3])
      isotopePattern <- cbind(c(0.0,1.99795),c(1.0,0.98))
      #if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = inputspectrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        #enter number (ID) of monoisotopic peak in isotope row:
        Br1[i] <- i
        Br1[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]*isotopePattern[2,2])) < peakliste[i,3]*inttol))] <- i
      }  
      
      #Brominated (2 times)
      #inputspectrum <- peakliste[abs(peakliste[,2]-peakliste[i,2])<RT_Tol,,drop = FALSE]
      #inputspectrum <- cbind(inputspectrum[,1],inputspectrum[,3])
      isotopePattern <- cbind(c(0.0,1.99795,3.9959),c(0.51,1,0.49))
      #if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = inputspectrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        #enter number (ID) of monoisotopic peak in row of both isotopes:
        Br2[i] <- i
        MainPeak <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]/isotopePattern[1,2])) < peakliste[,3]*inttol))
        Br2[MainPeak] <- i
        MainPeak <- MainPeak[which.min(abs(peakliste[i,4]-peakliste[MainPeak,4]))]
        MainPeak <- MainPeak[which.min(abs(peakliste[i,1]+isotopePattern[2,1]-peakliste[MainPeak,1]))]
        Br2[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[3,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[3,2])) < peakliste[MainPeak,3]*inttol))] <- i
      }  
      
      #Brominated (3 times)
      #inputspectrum <- peakliste[abs(peakliste[,2]-peakliste[i,2])<RT_Tol,,drop = FALSE]
      #inputspectrum <- cbind(inputspectrum[,1],inputspectrum[,3])
      isotopePattern <- cbind(c(0.0,1.99795,3.9959,5.99386),c(0.34,1,0.98,0.32))
      #if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = inputspectrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        #enter number (ID) of monoisotopic peak in row of all three isotopes:
        Br3[i] <- i
        MainPeak <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]/isotopePattern[1,2])) < peakliste[,3]*inttol))
        Br3[MainPeak] <- i
        MainPeak <- MainPeak[which.min(abs(peakliste[i,4]-peakliste[MainPeak,4]))]
        MainPeak <- MainPeak[which.min(abs(peakliste[i,1]+isotopePattern[2,1]-peakliste[MainPeak,1]))]
        Br3[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[3,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[3,2])) < peakliste[MainPeak,3]*inttol))] <- i
        Br3[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[4,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[4,2])) < peakliste[MainPeak,3]*inttol))] <- i
      }
      
      #Brominated (4 times)
      #inputspectrum <- peakliste[abs(peakliste[,2]-peakliste[i,2])<RT_Tol,,drop = FALSE]
      #inputspectrum <- cbind(inputspectrum[,1],inputspectrum[,3])
      isotopePattern <- cbind(c(0.0,1.99795,3.9959,5.99386,7.99181),c(0.17,0.68,1,0.65,0.16))
      #if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = inputspectrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
      if (hasSpecificPattern(masse = peakliste[i,1], inputspectrum = spektrum, isotopePattern = isotopePattern, mztol = ppm, inttol = inttol)) {
        #enter number (ID) of monoisotopic peak in row of all three isotopes:
        Br4[i] <- i
        MainPeak <- which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[3,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[i,3]/isotopePattern[1,2])) < peakliste[,3]*inttol))
        Br4[MainPeak] <- i
        MainPeak <- MainPeak[which.min(abs(peakliste[i,4]-peakliste[MainPeak,4]))]
        MainPeak <- MainPeak[which.min(abs(peakliste[i,1]+isotopePattern[3,1]-peakliste[MainPeak,1]))]
        Br4[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[2,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[2,2])) < peakliste[MainPeak,3]*inttol))] <- i
        Br4[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[4,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[4,2])) < peakliste[MainPeak,3]*inttol))] <- i
        Br4[which(abs(peakliste[,1]-(peakliste[i,1]+isotopePattern[5,1])) < (peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol) & (abs(peakliste[,3]-(peakliste[MainPeak,3]*isotopePattern[5,2])) < peakliste[MainPeak,3]*inttol))] <- i
      }
      
      #Na-Addukt
      #zwischenergebnis <- which((abs(peakliste[,1]-(peakliste[i,1]-21.98194)) < peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol))
      #if (length(zwischenergebnis) > 1) {
      # NaAddukt[i] <- zwischenergebnis[which.max(peakliste[zwischenergebnis,3])]
      #}
      #if (length(zwischenergebnis) == 1) {
      # NaAddukt[i] <- zwischenergebnis
      #}
      #if (length(zwischenergebnis) == 0) NaAddukt[i] <- 0
      NaAddukt[which((abs(peakliste[,1]-(peakliste[i,1]+21.98194)) < peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol))] <- i
      
      #K-Addukt
      #zwischenergebnis <- which((abs(peakliste[,1]-(peakliste[i,1]-37.95588)) < peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol))
      #if (length(zwischenergebnis) > 1) {
      #  KAddukt[i] <- zwischenergebnis[which.max(peakliste[zwischenergebnis,3])]
      #}
      #if (length(zwischenergebnis) == 1) {
      #  KAddukt[i] <- zwischenergebnis
      #}
      #if (length(zwischenergebnis) == 0) KAddukt[i] <- 0
      KAddukt[which((abs(peakliste[,1]-(peakliste[i,1]+37.95588)) < peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol))] <- i
      
      #NH4-Addukt
      #zwischenergebnis <- which((abs(peakliste[,1]-(peakliste[i,1]-17.02655)) < peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol))
      #if (length(zwischenergebnis) > 1) {
      # NH4Addukt[i] <- zwischenergebnis[which.max(peakliste[zwischenergebnis,3])]
      #}
      #if (length(zwischenergebnis) == 1) {
      # NH4Addukt[i] <- zwischenergebnis
      #}
      #if (length(zwischenergebnis) == 0) NH4Addukt[i] <- 0
      NH4Addukt[which((abs(peakliste[,1]-(peakliste[i,1]+17.02655)) < peakliste[i,1]/1000000*ppm) & (abs(peakliste[,2]-peakliste[i,2]) < RT_Tol))] <- i
      
      
    }
    
  #}) #end with progress
  
  ergebnis <- cbind(C13,NaAddukt,NH4Addukt,Cl1,Cl2,Cl3,Cl4,Br1,Br2,Br3,Br4,KAddukt,S1,S2)
  colnames(ergebnis) <- c("C13","NaAddukt","NH4Addukt","Cl1","Cl2","Cl3","Cl4","Br1","Br2","Br3","Br4","KAddukt","S1","S2")
  return(ergebnis)
}

# #' @export
# blankCorrection <- function(grouptable, intensityFactor = 10, deleteGrouped = TRUE) {
#   
#   blanks <- NULL
#   #Einbindung der Gruppenspalte erledigt
#   intensityCols <- grep("Int_", colnames(grouptable)) 
#   groupCols <- grep("gruppe_", colnames(grouptable))
#   ergebnis <- matrix(,ncol=ncol(grouptable),nrow=0)
#   
#   for (i in (1:length(headerList))) {
#     if (headerList[[i]]$sampleType == "Blank") 
#       blanks <- c(blanks,i)
#   }
#   #browser()
#   #Erstellt Matrix mit den intensivsten Peaks jeder Probe, jeder Gruppe in den "Unknown" Proben
#   grouptable_Int_max <- matrix(
#     1000000,
#     ncol=(2*length(groupCols[-blanks])),
#     nrow=max(grouptable[,groupCols[1:length(groupCols[-blanks])]]))
#   
#   for (i in 1:length(groupCols[-blanks])){
#     
#     for (j in 1:max(grouptable[,groupCols[i]])){
#       
#       if (0 != length(which(grouptable[,groupCols[i]] == j))) {
#         grouptable_Int_max[j,((2*i)-1)] <- max(
#           c(grouptable[which(grouptable[,groupCols[i]]==j),intensityCols[i]]))
#         }
#       grouptable_Int_max[j,(2*i)] <- j
#       
#     } #for j
#   } #for i
#   
#   
#   #Absuchen der grouptable nach Blanks
#   for (i in 1:nrow(grouptable)) {
#     
#     if (any(grouptable[i,intensityCols] > 
#             intensityFactor*max(grouptable[i,intensityCols[blanks]]))) {
#   
#       ergebnis <- rbind(ergebnis,grouptable[i,])
#       
#     } else if (deleteGrouped) {
#       
#       # Reihe ist ein Vector, welche anzeigt welcher Peak in den Proben 
#       #      1.) im Blank gefunden wurde und
#       #      2.) Gruppenintensivster Peak ist
#       #      mit grouptable[i,[intensityCols[Reihe]]] bekommt man sie angezeigt
#       
#       #Zähler, der durch die Anzahl an Proben(Col in der Grouptable)(ohne blank geht)
#       for (Probe in 1:length(intensityCols[-blanks])){ # Probe <- 1    
#         #Ist der Peak aus Probe auch ein Gruppenchef
#         Reihe <- which(grouptable[i,intensityCols[Probe]] == 
#                          grouptable_Int_max[,(2*Probe)-1])    
#         # Reihe ist die Gruppe des Peaks in Probe Probe
#         
#         if (length(Reihe)>0){
#           
#           t <- i   #Laufvariable t bekommt den Wert von i, da nicht die ganze grouptable durchsucht werden muss, sondern nur der noch kommende Teil
#           for (t in t:nrow(grouptable)){
#             if (grouptable[t, groupCols[Probe]] == Reihe) {
#               grouptable[t, (intensityCols[Probe]-3):(intensityCols[Probe]+2)] <- 0
#             }
#           }
#           
#           
#           if (length(ergebnis) > 0) {
#             for (f in 1:nrow(ergebnis)) {
#               if (ergebnis[f, groupCols[Probe]] == Reihe) 
#                 ergebnis[f,c((intensityCols[Probe]-3):(intensityCols[Probe]+2))] <- 0
#             }
#           }
#           
#           
#         } #if length(Reihe)
#         
#       } #for Probe
#       
#       #Sind die betroffenen Peaks die intensivsten in den betroffenen Gruppen?
#       # der logische Vektor gr_max gibt an, ob der untersuchte Peak der intensivste der Gruppe ist
#       
#       #Was passiert wenn der betroffene Peak das Int-max der Gruppe ist?
#       # Wenn gr_max TRUE ist wird die ganze Gruppe entfernt (peakID,mz,RT,Int,ms2scan, Gruppe <-0)
#       # Wenn gr_max FALSE passiert an dieser Stelle nichts und der Peak wird wie in Christians Orginal einfach nicht in
#       # die Ergebnisliste übertragen
#       
#     } #else
#     
#     
#   } # for
#   
#   if (0 < length(which(rowSums(ergebnis[, -c(1, 2, 3)]) == 0))) 
#     ergebnis <- ergebnis[-c(which(rowSums(ergebnis[,-c(1,2,3)])==0)),]
#   
#   # Eliminieren von ein-Peak-Gruppen  (Gruppe wird 0)
#   for (i in 1:max(ergebnis[,"Gruppe"])) {
#     if (length(which(ergebnis[,"Gruppe"]==i))==1) 
#       ergebnis[which(ergebnis[,"Gruppe"]==i),"Gruppe"] <- 0
#   }
#   
#   # Schließen der entstandenen Zwischenräume
#   V1 <- sort(unique(ergebnis[, "Gruppe"])[-(which(unique(ergebnis[, "Gruppe"]) == 0))])
#   V2 <- seq_along(V1) 
#   for (i in 1:length(V1)) {
#     ergebnis[which(ergebnis[, "Gruppe"] == V1[i]),"Gruppe"] <- V2[i]
#   }
#   
#   
#   return(ergebnis)   
#}

#' @export
checkAlignedIsotopesAndAdducts <- function(grouptable) {
  PeakIDcolumns <- seq(from = 4, to = ncol(grouptable), by = 5)
  toBeKept <- NULL
  MightBeKept <- NULL
  NumberOfTypes <- NULL
  Types <- matrix(,nrow = 0,ncol = length(PeakIDcolumns))
  
  for (i in (1:nrow(grouptable))) {
    
    IsMonoisotopicAndNoAdduct <- NULL
    IsotopeOrAdductType <- NULL
    
    for (IDcolumn in (1:length(PeakIDcolumns))) {
      
      if (grouptable[i,PeakIDcolumns[IDcolumn]] > 0) {
        IsMonoisotopicAndNoAdduct <- c(IsMonoisotopicAndNoAdduct,all(is.element(peaklist[[IDcolumn]][grouptable[i,PeakIDcolumns[IDcolumn]],18:ncol(peaklist[[IDcolumn]])],c(0,grouptable[i,PeakIDcolumns[IDcolumn]]))))
        type <- which((peaklist[[IDcolumn]][grouptable[i,PeakIDcolumns[IDcolumn]],18:ncol(peaklist[[IDcolumn]])] > 0) & (peaklist[[IDcolumn]][grouptable[i,PeakIDcolumns[IDcolumn]],18:ncol(peaklist[[IDcolumn]])] != grouptable[i,PeakIDcolumns[IDcolumn]]))
        #if (length(type) == 0) type <- 0
        IsotopeOrAdductType <- c(IsotopeOrAdductType,type)
      } 
      #else {
      #  IsotopeOrAdductType <- c(IsotopeOrAdductType,NA)
      #}
      
    }
    toBeKept <- c(toBeKept,all(IsMonoisotopicAndNoAdduct))
    MightBeKept <- c(MightBeKept,any(IsMonoisotopicAndNoAdduct))
    #NumberOfTypes <- c(NumberOfTypes, length(unique(IsotopeOrAdductType[IsotopeOrAdductType != 0])))
    NumberOfTypes <- c(NumberOfTypes, length(unique(IsotopeOrAdductType)))
    #Types <- rbind(Types,IsotopeOrAdductType)
  }
  ausgabe <- list()
  ausgabe$toBeKept <- toBeKept
  ausgabe$unsure <- which(MightBeKept != toBeKept)
  ausgabe$NumberOfTypes <- NumberOfTypes
  #ausgabe$Types <- Types
  return(ausgabe)
}

#' @export
correctAndRemoveAlignedIsotopesAndAdducts <- function(grouptable, int_threshold) {
  checkAnalyse <- checkAlignedIsotopesAndAdducts(grouptable)
  PeakIDcolumns <- seq(from = 4, to = ncol(grouptable), by = 5)
  vergleichsmatrix <- data.frame()
  
  for (unsureOnes in checkAnalyse$unsure) {
    #checkAnalyse$Types[unsureOnes,(checkAnalyse$Types[unsureOnes,] > 0) & (is.na(checkAnalyse$Types[unsureOnes,]) == FALSE)]
    
    for (IDcolumn in (1:length(PeakIDcolumns))) {
      if (grouptable[unsureOnes,PeakIDcolumns[IDcolumn]] > 0) {
        vergleichsmatrix <- rbind(vergleichsmatrix,cbind(IDcolumn,peaklist[[IDcolumn]][grouptable[unsureOnes,PeakIDcolumns[IDcolumn]],]))
      }
    }
    problematicColumns <- which(dim(table(vergleichsmatrix[,19:ncol(vergleichsmatrix)])) > 1)+18
    
    
    
    for (i in problematicColumns) {
      if (is.element(i,c(20,21,30))) { #was partly detected as adduct
        
        parentIonIntensity <- NULL
        for (sampleRow in 1:nrow(vergleichsmatrix)) {
          if (vergleichsmatrix[sampleRow,i] > 0) {
            parentIonIntensity <- c(parentIonIntensity,peaklist[[vergleichsmatrix[sampleRow,1]]]$Intensity[vergleichsmatrix[sampleRow,i]])
          } else {
            parentIonIntensity <- c(parentIonIntensity,0)
          }
        }
        intFactor <- mean(parentIonIntensity[parentIonIntensity > 0]/vergleichsmatrix$Intensity[which(parentIonIntensity > 0)])
        parentIonIntensity[parentIonIntensity == 0] <- vergleichsmatrix$Intensity[which(parentIonIntensity == 0)]*intFactor
        if (all(parentIonIntensity[which(vergleichsmatrix[,i] == 0)] < int_threshold)) {} #this one is an adduct, but the "parent" was too small in some cases 
      }
    }
    
  }
}

#' Calculates the optimum mz step based on mass differences. Takes a sample
#' of 100 spectra to make the calculation
#' @export
optimumMzStep <- function(daten, probs) {
  get_diff <- function(scan) {
    msspektrum <- xcms::getScan(daten, scan)
    -quantile(-diff(msspektrum[,1],lag=1),probs = probs)
  }
  diffs <- vapply(sample(seq_along(daten@scanindex), 100), get_diff, numeric(1))
  round(min(diffs), 2)
}


#' @export
findSingleTrend_LocalTemporary <- function(zeile, intensities, treatZeroIntensitiesAs = 0) {
  sn <- 10
  
  ausgabe <- matrix(,ncol=4,nrow=0)
  
  #x <- 1:((ncol(alignmentTable)-3)/4)
  #y <- alignmentTable[zeile,seq(6,ncol(grouped),by=4)]
  
  intensities[intensities == 0] <- treatZeroIntensitiesAs
  
  x <- seq(1,length(intensities), by = 1)
  trend <- cbind(x,intensities)
  deriv1 <- cbind(x[2:length(x)],diff(intensities))
  
  maxima <- deriv1[(deriv1[(2:(nrow(deriv1))-1),2] > 0) & (deriv1[(2:(nrow(deriv1))),2] <= 0),1]
  maxima <- maxima[(maxima > 1) & (maxima < nrow(trend))]
  
  left_end <- maxima
  right_end <- maxima
  
  if (length(maxima) > 0) {
    for (n in (1:length(maxima))) {
      left_end[n] <- max(trend[which(diff(trend[1:maxima[n],2]) <= 0),1])+1
      right_end[n] <- min(trend[(which(diff(trend[(maxima[n]+1):(nrow(trend)),2]) >= 0)+maxima[n]),1])
      if (is.infinite(left_end[n])) left_end[n] <- 1
      if (is.infinite(right_end[n])) right_end[n] <- nrow(trend)
      
      if ((n > 1) & (maxima[n] > 0)){
        if ((right_end[n-1] >= (left_end[n])) & (maxima[n-1] > 0))  {
          if (trend[left_end[n],2] > min(c(trend[maxima[n],2],trend[maxima[n-1],2]))/2) {  
            #in diesem Fall wird das Maximum dem intensiveren Peak zugeordnet:
            maxima[n] <- maxima[which.max(trend[c(maxima[n-1],maxima[n]),2])+n-2]
            #das Noiselevel (links des Peaks) und die Position des linken Rands wird vom vorherigen Maximum ?bernommen:
            left_end[n] <- left_end[n-1]
            #Das vorherige Maximum auf 0 setzen, damit es sp?ter gel?scht wird:
            maxima[n-1] <- 0
            left_end[n-1] <- 0
            right_end[n-1] <- 0
            
          }
        }
      }
      
    }
    
    left_end <- left_end[maxima[]>0]
    right_end <- right_end[maxima[]>0]
    maxima <- maxima[maxima[]>0]
    
  }
  
  noiselevel <- maxima
  
  
  if (length(maxima) > 0) {
    for (n in 1:length(maxima)) {
      noise <- cbind(trend[,1],trend[,2],0)
      for (m in which(((maxima < max(noise[,1])) & (maxima > min(noise[,1]))) | ((right_end < max(noise[,1])) & (right_end > min(noise[,1]))) | ((left_end < max(noise[,1])) & (left_end > min(noise[,1]))))) {
        noise[noise[,1] %in% trend[(left_end[m]:right_end[m]),1],3] <- 1
      }
      noise <- noise[1:left_end[n],,drop=FALSE]
      noise <- noise[noise[,3]==0,,drop=FALSE]
      
      
      if (nrow(noise) == 0) {noise <- trend[left_end[n],,drop=FALSE]}
      
      noiselevel[n] <- mean(noise[,2]) 
      
      if ((trend[maxima[n],2]) < noiselevel[n]*sn) {
        maxima[n] <- 0
        left_end[n] <- 0
        right_end[n] <- 0
        noiselevel[n] <- 0
      }
      
    }
    
    left_end <- left_end[maxima[]>0]
    right_end <- right_end[maxima[]>0]
    noiselevel <- noiselevel[maxima[]>0]
    maxima <- maxima[maxima[]>0]
  }
  if (length(maxima) > 0) {
    for (n in 1:length(maxima)) {
      ausgabe <- rbind(ausgabe,c(zeile,left_end[n],maxima[n],right_end[n]))
    }
  }
  return(ausgabe)
}

#' @export
findTrend <- function(samplingSite = "", samplingPosition = NULL, samplingDate = NULL, treatZeroIntensitiesAs = 0) {
  
  #samplingSite <- "Rhein"
  if (is.null(samplingPosition)) {
    samplesAtThisSite <- NULL
    samplingPositions <- NULL
    for (i in 1:((ncol(grouped)-3)/4)) {
      if ((headerList[[i]]$sampleType == "Unknown") & (headerList[[i]]$samplingSite == samplingSite) & (headerList[[i]]$samplingDate == samplingDate)) {
        samplesAtThisSite <- c(samplesAtThisSite,i)
        samplingPositions <- c(samplingPositions,headerList[[i]]$samplingPosition)
      }
    }
    
    
    #order by samplingPosition:
    samplesAtThisSite <- samplesAtThisSite[order(samplingPositions)]
    samplingPositions <- sort(samplingPositions)
    
    
    trendlist <- NULL
    #trendlist_list <- lapply(X = seq(1,nrow(grouped), by = 1), FUN = function(X)
    #{findSingleTrend(zeile = X)})
    trendlist_list <- lapply(X = seq(1,nrow(grouped), by = 1), FUN = function(X)
    {findSingleTrend_LocalTemporary(zeile = X, intensities = grouped[X,samplesAtThisSite*4+2], treatZeroIntensitiesAs = treatZeroIntensitiesAs)})
    trendlist <- rbind(trendlist, do.call(rbind, trendlist_list))
    return(trendlist)
  }
  
  if (is.null(samplingDate)) {
    samplesAtThisSite <- NULL
    samplingDates <- NULL
    for (i in 1:((ncol(grouped)-3)/4)) {
      if ((headerList[[i]]$sampleType == "Unknown") & (headerList[[i]]$samplingSite == samplingSite) & (headerList[[i]]$samplingPosition == samplingPosition)) {
        samplesAtThisSite <- c(samplesAtThisSite,i)
        samplingDates <- c(samplingDates, headerList[[i]]$samplingDate)
      }
    }
    
    
    #order by samplingPosition:
    samplesAtThisSite <- samplesAtThisSite[order(samplingDates)]
    samplingDates <- sort(samplingDates)
    
    
    trendlist <- NULL
    #trendlist_list <- lapply(X = seq(1,nrow(grouped), by = 1), FUN = function(X)
    #{findSingleTrend(zeile = X)})
    trendlist_list <- lapply(X = seq(1,nrow(grouped), by = 1), FUN = function(X)
    {findSingleTrend_LocalTemporary(zeile = X, intensities = grouped[X,samplesAtThisSite*4+2], treatZeroIntensitiesAs = treatZeroIntensitiesAs)})
    trendlist <- rbind(trendlist, do.call(rbind, trendlist_list))
    return(trendlist)
  }  
}

#' Peak picking function simplified for scripting
#' 
#' The original peak picking function can only be used together with the app.
#' This function can be used in a scripting environment. Results are the same.
#' This is not parallelized, since it assumed to calling function is already parallelized.
#'
#' @param daten 
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
#' @return
#' @export
#'
#' @examples
FindPeaks_BfG_scripting <- function(daten, 
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
  
  RT_Tol_scan <- 3
  mz_Tol_ppm <- 20
  
  rt_min_scan <- min(which(daten@scantime > rt_min))
  rt_max_scan <- max(which(daten@scantime < rt_max))
  if (is.infinite(rt_min_scan)) rt_min_scan <- 1
  if (is.infinite(rt_max_scan)) rt_max_scan <- length(daten@scantime)
  
  peaklist <- foreach(X = seq(mz_min, mz_max, by = mz_step*0.5), 
                      .packages = c("xcms", "ntsworkflow"), .combine = rbind) %do% {
                        peakpicking_BfG_cpp(i = X, rawData = daten, mz_step = mz_step,
                                            rt_min_scan = rt_min_scan,
                                            rt_max_scan = rt_max_scan, 
                                            sn=sn,
                                            int_threshold=int_threshold,
                                            NoiseScans=peak_NoiseScans,
                                            precursormzTol = precursormzTol,
                                            peakwidth_min=peakwidth_min,
                                            peakwidth_max=peakwidth_max,
                                            maxPeaksPerSignal=maxPeaksPerSignal
                        )
                      }
  
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
  if (is.null(peaklist)) 
    return(matrix(NA, 0, 18, dimnames = list(NULL, cn)))
  
  if (nrow(peaklist) >= 1) {
    
    peaklist <- peaklist[order(peaklist[,1]), , drop = F]
    
    #delete shoulder peaks which are inside the RT range of another peak:
    for (i in 1:nrow(peaklist)) {
      doppelte <- which((abs(peaklist[,1]-peaklist[i,1]) < peaklist[i,1]/1000000*mz_Tol_ppm) & 
                          (peaklist[,2] >= peaklist[i,5]) & (peaklist[,2] <= peaklist[i,6])) 
      doppelte <- doppelte[doppelte != i]
      peaklist[doppelte, 1] <- 0  
    }
    
    peaklist <- peaklist[peaklist[,1] > 0, , drop = F]
    
  }
  
  # calculation of SN
  if (nrow(peaklist) > 0)
    peaklist <- cbind(peaklist,peaklist[,3]/peaklist[,10],0)
  
  colnames(peaklist) <- cn
  
  # remove rows which are below intensity threshold and not within RT window
  peaklist <- peaklist[peaklist[, "Intensity"] >= int_threshold, , drop = F]
  peaklist <- peaklist[peaklist[, "RT"] >= rt_min, , drop = F]
  peaklist <- peaklist[peaklist[, "RT"] <= rt_max, , drop = F]
  peaklist
}
