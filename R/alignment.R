
#' Add additional peak and sample ID information to peaklists
#' 
#' @description takes a list of peaklists and adds sample ID and peak ID information
#'
#' @param peaklist `list` of peak tables
#' @param datenList `list` of `xcms::xcmsRAW` objects
#' 
#' @details This needs to be run before alignment 
#' @returns New peak table `list` in the same form as before but with additional columns
#' sample_id and peak_id_all. If these were already present, they are over-written
#' @export
new_peaklist <- function(peaklist, datenList) {
  peaklisttemp <- Map(function(pl, sid) transform(pl, sample_id = sid, peak_id_all = NA), 
                      peaklist, seq_along(datenList))
  new_pl <- do.call("rbind", peaklisttemp)
  new_pl$peak_id_all <- seq_len(nrow(new_pl))
  split(new_pl, new_pl$sample_id)
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


#' Alignment of peaks from difference samples
#' 
#' @description This a wrapper for the C++ function alignmentBfGC. It takes a list of peaklists
#' and produces an alignment table. 
#' 
#' @param peaklists list of peaklists
#' @param mz_dev m/z tolerance
#' @param DeltaRT rt tolerance in s
#' @param mz_dev_unit Units for m/z tolerance (ppm or mDa)
#' 
#' @returns Matrix of aligned peaks, one aligned feature per row
#' 
#' @export
alignment_BfG_cpp <- function(peaklists, mz_dev, DeltaRT, mz_dev_unit){
  
  # Convert peaklist to matrix for C++, only take the first 3 columns,
  # mz, rt, intensity
  peaklistC <- list()
  for (i in 1:length(peaklists)) {
    peaklistC[[i]] <- as.matrix(peaklists[[i]][,1:3])
  }
  
  # set mz_dev_unit to integer 1 or 2, 1 for ppm, 2 for mDa
  stopifnot(mz_dev_unit %in% c("ppm", "mDa"))
  mz_dev_unit <- switch(mz_dev_unit, ppm = 1, mDa = 2)
  
  # alignedtable has a column for each peaklist the number in the cell
  # represents the rownumber in that peaklist of the feature
  # so a one in column 1 is the first feature in the first peaklist.
  # A 0 means that no matching feature was found.
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
    
    # Goes throgh each row of alignedtable
    for (j in 1:(nrow(alignedtable))) {
      if (alignedtable[j,i] > 0) {
        grouped2[j,i*6-1] <- peaklists[[i]]$peak_id_all[alignedtable[j,i]]
        grouped2[j,i*6] <- peaklists[[i]]$mz[alignedtable[j,i]]
        grouped2[j,i*6+1] <- peaklists[[i]]$RT[alignedtable[j,i]]
        grouped2[j,i*6+2] <- peaklists[[i]]$Intensity[alignedtable[j,i]]
        grouped2[j,i*6+3] <- peaklists[[i]]$MS2scan[alignedtable[j,i]]
        grouped2[j,i*6+4] <- peaklists[[i]]$Gruppe[alignedtable[j,i]]
      }
      # Compute average m/z and RT for this row
      if (i == (ncol(alignedtable))) {
        grouped2[j,1] <- mean(grouped2[j,which(alignedtable[j,] > 0)*6])
        grouped2[j,2] <- mean(grouped2[j,which(alignedtable[j,] > 0)*6+1])
      }
    }
  }
  
  colnames(grouped2) <- spaltennamen
  
  grouped2
}

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow
