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




# Contains functions for filtering and reducing the alignment table



#' Keep only features which are the group leaders
#' 
#' This will go through the alignment table and check every feature if it is a group leader. Features
#' which are marked as group leader in ANY sample will kept. In other words, if the componentization is inconsistent
#' and a feature is sometimes marked as groupleader and sometimes not, the feature will be kept.
#'
#' @param alignment alignment table
#' @param pll peaklistlist
#' @param samples sample IDs of the samples which will be looked at, all other samples will be ignored 
#'
#' @return an alignment table with reduced number of rows
#' @export
#' @import foreach
onlyGroupLeaders <- function(alignment, pll, samples) {
  stopifnot(is.numeric(samples))
  stopifnot(is.list(pll))
  # loop through every row of alignment
  keep <- foreach(r = seq_len(nrow(alignment)), .combine = "c") %do% {
    # get a logical vector of samples, if the feature is group leader or not in each sample
    leader <- vapply(samples, function(sample) {
      ## get T/F for each sample by accessing id in peaklistlist
      currentPID <- alignment[r, paste0("PeakID_", sample)]
      if (currentPID == 0) {
        NA
      } else {
        subset(pll[[sample]], peak_id_all == currentPID, groupleader, drop = T)
      }
    }, logical(1))
    ## if any are groupleaders, keep that row
    any(leader, na.rm = TRUE)
  }
  alignment[keep, , drop = FALSE]
}

#' Get mean intensities for replicate injections
#' 
#' @description
#' This function will take an alignment table and average the intensities of replicate injections.
#' The first sample in each replicate will remain in the table with average int. Zeros will be
#' ignored (unless all intensities are zero).
#' 
#' @param alignment alignment table
#' @param samples numeric sample IDs to merge, length must be divisable by repNum
#' @param reps number of replicates
#'
#' @returns An alignment table with a reduced number of samples but the same number of rows.
#' @export
#'
averageReplicates <- function(alignment, samples, reps) {
  stopifnot(length(samples) %% reps == 0)
  # make list containing grouped samp ID which need to be merged
  sets <- split(samples, rep(1:(length(samples) / reps), each = reps))
  intColsSets <- lapply(sets, function(set) which(colnames(alignment) %in% paste0("Int_", set)))
  # loop through every row of alignment
  for (r in seq_len(nrow(alignment))) {
    # loop through the sets of replicates
    for (setNum in seq_along(sets)) {
      # get the intensity columns
      intCols <- intColsSets[[setNum]]
      allInt <- alignment[r, intCols]
      # if all intensities are 0, no change needed
      if (all(allInt == 0)) {
        next
      } else { # get average intensity (ignore 0s)
        allInt2 <- allInt[allInt != 0]
        meanInt <- mean(allInt2, na.rm = TRUE) 
      }
      stopifnot(!is.na(meanInt))
      stopifnot(meanInt != 0)
      # set average intensity to first sample in set
      alignment[r, intCols[1]] <- signif(meanInt)
    }
  }
  # keep only columns of first sample in each set, and all columns not in set
  # get all samples
  ic <- colnames(alignment)[grep("^Int_", colnames(alignment))]
  allSamp <- as.numeric(stringr::str_match(ic, "Int_(\\d+)$")[,2])
  otherSamp <- allSamp[Negate(is.element)(allSamp, samples)]
  firstOfSet <- append(otherSamp, sapply(sets, function(x) x[1]))
  keepCols <- c("mean_mz", "mean_RT", "MS2Fit", "Gruppe", 
                paste0("PeakID_", firstOfSet), 
                paste0("mz_", firstOfSet), 
                paste0("RT_", firstOfSet), 
                paste0("Int_", firstOfSet), 
                paste0("ms2scan_", firstOfSet), 
                paste0("gruppe_", firstOfSet),
                "alignmentID"
  )
  alignment[, which(colnames(alignment) %in% keepCols), drop = FALSE]
}


#' Only retain peaks found in all replicates
#' 
#' This function will check, for a given number of replicate injections, that a feature was found
#' in all replicates. Features will be set to zero if they are not found in replicates.
#' Aligned features not found in any samples (everything zero) will be removed. 
#' 
#' @param alignment alignment table
#' @param samples numeric vector sample IDs to merge, length must be divisable by reps
#' @param reps numer of replicates
#' @param least Peak found in at least this many replicates, must be <= reps
#'
#' @return an alignment table with the same number of columns but a reduced number of rows.
#' @export
keepReps <- function(alignment, samples, reps, least) {
  stopifnot(length(samples) %% reps == 0)
  stopifnot(least <= reps)
  stopifnot(ncol(alignment) >= length(samples) * 6 + 4) 
  
  # make list containing grouped samp ID which need to be merged
  sets <- split(samples, rep(1:(length(samples) / reps), each = reps))
  
  # determine the columns for each set once
  intColsSets <- lapply(sets, function(set) which(colnames(alignment) %in% paste0("Int_", set)))
  mzColsSets <- lapply(sets, function(set) which(colnames(alignment) %in% paste0("mz_", set)))
  rtColsSets <- lapply(sets, function(set) which(colnames(alignment) %in% paste0("RT_", set)))
  ms2ColsSets <- lapply(sets, function(set) which(colnames(alignment) %in% paste0("ms2scan_", set)))
  gruppeColsSets <- lapply(sets, function(set) which(colnames(alignment) %in% paste0("gruppe_", set)))
  pidColsSets <- lapply(sets, function(set) which(colnames(alignment) %in% paste0("PeakID_", set)))
  # loop through each row of alignment
  for (r in seq_len(nrow(alignment))) {
    # loop through each element of sets
    for (setNum in seq_along(sets)) {
      # get the intensity columns for this set
      
      intCols <- intColsSets[[setNum]]
      mzCols <- mzColsSets[[setNum]]
      rtCols <- rtColsSets[[setNum]]
      ms2Cols <- ms2ColsSets[[setNum]]
      gruppeCols <- gruppeColsSets[[setNum]]
      pidCols <- pidColsSets[[setNum]]
      
      # check which features are not zero
      found <- alignment[r, intCols] != 0
      # if all are zero already, skip to next
      if (all(!found))
        next
      # if number non-zero < least, then set all inten to zero in alignment 
      if (sum(found) < least) 
        alignment[r, c(intCols,mzCols,rtCols,ms2Cols,gruppeCols,pidCols)] <- 0
    }
  }
  # check if any rows of alignment are all zero intensity, if so remove.
  # get only intensity columns
  intensities <- alignment[, grep("^Int_", colnames(alignment))]
  # find which rows have any non-zero intensities
  keep <- apply(intensities, 1, function(x) any(x != 0))
  alignment[keep, , drop = FALSE]
}

#' Remove features found in less than x samples
#'
#' @param alignment alignment table
#' @param minimum number of samples in which a feature is found must be less than
#' total number of samples
#'
#' @return alignment table with rows removed
#' @export 
removeRare <- function(alignment, minimum) {
  # TODO this function should ignore blanks!
  stopifnot(minimum <= (ncol(alignment) - 4) / 6)  # check
  intensities <- alignment[, grep("^Int_", colnames(alignment))]
  keep <- apply(intensities, 1, function(x) sum(x != 0) >= minimum)
  alignment[keep, , drop = FALSE]
}

#' Remove features also found in the blanks
#' 
#' @return an alignment table with rows removed
#' @export
#' @import foreach
blankCorrection <- function(alignment, samplels, intensityFactor = 10, deleteGrouped = TRUE) {
  
  intensityCols <- grep("^Int_", colnames(alignment)) 
  groupCols <- grep("^gruppe_", colnames(alignment))
  
  # get columns in alignment table of blanks
  blanks <- samplels[samplels$sampleType == "Blank", "ID"]
  blankIntColsNames <- paste0("Int_", blanks)
  blankIntCols <- which(colnames(alignment) %in% blankIntColsNames)
  noBlankIntCols <- setdiff(intensityCols,blankIntCols)
  
  keep <- foreach(i = seq_len(nrow(alignment)), .combine = "c") %do% {
    any(alignment[i, noBlankIntCols] > intensityFactor * max(alignment[i, blankIntCols]))
  }
  
  ergebnis <- alignment[keep, , drop=F]
  
  
  if (deleteGrouped) {
    
    # When the Groupleader due to blankcorrection is deleted then the whole group
    # should be deleted in this sample.
    
    for (i in 1:length(groupCols[-blanks])){  # Loop over samples
      
      for (j in 1:max(alignment[,groupCols[i]])){  # Loop over groups
        
        # If the group exists the most intense feature in the group is determined
        
        if (0 != length(which(alignment[,groupCols[i]] == j))) {
          Int_max <- max(c(alignment[which(alignment[,groupCols[i]]==j),intensityCols[i]]))
          group_Int_max <- j
          
          # Questions:  A is the groupleader not in the result and
          #           B are other features of this group in the result.
          # If both are true the values of the columns from PeakID to gruppe 
          # for the other features of this group in this sample are set to 0
          
          if (!any(ergebnis[,sprintf("Int_%i",i)]==Int_max) &&        # A
              any(ergebnis[,sprintf("gruppe_%i",i)]==group_Int_max)){ # B
            
            ergebnis[ergebnis[,sprintf("gruppe_%i",i)]==group_Int_max, # rows
                     grep(sprintf("PeakID_%i",i), colnames(ergebnis)): # cols from PeakID_i
                       grep(sprintf("gruppe_%i",i), colnames(ergebnis))] <- 0}  # to gruppe_i
          
        } #if
      } #for j
    } #for i
    
    # Delete rows with all zero. Is this necessary? Yes!
    if (0 < length(which(rowSums(ergebnis[, -c(grep("^mean_mz$", colnames(ergebnis)),
                                               grep("^mean_RT$", colnames(ergebnis)),
                                               grep("^MS2Fit$", colnames(ergebnis)),
                                               grep("^Gruppe$", colnames(ergebnis)),
                                               grep("^alignmentID$", colnames(ergebnis)))]) == 0))){ 
      ergebnis <- ergebnis[-c(which(rowSums(ergebnis[,-c(grep("^mean_mz$", colnames(ergebnis)),
                                                         grep("^mean_RT$", colnames(ergebnis)),
                                                         grep("^MS2Fit$", colnames(ergebnis)),
                                                         grep("^Gruppe$", colnames(ergebnis)),
                                                         grep("^alignmentID$", colnames(ergebnis)))])==0)),]}
    
  } #if (deleteGrouped) 
  
  # Removal of one peak groups (Group is set to 0)
  for (i in 1:max(ergebnis[,"Gruppe"])) {
    if (length(which(ergebnis[,"Gruppe"]==i))==1) 
      ergebnis[which(ergebnis[,"Gruppe"]==i),"Gruppe"] <- 0
  }
  
  # Closing the gaps
  V1 <- sort(unique(ergebnis[, "Gruppe"])[-(which(unique(ergebnis[, "Gruppe"]) == 0))])
  V2 <- seq_along(V1) 
  for (i in 1:length(V1)) {
    ergebnis[which(ergebnis[, "Gruppe"] == V1[i]),"Gruppe"] <- V2[i]
  }
  
  ergebnis   
}


#' Keep only rows found in at least X consecutive samples
#'
#' @param alignment Alignment Table
#' @param samples samples to consider (everything else will be ignored), by ID
#' @param consecutive number of consecutive samples a peak should be present in
#'
#' @returns Alignment table with possibly some rows removed
#' @export
#'
#' @import foreach
keepConsecutive <- function(alignment, samples, consecutive) {
  stopifnot(is.numeric(samples))
  searchString <- paste0(rep("TRUE", consecutive), collapse = "")
  # loop through each row of alignment
  keep <- foreach(ar = seq_len(nrow(alignment)), .combine = "c") %do% {
    # get intensity column names
    intNames <- paste0("Int_", samples)
    # get intensities for this row
    ints <- alignment[ar, intNames, drop = T]
    peakPresent <- ints != 0
    # is peak present in x consecutive samples?
    grepl(searchString, paste0(peakPresent, collapse = ""))
  }
  alignment[keep, , drop = FALSE]
}


