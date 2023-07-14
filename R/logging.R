# Copyright 2016-2023 Bundesanstalt für Gewässerkunde
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


#' create_log_file
#'
#' @param no parameters needed
#'
#' @return log_file with used settings for peakpicking, alignment, blankcorrectio,normalization, filtering, annotation
#' @export
#'
#' @examples
create_log_file <- function () {
  
  if(!exists("log_file")){ 
    log_file <- list()
    log_file[6][[1]] <- list()}
  
  # Peakpicking Sampletype = "Unknown"
  if(any(sampleList$sampleType=="Unknown") & 
     any(ls(envir=parent.frame())=="adjust_tolerance")) {
      sample <- which(sampleList$sampleType=="Unknown")
      Settings_PP_sample <- list()
      Settings_PP_sample$Peakpicking_Settings <- "Unknown"
      Settings_PP_sample$Mass_Range <- peakPickSettings[[sample[1]]]$massrange
      Settings_PP_sample$mz_Step <- peakPickSettings[[sample[1]]]$mz_step
      Settings_PP_sample$RT_Range_min <- peakPickSettings[[sample[1]]]$rtrange/60
      Settings_PP_sample$Peakwidth_s <- peakPickSettings[[sample[1]]]$peakwidth
      Settings_PP_sample$Noise_Scans <- peakPickSettings[[sample[1]]]$NoiseScans
      Settings_PP_sample$SN <- peakPickSettings[[sample[1]]]$sn
      Settings_PP_sample$Int_Threshold <- peakPickSettings[[sample[1]]]$int_threshold
      Settings_PP_sample$MaxNumPeaks <- peakPickSettings[[sample[1]]]$maxNumPeaks
      
      adjust_tolerance <- get ("adjust_tolerance", parent.frame())
      if (adjust_tolerance) {Settings_PP_sample$Dynamic_Tolerance_componentization <- TRUE
      } else {Settings_PP_sample$Dynamic_Tolerance_componentization <- FALSE
              Settings_PP_sample$mz_margin_Componentization_ppm <- peakPickSettings[[sample[1]]]$ppm
              Settings_PP_sample$RT_Tolerance_Componentization <- peakPickSettings[[sample[1]]]$RT_Tol}
      
      log_file[[1]] <- Settings_PP_sample
  }
  
  # Peakpicking Sampletype = "Blank"
  if(any(sampleList$sampleType=="Blank") &
     any(ls(envir=parent.frame())=="adjust_tolerance")){
      blanksample <- which(sampleList$sampleType=="Blank")
      Settings_PP_blank <- list()
      Settings_PP_blank$Peakpicking_Settings <- "Blank"
      Settings_PP_blank$Mass_Range <- peakPickSettings[[blanksample[1]]]$massrange
      Settings_PP_blank$mz_Step <- peakPickSettings[[blanksample[1]]]$mz_step
      Settings_PP_blank$RT_Range_min <- peakPickSettings[[blanksample[1]]]$rtrange
      Settings_PP_blank$Peakwidth_s <- peakPickSettings[[blanksample[1]]]$peakwidth
      Settings_PP_blank$NoiseScans <- peakPickSettings[[blanksample[1]]]$NoiseScans
      Settings_PP_blank$SN <- peakPickSettings[[blanksample[1]]]$sn
      Settings_PP_blank$Int_Threshold <- peakPickSettings[[blanksample[1]]]$int_threshold
      Settings_PP_blank$MaxNumPeaks <- peakPickSettings[[blanksample[1]]]$maxNumPeaks
    
      adjust_tolerance <- get ("adjust_tolerance", parent.frame())
      if (adjust_tolerance) {Settings_PP_blank$Dynamic_Tolerance_componentization <- TRUE
       } else {Settings_PP_blank$Dynamic_Tolerance_componentization <- FALSE
               Settings_PP_blank$mz_margin_Componentization_ppm <- peakPickSettings[[blanksample[1]]]$ppm
               Settings_PP_blank$RT_Tolerance_Componentization <- peakPickSettings[[blanksample[1]]]$RT_Tol}
    
    log_file[[2]] <- Settings_PP_blank
  }
  
  # Alignment
  if(any(ls(envir=parent.frame())=="ppm_dev") &
     any(ls(envir=parent.frame())=="DeltaRT")){
      ppm_dev <- get ("ppm_dev", parent.frame())
      DeltaRT <- get ("DeltaRT", parent.frame())
      log_file[[3]] <- list(Alignment_Settings="",
                            Mass_Deviation_ppm=ppm_dev,
                            RT_Deviation_s=DeltaRT)
    # A new alignment table is created during a new alignment, therefore ...
    if(length(log_file)>3){ log_file[[4]] <- list()
                            log_file[[5]] <- list()
                            log_file[[6]] <- list()}
    if(length(log_file)>6) log_file[[7]][[1]] <- "Due to new data processing, the annotation may no longer be up to date and should be repeated."
    
    }
  
  # Blankcorrection
  if(any(ls(envir=parent.frame())=="intensityFactor") &
     any(ls(envir=parent.frame())=="deleteGrouped")){
      intensityFactor <- get ("intensityFactor", parent.frame())
      deleteGrouped <- get ("deleteGrouped", parent.frame())
      grouped_before <- get ("nrow_grouped_before_blankcorrection", parent.frame())
      number_of_removed_features <- grouped_before - nrow(grouped)
      log_file[[4]] <- list(Blankcorrection_Settings="",
                            IntensityFactor=intensityFactor,
                            DeleteGrouped=deleteGrouped,
                            NumRemoveFeat=number_of_removed_features)
    # When blank correction is repeated, any previously applied filter steps are overwritten, therefore...
      if(length(log_file)>4){log_file[[5]] <- list()
                             log_file[[6]] <- list()}
      if(length(log_file)>6) log_file[[7]][[1]] <- "Due to new data processing, the annotation may no longer be up to date and should be repeated."
    }
  
  # Normalize
  if(any(ls(envir=parent.frame())=="intCols_normalize")){
    mean_mz_chosen_feature <- get ("mean_mz_chosen_feature", parent.frame())
    mean_RT_chosen_feature <- get ("mean_RT_chosen_feature", parent.frame())/60
    row_number <- get ("row_number", parent.frame())
    log_file[[5]] <- list(normalization_settings="",
                          Mean_mz_chosen_Feature=mean_mz_chosen_feature,
                          Mean_RT_chosen_Feature_min=mean_RT_chosen_feature,
                          RowNumber_chosen_Feature=row_number)
  } 
  
  # Filter
  # Remove rows with fewer than X rows
  if(any(ls(envir=parent.frame())=="minDetections")){
    minDetections <- get ("minDetections", parent.frame())
    grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
    number_of_removed_features <- grouped_before - nrow(grouped)
    log_file[[6]][[1]] <- list(Filter_Settings="Remove rows with fewer than X detections",
                               MinSampleDetect=minDetections,
                               NumRemoveFeat=number_of_removed_features)
    if(length(log_file)>6) log_file[[7]][[1]] <- "Due to new data processing, the annotation may no longer be up to date and should be repeated."
    }    
  
  # Remove features which are not the highest feature in their component in files x:y
  if(any(ls(envir=parent.frame())=="leader_first_sample") &
    any(ls(envir=parent.frame())=="leader_last_sample")){
      leader_first_sample <- get ("leader_first_sample", parent.frame())
      leader_last_sample <- get ("leader_last_sample", parent.frame())
      grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
      number_of_removed_features <- grouped_before - nrow(grouped)
      log_file[[6]][[2]] <- list(Filter_Settings="Remove features which are not the highest feature in their component",
                                 First_sample=leader_first_sample,
                                 Last_sample=leader_last_sample,
                                 NumRemoveFeat=number_of_removed_features)
      if(length(log_file)>6) log_file[[7]][[1]] <- "Due to new data processing, the annotation may no longer be up to date and should be repeated."
  }
  
  # Remove features which are not found in replicate injections
  if(any(ls(envir=parent.frame())=="rep_first_sample") &
     any(ls(envir=parent.frame())=="rep_last_sample") &
     any(ls(envir=parent.frame())=="rep_number_replicates") &
     any(ls(envir=parent.frame())=="in_at_least_samples")){
      rep_first_sample <- get ("rep_first_sample", parent.frame())
      rep_last_sample <- get ("rep_last_sample", parent.frame())
      number_replicates <- get ("rep_number_replicates", parent.frame())
      in_at_least_samples <- get ("in_at_least_samples", parent.frame())
      grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
      number_of_removed_features <- grouped_before - nrow(grouped)
      log_file[[6]][[3]] <- list(Filter_Settings="Remove features which are not found in replicate injections",
                                Replicate_first_sample=rep_first_sample,
                                Replicate_last_sample=rep_last_sample,
                                Number_Replicates=number_replicates,
                                In_at_least_samples=in_at_least_samples,
                                NumRemoveFeat=number_of_removed_features)
      if(length(log_file)>6) log_file[[7]][[1]] <- "Due to new data processing, the annotation may no longer be up to date and should be repeated."
      }
  
  # Get average intensity of features in replicates
  if(any(ls(envir=parent.frame())=="repAve_first_sample") &
     any(ls(envir=parent.frame())=="repAve_last_sample") &
     any(ls(envir=parent.frame())=="repAve_number_replicates")){
      repAve_first_sample <- get ("repAve_first_sample", parent.frame())
      repAve_last_sample <- get ("repAve_last_sample", parent.frame())
      number_replicates <- get ("repAve_number_replicates", parent.frame())
      log_file[[6]][[4]] <- list(Filter_Settings="Get average intensity of features in replicates",
                                RepAve_first_sample=repAve_first_sample,
                                RepAve_last_sample=repAve_last_sample,
                                Number_Replicates=number_replicates)
      if(length(log_file)>6) log_file[[7]][[1]] <- "Due to new data processing, the annotation may no longer be up to date and should be repeated."
      }
  
  # Remove rows where features are found in fewer than X consecutive files
  #browser()
  if(any(ls(envir=parent.frame())=="consec_first_sample") &
     any(ls(envir=parent.frame())=="consec_last_sample") &
     any(ls(envir=parent.frame())=="number_consecutives")){
      consec_first_sample <- get ("consec_first_sample", parent.frame())
      consec_last_sample <- get ("consec_last_sample", parent.frame())
      number_consecutives <- get ("number_consecutives", parent.frame())
      grouped_before <- get ("nrow_grouped_before_filter", parent.frame())
      number_of_removed_features <- grouped_before - nrow(grouped)
      log_file[[6]][[5]] <- list(Filter_Settings="Remove rows where features are found in fewer than X consecutive files",
                                Consec_first_sample=consec_first_sample,
                                Consec_last_sample=consec_last_sample,
                                Number_consecutives=number_consecutives,
                                NumRemoveFeat=number_of_removed_features)
      if(length(log_file)>6) log_file[[7]][[1]] <- "Due to new data processing, the annotation may no longer be up to date and should be repeated."
      }
  
  # Annotation Settings
  if(any(ls(envir=parent.frame())=="threshold_score") &
     any(ls(envir=parent.frame())=="mztolu") &
     any(ls(envir=parent.frame())=="rttol") &
     any(ls(envir=parent.frame())=="polarity") &
     any(ls(envir=parent.frame())=="CE") &
     any(ls(envir=parent.frame())=="CES") &
     any(ls(envir=parent.frame())=="mztolu_ms2") &
     any(ls(envir=parent.frame())=="rtoffset")){
      threshold_score <- get ("threshold_score", parent.frame())
      mztolu <- get ("mztolu", parent.frame())
      rttol <- get ("rttol", parent.frame())
      polarity <- get ("polarity", parent.frame())
      CE <- get ("CE", parent.frame())
      CES <- get ("CES", parent.frame())
      mztolu_ms2 <- get ("mztolu_ms2", parent.frame())
      rtoffset <- get ("rtoffset", parent.frame())
      groupedWithAnnotation <- get ("groupedWithAnnotation", parent.frame())
      log_file[[7]] <- list(Annotation_Settings="",
                            Threshold_Score=threshold_score,
                            mz_Tolerance_mDa=mztolu*1000,
                            RT_Tolerance_min=rttol,
                            Polarity=polarity,
                            CE=CE,
                            CES=CES,
                            mz_Tolerance_MS2_mDa=mztolu_ms2*1000,
                            RT_Offset_min=rtoffset,
                            NumAnnoCompounds=length(unique(groupedWithAnnotation$name))-1)
  }
  
  return(log_file)
  
}


