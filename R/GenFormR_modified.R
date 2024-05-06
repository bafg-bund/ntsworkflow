# Sending data to GenForm and getting results
# E. Schymanski, 9/6/2015
# built into function form 26/11/2015
# modified Jan&Feb 2017.
# modified for upload to GitHub April 2018.
# 
# K. Jewell downloaded from GitHub and modified (KJ) August 2020 for use on CentOS

#' Run GenForm with various Commands
#'
#' @description This is a wrapper function for GenForm, requiring MS, MS/MS, ion settings
#' and directories. Original function Copyright Emma Schymanski (schymane/GenFormR). 
#' Modified to work on CentOS.
#'
#' @usage RunGenForm(ms_file, msms_file, ion, GenFormDir, ResultsDir,ppm=5, acc=15,
#' elements="CHNOPS",ff="", dbe=TRUE, oei=TRUE, mz=0, exist=TRUE, an_loss=TRUE,an_loss_limit=500,
#' sort=FALSE, cleanMSMS=FALSE)
#'
#' @param ms_file File name of a plain text, two column file (e.g. \code{MS.txt}) containing
#' space-separated mz and intensity values of the MS1 spectrum (full scan). The most intense
#' peak is taken as "m" unless this is defined with \code{"mz"}.
#' @param msms_file A plain text, two column file (e.g. \code{MSMS.txt}) containing
#' space-separated mz and intensity values of the tandem (MS2) spectrum, i.e. the fragments.
#' While GenForm can be run without MSMS, this function requires its use at this stage.
#' @param ion Ion setting to use. Options: \code{c("-e","+e","-H","+H","+Na")}.
#' @param GenFormDir Directory containing \code{GenForm.exe}.
#' @param ResultsDir Directory to save the GenForm results. Existing files will be overwritten.
#' @param ppm Accuracy setting for MS1 (default \code{5} ppm).
#' @param acc Accuracy setting for MS2 (default \code{15} ppm).
#' @param elements Selection of elements to consider. Default \code{"CHNOPS"}, implemented
#' are \code{"CHBrClFINOPSSi"}.
#' @param ff A "fuzzy formula" input to offer more control beyond \code{elements}. This should be
#' in format \code{C0-6H0-20O1-3}.
#' @param dbe Default \code{TRUE} writes double bond equivalents to output; \code{FALSE} suppresses
#' this output.
#' @param oei Default \code{TRUE} allows odd electron ions for explaining MS/MS peaks, \code{FALSE}
#' means only even electron ions are considered.
#' @param mz If mz is >0, this is the mass GenForm will use to calculate the formulas, considering
#' the \code{ion} setting. Default is the most intense peak in \code{ms_file}.
#' @param exist Default \code{TRUE} activates the existance filter (allowing only formulae for which
#' at least one structural formula exists). \code{FALSE} is not recommended for most use cases.
#' @param an_loss Activate the annotation of MSMS fragments with subformulas and losses.
#' @param an_loss_limit Define the cut-off limit for annotating subformulas and losses. For
#' many peaks and formulas, this becomes time consuming and produces massive files of limited use.
#' @param sort Default \code{FALSE} will produce default output. If \code{TRUE}, generates an additional
#' sorted file sorted by combined match value.
#' @param cleanMSMS Default \code{FALSE}. If \code{TRUE}, produces an additional cleaned MSMS file
#' containing only peaks with subformula annotation.
#'
#' @return Several files containing various GenForm outputs.
#'
#' @author Emma Schymanski (R wrapper, \code{\link{emma.schymanski@@uni.lu}}) and Markus Meringer (GenForm)
#'
#' @references GenForm on Source Forge: https://sourceforge.net/projects/GenForm/
#'
#' @export
#'
#'
RunGenForm <- function(ms_file, msms_file, ion, GenFormDir, ResultsDir,
                       ppm=5, acc=15, elements="CHNOPS",ff="", dbe=TRUE, oei=TRUE, mz=0,
                       exist=TRUE, an_loss=TRUE,an_loss_limit=500, sort=FALSE, cleanMSMS=TRUE) {
  curr_dir <- getwd()
  # check files exist
  if(!file.exists(ms_file)) {
    stop(paste(ms_file, "doesn't exist, please try again"),sep=" ")
  }
  if(!file.exists(msms_file)) {
    stop(paste(msms_file, "doesn't exist, please try again"),sep=" ")
  }
  if(!file.exists(GenFormDir)) {
    stop(paste(GenFormDir, "doesn't exist, please try again"),sep=" ")
  }
  if(!file.exists(ResultsDir)) {
    dir.create(ResultsDir)
  }
  possible_ions <- c("+H","-H","-e","+e")
  if(!(ion %in% possible_ions)) {
    stop(paste("Incorrect ion state defined, try one of", paste(possible_ions,collapse=", ")))
  }
#   possible_sort <- c("ppm","msmv","msmsmv","combmv","none")
#   if(!(sort %in% possible_sort)) {
#     stop(paste("Incorrect sort criteria defined, try one of", paste(possible_sort,collapse=", ")))
#   }

  # build commands
  ms_file <- ms_file
  msms_file <- msms_file
  msms_file_name <- basename(msms_file)
  out_file <- paste(ResultsDir,"/",sub(".txt","_GenForm.txt",msms_file_name),sep="")
  log_file <- paste(ResultsDir,"/",sub(".txt","_GenForm_log.txt",msms_file_name),sep="")
  other_options <- ""
  if(dbe) other_options <- paste(other_options,"dbe")
  if(oei) other_options <- paste(other_options,"oei")
  if(exist) other_options <- paste(other_options,"exist")
#   if(!grepl("none",sort)) other_options <- paste(other_options," sort=",sort,sep="")
  # check if mz is way too high
  ms <- read.table(ms_file)
  #ms <- suppressWarnings(read.table("C:/DATA/Workflow/GenForm/121m0291a_MS.txt"))
  # build the command
  if (nchar(ff)>0 && mz>0) {
    GenFormBaseCommand <- paste("genform ms=",ms_file," msms=",msms_file," ion=",ion," ppm=",ppm," acc=",acc,
                            " m=",mz," ff=",ff,other_options,sep="")
  } else if (nchar(ff)>0 && mz==0) {
    GenFormBaseCommand <- paste("genform ms=",ms_file," msms=",msms_file," ion=",ion," ppm=",ppm," acc=",acc,
                            " ff=",ff,other_options,sep="")
  } else if (nchar(ff)==0 && mz>0) {
    GenFormBaseCommand <- paste("genform ms=",ms_file," msms=",msms_file," ion=",ion," ppm=",ppm," acc=",acc,
                            " m=",mz," el=",elements,other_options,sep="")
  } else {
    GenFormBaseCommand <- paste("genform ms=",ms_file," msms=",msms_file," ion=",ion," ppm=",ppm," acc=",acc,
                            " el=",elements,other_options,sep="")
  }
  # KJ added 
  GenFormBaseCommand <- paste0(GenFormDir, GenFormBaseCommand)
  # browser()
  GenFormCommand <- paste(GenFormBaseCommand," out=",out_file,sep="")
  system(GenFormCommand)
  #write(logmsg1, log_file,append=TRUE)
  write(GenFormCommand, log_file,append=TRUE)
  # test to see whether to run analyze/loss
  mz_to_test <- ms[which.max(ms[,2]),1]
  if(an_loss && (mz_to_test<an_loss_limit)) {
    out_file_an <- paste(ResultsDir,"/",sub(".txt","_GenForm_an.txt",msms_file_name),sep="")
    out_file_loss <- paste(ResultsDir,"/",sub(".txt","_GenForm_loss.txt",msms_file_name),sep="")
    GenFormCommand_an <- paste(GenFormBaseCommand," analyze out=",out_file_an,sep="")
    GenFormCommand_loss <- paste(GenFormBaseCommand," analyze loss out=",out_file_loss,sep="")
    system(GenFormCommand_an)
#    write(logmsg2, log_file,append=TRUE)
    system(GenFormCommand_loss)
    write(GenFormCommand_an, log_file,append=TRUE)
    write(GenFormCommand_loss, log_file,append=TRUE)
    #    write(logmsg3, log_file,append=TRUE)
  } else {
    print("No analyze/loss output performed, either mz is too large or option is off.")
  }
  #create a sorted output file too
  if (sort) {
    sorted_out_file <- paste(ResultsDir,"/",sub(".txt","_GenForm_sorted.txt",msms_file_name),sep="")
    output <- read.table(out_file)
    o <- order(output$V6,decreasing=TRUE)
    ordered_out <- output[o,]
    write.table(ordered_out,sorted_out_file,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
  }
  if (cleanMSMS) {
    cleanedMSMS_file <- paste(ResultsDir,"/",sub(".txt","_GenForm_cleanMSMS.txt",msms_file_name),sep="")
    GenFormCommand_clean <- paste(GenFormBaseCommand," oclean=",cleanedMSMS_file, " out",sep="")
    system(GenFormCommand_clean)
    write(GenFormCommand_clean, log_file,append=TRUE)
  }
  setwd(curr_dir)
}



#' Annotate Cleaned MSMS
#'
#' @description This takes the output from \code{\link{RunGenForm}} and creates a peaklist with the
#' annotations for further use.
#'
#' @usage annotateCleanMSMS(cleanMSMS, GenForm_an_out, GenForm_loss_out, lowestPpm=TRUE, round=TRUE,
#' removeAnnoSpace=FALSE)
#'
#' @param cleanMSMS File name for Cleaned MSMS text file output from \code{\link{RunGenForm}}
#' @param GenForm_an_out File name for annotated MSMS text file output from \code{\link{RunGenForm}}
#' @param GenForm_loss_out File name for annotated "loss" MSMS text file output from \code{\link{RunGenForm}}
#' @param lowestPpm If \code{TRUE}, choses the lowest ppm subformula for assignment of formula
#' to fragments. If \code{FALSE}, takes the first formula.
#' @param round If \code{TRUE}, rounds the output to 4 (mz) and 1 (Int) decimal places for aesthetics.
#' If \code{FALSE}, no rounding is done.
#' @param removeAnnoSpace If \code{TRUE}, removes the space in the \code{"no loss"} annotation.
#'
#' @return Returns annotated peaks for further use
#'
#' @author Emma Schymanski (\code{\link{emma.schymanski@@uni.lu}})
#'
#' @export
#'
#'
annotateCleanMSMS <- function(cleanMSMS, GenForm_an_out, GenForm_loss_out, lowestPpm=TRUE,
                              round=TRUE, removeAnnoSpace=FALSe) {
  peaks <- read.table(cleanMSMS)
  colnames(peaks) <- c("mz", "Int")
  peaks$Int <- 999*peaks$Int
  subform_an <- read.table(GenForm_an_out, header=F, sep="\t", skip=1)
  colnames(subform_an) <- c("mz", "subform_an", "ppm")
  subform_loss <- read.table(GenForm_loss_out, header=F, sep="\t", skip=1)
  colnames(subform_loss) <- c("mz", "subform_loss", "ppm")
  subform_an[,4] <- subform_loss$subform_loss
  subform_an[,5] <- 0
  colnames(subform_an) <- c("mz", "subform_an", "ppm", "subform_loss", "Int")
  # replace missing mzs with numbers, add intensity
  for (i in 1:length(subform_an$mz)) {
    if (is.na(subform_an$mz[i])) {
      subform_an$mz[i] <- subform_an$mz[(i-1)]
      #subform_loss$mz[i] <- subform_loss$mz[(i-1)]
    }
    subform_an$Int[i] <- peaks$Int[which(peaks$mz==subform_an$mz[i])]
  }
  # now deal with duplicate entries ...
  peaks$subform_an <- ""
  peaks$subform_loss <- ""
  peaks$ppm <- ""
  peaks$n_formulas <- 0
  for (n in 1:length(peaks$mz)) {
    subform_an_entries <- subform_an[which(subform_an$mz==peaks$mz[n]),]
    if (length(subform_an_entries$mz)==1) {
      peaks$subform_an[n] <- as.character(subform_an_entries$subform_an[1])
      peaks$subform_loss[n] <- as.character(subform_an_entries$subform_loss[1])
      peaks$ppm[n] <- subform_an_entries$ppm[1]
      peaks$n_formulas[n] <- 1
    } else if (length(subform_an_entries$mz)>1) {
      # choose either first entry or smallest ppm value
      if (lowestPpm) {
        ppm_ind <- which(abs(subform_an_entries$ppm)==min(abs(subform_an_entries$ppm)))
        peaks$subform_an[n] <- as.character(subform_an_entries$subform_an[ppm_ind[1]])
        peaks$subform_loss[n] <- as.character(subform_an_entries$subform_loss[ppm_ind[1]])
        peaks$ppm[n] <- subform_an_entries$ppm[ppm_ind[1]]
        peaks$n_formulas[n] <- length(subform_an_entries$mz)
      } else { #take first entry
        peaks$subform_an[n] <- as.character(subform_an_entries$subform_an[1])
        peaks$subform_loss[n] <- as.character(subform_an_entries$subform_loss[1])
        peaks$ppm[n] <- subform_an_entries$ppm[1]
        peaks$n_formulas[n] <- length(subform_an_entries$mz)
      }
    }
  }
  # adjust the rounding ... to make output look nicer
  if (round) {
    peaks$Int <- round(peaks$Int,1)
    peaks$mz <- round(peaks$mz,4)
  }

  return(peaks)

}

#' Fix Whitespace Error in GenForm Loss Annotation
#'
#' @description This takes a filename with e.g. the output from \code{\link{annotateCleanMSMS}} and
#' removes a space in the \code{"no loss"} annotation to allow for easier handling.
#'
#' @usage fixNoLossSpace(anno_MS_file, replace=TRUE)
#'
#' @param anno_MS_file File name for annotated cleaned MSMS text file output from \code{\link{annotateCleanMSMS}}
#' @param replace If \code{TRUE}, replace the input file. If \code{FALSE}, saved under the new name
#' \code{"*_ed.txt"}.
#'
#' @return File name of the "fixed" file.
#'
#' @author Emma Schymanski (\code{\link{emma.schymanski@@uni.lu}})
#'
#' @seealso \code{\link{annotateCleanMSMS}}
#'
#' @export
#'
fixNoLossSpace <- function(anno_MS_file, replace=TRUE) {
  fileConnection <- file(anno_MS_file)
  anno_MS_lines <- readLines(fileConnection)
  close(fileConnection)
  anno_MS_lines[length(anno_MS_lines)] <- sub("no loss", "noLoss",
                                              anno_MS_lines[length(anno_MS_lines)], fixed=TRUE)
  if (replace) {
    file_name <- anno_MS_file
  } else {
    file_name <- sub(".txt", "_ed.txt",anno_MS_file, fixed=TRUE)
  }
  write.table(anno_MS_lines, file_name, row.names=F, col.names = F, quote=F)
  return(file_name)
}

