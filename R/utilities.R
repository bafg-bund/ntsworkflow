# Utilities ####





#' Calculate the mass/charge for proton adduct from molecular formula
#'
#' If adduct is not given or [M] assumes a [M+H]+ or [M-H]- adduct if charge is not 0.
#'
#'
#' @param formula String of molecular formula, insert space before  atomic
#'   numbers
#' @param charge Integer
#'
#' @return Mass/charge ratio
#' @export
#'
#' @examples
#' # [M+H]+ Mass of Carbamazepine-13C15N
#' get_mass("C14 13CH12 15NNO", 1)
get_mass <- function(formula, charge, adduct = "[M]") {

  stopifnot(adduct %in% c("[M+H]+", "[M-H]-", "[M]+", "[M]-", "[M]"))
  chargeGiven <- charge
  # Replace isotopes with letter codes
  # 2H = D
  # 13C = X
  # 15N = L
  # 37Cl = Q
  # etc to add more

  formula <- gsub(" 37Cl", "Q", formula)
  formula <- gsub(" 2H", "D", formula)
  formula <- gsub(" 13C", "X", formula)
  formula <- gsub(" 15N", "L", formula)

  #load edited method from OrgMassSpecR
  ListFormulaIsot <- function (elemental.formula)
  {
    chr <- gregexpr("[[:upper:]][[:lower:]]{0,1}", elemental.formula)
    for (i in 1:length(chr[[1]])) {
      y <- attr(chr[[1]], which = "match.length")[i]
      z <- substr(elemental.formula, chr[[1]][i], chr[[1]][i] +
                    y - 1)
      if (!(z == "C" | z == "H" | z == "N" | z == "O" | z ==
            "S" | z == "P" | z == "Br" | z == "Cl" | z == "F" |
            z == "Si" | z == "Sn" | z == "D" | z == "X" | z == "Q" | z == "L" | z == "I"))
        stop(paste("Elemental formula", elemental.formula,
                      "contains element not of C,H,N,O,S,P,Br,Cl,F,Si,Sn."))
    }
    GetAtoms <- function(elemental.formula, element) {
      reg.exp <- paste(element, "[[:digit:]]*(?![[:lower:]])",
                       sep = "")
      x <- gregexpr(reg.exp, elemental.formula, perl = TRUE)
      if (x[[1]][1] != -1) {
        n <- vector(mode = "numeric", length = length(x[[1]]))
        for (i in 1:length(x[[1]])) {
          y <- attr(x[[1]], which = "match.length")[i]
          z <- substr(elemental.formula, x[[1]][i], x[[1]][i] +
                        y - 1)
          number <- as.numeric(strsplit(z, split = element)[[1]][2])
          if (is.na(number)) {
            n[i] <- 1
          }
          else {
            n[i] <- number
          }
          atoms <- sum(n)
        }
      }
      else {
        atoms <- 0
      }
      return(atoms)
    }
    elements <- c("C", "H", "N", "O", "S", "P", "Br", "Cl", "F",
                  "Si", "Sn", "D", "X", "Q", "L", "I")
    result <- as.list(sapply(elements, function(x) {
      GetAtoms(elemental.formula, x)
    }))
    return(result)
  }

  #load edited method OrgMassSpecR
  MonoisotopicMassIsot <- function (formula = list(), isotopes = list(), charge = 0)
  {
    defaultFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0,
                           P = 0, Br = 0, Cl = 0, F = 0, Si = 0, D = 0, X = 0, Q = 0, L = 0, I = 0)
    defaultFormula[names(formula)] <- formula
    defaultIsotopes <- list(C = 12, H = 1.0078250321, N = 14.0030740052,
                            O = 15.9949146221, S = 31.97207069, P = 30.97376151,
                            Br = 78.9183376, Cl = 34.96885271, F = 18.9984032, Si = 27.9769265327, D = 2.0141017780, X = 13.0033548378, L = 15.0001088984, Q = 36.9665, I = 126.9050)
    defaultIsotopes[names(isotopes)] <- isotopes
    if (charge < 0 & abs(charge) > defaultFormula$H + defaultFormula$D)
      stop("the number of negative charges exceeds the number of hydrogens in the formula list")
    mass <- (defaultFormula$C * defaultIsotopes$C + defaultFormula$H *
               defaultIsotopes$H + defaultFormula$N * defaultIsotopes$N +
               defaultFormula$O * defaultIsotopes$O + defaultFormula$S *
               defaultIsotopes$S + defaultFormula$P * defaultIsotopes$P +
               defaultFormula$Br * defaultIsotopes$Br + defaultFormula$Cl *
               defaultIsotopes$Cl + defaultFormula$F * defaultIsotopes$F +
               defaultFormula$Si * defaultIsotopes$Si + defaultFormula$D *
               defaultIsotopes$D + defaultFormula$X *
               defaultIsotopes$X + defaultFormula$L *
               defaultIsotopes$L + defaultFormula$Q *
               defaultIsotopes$Q + defaultFormula$I * defaultIsotopes$I)
    if (adduct == "[M+H]+") {
      stopifnot(charge == 1)
      mass <- abs((mass + charge * 1.007276466)/charge)
    } else if (adduct == "[M-H]-") {
      stopifnot(charge == -1)
      mass <- abs((mass + charge * 1.007276466)/charge)
    } else if (adduct == "[M]+") {
      stopifnot(charge == 1)
      mass <- mass - 0.00054858  # mass of electron
    } else if (adduct == "[M]-") {
      stopifnot(charge == -1)
      mass <- mass + 0.00054858
    } else if (charge != 0 && adduct == "[M]"){
      mass <- abs((mass + charge * 1.007276466)/charge)
    } else if (charge == 0 && adduct == "[M]") {
      mass <- mass
    } else {
      stop(paste("cannot calculate mass for", formula, "with charge", charge, "and adduct", adduct))
    }

    return(mass)
  }

  # convert formula to list
  formula.list <- ListFormulaIsot(formula)

  # calculation of mass
  MonoisotopicMassIsot(formula = formula.list, charge = chargeGiven)
}

# Infix function to return default value ####
`%||%` <- function(a, b) if (!is.null(a)) a else b
`%notin%` <- function(x, y) !(x %in% y)

dot_every <- function(n, f) {
  i <- 1
  function(...) {
    if (i %% n == 0) message(".", appendLF = FALSE)
    i <<- i + 1
    f(...)
  }
}


#' Simple function to create window tolerance
#' 
#' Internal package use only
#' 
#' @param value 
#' @param tol 
#'
#' @return two numbers
wind <- function(value, tol) {
  c(value - tol, value + tol)
} 


#' Adjust a formula to an adduct
#'
#' @param formula
#' @param adduct
#'
#' @return
#' @export
#'
#' @examples
correct_formula <- function(formula, adduct) {

  if (adduct == "[M+H]+") {
    formula <- paste0(formula, "H")
    thisCharge <- 1
  } else if (adduct == "[M+Na]+") {
    formula <- paste0(formula, "Na")
    thisCharge <- 1
  } else if (adduct == "[M]+") {
    thisCharge <- 1
  } else if (adduct == "[M-H]-") {
    thisCharge <- -1
    # adjust number of protons
    #browser()
    noProt <- as.numeric(stringr::str_match(formula, "\\SH(\\d+)")[, 2]) - 1
    formula <- if (is.na(noProt) && grepl("H", formula)) {
      stringr::str_replace(formula, "H([A-Z])", "\\1")
    } else if (is.numeric(noProt) && noProt != 0) {
      stringr::str_replace(formula, "(\\SH)\\d+", paste0("\\1", noProt))
    } else {
      formula
    }
  } else if (adduct == "[M+NH4]+") {
    thisCharge <- 1
    formula <- paste0(formula, "NH4")
  } else if (adduct == "[M+H2CO2-H]-") {
    thisCharge <- -1
    formula <- paste0(formula, "HCO2")
  } else if (adduct == "[M-H2O+H]+") {
    thisCharge <- 1
    # adjust number of protons
    noProt <- as.numeric(stringr::str_match(formula, "H(\\d+)")[, 2]) - 1
    formula <- paste0(stringr::str_match(formula, "(.*H)\\d+")[, 2],
                      noProt, stringr::str_match(formula, ".*H\\d+(.*)")[, 2])
    #check number oxygen atoms
    oxform <- as.numeric(stringr::str_match(formula, "O(\\d+)")[, 2])
    if (is.na(oxform)) {
      #only one oxygen atom delete O from formular
      formula <- paste0(stringr::str_replace(formula,"O",""))
    } else if(oxform1 > 1) {
      #more oxygen atoms reduce number of atoms
      oxform <- as.numeric(stringr::str_match(formula, "O(\\d+)")[, 2]) - 1
      formula <- paste0(stringr::str_match(formula, "(.*O)\\d+")[, 2],
                        oxform, stringr::str_match(formula, ".*O\\d+(.*)")[, 2])
    } else if (oxform <= 1) {
      formula <- paste0(stringr::str_replace(formula,"O1",""))
    } else {
      stop("error in no. of oxygens")
    }
  } else {
    formula <- NA
    thisCharge <- -1
  }

  list(form = formula, charge = thisCharge)
}


compact <- function(x) {
  Filter(Negate(is.null), x)
  }


#' Export Function to For-Ident.
#' 
#' For batch processing in For-Ident, txt file is created for import on For-Ident platform. MS2 
#' spectra are included for Met-frag or Massbank search. The spectrum for each row is taken from the
#' sample with the highest intensity in this row. The function takes approx. 1 min to process 100 
#' rows, so limit the size of the alignment matrix passed to the app. If no MS2 is available, row is
#' skipped.
#' 
#' 
#' @param align_matrix Alignment \code{matrix} from Christian's App
#' @param sampleList sampleList \code{matrix} from App
#' @param rawDataList datenList \code{list} from App
#' @param intThreshRel Intensity threshold for including peaks in MS2 spectra
#' @param report_nm Name of output file
#'   
#' @return A txt file with each found mass, retention time and MS2 spectrum (if available).
#'   
#' @return
#' @export
#' 
#' @examples
forIdentExport <- function(align_matrix, sampleList, rawDataList, intThreshRel = 0.01, report_nm) {
  
  # change sample_list_save to dataframe with real variables
  sl <- as.data.frame(sampleList, stringsAsFactors = FALSE)
  sl$ID <- as.numeric(sl$ID)
  sl$File <- as.character(sl$File)
  sl$RAM <- as.logical(sl$RAM)
  
  al <- as.data.frame(align_matrix)
  
  # loop through each row in alignment matrix,
  # assume first that all samples are loaded into the RAM
  
  for (i in 1:nrow(al)) { #i <- 3
    if (i %% 100 == 0) message(paste(i, "completed."), appendLF = TRUE)
    # get sample ID of most intense peak in the row
    intens <- dplyr::select(al[i, ], dplyr::matches("Int_"))
    sampId <- as.numeric(stringr::str_match(colnames(intens)[which.max(intens)], "Int_(\\d+)$")[, 2])
    
    filenm <- stringr::str_match(basename(sl[sl$ID == sampId, "File"]), "(.*)\\.mzX?ML$")[, 2]
    scan_i <- al[i, sprintf("ms2scan_%i", sampId)]
    id_i <- al[i, sprintf("PeakID_%i", sampId)]
    
    if (scan_i == 0 || id_i == 0)
      next
    
    if (!sl[sl$ID == sampId, "RAM"]) {
      rawDataList[[sampId]] <- xcms::xcmsRaw(rawDataList[[sampId]]@filepath@.Data, includeMSn = TRUE)
      sl[sl$ID == sampId, "RAM"] <- TRUE
    }
    
    nm <- paste0(filenm, "_unknown", id_i)
    mz <- al[i, sprintf("mz_%i", sampId)]
    rt <- al[i, sprintf("RT_%i", sampId)] / 60
    form <- "No data for Formula Finder Composition Name xcm"
    spec <- getMsnScan(rawDataList[[sampId]], scan_i)
    # subset spec to mz and int threshold
    intThresh <- max(spec[, "intensity"]) * intThreshRel
    fn <- paste0(report_nm, ".txt")
    spec <- spec[spec[, "mz"] < mz + 1 & spec[, "intensity"] > intThresh, , drop = FALSE]
    cat("\n\nNAME: ", nm, "\n", file = fn, append = TRUE, sep = "")
    cat("RETENTIONTIME: ", rt, "\n", file = fn, append = TRUE, sep = "")
    cat("PRECURSORMZ: ", mz, "\n", file = fn, append = TRUE, sep = "")
    cat("Formula: ", form, "\n", file = fn, append = TRUE, sep = "")
    
    for (r in 1:nrow(spec)) {
      cat(sprintf("%.4f", spec[r, 1]), sprintf("%.4f ", spec[r, 2],5), file = fn, append = TRUE)
    }
    cat("\n//", file = fn, append = TRUE)
    
  }
  
}

#' Search for trues in a row
#'
#' @param set logical 
#' @param in_a_row 
#'
#' @return Logical of length 1 
#' 
#' @export
#'
#' @examples
trueInaRow<- function(set, in_a_row = 2) { # set <- c(T,F,T,T,T,F,T,F,T,T)
  # need to make the set multiple of pattern
  #browser()
  pattern = rep(T, in_a_row)
  left_over <- length(set) %% length(pattern)
  set <- append(set, rep(F, left_over))
  # split set at every position from 1 to pattern length - 1 (the last subset goes to end of set) 
  step_size <- length(pattern) - 1
  start_pos <- 1:(length(set) - step_size)
  subsets <- list()
  for (i in start_pos)
    subsets[[i]] <- unname(set[i:(i + step_size)])
  any(vapply(subsets, identical, logical(1), y = pattern))
}

#' Load a report file using a dialog box or by giving the path to the file
#' 
#' @export
loadReport <- function(dialog = TRUE, path = NULL) {
  if (dialog)
    path <- rstudioapi::selectFile("Select Report", filter = "DBscreening report file (*.report)")
  stopifnot(file.exists(path))
  readRDS(path)
}

#' Merge two Report objects
#' @return a Report object
#' @export
mergeReport <- function(report1, report2) {
  # copy to avoid reference semantics
  x <- report1$copy()
  y <- report2$copy()
  # check that no samples are the same
  stopifnot(length(intersect(basename(x$rawFiles), basename(y$rawFiles))) == 0)
  # check that IS is the same
  stopifnot(identical(x$IS, y$IS))
  # check that rawData is clear
  stopifnot(length(x$rawData) == 0 && length(y$rawData) == 0)
  #browser()
  xPeakID <- x$currentPeakID - 1
  xISpeakID <- x$currentISpeakID - 1
  # change IDs of report y
  y$peakList$peakID <- y$peakList$peakID + xPeakID
  y$EIC$peakID <- y$EIC$peakID + xPeakID
  y$MS1$peakID <- y$MS1$peakID + xPeakID
  y$MS2$peakID <- y$MS2$peakID + xPeakID
  y$ISresults$ISpeakID <- y$ISresults$ISpeakID + xISpeakID
  
  # append data
  x$currentPeakID <- max(y$peakList$peakID) + 1
  x$currentISpeakID <- max(y$ISresults$ISpeakID) + 1
  x$rawFiles <- append(x$rawFiles, y$rawFiles) 
  x$rawFilesCompl <- rbind(x$rawFilesCompl, y$rawFilesCompl)
  x$peakList <- rbind(x$peakList, y$peakList)
  x$EIC <- rbind(x$EIC, y$EIC)
  x$MS1 <- rbind(x$MS1, y$MS1)
  x$MS2 <- rbind(x$MS2, y$MS2)
  x$ISresults <- rbind(x$ISresults, y$ISresults)
  
  x
}


