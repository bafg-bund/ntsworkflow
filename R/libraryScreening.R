# Library screening ####

# Includes functions for library screening used in conjunction with Report class
# (dbas) and the non-target app (peak picking).



#' Library screening without prior peak-picking
#'
#' @description
#' A direct library screening approach using measurement files and the CSL
#' spectral library in SQLite format.
#'
#'
#' @param data_path Path to rawfiles, can be a vector of file locations or a list of xcmsRaw objects
#' @param csl_path Path to the Collective Spectral Library (CSL)
#' @param rttolm Retention time tolerance in minutes
#' @param mztolu m/z tolerance in Da (MS2 precursor mass)
#' @param mztolu_fine m/z tolerance in Da (second stage)
#' @param chromatography Chromatography name to choose the correct retention times from the library
#' @param pol Polarity
#' @param CE_s Collision energy
#' @param CES_s Collision energy spread
#' @param instr Instruments allowed, must match the instruments listed in the library
#' @param ceunit Collision energy units
#' @param comparison How the spectra should be compared, at the moment only "dot_product" is allowed
#' @param threshold Threshold score for spectral comparison (for dot_product 1000 is max)
#' @param rt_res Resolution of two peaks in chromatography (min)
#' @param rtoffset Retention time offset (will be added to library RTs)
#' @param ndp_m Peak intensity weighting factor for dot-product
#' @param ndp_n m/z weighting factor for dot-product
#' @param mztolu_ms2 m/z window for dot-product
#' @param compounds character vector of compounds which should be processed
#'
#' @details
#' The function is used in the Report RC to obtain the peak list. No prior peak
#' picking is needed for this function. The function starts from the list of
#' masses provided in the spectral library and searches for these in the
#' measurement data. The measurement data is in mz(X)ML format and contains DDA
#' MS² spectra.
#' \strong{Multiple peaks for one compound}
#' Will be sorted by intensity and labelled A-Z  # TODO: Check if this is already the case
#'
#' ## References
#'   Jewell, K. S., Kunkel, U., Ehlig, B., Thron, F., Schlüsener, M., Dietrich,
#'   C., . . . Ternes, T. A. (2019). Comparing mass, retention time and MS2
#'   spectra as criteria for the automated screening of small molecules in
#'   aqueous environmental samples analyzed by LC-QToF-MS/MS. Rapid
#'   Communications in Mass Spectrometry, 34, e8541. doi:10.1002/rcm.8541
#'
#' @returns a `data.frame` with hits from the CSL
#'
#' @seealso [annotate_grouped()]
#'
#' @export
#' @import dplyr
ms2_search <- function(data_path, csl_path,
                       chromatography = "bfg_nts_rp1",
                       pol = "pos",
                       CE_s = 30:40,
                       CES_s = 0:15,
                       instr = "LC-ESI-QTOF TripleTOF 5600 SCIEX",
                       rtoffset = 0,
                       rt_res = 1.5,
                       rttolm = 1,
                       mztolu = 0.5,
                       mztolu_fine = 0.005,
                       mztolu_ms2 = 0.015,
                       threshold = 400,
                       ndp_m = 2,
                       ndp_n = 1,
                       compounds = NULL,
                       ceunit = c("V", "eV"),  # Todo: Is not used
                       comparison = "dot_product"  # Todo: Is not used
) {
  
  # Validate and normalize data file paths
  stopifnot(inherits(data_path, "character") || inherits(data_path, "list"))
  if (inherits(data_path, "character")) {
    data_path <- normalizePath(data_path)
  }
  
  # Collect parameter settings
  par <- list(
    rtoffset = rtoffset,
    rt_res = rt_res,
    rttolm = rttolm,
    mztolu = mztolu,
    mztolu_fine = mztolu_fine,
    mztolu_ms2 = mztolu_ms2,
    ndp_m = ndp_m,
    ndp_n = ndp_n,
    threshold = threshold
  )
  
  # Connect to the CSL database and filter relevant data
  csl_data <- filter_csl_data(csl_path, chromatography, pol, CE_s, CES_s, instr)
  
  # Filter suspects by specified compound names
  csl_data$suspects <- filter_suspects_by_name(csl_data$suspects, compounds)
  
  # Process all measurement files
  all_samps <- lapply(data_path, function(pth) {
    eval_samp(pth, csl_data, par)
  })
  all_samps <- compact(all_samps)
  
  # Combine results
  if (length(all_samps) == 0) return(NULL)
  all_samps <- do.call("rbind", all_samps)
  attr(all_samps, "processing_date") <- date()
  attr(all_samps, "csl_path") <- csl_path
  
  return(all_samps)
}


# Subfunctions ####

filter_csl_data <- function(csl_path, chromatography, pol, CE_s, CES_s, instr){
  
  # Connect to the CSL database
  db <- DBI::dbConnect(RSQLite::SQLite(), csl_path)
  
  # Get tables from the CSL database
  expTable <- tbl(db, "experiment")
  paraTable <- tbl(db, "parameter")
  rtTable <- tbl(db, "retentionTime")
  compTable <- tbl(db, "compound")
  
  # Filter suspects using parameters
  suspects <- compTable %>%
    left_join(rtTable, by = "compound_id") %>%
    filter(chrom_method == chromatography || is.na(chrom_method)) %>%
    left_join(expTable, by = "compound_id") %>%
    left_join(paraTable, by = "parameter_id") %>%
    filter(polarity == pol) %>%
    filter(CE %in% CE_s) %>%
    filter(CES %in% CES_s) %>%
    filter(instrument %in% instr) %>%
    select(name, CAS, mz, rt, adduct, isotope, experiment_id, compound_id) %>%
    dplyr::collect() %>%
    distinct()
  
  # Collect other tables
  exptbl <- tbl(db, "experiment") %>% collect()
  comptbl <- tbl(db, "compound") %>%
    select(-chem_list_id) %>%
    collect()
  fragtbl <- tbl(db, "fragment") %>% collect()
  
  # Disconnect the CSL database
  DBI::dbDisconnect(db)
  
  # Bundle results into a list
  csl_data <- list(
    suspects = suspects,
    exptbl = exptbl,
    comptbl = comptbl,
    fragtbl = fragtbl
  )
  
  return(csl_data)
}


filter_suspects_by_name <- function(suspects, compounds){
  
  # Filter suspects for compounds
  if (!is.null(compounds)) {
    filtered_suspects <- filter(suspects, name %in% compounds)
    if (nrow(suspects) == 0) {
      warning("Chosen compounds and settings not available in DB")
      return(NULL)
    } else {
      message("Searching compounds ", paste(unique(suspects$name), collapse = "; "))
    }
  } else {
    filtered_suspects <- suspects
  }
  
  return(filtered_suspects)
}


eval_samp <- function(pth, csl_data, par) {
  
  # Create xcmsRaw class
  # pth = data_path[1]  # TODO: testing
  if (inherits(pth, "character")) {
    raw_data <- xcms::xcmsRaw(pth, includeMSn = TRUE)
  } else if (inherits(pth, "xcmsRaw")) {
    raw_data <- pth
  } else {
    stop("data_path is not character or list of xcmsRaw")
  }
  
  # Separates available CSL compounds by unique name-adduct-isotope key   
  csl_data$suspects$comp_adduct_isot <- paste(csl_data$suspects$name, csl_data$suspects$adduct, csl_data$suspects$isotope, sep = "_")
  splitByComp <- split(csl_data$suspects, csl_data$suspects$comp_adduct_isot)
  
  # Run evaluation for each available CSL compounds
  # comp = splitByComp[1]  # TODO: testing
  all_csl_comps <- lapply(splitByComp, function(csl_comp) {
    eval_comp(csl_comp, raw_data, csl_data, par)
  })
  
  all_csl_comps <- Filter(function(x) !inherits(x, "try-error"), all_csl_comps)
  
  all_csl_comps <- compact(all_csl_comps)
  
  if (inherits(pth, "character")) {
    rm(raw_data)
  }
  
  if (length(all_csl_comps) == 0) {
    return(NULL)
  }
  
  # Join all compounds
  all_csl_comps <- do.call("rbind", all_csl_comps)
  return(all_csl_comps)
}


# For each compound+adduct+isotopologue in the CSL, search raw data to find 
# MS2 spectra with matching precursor mass and RT
eval_comp <- function(csl_comp, raw_data, csl_data, par) {
  
  # Create data frame with precursor mass and RT data
  df <- data.frame(
    mz = raw_data@msnPrecursorMz,
    rt = raw_data@msnRt,
    index = seq_along(raw_data@msnPrecursorMz)
  )
  
  # Filter matching precursor mass and RT of sample and CSL compound within the set tolerances
  if (is.na(csl_comp$rt[1])) {
    # If RT is unavailable, filter only by precursor mass
    df <- df[abs(df[, "mz"] - csl_comp$mz[1]) <= mztolu, , drop = FALSE]
  } 
  else {
    # Filter by precursor mass and RT with added RT offset
    csl_comp$rt <- csl_comp$rt + par$rtoffset
    combination <- cbind(
      abs(df[, "mz"] - csl_comp$mz[1]) <= par$mztolu,
      abs(df[, "rt"] - csl_comp$rt[1] * 60) <= par$rttolm * 60
    )
    df <- df[apply(combination, 1, all), , drop = FALSE]
  }
  
  # If no entries (matches) remain, return NULL
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # Retrieves MS2 spectra for each match
  data_specs <- lapply(df[, "index"], xcms::getMsnScan, object = raw_data)
  data_specs <- lapply(data_specs, as.data.frame)
  
  # Add the matched compound name to each spectrum
  data_specs <- lapply(data_specs, function(x) {
    attr(x, "comp_name") <- csl_comp$name[1]
    x
  })
  
  # Add the precursor mass to each spectrum
  data_specs <- Map(set_mz_attribute, data_specs, df[, "mz"])
  
  # Normalize each MS2 data spectrum 
  data_specs <- lapply(data_specs, filter_and_normalize_ms2)
  
  # Remove df entries with NULL spectrum
  stopifnot(nrow(df) == length(data_specs))
  keep <- !vapply(data_specs, is.null, logical(1))
  df <- df[keep, ]
  data_specs <- data_specs[keep]
  
  # Check if any spectrum remains
  if (is.null(data_specs) || length(data_specs) == 0) {
    return(NULL)
  }
  stopifnot(nrow(df) == length(data_specs))
  
  # Collect spectra from the CSL database
  # Note: parallel cores can not read the .db file
  exptt <- csl_data$exptbl %>% select(experiment_id, compound_id, mz)
  compt <- csl_data$comptbl %>% select(compound_id, CAS, name)
  spectra <- filter(csl_data$fragtbl, experiment_id %in% !!csl_comp$experiment_id) %>%
    select(mz, int, experiment_id) %>%
    left_join(exptt, by = "experiment_id") %>%
    left_join(compt, by = "compound_id")
  
  # Group spectra by unique experiment ID
  csl_specs <- split(spectra, spectra$experiment_id)
  
  # Add attributes to each CSL spectrum and rename columns
  csl_specs <- lapply(csl_specs, function(spec) {
    res <- spec[, c("mz.x", "int")]
    colnames(res) <- c("mz", "int")
    attr(res, "precursor_mz") <- spec$mz.y[1]
    attr(res, "comp_name") <- spec$name[1]
    attr(res, "CAS") <- spec$CAS[1]
    attr(res, "db_exp_ID") <- spec$experiment_id[1]
    res
  })
  
  # Normalize each MS2 CSL spectrum 
  csl_specs <- lapply(csl_specs, filter_and_normalize_ms2)
  
  # Remove NULL entries and check if any spectrum remains
  csl_specs <- compact(csl_specs)
  if (is.null(csl_specs) || length(csl_specs) == 0) {
    return(NULL)
  }
  
  # Compare MS2 spectra and add best match and similarity score.
  df <- filter_matching_ms2(df, data_specs, csl_specs, par)
  
  
  # Keep proposed CSL spectra with a similarity score over the set threshold
  df <- df[df$score >= par$threshold, , drop = FALSE]
  
  # Check if any spectrum remains
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # Search MS1 scans near the MS2 RT for a precursor ion matching the expected m/z
  # and add the m/z and intensity of the closest matching peak
  df <- search_ms1_for_precursor(df, raw_data, par, csl_comp)
  
  # Group and label separate peaks based on RT
  df <- label_separate_peaks(df, par, csl_comp)
  
  # Add metadata
  df$rt_min <- round(df$rt / 60, 2)
  df$comp_id <- csl_comp$compound_id[1]
  df$samp <- basename(raw_data@filepath)
  df$comp_name <- csl_comp$name[1]
  df$comp_CAS <- csl_comp$CAS[1]
  df$mz_error_mDa <- round(abs(csl_comp$mz[1] - df$real_mz) * 1000, 3)
  df$adduct <- csl_comp$adduct[1]
  df$isotopologue <- csl_comp$isotope[1]
  if (is.na(csl_comp$rt[1])) {
    df$rt_error_min <- NA
  } else {
    df$rt_error_min <- round(abs(csl_comp$rt[1] - df$rt / 60), 2)
  }
  if (inherits(df, "try-error")) {
    return(NULL)
  }
  return(df)
}


# Helper functions eval_comp ####

set_mz_attribute <- function(spec, mz) {
  # Adds attribute to object
  attr(spec, "precursor_mz") <- mz
  return(spec)
}


#' Filter and normalize MS2 spectra by removing precursor peaks and noise.
#'
#' This function filters out peaks close to the precursor m/z and low-intensity 
#' noise, normalizes intensities to a maximum of 1, and retains relevant metadata.
#' Returns NULL if no valid peaks remain.
filter_and_normalize_ms2 <- function(x) {  #Todo: x bitte umbenennen in csl_spec
  # Retrieve attributes
  compname <- attr(spec, "comp_name")
  compmz <- attr(spec, "precursor_mz")
  expID <- try(attr(spec, "db_exp_ID"))  # Todo use tryCatch
  
  # Remove precursor otherwise essentially MS1 checked again # TODO: What does the last part mean?
  # TODO: Remove peaks within 0.02 units of the precursor m/z
  spec <- spec[spec$mz < attr(spec, "precursor_mz") - 0.02, ]
  
  # Remove low-intensity noise
  spec <- spec[spec$int >= 0.1, ]
  
  # Check if any peak remains
  if (length(spec) == 0 || nrow(spec) == 0) {
    return(NULL)
  }
  
  # Normalize intensity
  spec_clean <- data.frame(mz = spec$mz, int = spec$int / max(spec$int))
  
  # Add attributes
  attr(spec_clean, "comp_name") <- compname
  attr(spec_clean, "precursor_mz") <- compmz
  if (is.numeric(expID)) {
    attr(spec_clean, "db_exp_ID") <- expID
  }
  return(spec_clean)
}


#' Compare MS2 spectra to a reference database and record the best match and similarity score.
filter_matching_ms2 <- function(df, data_specs, csl_specs, par){

  # Define function to compute a similarity score between two MS2 spectra
  custom_calc_ndp <- function(d_spec, csl_spec) {
    calc_ndp(d_spec, csl_spec, par$ndp_m, par$ndp_n, par$mztolu_ms2)
  }
  # Compute similarity matrix (pairwise comparison of all data and CSL spectra)
  scorA <- outer(data_specs, csl_specs, Vectorize(custom_calc_ndp))
  
  # Record maximum similarity score for each data spectrum
  df$score <- apply(scorA, 1, max, na.rm = TRUE)
  
  # Record experiment ID of the best match
  allExpIds <- sapply(csl_specs, attr, which = "db_exp_ID")
  df$expID <- allExpIds[apply(scorA, 1, which.max)]
  
  return(df)
}


#' Calculate Normalized Dot Product Similarity Score Between MS2 Spectra
#'
#' @param data_spec Spectrum from data
#' @param csl_spec Spectrum from CSL database
#' @param ndp_m Peak intensity weighting factor
#' @param ndp_n m/z weighting factor
#' @param mztolu m/z tolerance for binning peaks (u)
#'
#' @details See MetCirc::NDP for more information
#'
#' @returns Similarity score from 0 to 1000, where 1000 indicates a perfect match
#' @export
#'
calc_ndp <- function(data_spec, csl_spec, ndp_m = 2, ndp_n = 1, mztolu_ms2 = 0.015) {
  
  # Find all pairs of peaks that match within the m/z tolerance
  ar <- which(abs(outer(data_spec[, 1], csl_spec[, 1], "-")) <= mztolu_ms2,
              arr.ind = TRUE
  )
  
  # If no matching peaks are found, return a similarity score of 0
  if (nrow(ar) == 0) {
    return(0)
  }
  
  # Extract m/z values of matched peak pairs from both spectra
  m <- cbind(data_spec[, 1][ar[, 1]], csl_spec[, 1][ar[, 2]])
  
  # Extract intensities of matched peaks
  data_int <- data_spec[, 2][ar[, 1]]
  csl_int <- csl_spec[, 2][ar[, 2]]
  
  # Compute mean m/z for each matching pair
  masses <- rowMeans(m)
  
  # Extract unmatched peaks from both spectra
  data_nonmatching <- data_spec[-ar[, 1], ]
  csl_nonmatching <- csl_spec[-ar[, 2], ]
  
  # Add unmatched sample peaks to vector
  data_int <- append(data_int, data_nonmatching[, 2])
  masses <- append(masses, data_nonmatching[, 1])
  csl_int <- append(csl_int, rep(0, nrow(data_nonmatching)))
  
  # Add unmatched reference peaks to vector
  csl_int <- append(csl_int, csl_nonmatching[, 2])
  masses <- append(masses, csl_nonmatching[, 1])
  data_int <- append(data_int, rep(0, nrow(csl_nonmatching)))
  
  # Apply weighting to intensity and m/z to get weighted spectra
  WS1 <- data_int^ndp_m * masses^ndp_n
  WS2 <- csl_int^ndp_m * masses^ndp_n
  
  # Calculate normalized dot product between the two weighted spectra
  r_ndp <- (sum(WS1 * WS2))^2 / (sum(WS1^2) * sum(WS2^2))
  
  # Special case: if either spectrum has only one fragment, then the most intense
  # fragment in both spectra must be the same mass, otherwise return 0
  if (nrow(data_spec) == 1 || nrow(csl_spec) == 1) {
    mzD <- data_spec[which.max(data_spec[, 2]), 1]
    mzS <- csl_spec[which.max(csl_spec[, 2]), 1]
    if (abs(mzD - mzS) > mztolu_ms2) {
      r_ndp <- 0
    }
  }
  
  return(r_ndp * 1000)
}


#' Search MS1 scans near the MS2 retention time for a precursor ion matching the expected m/z.
search_ms1_for_precursor <- function(df, raw_data, par, csl_comp){
  # TODO: Include an isotope check (calculate match to theoretical)
  
  # Initialize output columns
  df$mz_ok <- NA
  df$real_mz <- NA
  df$int_h <- NA
  
  # Loop over all MS2 spectra
  for (i in 1:nrow(df)) {
    
    # Find the MS1 scan closest in retention time
    ind <- which.min(abs(raw_data@scantime - df$rt[i]))
    
    # Extract peaks within ±0.05 Da of the target m/z
    ms1 <- xcms::getScan(
      raw_data, ind,
      c(df$mz[i] - 0.05, df$mz[i] + 0.05)
    )
    
    # If no peaks were found, try a wider window (±0.5 Da)
    if (nrow(ms1) == 0) {
      ms1 <- xcms::getScan(
        raw_data, ind,
        c(df$mz[i] - 0.5, df$mz[i] + 0.5)
      )
    }
    
    # If still no peaks were found, then log a warning and skip to next spectrum
    if (nrow(ms1) == 0) {
      warning(sprintf(
        "Error: no precursor found for mz %.4f @ %.2f min",
        df$mz[i], df$rt[i] / 60
      ))
      df$mz_ok[i] <- FALSE
      next
    }
    
    # Check if any peaks in MS1 are within fine m/z tolerance of expected precursor
    mz_error <- abs(ms1[, "mz"] - csl_comp$mz[1])
    if (!any(mz_error <= par$mztolu_fine)) {
      df$mz_ok[i] <- FALSE
      next
    }
    
    # Store the m/z and intensity of the closest matching peak (lowest mz-error)
    df$real_mz[i] <- ms1[which.min(mz_error), "mz"]
    df$int_h[i] <- ms1[which.min(mz_error), "intensity"]
    df$mz_ok[i] <- TRUE
  }
  
  # If no spectra had matching precursor peaks, return NULL
  if (!any(df$mz_ok)) {
    return(NULL)
  }
  
  # Keep only rows with valid precursor matches
  df <- df[df$mz_ok, , drop = FALSE]
  
  return(df)
}


#' Group and label separate chromatographic peaks (A-Z) based on retention time, 
#' keeping only the most intense from each group.
label_separate_peaks <- function(df, par, csl_comp){

  # Case: more than one MS2 match — cluster peak groups by retention time
  if (nrow(df) > 1) {
    df <- df[order(df$rt), ]     # Sort by retention time
    peakNum <- 1                 # Start with first peak group
    df$peak <- LETTERS[peakNum]  # Label first row as "A"
    
    # Loop through rows (MS2 spectra) to separate into distinct peaks
    for (i in 2:nrow(df)) {
      diff <- df$rt[i] - df$rt[i - 1]
      
      # If RT difference exceeds threshold (par$rt_res), start a new peak group
      if (diff > par$rt_res * 60) {
        peakNum <- peakNum + 1
        df[i:nrow(df), "peak"] <- LETTERS[peakNum]
      }
    }
    
    # From each peak group, keep only the row with highest intensity
    good <- by(df, df$peak,
               function(x) x[which.max(x$int_h), ],
               simplify = FALSE)
    
    # Combine selected rows back into one data.frame
    if (length(good) == 1) {
      df <- good[[1]]
    } else {
      df <- do.call("rbind", good)
    }
    # TODO: sort peaks by intensity (highest is peak A)
  } 
  
  # Case: only one peak found — label it "A" 
  else if (nrow(df) == 1) {
    df$peak <- "A"
  } 
  
  # Case: no valid peaks — stop and print error
  else {
    stop(sprintf(
      "error in CSL compound %s during collection of plausible peaks",
      csl_comp$name[1]
    ))
  }
  
  return(df)
}


# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow