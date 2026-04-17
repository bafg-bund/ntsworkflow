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
#' @param ce_s Collision energy
#' @param ces_s Collision energy spread
#' @param instr Instruments allowed, must match the instruments listed in the library
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
ms2_search <- function(data_path, 
                       csl_path,
                       chromatography = "bfg_nts_rp1",
                       pol = "pos",
                       ce_s = 30:40,
                       ces_s = 0:15,
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
                       compounds = NULL
) {
  
  # Validate and normalize data file paths
  stopifnot(inherits(data_path, "character") || inherits(data_path, "list"))
  if (inherits(data_path, "character")) {
    data_path <- normalizePath(data_path)
  }
  
  # Collect parameter settings
  params <- list(
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
  
  # Connect to the CSL, get relevant data, and filter by specified compounds
  csl_data <- get_csl_data(csl_path, chromatography, pol, ce_s, ces_s, instr, 
                           compounds)
  
  # Process all measurement files
  all_samps <- lapply(data_path, function(pth) {
    eval_samp(pth, csl_data, params)
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

get_csl_data <- function(csl_path, chromatography, pol, ce_s, ces_s, instr, 
                         compounds){
  
  # Connect to the CSL database
  db <- DBI::dbConnect(RSQLite::SQLite(), csl_path)
  
  # Get tables from the CSL database
  expTable <- tbl(db, "experiment")
  paraTable <- tbl(db, "parameter")
  rtTable <- tbl(db, "retention_time")
  compTable <- tbl(db, "compound")
  
  # Filter suspects using parameters
  suspects <- compTable %>%
    left_join(rtTable, by = "compound_id") %>%
    filter(chrom_method == chromatography || is.na(chrom_method)) %>%
    left_join(expTable, by = "compound_id") %>%
    left_join(paraTable, by = "parameter_id") %>%
    filter(polarity == pol) %>%
    filter(ce %in% ce_s) %>%
    filter(ces %in% ces_s) %>%
    filter(instrument %in% instr) %>%
    select(name, cas, mz, rt, adduct, isotope, experiment_id, compound_id) %>%
    dplyr::collect() %>%
    distinct()
  
  # Collect other tables
  exptbl <- tbl(db, "experiment") %>% collect()
  comptbl <- tbl(db, "compound") %>% collect()
  fragtbl <- tbl(db, "fragment") %>% collect()
  
  # Disconnect the CSL database
  DBI::dbDisconnect(db)
  
  # Filter suspects for compounds
  if (!is.null(compounds)) {
    suspects <- filter(suspects, name %in% compounds)
    if (nrow(suspects) == 0) {
      warning("Chosen compounds and settings not available in DB")
      return(NULL)
    } else {
      message("Searching compounds ", paste(unique(suspects$name), 
                                            collapse = "; "))
    }
  } 
  
  # Bundle results into a list
  csl_data <- list(
    suspects = suspects,
    exptbl = exptbl,
    comptbl = comptbl,
    fragtbl = fragtbl
  )
  
  return(csl_data)
}


# Process one measurement file with the selected library data
eval_samp <- function(pth, csl_data, params) {
  
  # Create xcmsRaw class
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
  all_csl_comps <- lapply(splitByComp, function(csl_comp) {
    eval_comp(csl_comp, raw_data, csl_data, params)
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
  row.names(all_csl_comps) <- NULL
  return(all_csl_comps)
}

# For each compound+adduct+isotopologue in the CSL, search raw data to find 
# MS2 spectra with matching precursor mass and RT
eval_comp <- function(csl_comp, raw_data, csl_data, params) {
  
  # Create data frame with precursor mass and RT data
  df <- data.frame(
    mz = raw_data@msnPrecursorMz,  # Precursor mass (Dalton)
    rt = raw_data@msnRt,            # Retention time (s)
    index = seq_along(raw_data@msnPrecursorMz)
  )

  # Match sample data by precursor mass and RT of csl compound within tolerances
  df <- filter_matching_rt_mz(df, csl_comp$rt, csl_comp$mz, params$mztolu, 
                              params$rtoffset, params$rttolm)
  
  if (nrow(df) == 0) 
    return(NULL)
  
  data_specs <- collect_spectra_from_data(df, raw_data, comp_name = csl_comp$name[1])
  
  # Remove df entries with NULL spectrum
  stopifnot(nrow(df) == length(data_specs))
  keep_data_spec <- !vapply(data_specs, is.null, logical(1))
  df <- df[keep_data_spec, ]
  data_specs <- data_specs[keep_data_spec]
  
  # Check if any spectrum remains
  if (is.null(data_specs) || length(data_specs) == 0) 
    return(NULL)
  
  csl_specs <- collect_spectra_from_csl(csl_data, csl_comp$experiment_id)
  
  if (is.null(csl_specs) || length(csl_specs) == 0) 
    return(NULL)
  
  # Compare MS2 spectra and add best match and similarity score.
  df <- filter_matching_ms2(df, data_specs, csl_specs, params)
  # Keep proposed CSL spectra with a similarity score over the set threshold
  df <- df[df$score >= params$threshold, , drop = FALSE]
  
  if (nrow(df) == 0) 
    return(NULL)
  
  # Search MS1 scans near the MS2 RT for a precursor ion matching the expected m/z
  # and add the m/z and intensity of the closest matching peak
  df <- search_ms1_for_precursor(df, raw_data, params, csl_comp)

  if (nrow(df) == 0) 
    return(NULL)
  
  df <- label_separate_peaks_by_rt(df, params, csl_comp)
  df <- add_comp_metadata(df, csl_comp, raw_data)
  row.names(df) <- NULL

  return(df)
}


# Helper functions eval_comp ####

#' Match sample data by precursor mass and RT of csl compound within tolerances
filter_matching_rt_mz <- function(df, rt, mz, mztolu, rtoffset, rttolm) {
  # Check if RT exists
  if (is.na(rt[1])) {
    # If RT is unavailable, filter only by precursor mass
    df <- df[abs(df[, "mz"] - mz[1]) <= mztolu, , drop = FALSE]
  } 
  else {
    # Filter by precursor mass and RT with added RT offset
    rt <- rt + rtoffset
    combination <- cbind(
      abs(df[, "mz"] - mz[1]) <= mztolu,
      abs(df[, "rt"] - rt[1] * 60) <= rttolm * 60 # RT expected to be in seconds
    )
    df <- df[apply(combination, 1, all), , drop = FALSE]
  }
  
  return(df)
}

collect_spectra_from_data <- function(df, raw_data, comp_name) {
  data_specs <- lapply(df[, "index"], xcms::getMsnScan, object = raw_data)
  data_specs <- lapply(data_specs, as.data.frame)
  
  # Add the matched compound name to each spectrum
  data_specs <- lapply(data_specs, function(x) {attr(x, "comp_name") <- comp_name; x})
  
  # Add the precursor mass to each spectrum
  data_specs <- Map(set_mz_attribute, data_specs, df[, "mz"])
  
  # Normalize each MS2 data spectrum 
  lapply(data_specs, clean_and_normalize_ms2)
}

set_mz_attribute <- function(spec, mz) {
  # Adds attribute to object
  attr(spec, "precursor_mz") <- mz
  return(spec)
}

collect_spectra_from_csl <- function(csl_data, experiment_ids) {
  exptt <- csl_data$exptbl %>% select(experiment_id, compound_id, mz)
  compt <- csl_data$comptbl %>% select(compound_id, cas, name)
  spectra <- filter(csl_data$fragtbl, experiment_id %in% experiment_ids) %>%
    select(mz, int, experiment_id) %>%
    left_join(exptt, by = "experiment_id") %>%
    left_join(compt, by = "compound_id")
  
  # Group spectra by unique experiment ID
  csl_specs <- split(spectra, spectra$experiment_id)
  
  # Add attributes to each CSL spectrum and rename columns
  csl_specs <- lapply(csl_specs, comp_metadata_as_attribute)
  
  # Normalize each MS2 CSL spectrum 
  csl_specs <- lapply(csl_specs, clean_and_normalize_ms2)
  
  # Remove NULL entries and check if any spectrum remains
  compact(csl_specs)
}

comp_metadata_as_attribute <- function(spec) {
  res <- spec[, c("mz.x", "int")]
  colnames(res) <- c("mz", "int")
  attr(res, "precursor_mz") <- spec$mz.y[1]
  attr(res, "comp_name") <- spec$name[1]
  attr(res, "cas") <- spec$cas[1]
  attr(res, "db_exp_ID") <- spec$experiment_id[1]
  res
}

#' Extraction of a MS2-Spectrum for a given ExperimentID from spectral library
#'
#' @param db Either a DBI connection object OR a csl_data list as returned by
#'   get_csl_data() (containing $exptbl, $comptbl, $fragtbl)
#' @param Exp.ID integer, a single experiment ID
#'
#' @return MS2 spectrum as a data.frame (mz, int)
#' @export
dbGetSpectrum <- function(db, Exp.ID) {
  stopifnot(length(Exp.ID) == 1)
  
  # Resolve tables from either a db connection or a csl_data list
  if (inherits(db, "DBIConnection")) {
    exptbl  <- DBI::dbReadTable(db, "experiment")
    comptbl <- DBI::dbReadTable(db, "compound")
    fragtbl <- DBI::dbReadTable(db, "fragment")
  } else if (is.list(db) && all(c("exptbl", "comptbl", "fragtbl") %in% names(db))) {
    exptbl  <- db$exptbl
    comptbl <- db$comptbl
    fragtbl <- db$fragtbl
  } else {
    stop("`db` must be a DBI connection or a csl_data list with exptbl/comptbl/fragtbl")
  }
  
  # Join tables
  exptt <- exptbl[, c("experiment_id", "compound_id", "mz")]
  compt <- comptbl[, c("compound_id", "cas", "name")]

  frags <- fragtbl[fragtbl$experiment_id == Exp.ID, c("mz", "int", "experiment_id")]
  
  # Merge experiment and compound metadata
  spectrum <- merge(frags, exptt, by = "experiment_id", all.x = TRUE)
  spectrum <- merge(spectrum, compt, by = "compound_id", all.x = TRUE)
  
  # Build the spectrum df
  spec <- comp_metadata_as_attribute(spectrum)
  
  return(spec)
}

#' Filter and normalize MS2 spectra by removing precursor peaks and noise.
#'
#' This function filters out peaks close to the precursor m/z and low-intensity 
#' noise, normalizes intensities to a maximum of 1, and retains relevant metadata.
#' Returns NULL if no valid peaks remain.
clean_and_normalize_ms2 <- function(spec) {
  # Retrieve attributes
  compname <- attr(spec, "comp_name")
  compmz <- attr(spec, "precursor_mz")
  expID <- expID <- tryCatch(attr(spec, "db_exp_ID"), error = function(e) NULL)
  
  # Remove peaks within 0.02 units of the precursor m/z 
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
filter_matching_ms2 <- function(df, data_specs, csl_specs, params){

  # Define function to compute a similarity score between two MS2 spectra
  custom_calc_ndp <- function(d_spec, csl_spec) {
    calc_ndp(d_spec, csl_spec, params$ndp_m, params$ndp_n, params$mztolu_ms2)
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
search_ms1_for_precursor <- function(df, raw_data, params, csl_comp){
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
    if (!any(mz_error <= params$mztolu_fine)) {
      df$mz_ok[i] <- FALSE
      next
    }
    
    # Store the m/z and intensity of the closest matching peak (lowest mz-error)
    df$real_mz[i] <- ms1[which.min(mz_error), "mz"]
    df$int_h[i] <- ms1[which.min(mz_error), "intensity"]
    df$mz_ok[i] <- TRUE
  }
  # Keep only rows with valid precursor matches
  return(df[df$mz_ok, , drop = FALSE])
}


#' Group and label separate chromatographic peaks (A-Z) based on retention time, 
#' keeping only the most intense from each group.
label_separate_peaks_by_rt <- function(df, params, csl_comp){
  # Case: more than one MS2 match — cluster peak groups by retention time
  if (nrow(df) > 1) {
    df <- df[order(df$rt), ]     # Sort by retention time
    peakNum <- 1                 # Start with first peak group
    df$peak <- LETTERS[peakNum]  # Label first row as "A"
    
    # Loop through rows (MS2 spectra) to separate into distinct peaks
    for (i in 2:nrow(df)) {
      diff <- df$rt[i] - df$rt[i - 1]
      
      # If RT difference exceeds threshold (par$rt_res), start a new peak group
      if (diff > params$rt_res * 60) {
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
  } else if (nrow(df) == 1) { # Case: only one peak found — label it "A" 
    df$peak <- "A"
  } 
  return(df)
}

add_comp_metadata <- function(df, csl_comp, raw_data) {
  df$rt_min <- round(df$rt / 60, 2)
  df$comp_id <- csl_comp$compound_id[1]
  df$samp <- basename(raw_data@filepath)
  df$comp_name <- csl_comp$name[1]
  df$comp_cas <- csl_comp$cas[1]
  df$mz_error_mDa <- round(abs(csl_comp$mz[1] - df$real_mz) * 1000, 3)
  df$adduct <- csl_comp$adduct[1]
  df$isotopologue <- csl_comp$isotope[1]
  if (is.na(csl_comp$rt[1])) {
    df$rt_error_min <- NA
  } else {
    df$rt_error_min <- round(abs(csl_comp$rt[1] - df$rt / 60), 2)
  }
  df
}


# Helper functions shiny app ####

#' @export
load_csl_data <- function(csl_path) {
  
  # Connection to the CSL database
  con <- DBI::dbConnect(RSQLite::SQLite(), csl_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  tables <- DBI::dbListTables(con)
  
  res <- list()
  
  if ("data_source" %in% tables) {
    res$data_source <- DBI::dbReadTable(con, "data_source")
  }
  
  if ("retention_time" %in% tables) {
    res$retention_time <- DBI::dbReadTable(con, "retention_time")
  }
  
  res
}


# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow