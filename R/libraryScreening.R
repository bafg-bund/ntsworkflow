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
#' @param db_path Path to spectral library
#' @param rttolm Retention time tolerance in minutes
#' @param mztolu m/z tolerance in Da (MS2 precursor mass)
#' @param mztolu_fine m/z tolerance in Da (second stage)
#' @param chromatography Chromatography name to choose the correct retention times from the library
#' @param pol Polarity used
#' @param CE_s Collision energy
#' @param CES_s Collision energy spread
#' @param instr Instruments allowed, must match the instruments listed in the library
#' @param ceunit Collision energy units
#' @param comparison How the spectra should be compared, at the moment only "dot_product" is allowed
#' @param threshold Threshold score for spectral comparison (for dot_product 1000 is max)
#' @param rt_res resolution of two peaks in chromatography
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
#'
#' ## References
#'   Jewell, K. S., Kunkel, U., Ehlig, B., Thron, F., Schlüsener, M., Dietrich,
#'   C., . . . Ternes, T. A. (2019). Comparing mass, retention time and MS2
#'   spectra as criteria for the automated screening of small molecules in
#'   aqueous environmental samples analyzed by LC-QToF-MS/MS. Rapid
#'   Communications in Mass Spectrometry, 34, e8541. doi:10.1002/rcm.8541
#'
#' @returns a `data.frame` with hits from the spectral library
#'
#' @seealso [annotate_grouped()]
#'
#' @export
#' @import dplyr
ms2_search <- function(data_path, db_path,
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
  
  # Check measurement file path
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
  
  # Load CSL database and filter relevant data
  csl_data <- filter_csl_data(db_path, chromatography, pol, CE_s, CES_s, instr)
  
  # Filter suspects by specified compound names
  csl_data$suspects <- filter_suspects_by_name(csl_data$suspects, compounds)
  
  # Process all measurement files
  all_samps <- lapply(data_path, function(pth) {
    eval_samp(pth, csl_data, par)
  })
  all_samps <- compact(all_samps)
  
  # Combine results
  if (length(all_samps) == 0) return(NULL)
  all_samps <- Reduce(rbind, all_samps)
  class(all_samps) <- c("sus_search_ms2", "sus_search", "data.frame")
  attr(all_samps, "processing_date") <- date()
  attr(all_samps, "db_path") <- db_path
  
  return(all_samps)
}


# Subfunctions ####

filter_csl_data <- function(db_path, chromatography, pol, CE_s, CES_s, instr){
  
  # Connect to CSL database
  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  
  # Get tables from CSL database
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
  
  # Disconnect CSL database
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
  
  # Check provided path
  if (inherits(pth, "character")) {
    raw_data <- xcms::xcmsRaw(pth, includeMSn = TRUE)
  } else if (inherits(pth, "xcmsRaw")) {
    raw_data <- pth
  } else {
    stop("data_path is not character or list of xcmsRaw")
  }
  
  # Add adduct and isotope information
  csl_data$suspects$comp_adduct_isot <- paste(csl_data$suspects$name, csl_data$suspects$adduct, csl_data$suspects$isotope, sep = "_")
  splitByComp <- split(csl_data$suspects, csl_data$suspects$comp_adduct_isot)
  
  # Do ?? for all available CSL compounds  # Todo: add description for "??"
  all_comps <- lapply(splitByComp, function(comp) {
    eval_comp(comp, raw_data, csl_data, par)
  })
  
  all_comps <- Filter(function(x) !inherits(x, "try-error"), all_comps)
  
  all_comps <- compact(all_comps)
  
  if (inherits(pth, "character")) {
    rm(raw_data)
  }
  
  if (length(all_comps) == 0) {
    return(NULL)
  }
  
  # Join all compounds
  all_comps <- do.call("rbind", all_comps)
  return(all_comps)
}


# For each compound+adduct+isotopologue in db, search raw data to find ms2 spectra with the
# correct precursor mass and correct rt
eval_comp <- function(compound, raw_data, csl_data, par) {
  inf <- data.frame(     # needs a new name
    mz = raw_data@msnPrecursorMz,
    rt = raw_data@msnRt,
    index = seq_along(raw_data@msnPrecursorMz)
  )
  
  # in case there is no rt, filter only by mass, otherwise, filter by mass
  # and rt
  if (is.na(compound$rt[1])) {
    inf <- inf[abs(inf[, "mz"] - compound$mz[1]) <= mztolu, , drop = FALSE]
  } else {
    # add rtoffset
    compound$rt <- compound$rt + par$rtoffset
    combination <- cbind(
      abs(inf[, "mz"] - compound$mz[1]) <= par$mztolu,
      abs(inf[, "rt"] - compound$rt[1] * 60) <= par$rttolm * 60
    )
    inf <- inf[apply(combination, 1, all), , drop = FALSE]
  }
  
  
  if (nrow(inf) == 0) {
    return(NULL)
  }
  data_specs <- lapply(inf[, "index"], xcms::getMsnScan, object = raw_data)
  # run dot product comparison on all spectra
  data_specs <- lapply(data_specs, as.data.frame)
  data_specs <- lapply(data_specs, function(x) {
    attr(x, "comp_name") <- compound$name[1]
    x
  })
  
  data_specs <- Map(set_mz_attribute, data_specs, inf[, "mz"])
  
  data_specs <- lapply(data_specs, normalize_ms2)
  # remove rows with NULL spectrum
  stopifnot(nrow(inf) == length(data_specs))
  inf$keep <- TRUE
  for (i in seq_along(data_specs)) {
    if (is.null(data_specs[[i]])) {
      inf$keep[i] <- FALSE
    }
  }
  inf <- inf[inf$keep, ]
  inf$keep <- NULL
  
  data_specs <- compact(data_specs)
  if (is.null(data_specs) || length(data_specs) == 0) {
    return(NULL)
  }
  stopifnot(nrow(inf) == length(data_specs))
  
  # Collect spectra from database
  # parallel cores can not read db
  exptt <- csl_data$exptbl %>% select(experiment_id, compound_id, mz)
  compt <- csl_data$comptbl %>% select(compound_id, CAS, name)
  spectra <- filter(csl_data$fragtbl, experiment_id %in% !!compound$experiment_id) %>%
    select(mz, int, experiment_id) %>%
    left_join(exptt, by = "experiment_id") %>%
    left_join(compt, by = "compound_id")
  
  db_specs <- split(spectra, spectra$experiment_id)
  db_specs <- lapply(db_specs, function(sp) {
    res <- sp[, c("mz.x", "int")]
    colnames(res) <- c("mz", "int")
    attr(res, "precursor_mz") <- sp$mz.y[1]
    attr(res, "comp_name") <- sp$name[1]
    attr(res, "CAS") <- sp$CAS[1]
    attr(res, "db_exp_ID") <- sp$experiment_id[1]
    res
  })
  
  db_specs <- lapply(db_specs, normalize_ms2)
  db_specs <- compact(db_specs)
  if (is.null(db_specs) || length(db_specs) == 0) {
    return(NULL)
  }
  
  # TODO: Add description (it extends "inf" by score and expID)
  inf <- filter_matching_ms2_spectra(inf, data_specs, db_specs, par)
  
  
  # if any of the proposed db spectra give score over threshold, keep this
  # feature
  inf <- inf[inf$score >= par$threshold, , drop = FALSE]
  # inf <- inf[apply(scorA, 1, function(x) any(x >= threshold)), , drop = FALSE]
  if (nrow(inf) == 0) {
    return(NULL)
  }
  
  # TODO: Add description (it extends "inf" by mz_ok, real_mz, int_h)
  inf <- search_ms1_for_precursor(inf, raw_data, par, compound)
  
  # TODO: Add description (it extends "inf" by peak)
  inf <- label_separate_peaks(inf, par, compound)
  
  # Add metadata
  inf$rt_min <- round(inf$rt / 60, 2)
  inf$comp_id <- compound$compound_id[1]
  inf$samp <- basename(raw_data@filepath)
  inf$comp_name <- compound$name[1]
  inf$comp_CAS <- compound$CAS[1]
  inf$mz_error_mDa <- round(abs(compound$mz[1] - inf$real_mz) * 1000, 3)
  inf$adduct <- compound$adduct[1]
  inf$isotopologue <- compound$isotope[1]
  if (is.na(compound$rt[1])) {
    inf$rt_error_min <- NA
  } else {
    inf$rt_error_min <- round(abs(compound$rt[1] - inf$rt / 60), 2)
  }
  if (inherits(inf, "try-error")) {
    return(NULL)
  }
  return(inf)
}


# Helper functions eval_comp ####

set_mz_attribute <- function(spec, mz) {
  # Adds attribute to object
  attr(spec, "precursor_mz") <- mz
  spec
}


normalize_ms2 <- function(x) {
  # remove precursor and noise
  compname <- attr(x, "comp_name")
  compmz <- attr(x, "precursor_mz")
  expID <- try(attr(x, "db_exp_ID"))
  # Remove precursor otherwise essentially MS1 checked again
  x <- x[x$mz < attr(x, "precursor_mz") - 0.02, ]
  x <- x[x$int >= 0.1, ]
  if (length(x) == 0 || nrow(x) == 0) {
    return(NULL)
  }
  spec <- data.frame(mz = x$mz, int = x$int / max(x$int))
  attr(spec, "precursor_mz") <- compmz
  attr(spec, "comp_name") <- compname
  if (is.numeric(expID)) {
    attr(spec, "db_exp_ID") <- expID
  }
  spec
}


filter_matching_ms2_spectra <- function(inf, data_specs, db_specs, par){
  # Todo: add function description (one sentence)
  
  
  # Run the function calc_ndp with current settings
  custom_calc_ndp <- function(d_spec, db_spec) {
    calc_ndp(d_spec, db_spec, par$ndp_m, par$ndp_n, par$mztolu_ms2)
  }
  # perform comparison, each db_spec with each data_spec, result as matrix
  scorA <- outer(data_specs, db_specs, Vectorize(custom_calc_ndp))
  
  # record maximum score for each data_spec
  inf$score <- apply(scorA, 1, max, na.rm = TRUE)
  
  # record exp ID used in each case
  allExpIds <- sapply(db_specs, attr, which = "db_exp_ID")
  inf$expID <- allExpIds[apply(scorA, 1, which.max)]
  
  return(inf)
}


#' Calculate Normalized Dot-Product
#'
#' @param data_spec Spectrum from data
#' @param stan_spec Spectrum from database
#' @param ndp_m Peak intensity weighting factor
#' @param ndp_n m/z weighting factor
#' @param mztolu m/z tolerance for binning peaks (u)
#'
#' @details See MetCirc::NDP for more information
#'
#' @returns A numeric from 0 to 1000 with 1000 a perfect match
#' @export
#'
calc_ndp <- function(d_spec, db_spec, ndp_m = 2, ndp_n = 1, mztolu_ms2 = 0.015) {
  ar <- which(abs(outer(d_spec[, 1], db_spec[, 1], "-")) <= mztolu_ms2,
              arr.ind = TRUE
  ) # find matching masses
  if (nrow(ar) == 0) {
    return(0)
  }
  
  m <- cbind(d_spec[, 1][ar[, 1]], db_spec[, 1][ar[, 2]]) # extract matching
  
  d_int <- d_spec[, 2][ar[, 1]]
  db_int <- db_spec[, 2][ar[, 2]]
  masses <- rowMeans(m)
  # extract non-matching
  d_nonmatching <- d_spec[-ar[, 1], ]
  db_nonmatching <- db_spec[-ar[, 2], ]
  
  d_int <- append(d_int, d_nonmatching[, 2])
  masses <- append(masses, d_nonmatching[, 1])
  db_int <- append(db_int, rep(0, nrow(d_nonmatching)))
  
  db_int <- append(db_int, db_nonmatching[, 2])
  masses <- append(masses, db_nonmatching[, 1])
  d_int <- append(d_int, rep(0, nrow(db_nonmatching)))
  
  WS1 <- d_int^ndp_m * masses^ndp_n
  WS2 <- db_int^ndp_m * masses^ndp_n
  
  r_ndp <- (sum(WS1 * WS2))^2 / (sum(WS1^2) * sum(WS2^2))
  
  # if either of the two spectra only have 1 fragment, then the most intense
  # fragment in both spectra must be the same mass, otherwise return 0
  if (nrow(d_spec) == 1 || nrow(db_spec) == 1) {
    mzD <- d_spec[which.max(d_spec[, 2]), 1]
    mzS <- db_spec[which.max(db_spec[, 2]), 1]
    if (abs(mzD - mzS) > mztolu_ms2) {
      r_ndp <- 0
    }
  }
  
  r_ndp * 1000
}



search_ms1_for_precursor <- function(inf, raw_data, par, compound){
  # TODO: Add function description (one sentence)
  
  # get the real mz and intensity from raw data, check that mz matches 
  # TODO:
  # here you can also include an isotope check (calculate match to
  # theoretical)
  inf$mz_ok <- NA
  inf$real_mz <- NA
  inf$int_h <- NA
  for (i in 1:nrow(inf)) {
    ind <- which.min(abs(raw_data@scantime - inf$rt[i]))
    ms1 <- xcms::getScan(
      raw_data, ind,
      c(inf$mz[i] - 0.05, inf$mz[i] + 0.05)
    )
    # widen search if nothing found within 0.1 Da
    if (nrow(ms1) == 0) {
      ms1 <- xcms::getScan(
        raw_data, ind,
        c(inf$mz[i] - 0.5, inf$mz[i] + 0.5)
      )
    }
    
    # if still nothing found then error, skip to next
    if (nrow(ms1) == 0) {
      warning(sprintf(
        "Error: no precursor found for mz %.4f @ %.2f min",
        inf$mz[i], inf$rt[i] / 60
      ))
      inf$mz_ok[i] <- FALSE
      next
    }
    
    # first check that mz fits, if no mz fits, skip to next
    mz_error <- abs(ms1[, "mz"] - compound$mz[1])
    if (!any(mz_error <= par$mztolu_fine)) {
      inf$mz_ok[i] <- FALSE
      next
    }
    
    # Choose intensity with the lowest mz-error
    inf$real_mz[i] <- ms1[which.min(mz_error), "mz"]
    inf$int_h[i] <- ms1[which.min(mz_error), "intensity"]
    inf$mz_ok[i] <- TRUE
  }
  
  # if none of the mz matches, return NULL
  if (!any(inf$mz_ok)) {
    return(NULL)
  }
  # keep only matching masses
  inf <- inf[inf$mz_ok, , drop = FALSE]
  
  
  return(inf)
}


label_separate_peaks <- function(inf, par, compound){
  # TODO: Add function description (one sentence)
  
  
  # Of all peaks which match, cluster according to retention time, find
  # the most intense in each cluster and keep these only
  # collect all matching peaks and name
  # these A, B etc. if there is only one match, this becomes peak "A"
  if (nrow(inf) > 1) {
    inf <- inf[order(inf$rt), ]
    peakNum <- 1
    inf$peak <- LETTERS[peakNum]
    for (i in 2:nrow(inf)) {
      diff <- inf$rt[i] - inf$rt[i - 1]
      if (diff > par$rt_res * 60) {
        peakNum <- peakNum + 1
        inf[i:nrow(inf), "peak"] <- LETTERS[peakNum]
      }
    }
    good <- by(inf, inf$peak,
               function(x) x[which.max(x$int_h), ],
               simplify = FALSE
    )
    if (length(good) == 1) {
      inf <- good[[1]]
    } else {
      inf <- Reduce(rbind, good)
    }
  } else if (nrow(inf) == 1) {
    inf$peak <- "A"
  } else {
    stop(sprintf(
      "error in compound %s during collection of plausible peaks",
      compound$name[1]
    ))
  }
  
  return(inf)
}


# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow