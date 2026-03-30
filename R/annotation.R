
#' Annotate alignment table based only on m/z and RT comparison
#'
#' Enables also annotation based only on m/z. Uses a csv file substance
#' library.
#'
#' @param alig alignment table from shiny app
#' @param compLibPath path to csv file
#' @param mztol tolerance in u
#' @param rttol tolerance in minutes
#'
#' @returns Annotation table (`data.frame`) which gives the annotations of the
#'   alignment table linked by `alignmentID`
#' @export
annotate_grouped_mz_rt <- function(alig, compLibPath, mztol, rttol) {
  intCols <- alig[, c(
    grep("^Int_", colnames(alig)),
    grep("alignmentID", colnames(alig))
  )]
  intCols <- as.data.frame(intCols)
  alig <- alig[, c("mean_mz", "mean_RT", "alignmentID")]
  alig <- as.data.frame(alig)
  alig$mean_RT <- round(alig$mean_RT / 60, 2)
  compLib <- read.csv(compLibPath)
  if (any(sapply(compLib$rt, is.na)) || any(sapply(compLib$mz, is.na))) {
    warning("Compound library format not recognized, check for NAs")
    return(NULL)
  }
  if (!all(c("mz", "name") %in% colnames(compLib)) ||
      nrow(compLib) == 0 || !inherits(compLib$mz, "numeric") ||
      !inherits(compLib$name, "character")) {
    warning("Compound library format not recognized")
    return(NULL)
  }
  
  if (any(duplicated(compLib$name))) {
    warning("There are duplicated compound names")
    return(NULL)
  }
  
  # Filter alig
  filter_alig <- function(compLibRow) {
    stopifnot(nrow(compLibRow) == 1)
    aligFilt <- alig[abs(alig$mean_mz - compLibRow$mz) <= mztol, ]
    if (nrow(aligFilt) > 0) {
      aligFilt$mzDB <- compLibRow$mz
    }
    
    if (nrow(aligFilt) > 0 && is.element("rt", colnames(compLibRow))) {
      stopifnot(inherits(compLibRow$rt, "numeric"))
      aligFilt <- aligFilt[abs(aligFilt$mean_RT - compLibRow$rt) <= rttol, ]
      if (nrow(aligFilt) > 0) {
        aligFilt$rtDB <- compLibRow$rt
      }
    }
    
    if (nrow(aligFilt) > 0) {
      aligFilt$name <- compLibRow$name
      # highest intensity of each row
      aligFilt$sample <- numeric(nrow(aligFilt))
      for (i in seq_len(nrow(aligFilt))) {
        intens <- intCols[
          intCols$alignmentID == aligFilt$alignmentID[i],
          grep("^Int_", colnames(intCols)),
          drop = T
        ]
        intensVec <- as.numeric(intens)
        highSamp <- as.numeric(
          stringr::str_match(
            names(intens[which.max(intensVec)]),
            "^Int_(\\d+)$"
          )[, 2]
        )
        aligFilt[i, "sample"] <- highSamp
      }
    }
    aligFilt
  }
  
  re <- do.call("rbind", by(compLib, compLib$name, filter_alig, simplify = FALSE))
  rownames(re) <- NULL
  
  if (nrow(re) > 0) {
    colnames(re) <- sub("mean_mz", "mzData", colnames(re))
    colnames(re) <- sub("mean_RT", "rtData", colnames(re))
    if (!is.element("rtDB", colnames(re))) {
      re$rtDB <- NA
    }
    re$CAS <- NA
    re$expID <- NA
    re$formula <- NA
    re$SMILES <- NA
    re$adduct <- NA
    re$isotope <- NA
    re$score <- NA
    re$db_available <- TRUE
    re$CE <- NA
    re$CES <- NA
  } else {
    message("Nothing found")
    return(NULL)
  }
  re
}


#' Annotate alignment table with compound names by spectral library screening
#'
#' @description
#' Using the SQLite spectral library (CSL), this function will pass over the alignment
#' table, open the measurement files stored in datenList and compare m/z,
#' retention time and MS² spectra found in the data and the library.
#' 
#'
#' @param sampleList data.frame of filenames and paths
#' @param peakListList list of data.frames of the peak-lists
#' @param alignmentTable Knonw as "grouped" in the app
#' @param db_path Path to the spectral library (CSL in SQLite format)
#' @param threshold_score dot product threshold score
#' @param mztolu m/z tolerance in Da
#' @param rttol RT tolerance in min.
#' @param polarity Polarity of the measurement
#' @param CE Collision energy
#' @param CES Collision energy spread
#' @param instrument Instruments allowed, must match the instruments listed in
#'   the library
#' @param chrom_method Chromatography name to choose the correct retention times
#'   from the library
#' @param ndp_m Peak intensity weighting factor for dot-product
#' @param ndp_n m/z weighting factor for dot-product
#' @param mztolu_ms2 m/z window for dot-product
#' @param rtoffset Offset the database to make it match your values (in min)
#' @param intCutData Cut off intensity under which fragments are ignored (in
#'   data spectra)
#' @param numcores Number of cores for parallelization (currently not used)
#'
#' @details This function operates either in conjunction with the non-target app
#'   but can also work in stand-alone scripts, see examples in ntsportal. If
#'   data files are not in memory (xcmsRaw objects), they will be temporary
#'   loaded. `datenListLocal` is a shallow copy. Therefore, changes to the
#'   underlying xcmsRaw objects affect the calling environment.
#'    
#'   ## Spectral library
#'   The spectral library is an SQLite file (see Jewell et al. 2019). 
#'   
#'   ## References
#'   Jewell, K. S., Kunkel, U., Ehlig, B., Thron, F., Schlüsener, M., Dietrich,
#'   C., . . . Ternes, T. A. (2019). Comparing mass, retention time and MS2
#'   spectra as criteria for the automated screening of small molecules in
#'   aqueous environmental samples analyzed by LC-QToF-MS/MS. Rapid
#'   Communications in Mass Spectrometry, 34, e8541. doi:10.1002/rcm.8541
#' 
#' @returns `data.frame` with the annotations for the alignment table referenced
#'   by the column "alignementID" (compound identifiers and results of the
#'   spectral comparison)
#'
#' @export
annotate_grouped <- function(sampleListLocal,
                             peakListList,
                             alignmentTable,
                             db_path, # can be .db or .yaml
                             threshold_score = 500,
                             mztolu = 0.005,
                             rttol = 1,
                             polarity = "pos",
                             CE = c(30, 40),
                             CES = c(0, 15),
                             instrument = c(
                               "LC-ESI-QTOF TripleTOF 5600 SCIEX",
                               "LC-ESI-QTOF TripleTOF 6600 SCIEX",
                               "LC-ESI-Orbitrap QExactive Thermo",
                               "LC-ESI-QTOF Agilent 6500 Series Q-TOF",
                               "LC-ESI-QTOF TripleTOF X500R SCIEX",
                               "LC-ESI-QTOF ZenoTOF 7600 SCIEX",
                               "LC-ESI-Orbitrap Exploris 240 Thermo"
                             ),
                             chrom_method =
                               "bfg_nts_rp1",
                             ndp_m = 2,
                             ndp_n = 1,
                             mztolu_ms2 = 0.015,
                             rtoffset = 0,
                             intCutData = 0,
                             numcores = 1,
                             datenListLocal,
                             expGroups = "bfg") {
  # open database connection or open custom database (.yaml)
  if (grepl("\\.yaml$", db_path)) {
    useCustom <- TRUE
    customdb <- yaml::read_yaml(db_path)
    getValue <- function(x, n) {
      if (!is.null(x[[n]])) x[[n]] else NA
    } # if values are missing, return NA
    allExps <- data.frame(
      experiment_id = NA,
      mz = sapply(customdb, getValue, n = "mz"),
      adduct = sapply(customdb, getValue, n = "adduct"),
      isotope = NA,
      CE = NA,
      CES = NA,
      CAS = NA,
      formula = sapply(customdb, getValue, n = "formula"),
      SMILES = NA,
      name = names(customdb),
      rt = sapply(customdb, getValue, n = "rt"),
      stringsAsFactors = FALSE
    )
  } else if (grepl("\\.db$", db_path)) {
    useCustom <- FALSE
    dbi <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    # produce table of all possible experiments
    # change all variable names to allow NSE. These names have _i after them:
    oldNames <- c("polarity", "CE", "CES", "instrument", "chrom_method")
    for (name in oldNames) {
      assign(paste0(name, "_i"), get(eval(name)))
    }
    rm(list = oldNames)
    
    expTable <- tbl(dbi, "experiment")
    allCe <- CE_i[1]:CE_i[2]
    allces <- CES_i[1]:CES_i[2]
    compTable <- tbl(dbi, "compound")
    
    
    # Spectra_DB or Label_DB, decided by the existence of the table "sample"
    if (DBI::dbExistsTable(dbi, "sample")) {
      sampTable <- tbl(dbi, "sample")
      paraTable <- tbl(dbi, "parameter") %>%
        filter(
          chrom_method == chrom_method_i, polarity == polarity_i, CE %in% allCe, 
          CES %in% allces, instrument %in% instrument_i
        )
      
      allExps <- expTable %>%
        inner_join(paraTable, by = "parameter_id") %>%
        select(experiment_id, compound_id, sample_id, mz, rt, adduct, isotope, CE, CES) %>%
        left_join(compTable, by = "compound_id") %>%
        left_join(sampTable, by = "sample_id") %>%
        select(
          -formula, -SMILES, -chem_list_id, -compound_id, -sample_id,
          -data_location, -date, -leaching_id, -extraction_id
        ) %>%
        dplyr::collect()
    } else {
      paraTable <- tbl(dbi, "parameter") %>%
        filter(
          polarity == polarity_i, ce %in% allCe, ces %in% allces,
          instrument %in% instrument_i
        )
      rtTable <- tbl(dbi, "retention_time") %>% filter(chrom_method == chrom_method_i)
      
      allExps <- expTable %>%
        inner_join(paraTable, by = "parameter_id") %>%
        select(experiment_id, compound_id, mz, adduct, isotope, ce, ces) %>%
        left_join(compTable, by = "compound_id") %>%
        left_join(rtTable, by = "compound_id") %>%
        select(-chrom_method, -retention_time_id, -compound_id) %>%
        dplyr::collect()
      
      # filter by experiment group (user source selection)
      expIdsSource <- tbl(dbi, "data_source") %>%
        filter(name %in% expGroups) %>%
        left_join(tbl(dbi, "experiment"), by = "data_source_id") %>%
        select(experiment_id) %>%
        collect() %>%
        .$experiment_id
      
      allExps <- allExps[allExps$experiment_id %in% expIdsSource, ]
    }
    
    allExps$rt <- allExps$rt + rtoffset
  } else {
    stop("Database must be .yaml or .db")
  }
  
  
  # first group the rows of the alignment table so that you check all the
  # features from one file at a time. that way you only need to load each file
  # once. Get highest intensity row for each row of alignmentTable
  
  intCols <- alignmentTable[, grep("^Int", colnames(alignmentTable))]
  ms2Cols <- alignmentTable[, grep("^ms2scan", colnames(alignmentTable))]
  getBestSamp <- function(rowInt, rowMs2) {
    if (all(rowMs2 == 0)) {
      return(NA)
    }
    x <- cbind(seq_along(rowInt), rowInt, rowMs2)
    x <- x[order(x[, 2], decreasing = TRUE), ]
    x <- x[x[, 3] != 0, , drop = FALSE]
    if (nrow(x) == 0) {
      return(NA)
    }
    x[1, 1]
  }
  intMax <- mapply(getBestSamp,
                   split(intCols, seq_len(nrow(alignmentTable))),
                   split(ms2Cols, seq_len(nrow(alignmentTable))),
                   SIMPLIFY = TRUE, USE.NAMES = FALSE
  )
  
  
  # reorder alignment table by max samples
  alignmentTable <- cbind(alignmentTable, intMax)
  alignmentTableDf <- as.data.frame(alignmentTable)
  
  # NAs are ignored, but that doesnt matter anyway
  alignmentTableDfList <- split(alignmentTableDf, alignmentTableDf$intMax)
  
  # preload tables for speed
  if (!useCustom) {
    fragTable <- tbl(dbi, "fragment") %>% collect()
    expTable <- tbl(dbi, "experiment") %>% collect()
    compTable <- tbl(dbi, "compound") %>% collect()
  }
  
  # Process list of alignmentTables
  hitsBySample <- parallel::mclapply(alignmentTableDfList, function(aligSamp) {
    sample_highest <- aligSamp[, "intMax"]
    sample_highest <- unique(sample_highest)
    stopifnot(length(sample_highest) == 1, !is.na(sample_highest), is.numeric(sample_highest))
    
    # check that file is loaded, if not, load file
    
    if (!sampleListLocal[sampleListLocal$ID == sample_highest, "RAM"]) {
      currentDataFile <- suppressMessages(xcms::xcmsRaw(
        datenListLocal[[sample_highest]]@filepath@.Data,
        includeMSn = TRUE
      ))
    } else {
      currentDataFile <- datenListLocal[[sample_highest]]
    }
    
    
    
    hitsForSamp <- lapply(seq_len(nrow(aligSamp)), function(row) {
      # get mz and RT of most intense peak, get DB experiments
      # remove cross-referencing in table, turn into long form for easier data manipulation
      # ignore unnecessary columns
      value <- subset(aligSamp, , -c(mean_mz, mean_RT, MS2Fit, Gruppe, alignmentID, intMax))
      value <- value[row, ]
      file <- as.numeric(stringr::str_match(names(value), "_(\\d+)$")[, 2])
      type <- stringr::str_match(names(value), "^(\\w+)_")[, 2]
      feature <- data.frame(type, file, value = as.numeric(value), stringsAsFactors = FALSE)
      
      # get sample MS2 scan number from most intense peak with an MS2
      hasMS2 <- feature[feature$type == "ms2scan" & feature$value != 0, "file"]
      feature <- feature[feature$file %in% hasMS2, ]
      ms2_scan_highest <- round(feature[feature$type == "ms2scan" &
                                          feature$file == sample_highest, "value"])
      # mz and rt of peak
      peakMz <- feature[feature$type == "mz" & feature$file == sample_highest, "value"]
      peakRt <- feature[feature$type == "RT" & feature$file == sample_highest, "value"] / 60
      
      # get available DB or custom spectra
      boolMz <- is.na(allExps$mz) | abs(peakMz - allExps$mz) <= mztolu # some compounds no mz
      boolRt <- is.na(allExps$rt) | abs(peakRt - allExps$rt) <= rttol # some compounds no RT
      viabExp <- allExps[boolMz & boolRt, ]
      if (nrow(viabExp) == 0) {
        return(NULL)
      } # if no db entry, no result
      
      if (useCustom) {
        # extract spectra from yaml
        db_spectra <- lapply(viabExp$name, function(x) {
          # three cases: only fragments, only neutral losses or both
          av <- names(customdb[[x]])
          if ("fragments" %in% av && "neutral_losses" %in% av) {
            spec <- customdb[[x]]$fragments
            nl <- customdb[[x]]$neutral_losses
            return(list(frags = spec, losses = nl))
          } else if ("fragments" %in% av) {
            spec <- customdb[[x]]$fragments
            if (is.list(spec)) {
              spec <- do.call("rbind", spec)
              colnames(spec) <- c("mz", "int")
            }
            return(spec)
          } else if ("neutral_losses" %in% av) {
            nl <- customdb[[x]]$neutral_losses
            return(nl)
          } else {
            stop("need fragments or neutral losses")
          }
        })
        
        viabExp$db_available <- TRUE
        stopifnot(all(!is.null(db_spectra)))
      } else {
        db_spectra <- lapply(
          viabExp$experiment_id,
          get_spectrum_preloaded,
          fragTable = fragTable,
          expTable = expTable,
          compTable = compTable
        )
        db_spectra <- lapply(db_spectra, clean_and_normalize_ms2)
        viabExp$db_available <- vapply(db_spectra, Negate(is.null), logical(1))
        viabExp <- viabExp[viabExp$db_available, ]
        if (nrow(viabExp) == 0) {
          return(NULL)
        } # check that spectra are available from DB
        compact <- function(x) {
          Filter(Negate(is.null), x)
        }
        db_spectra <- compact(db_spectra)
      }
      
      # get MS2 from file
      ms2spektrum <- xcms::getMsnScan(currentDataFile, ms2_scan_highest)
      ms2spektrum <- as.data.frame(ms2spektrum)
      attr(ms2spektrum, "comp_name") <- sprintf("Sample%i_%.4f_%.2f", sample_highest, 
                                                peakMz, peakRt)
      attr(ms2spektrum, "precursor_mz") <- peakMz
      colnames(ms2spektrum) <- c("mz", "int")
      ms2spektrumNorm <- clean_and_normalize_ms2(ms2spektrum)
      # cut off low intensity fragments if desired
      if (intCutData != 0) {
        ms2spektrumNorm <- ms2spektrumNorm[ms2spektrumNorm$int >= intCutData, ]
      }
      # if the db_spectra are just a vector of masses, do a simple search for the fragments (yes/no)
      
      # simple search of fragments
      simpleSearch <- function(d_spec, db_spec) {
        all(vapply(db_spec, function(x) any(abs(d_spec$mz - x) <= mztolu_ms2), logical(1)))
      }
      
      # compare with data spectrum with DB spectra
      custom_calc_ndp <- function(d_spec, db_spec) {
        ntsworkflow::calc_ndp(d_spec, db_spec, ndp_m = ndp_m, ndp_n = ndp_n, 
                              mztolu_ms2 = mztolu_ms2)
      }
      # prepare spec for neutral loss search, returns vector of all differences
      # spec must be normalized already!!
      # only get mass area within 200 Da of precursor
      get_nl <- function(spec) {
        mz <- spec[spec$mz > peakMz - 200, "mz"]
        mz <- append(mz, peakMz) # append precursor
        all_diff <- abs(outer(mz, mz, "-"))
        all_diff <- as.vector(all_diff[upper.tri(all_diff)])
        if (length(all_diff) == 0) {
          data.frame(mz = 0, int = 0)
        } else {
          data.frame(mz = all_diff, int = 1)
        }
      }
      
      # the type of search done depends on what was provided
      
      if (useCustom) {
        # which compounds are "normal" fragments, simple fragments, simple fragments and nl and just nl
        # normal comparison
        hasSpec <- sapply(db_spectra, is.matrix)
        normal <- which(hasSpec)
        viabExp[normal, "score"] <- vapply(db_spectra[normal], custom_calc_ndp,
                                           numeric(1), d_spec = ms2spektrumNorm)
        
        # okay now the complicated cases
        provided <- lapply(customdb[viabExp$name], names)
        hasFragments <- vapply(provided, is.element, logical(1), el = "fragments")
        hasNl <- vapply(provided, is.element, logical(1), el = "neutral_losses")
        # just simple frag
        justFrag <- which(hasFragments & !hasNl & !hasSpec)
        if (length(justFrag) != 0) {
          found <- vapply(db_spectra[justFrag], simpleSearch, logical(1), 
                          d_spec = ms2spektrumNorm)
          viabExp[justFrag, "score"] <- ifelse(found, 1000, 0)
        }
        # just nl
        justNl <- which(!hasFragments & hasNl & !hasSpec)
        if (length(justNl) != 0) {
          found <- vapply(db_spectra[justNl], simpleSearch, logical(1), 
                          d_spec = get_nl(ms2spektrumNorm))
          viabExp[justNl, "score"] <- ifelse(found, 1000, 0)
        }
        # combined search, frag and nl
        combi <- which(hasFragments & hasNl & !hasSpec)
        if (length(combi) != 0) {
          found_frags <- vapply(lapply(db_spectra[combi], "[[", "frags"), simpleSearch, 
                                logical(1), d_spec = ms2spektrumNorm)
          found_losses <- vapply(lapply(db_spectra[combi], "[[", "losses"), simpleSearch, 
                                 logical(1), d_spec = get_nl(ms2spektrumNorm))
          viabExp[combi, "score"] <- ifelse(found_frags & found_losses, 1000, 0)
        }
      } else { # for normal db search
        viabExp$score <- vapply(db_spectra, custom_calc_ndp, numeric(1), d_spec = ms2spektrumNorm)
      }
      
      stopifnot(all(!is.na(viabExp$score)))
      
      # are any above threshold?
      if (all(viabExp$score < threshold_score, na.rm = TRUE)) {
        return(NULL)
      }
      
      # get matching compounds from DB
      # for those above threshold save names and CAS in a table, link this with peaklist table
      # by "Zeile"
      viabExp <- viabExp[viabExp$score >= threshold_score, ]
      # for each compound found, select the highest score
      if (any(is.na(viabExp$name))) {
        viabExp$name[is.na(viabExp$name)] <- "unknown compound"
      }
      uniques <- by(viabExp, viabExp$name, function(r) r[which.max(r$score), ], simplify = TRUE)
      viabExp <- Reduce(rbind, uniques)
      viabExp <- viabExp[order(viabExp$score, decreasing = TRUE), ]
      viabExp$alignmentID <- aligSamp[row, "alignmentID"]
      viabExp$datenListVerwendet <- sample_highest
      viabExp$mzData <- peakMz
      viabExp$rtData <- round(peakRt, 2)
      
      viabExp
    })
    
    hitsForSamp <- do.call("rbind", hitsForSamp)
    hitsForSamp
  }, mc.preschedule = FALSE, mc.cores = 1) # at the moment does not work multicore, do not know why
  
  hits <- do.call("rbind", hitsBySample)
  
  # reformat the hits table
  if (is.null(hits)) {
    return(NULL)
  }
  
  if (!useCustom && DBI::dbExistsTable(dbi, "sample")) {
    hits <- hits[, c(
      "alignmentID", "mzData", "rtData", "samplename", "sample_type",
      "name", "CAS", "mz", "rt", "experiment_id", "adduct",
      "isotope", "score", "db_available", "CE", "CES", "enrichment", "enrichment_factor",
      "contact", "project", "datenListVerwendet"
    )]
    colnames(hits) <- c(
      "alignmentID", "mzData", "rtData", "samplename", "sample_type",
      "name", "CAS", "mzDB", "rtDB", "expID", "adduct",
      "isotope", "score", "db_available", "CE", "CES", "enrichment", "enrichment_factor",
      "contact", "project", "sample"
    )
  } else {
    hits <- hits[, c(
      "alignmentID", "mzData", "rtData", "name",
      "cas", "mz", "rt", "experiment_id", "formula", "smiles", "adduct",
      "isotope", "score", "db_available", "ce", "ces", "datenListVerwendet"
    )]
    colnames(hits) <- c(
      "alignmentID", "mzData", "rtData", "name",
      "CAS", "mzDB", "rtDB", "expID", "formula", "SMILES", "adduct", "isotope",
      "score", "db_available", "CE", "CES", "sample"
    )
  }
  hits$rtDB <- round(hits$rtDB, 2)
  hits$score <- round(hits$score)
  # reorder by alignemntID
  hits <- hits[order(hits$alignmentID), ]
  
  # close db connection
  if (!useCustom) {
    DBI::dbDisconnect(dbi)
  }
  
  hits
}


#' @export
#' @import dplyr
get_spectrum_preloaded <- function(fragTable, expTable, compTable, Exp.ID) {
  stopifnot(length(Exp.ID) == 1)
  
  # Selection of Spectrum by Exp.ID
  dbSpectrum <- fragTable %>% # Selection of database
    filter(experiment_id == !!Exp.ID) %>% # Selection of one Experiment
    select(mz, int)
  
  mz_i <- expTable %>%
    filter(experiment_id == !!Exp.ID) %>%
    .$mz
  
  meta_data <- id_to_name_preloaded(Exp.ID, expTable, compTable)
  
  attr(dbSpectrum, "precursor_mz") <- mz_i
  attr(dbSpectrum, "comp_name") <- meta_data$name
  attr(dbSpectrum, "CAS") <- meta_data$cas
  attr(dbSpectrum, "db_exp_ID") <- Exp.ID
  dbSpectrum
}

#' @export
#' @import dplyr
id_to_name_preloaded <- function(id, expTable, compTable) {
  name <- expTable %>%
    filter(experiment_id == !!id) %>%
    inner_join(compTable, by = "compound_id") %>%
    select(name, cas)
  name$exp_ID <- id
  name
}
