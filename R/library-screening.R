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


# Library screening ####

# Includes functions for library screening used in conjunction with Report class
# (dbas) and the non-target app (peak picking).


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
                               "LC-ESI-Orbitrap QExactive"
                             ),
                             chrom_method =
                               "dx.doi.org/10.1016/j.chroma.2015.11.014",
                             ndp_m = 2,
                             ndp_n = 1,
                             mztolu_ms2 = 0.015,
                             rtoffset = 0,
                             intCutData = 0,
                             numcores = 1,
                             datenListLocal,
                             expGroups = "BfG") {
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
          polarity == polarity_i, CE %in% allCe, CES %in% allces,
          instrument %in% instrument_i
        )
      rtTable <- tbl(dbi, "retention_time") %>% filter(chrom_method == chrom_method_i)

      allExps <- expTable %>%
        inner_join(paraTable, by = "parameter_id") %>%
        select(experiment_id, compound_id, mz, adduct, isotope, CE, CES) %>%
        left_join(compTable, by = "compound_id") %>%
        left_join(rtTable, by = "compound_id") %>%
        select(-chrom_method, -ret_time_id, -chem_list_id, -compound_id) %>%
        dplyr::collect()

      # filter by experiment group (user source selection)
      expIdsSource <- tbl(dbi, "experimentGroup") %>%
        filter(name %in% expGroups) %>%
        left_join(tbl(dbi, "expGroupExp"), by = "experimentGroup_id") %>%
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
    compTable <- tbl(dbi, "compound") %>%
      select(-chem_list_id) %>%
      collect()
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
          ntsworkflow::get_spectrum_preloaded,
          fragTable = fragTable,
          expTable = expTable,
          compTable = compTable
        )
        db_spectra <- lapply(db_spectra, ntsworkflow::normalizeMs2)
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
      ms2spektrumNorm <- ntsworkflow::normalizeMs2(ms2spektrum)
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
      "CAS", "mz", "rt", "experiment_id", "formula", "SMILES", "adduct",
      "isotope", "score", "db_available", "CE", "CES", "datenListVerwendet"
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


#' Library screening without prior peak-picking
#' 
#' @description
#' A direct library screening approach using measurement files and the CSL 
#' spectral libary in SQLite format.
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
                       rttolm = 1, mztolu = 0.5,
                       mztolu_fine = 0.005,
                       chromatography =
                         "dx.doi.org/10.1016/j.chroma.2015.11.014",
                       pol = "pos", CE_s = 30:40, CES_s = 0:15,
                       instr = "LC-ESI-QTOF TripleTOF 5600 SCIEX",
                       ceunit = c("V", "eV"),
                       comparison = "dot_product", threshold = 400,
                       rt_res = 1.5, rtoffset = 0,
                       ndp_m = 2, ndp_n = 1, mztolu_ms2 = 0.015,
                       compounds = NULL) {
  stopifnot(inherits(data_path, "character") || inherits(data_path, "list"))
  if (inherits(data_path, "character")) {
    data_path <- normalizePath(data_path)
  }
  #browser()
  if (db_path == "Z:\\G\\G2\\HRMS\\Spektrendatenbank\\sqlite\\MS2_db_v7.db") {
    stop("Copy database to the local harddrive, do not use copy on Z")
  }

  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # get list of all compounds in db with rt
  expTable <- tbl(db, "experiment")
  paraTable <- tbl(db, "parameter")
  rtTable <- tbl(db, "retention_time")
  compTable <- tbl(db, "compound")
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

  exptbl <- tbl(db, "experiment") %>% collect()
  comptbl <- tbl(db, "compound") %>%
    select(-chem_list_id) %>%
    collect()
  fragtbl <- tbl(db, "fragment") %>% collect()

  DBI::dbDisconnect(db)

  # if searching for specific compounds, check that these are present
  if (!is.null(compounds)) {
    suspects <- filter(suspects, name %in% compounds)
    if (nrow(suspects) == 0) {
      warning("Chosen compounds and settings not available in DB")
      return(NULL)
    } else {
      message("Searching compounds ", paste(unique(suspects$name), collapse = "; "))
    }
  }

  # evaluate each sample
  eval_samp <- function(pth, exptbl, comptbl, fragtbl) {
    if (inherits(pth, "character")) {
      raw_data <- xcms::xcmsRaw(pth, includeMSn = TRUE)
    } else if (inherits(pth, "xcmsRaw")) {
      raw_data <- pth
    } else {
      stop("data_path is not character or list of xcmsRaw")
    }

    # For each compound+adduct+isotopologue in db, search raw data to find ms2 spectra with the
    # correct precursor mass and correct rt
    eval_comp <- function(compound, exptbl, comptbl, fragtbl) {
      inf <- data.frame(
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
        compound$rt <- compound$rt + rtoffset
        combination <- cbind(
          abs(inf[, "mz"] - compound$mz[1]) <= mztolu,
          abs(inf[, "rt"] - compound$rt[1] * 60) <= rttolm * 60
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

      set_att <- function(spec, mz) {
        attr(spec, "precursor_mz") <- mz
        spec
      }
      data_specs <- Map(set_att, data_specs, inf[, "mz"])

      data_specs <- lapply(data_specs, normalizeMs2)
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
      exptt <- exptbl %>% select(experiment_id, compound_id, mz)
      compt <- comptbl %>% select(compound_id, CAS, name)
      spectra <- filter(fragtbl, experiment_id %in% !!compound$experiment_id) %>%
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

      db_specs <- lapply(db_specs, normalizeMs2)
      db_specs <- compact(db_specs)
      if (is.null(db_specs) || length(db_specs) == 0) {
        return(NULL)
      }

      # prepare calc_ndp with current settings

      custom_calc_ndp <- function(d_spec, db_spec) {
        calc_ndp(d_spec, db_spec, ndp_m = ndp_m, ndp_n = ndp_n, mztolu_ms2 = mztolu_ms2)
      }
      # perform comparison, each db_spec with each data_spec, result as matrix
      scorA <- outer(data_specs, db_specs, Vectorize(custom_calc_ndp))


      # record maximum score for each data_spec
      inf$score <- apply(scorA, 1, max, na.rm = TRUE)

      # record exp ID used in each case
      allExpIds <- sapply(db_specs, attr, which = "db_exp_ID")
      inf$expID <- allExpIds[apply(scorA, 1, which.max)]

      # if any of the proposed db spectra give score over treshold, keep this
      # feature
      inf <- inf[inf$score >= threshold, , drop = FALSE]
      # inf <- inf[apply(scorA, 1, function(x) any(x >= threshold)), , drop = FALSE]
      if (nrow(inf) == 0) {
        return(NULL)
      }

      # get the real mz and intensity from raw data, check that mz matches TODO:
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
        if (!any(mz_error <= mztolu_fine)) {
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
          if (diff > rt_res * 60) {
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
      inf
    }
    suspects$comp_adduct_isot <- paste(suspects$name, suspects$adduct, suspects$isotope, sep = "_")
    splitByComp <- split(suspects, suspects$comp_adduct_isot)

    all_comps <- lapply(
      splitByComp,
      eval_comp,
      exptbl = exptbl,
      comptbl = comptbl,
      fragtbl = fragtbl
    )


    all_comps <- Filter(function(x) !inherits(x, "try-error"), all_comps)

    all_comps <- compact(all_comps)

    if (inherits(pth, "character")) {
      rm(raw_data)
    }

    if (length(all_comps) == 0) {
      return(NULL)
    }
    # join all compounds together
    # all_comps <- data.table::rbindlist(all_comps)
    all_comps <- do.call("rbind", all_comps)
    all_comps
  }
  all_samps <- lapply(data_path, eval_samp, exptbl = exptbl, comptbl = comptbl, 
                      fragtbl = fragtbl)
  all_samps <- compact(all_samps)
  if (length(all_samps) == 0) {
    message("no compounds found in sample(s)")
    return(NULL)
  }
  # join all samples together
  all_samps <- Reduce(rbind, all_samps)
  class(all_samps) <- c("sus_search_ms2", "sus_search", "data.frame")
  attr(all_samps, "processing_date") <- date()
  attr(all_samps, "db_path") <- db_path


  all_samps
}


# Raw data access functions ####

#' Extracting Spetra from mzML, mzXML data
#'
#' @param data.mzR object resulting from mzR::openMSfile
#' @param mass m/z of MS1 peak
#' @param rt retention time of MS1 peak in min
#' @param p_mztolu mz tol in u
#' @param p_rttol rt tol in min
#'
#' @return spectrum as data.frame (headers: mz, int), if no MS2 is found, returns NULL
#' @export
#' @import tidyr
#' @import dplyr
#' @import mzR
#'
get_MS2spectrum <- function(data.mzR, mass, rt, p_mztolu = 0.5, p_rttol = 0.2) {
  delta.mass2 <- p_mztolu
  rt.mzML.tol <- p_rttol

  # extract hrms from mzR data
  df.hrms <- mzR::header(data.mzR) %>% # header()
    tbl_df() %>% #
    mutate(., retentionTime = retentionTime / 60) # RT [s] <-> [min]

  # isolate the right ms2 experiment
  ms2.exp <- df.hrms %>%
    filter(., msLevel == 2 & precursorMZ != 0) %>%
    # only MS2-experiments & no empty MS2-experiments
    filter(., abs(precursorMZ - mass) <= delta.mass2) # suitable mass range

  if (nrow(ms2.exp) != 0) {
    ms2.exp2 <- ms2.exp %>%
      filter(abs(retentionTime - rt) <= rt.mzML.tol) # suitable rt range
    if (nrow(ms2.exp2) != 0) {
      ms2.exp3 <- ms2.exp2 %>%
        mutate(., delta.RT = abs(retentionTime - rt)) %>% # min(delta rt)
        arrange(., delta.RT) %>%
        slice(., 1) %>% # take the spectrum closest to rt
        select(., acquisitionNum) %>% # Spectrum ID
        as.numeric()

      # extract the spectrum
      spectrum2 <- data.mzR %>%
        mzR::peaks(., ms2.exp3) %>% # peaks() + Spectrum ID info -> spectrum
        as.data.frame() %>%
        tbl_df() %>%
        rename(., mz = V1, int = V2)
      if (nrow(spectrum2) == 0) {
        spectrum <- t(as.matrix(c(0, 0)))
        colnames(spectrum) <- c("mz", "int")
      } else {
        spectrum <- spectrum2
      }
    } else {
      spectrum <- t(as.matrix(c(0, 0)))
      colnames(spectrum) <- c("mz", "int")
    }
  } else {
    spectrum <- t(as.matrix(c(0, 0)))
    colnames(spectrum) <- c("mz", "int")
  }

  spectrum <- as.data.frame(spectrum)

  if (spectrum$mz[1] == 0) {
    return(NULL)
  }
  attr(spectrum, "precursor_mz") <- mass
  attr(spectrum, "rt") <- rt
  spectrum
}


# Database access functions ####

id_to_name <- function(id, db) {
  compTbl <- tbl(db, "compound")
  name <- tbl(db, "experiment") %>%
    filter(experiment_id == id) %>%
    inner_join(compTbl, by = "compound_id") %>%
    select(name, CAS) %>%
    dplyr::collect()
  name$exp_ID <- id
  name
}

#' Extraction of a MS2-Spectrum for a given ExperimentID from spectral library
#'
#' @param db database connection object using dplyr::src_sqlite
#' @param Exp.ID integer of experiment ID (only 1)
#'
#' @return MS2-Spectrum as data.frame (mz, int)
#' @export
#' @import dplyr
dbGetSpectrum <- function(db, Exp.ID) {
  stopifnot(length(Exp.ID) == 1)

  # Selection of Spectrum by Exp.ID
  dbSpectrum <- db %>% # Selection of database
    tbl(., "fragment") %>% # Selection of Fragment-Table (id, mz, int, experiment_id)
    filter(., experiment_id == Exp.ID) %>% # Selection of one Experiment
    select(., mz, int) %>% # mz+int- Table
    dplyr::collect()

  mz_i <- db %>%
    tbl("experiment") %>%
    filter(experiment_id == Exp.ID) %>%
    dplyr::collect() %>%
    .$mz

  meta_data <- id_to_name(Exp.ID, db)

  attr(dbSpectrum, "precursor_mz") <- mz_i
  attr(dbSpectrum, "comp_name") <- meta_data$name
  attr(dbSpectrum, "CAS") <- meta_data$CAS
  attr(dbSpectrum, "db_exp_ID") <- Exp.ID
  return(dbSpectrum)
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

  meta_data <- ntsworkflow::id_to_name_preloaded(Exp.ID, expTable, compTable)

  attr(dbSpectrum, "precursor_mz") <- mz_i
  attr(dbSpectrum, "comp_name") <- meta_data$name
  attr(dbSpectrum, "CAS") <- meta_data$CAS
  attr(dbSpectrum, "db_exp_ID") <- Exp.ID
  dbSpectrum
}

#' @export
#' @import dplyr
id_to_name_preloaded <- function(id, expTable, compTable) {
  name <- expTable %>%
    filter(experiment_id == !!id) %>%
    inner_join(compTable, by = "compound_id") %>%
    select(name, CAS)
  name$exp_ID <- id
  name
}


# Spectral comparison functions ####

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
  # fragement in both spectra must be the same mass, otherwise return 0
  if (nrow(d_spec) == 1 || nrow(db_spec) == 1) {
    mzD <- d_spec[which.max(d_spec[, 2]), 1]
    mzS <- db_spec[which.max(db_spec[, 2]), 1]
    if (abs(mzD - mzS) > mztolu_ms2) {
      r_ndp <- 0
    }
  }

  r_ndp * 1000
}


# Misc and helper functions ####


normalizeMs2 <- function(x) {
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
