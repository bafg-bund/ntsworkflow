
#' A reference class for database screening.
#'
#' This is a class for conducting non-target screening data evaluation based on database searching by m/z, retention time and MS2 spectral comparison. View the vignette "Database screening with Report in ntsworkflow" for a step-by-step guide on how to use this function.
#'
#'
#' @details Initialize a Report by calling Report$new() this starts a new Report with the default settings
#' (see settings field, below). Settings can be changed by calling the \code{changeSettings} method.
#' Load files by calling the \code{addRawFiles} method and process files with the \code{process_all}
#' method. After processing the Report can be saved for viewing in the visuallization tool (susS app).
#'
#' @field rawFiles A character vector with file paths in the order in which they are to be displayed.
#'
#' @field settings A list with settings for all processing steps.
#' The settings list contains the following fields and defaults:
#' \code{list(db_path = NULL, rttolm = 1, mztolu = 0.05, mztolu_fine = 0.005, chromatography =
#' "dx.doi.org/10.1016/j.chroma.2015.11.014", pol = "pos", CE_s = 30:40, CES_s = 0:15, instr =
#' "LC-ESI-QTOF TripleTOF 5600 SCIEX", ceunit = c("V", "eV"), comparison = "dot_product", threshold
#' = 400, rt_res = 1.5, EIC_extraction = 0.05, baseline_noise_MS1 = 0.5, sn = 3, rtoffset = 0,
#' IS_rtoffset = 0, ISrttolm = 1, blank_int_factor = 3, rtTolReinteg = 1, mzTolReinteg = 0.005,
#' ndp_m = 0.4, ndp_n = 1, mztolu_ms2 = 0.015, area_threshold = 1, height_threshold = 1,
#' use_int_threshold = "area"}.
#' \code{db_path}: path to sqlite spectral-database (SDB).
#' \code{rttolm}: retention time tolerance for the suspect search.
#' \code{mztolu}: m/z tolerance for the spectrum extraction (MS2 precursor mass).
#' \code{mztolu_fine}: m/z tolerance for precursor mass in MS1.
#' \code{chromatography} Is the chromatographic method allowed from the SDB, use DOI of paper in
#' which method is described.
#' \code{pol}: polarity.
#' \code{CE_s}: allowed collision energies from SDB.
#' \code{CES_s}: allowed collision energy spreads from SDB.
#' \code{instr}: allowed instruments from SDB.
#' \code{ceunit}: allowed collision energy unit from SDB.
#' \code{comparison}: algorithm to use for spectral comparison, currently only the default is allowed.
#' \code{threshold}: score of spectral comparison under which results are rejected automatically,
#' 1000=perfect match.
#' \code{rt_res}: Resolution of chrom. peaks (see \code{\link{ms2_search}}) in min.
#' Peaks with RT difference > \code{rt_res} are considered to be from different substances.
#' \code{EIC_extraction}: Extraction width to produce EIC, affects integration.
#' \code{baseline_noise_MS1}: Signals under this intensity are ignored and also not
#' included in the saved spectra.
#' \code{sn}: signal-to-noise ratio limit for peak integration.
#' \code{rtoffset}: Retention time offset between samples and database.
#' \code{ISrttolm}: RT tolerance for IS peak finding and integration.
#' \code{blank_int_factor}: Intensity factor to remove compounds found in blank. Peak in sample
#' must be within \code{blank_int_factor}: * intensity in blank.
#' \code{rtTolReinteg}: is the retention time tolerance for reintegration.
#' \code{mzTolReinteg}: is the m/z tolerance for reintegration.
#' \code{ndp_m}: Peak intensity weighting factor for dot-product,
#' \code{ndp_n}: m/z weighting factor for dot-product,
#' \code{mztolu_ms2}: m/z tolerance for dot-product fragment mass binning,
#' \code{area_threshold}: Peak area intensity threshold.
#' \code{height_threshold}: Peak height intensity threshold.
#' \code{use_int_threshold}: can be either "area", "height", or "none"
#' \code{peaksPerPeak}:
#' \code{mustFindChromPeak}: default FALSE, if TRUE, peaks without area are deleted
#'
#' @field rawFilesCompl \code{data.frame} with processed files and date of processing.
#' @field peakList \code{data.frame} with all suspect search results.
#' @field currentPeakID numeric to keep track of peakIDs, do not change.
#' @field MS1 \code{data.frame} holding all MS1 spectra.
#' @field EIC \code{data.frame} holding all chromatograms.
#' @field MS2 \code{data.frame} holding all MS2 spectra.
#' @field IS \code{data.frame} with list of IS to process, this has to be imported from a csv file
#' (comma sep) using this method \code{addIS()}. The csv file should have the columns name, formula,
#' rt, adduct. Where formula is in the form e.g. for CBZ: C14 13CH12N 15NO, and adduct is in the
#' form [M+H]+ or [M]+.
#' @field ISresults \code{data.frame} with results of IS processings.
#' @field falsePos \code{data.frame} recording which false positives should be deleted in which files
#' @field integRes
#'
#' @import parallel
#' @import shiny
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @export Report
#' @exportClass Report
Report <- setRefClass(
  "Report",
  fields = list(
    rawFiles = "character",
    rawData = "list",
    settings = "list",
    rawFilesCompl = "data.frame",
    peakList = "data.frame",
    currentPeakID = "numeric",
    currentISpeakID = "numeric",
    MS1 = "data.frame",
    EIC = "data.frame",
    MS2 = "data.frame",
    IS = "data.frame", # table to keep track of internal standards to be evaluated.
    ISresults = "data.frame",
    falsePos = "data.frame", # list of compounds which will always be deleted
    integRes = "data.frame",
    numcores = "numeric"
  ),
  methods = list(
    initialize = function(...) {
      # Standard settings ####
      settings <<- list(
        db_path = NULL,
        rttolm = 1, mztolu = 0.05,
        mztolu_fine = 0.005,
        chromatography =
          "dx.doi.org/10.1016/j.chroma.2015.11.014",
        pol = "pos", CE_s = 30:40, CES_s = 0:15,
        instr = c("LC-ESI-QTOF TripleTOF 5600 SCIEX", "LC-ESI-QTOF TripleTOF 6600 SCIEX"),
        ceunit = c("V", "eV"),
        comparison = "dot_product", threshold = 400,
        rt_res = 0.333,
        EIC_extraction = 0.02,
        EIC_extraction_range = c(0.01, 0.1),
        baseline_noise_MS1 = 0.6,
        sn = 2, rtoffset = 0, IS_rtoffset = 0, ISrttolm = 1, ISmztol = 0.005,
        blank_int_factor = 5,
        rtTolReinteg = 1, # min rt tolerance for reintegrating peaks
        mzTolReinteg = 0.005,
        ndp_m = 2,
        ndp_n = 1,
        mztolu_ms2 = 0.015,
        area_threshold = 1, # peak area intensity threshold
        height_threshold = 1, # peak height intensity threshold
        use_int_threshold = "area",
        numcores = 1, # number of cores to use for processing
        peaksPerPeak = 10,
        mustFindChromPeak = TRUE # if true, peaks without area are deleted
      )

      rawFilesCompl <<- data.frame(
        path = character(),
        date = character(),
        stringsAsFactors = FALSE
      )
      peakList <<- data.frame(
        peakID = integer(),
        mz = numeric(),
        rt = numeric(),
        index = integer(),
        score = integer(),
        expID = integer(),
        mz_ok = logical(),
        real_mz = numeric(),
        int_h = integer(),
        int_a = integer(),
        peak = character(),
        adduct = character(),
        isotopologue = character(),
        rt_min = numeric(),
        real_rt_min = numeric(),
        comp_id = integer(),
        samp = character(),
        comp_name = character(),
        comp_CAS = character(),
        mz_error_mDa = numeric(),
        rt_error_min = numeric(),
        peak_start = numeric(),
        peak_end = numeric(),
        s_to_n = numeric(),
        duplicate = numeric(),
        eic_extraction_width = numeric(),
        stringsAsFactors = FALSE
      )
      # start counting for peakID
      currentPeakID <<- 1L
      currentISpeakID <<- 1L
      numcores <<- if (Sys.info()["sysname"] == "Linux") detectCores() / 2 else 1

      MS1 <<- data.frame(
        peakID = integer(),
        mz = numeric(),
        int = numeric()
      )
      MS2 <<- data.frame(
        peakID = integer(),
        mz = numeric(),
        int = numeric()
      )
      EIC <<- data.frame(
        peakID = integer(),
        scan = integer(),
        time = numeric(),
        int = numeric()
      )
      IS <<- data.frame(
        name = character(), formula = character(), rt = numeric(),
        default = logical(), adduct = character(),
        stringsAsFactors = FALSE
      )
      ISresults <<- data.frame(
        samp = character(), IS = character(),
        mz = numeric(), rt = numeric(), int_h = integer(),
        int_a = integer(), peak_start = numeric(),
        peak_end = numeric(),
        ISpeakID = integer(),
        eic_extraction_width = numeric(),
        stringsAsFactors = FALSE
      )
      integRes <<- data.frame(
        samp = character(), 
        comp_name = character(), 
        adduct = character(),
        isotopologue = character(),
        int_h = integer(),
        int_a = integer(), s_to_n = numeric(),
        rt_error_min = numeric(), # rt error to the average found from MS2 search
        eic_extraction_width = numeric(),
        real_mz = numeric(), real_rt_min = numeric(),
        stringsAsFactors = FALSE
      )
      falsePos <<- data.frame(name = character(), sampNum = numeric())
      callSuper(...)
    },
    addIS = function(dialog = TRUE, fileName = NULL) {
      "Include a list of internal standards in the report. Dialog indicates use of interactive file
      choosing dialog. The file should be a csv with 4 columns: 'name', 'formula', 'rt', 'adduct'."
      if (dialog) {
        fileName <- rstudioapi::selectFile(path = getwd(), filter = "*.csv")
      }
      df <- read.csv(fileName, stringsAsFactors = FALSE)
      # check that the table is ok
      stopifnot(all.equal(colnames(df), c("name", "formula", "rt", "adduct")))
      stopifnot(inherits(df$rt, "numeric"))

      test <- rbind(IS, df)

      if (any(duplicated(test$name))) {
        stop("There are duplicated IS")
      }

      IS <<- test
    },
    addRawFilesDir = function(dialog = TRUE, dir_path = NULL) {
      "Add all raw files in the chosen directory"
      if (dialog) {
        dir_path <- rstudioapi::selectDirectory(path = "~", caption = "Select directory containing mzXML files")
      }

      newFiles <- list.files(dir_path, pattern = "*\\.mzXML$", full.names = TRUE)
      testMe <- append(basename(rawFiles), basename(newFiles))
      if (anyDuplicated(testMe)) {
        stop("All file names must be unique")
      }
      # check that files exist
      rawFiles <<- append(rawFiles, newFiles)
    },
    addRawFiles = function(dialog = TRUE, file_list = NULL) {
      "Add raw files, (mzML or mzXML) Dialog indicates use of interactive file choosing dialog.
      Does not load files to RAM, only marks file path."
      if (dialog && .Platform$OS.type == "windows") {
        newFiles <- choose.files(
          default = getwd(),
          filters = matrix(
            data = c("mzML or mzXML files", "*.mzML;*.mzXML"),
            ncol = 2
          )
        )
      } else if (dialog && .Platform$OS.type == "unix") {
        newFiles <- rstudioapi::selectFile(path = "~", filter = "mzXML file (*.mzXML)")
      } else if (!dialog) {
        newFiles <- normalizePath(file_list)
      } else {
        stop("No files selected")
      }
      testMe <- append(basename(rawFiles), basename(newFiles))
      if (anyDuplicated(testMe)) {
        stop("All file names must be unique")
      }
      # check that files exist
      rawFiles <<- append(rawFiles, newFiles)
    },
    addDB = function(dialog = TRUE, dbPath = NULL) {
      if (dialog) {
        dbPath <- rstudioapi::selectFile(path = "~", filter = "sqlite file (*.db)")
      }
      if (file.exists(dbPath)) {
        .self$changeSettings("db_path", dbPath)
      } else {
        warning("Db not found.")
      }
    },
    remRawFiles = function(indices) {
      "Delete files based on their indices"
      delIDs <- peakList[peakList$samp %in% basename(rawFiles[indices]), "peakID"]

      rawFilesCompl <<- rawFilesCompl[rawFilesCompl$path %notin% rawFiles[indices], ]
      ISresults <<- ISresults[!(ISresults$samp %in% basename(rawFiles[indices])), ]
      integRes <<- integRes[integRes$samp %notin% basename(rawFiles[indices]), ]
      peakList <<- peakList[!(peakList$peakID %in% delIDs), ]
      MS1 <<- MS1[!(MS1$peakID %in% delIDs), ]
      MS2 <<- MS2[!(MS2$peakID %in% delIDs), ]
      EIC <<- EIC[!(EIC$peakID %in% delIDs), ]
      rawFiles <<- rawFiles[-indices]
    },
    moveRawFile = function(index, direction) {
      "Move a raw file in direction 'up' or 'down'"
      if (length(rawFiles) == 1) {
        return(NULL)
      }
      stopifnot(index %in% seq_along(rawFiles))
      stopifnot(direction %in% c("up", "down"))
      indices <- seq_along(rawFiles)
      endInd <- length(rawFiles)
      if (index == 1 && direction == "up") {
        return(NULL)
      }
      if (index == endInd && direction == "down") {
        return(NULL)
      }

      newIndices <- switch(direction,
        up = replace(indices, (index - 1):index, index:(index - 1)),
        down = replace(indices, index:(index + 1), (index + 1):index)
      )

      rawFiles <<- rawFiles[newIndices]
    },
    changeSettings = function(parameter, value) {
      "Change any setting by name and then value. See settings field for more details."
      stopifnot(parameter %in% names(settings))
      # run various tests
      if (parameter == "EIC_extraction_range" && (length(value) != 2 || !is.numeric(value))) {
        stop("EIC_extraction_range must be a length 2 numeric vector")
      }
      if (parameter == "pol" && (length(value) != 1 || !(value %in% c("pos", "neg")))) {
        stop("Polarity must 'pos' or 'neg'")
      }
      if (parameter == "use_int_threshold" && (length(value) != 1 || !(value %in% c("area", "height")))) {
        stop("Polarity must 'area' or 'height'")
      }
      settings[[parameter]] <<- value
    },
    saveSettings = function(path = getwd()) {
      "Save settings to a json file. Name of file is generated from current time."
      jsonStr <- jsonlite::toJSON(settings, pretty = TRUE)
      write(jsonStr, file = file.path(
        path, paste0(format(Sys.time(), "%y%m%d-%H%M"), "_Settings.json")
      ))
    },
    loadSettings = function() {
      "Load previously settings file from current working directory. This will fail if there is more
      than one file present in the directory."
      path <- choose.files(
        default = getwd(),
        filters = matrix(
          data = c("JSON settings files", "*.json"),
          ncol = 2
        )
      )
      settings <<- jsonlite::fromJSON(path)
    },
    loadData = function(all = FALSE, indices = NULL) {
      "Load files that still need to be processed into RAM for fast access, if all = TRUE then
      all files are loaded regardless, if indices is an integer vector, only these samps will
      be loaded"
      to_process <- if (all) {
        rawFiles
      } else if (!is.null(indices)) {
        rawFiles[indices]
      } else {
        setdiff(rawFiles, rawFilesCompl$path)
      }

      to_process <- setdiff(to_process, names(rawData)) # only those which are not already loaded

      if (length(to_process) == 0) {
        return(NULL)
      }

      rawData_temp <- mclapply(to_process,
        function(x) suppressMessages(xcms::xcmsRaw(x, includeMSn = TRUE)),
        mc.cores = settings$numcores, mc.preschedule = FALSE
      )
      names(rawData_temp) <- to_process
      if (!all(vapply(rawData_temp, inherits, what = "xcmsRaw", logical(1)))) {
        stop("Error in loading samples")
      }
      rawData <<- append(rawData, rawData_temp)

      message(paste0(basename(to_process), collape = ", "))
    },
    clearData = function(indices = NULL) {
      "Remove data from RAM to clear memory, use indices of raw files"
      if (is.null(indices)) {
        rawData <<- vector("list", 0)
      } else {
        paths_to_clear <- rawFiles[indices]
        rawData[paths_to_clear] <<- NULL
      }
      gc()
    },
    getPeak = function(rawLinki, comp_mzi, comp_rti, minIndi, maxIndi, width,
                       mztoli, rttoli) {
      "Internal function to integrate peaks during processing"

      getPeakAtWidth <- function(widthi) {
        resul <- pickPeaksOneEic(
          i = comp_mzi - widthi / 2,
          rawData = rawLinki, mz_step = widthi,
          rt_min_scan = minIndi,
          rt_max_scan = maxIndi,
          sn = settings$sn, int_threshold = settings$baseline_noise_MS1,
          NoiseScans = 60, peakwidth_min = 4,
          peakwidth_max = 100, # large values chosen
          maxPeaksPerSignal = settings$peaksPerPeak,
          precursormzTol = 5
        )
        if (is.null(resul)) {
          resul <- matrix(nrow = 0, ncol = 16)
        }
        # remove unknown extra column after switching to cpp
        resul <- resul[, -4, drop = FALSE]
        colnames(resul) <- c(
          "exactmass", "scantime", "peak_intens",
          "maxima", "scantimeleft_end", "scantimeright_end",
          "left_end", "right_end", "noisedeviation",
          "peakArea", "FWHM_left", "FWHM_right", "noiselevel",
          "i", "ms2scan"
        )
        resul <- as.data.frame(resul)
        resul$scantime_min <- resul[, "scantime"] / 60
        # filter according to known rt from MS2 spectrum
        resul <- resul[(abs(resul$scantime_min - comp_rti) < rttoli) &
          (abs(resul$exactmass - comp_mzi) < mztoli), ]
        if (nrow(resul) != 0) {
          resul$e_width <- widthi
        }
        resul
      }
      resu <- getPeakAtWidth(width)
      # if no peak was found, try the whole range and take the result
      # closest to the target width

      if (nrow(resu) == 0 && diff(settings$EIC_extraction_range) != 0) {
        ra <- settings$EIC_extraction_range
        widths <- seq(ra[1], ra[2], length.out = 10)
        widths <- widths[-which(sapply(widths, all.equal, target = width) == "TRUE")]
        widths <- widths[order(abs(widths - width))]
        for (w in widths) {
          resu <- getPeakAtWidth(w)
          if (nrow(resu) != 0) {
            resu$e_width <- w
            break
          }
        }
      }
      resu
    },
    process_all = function(comp_names = NULL, alsoDeleteFP = TRUE) {
      "Process all currently unprocessed files in the report object. Previously recorded false
      positives (without specified sample) are deleted by default."
      # find out which are left to process
      
      to_process <- setdiff(rawFiles, rawFilesCompl$path)
      .self$loadData(indices = which(rawFiles %in% to_process))

      if (nrow(IS) == 0) {
        warning("No internal standards loaded")
      }
      
      if (is.null(settings$db_path)) {
        stop("There is no DB")
      }

      # loop through each file and peak, collect MS1, EIC and MS2, store in corresponding tables
      for (datFile in to_process) {
        message(paste("Processing", basename(datFile)))
        stopifnot(inherits(rawData[[datFile]], "xcmsRaw"))
        rawLink <- rawData[[datFile]]

        results <- ntsworkflow::ms2_search(
          data_path = list(rawLink),
          db_path = settings$db_path, rttolm = settings$rttolm, mztolu = settings$mztolu,
          mztolu_fine = settings$mztolu_fine,
          chromatography = settings$chromatography,
          pol = settings$pol, CE_s = settings$CE_s, CES_s = settings$CES_s,
          instr = settings$instr,
          ceunit = settings$ceunit,
          comparison = settings$comparison, threshold = settings$threshold,
          rt_res = settings$rt_res, rtoffset = settings$rtoffset,
          ndp_m = settings$ndp_m, ndp_n = settings$ndp_n, mztolu_ms2 = settings$mztolu_ms2,
          compounds = comp_names
        )
        # if no ms2 comparison is desired i.e. threshold score = 0, for compounds not found
        # previously, search again based only on mz and rt
        if (settings$threshold == 0) {
          stop("Currently it is not possible to turn off the MS2 comparison")
          # remove compounds that have MS2 score below 500 from previous list (this is leading to errors)
          results <- results[results$score >= 500, ]

          # get list of remaining compounds
          db <- DBI::dbConnect(RSQLite::SQLite(), settings$db_path)

          # get list of all compounds in db with rt
          expTable <- tbl(db, "experiment")
          paraTable <- tbl(db, "parameter")
          rtTable <- tbl(db, "retention_time")
          compTable <- tbl(db, "compound")
          suspects <- compTable %>%
            left_join(rtTable, by = "compound_id") %>%
            filter(chrom_method == settings$chromatography || is.na(chrom_method)) %>%
            left_join(expTable, by = "compound_id") %>%
            left_join(paraTable, by = "parameter_id") %>%
            filter(polarity == settings$pol) %>%
            filter(CE %in% settings$CE_s) %>%
            filter(CES %in% settings$CES_s) %>%
            filter(instrument %in% settings$instr) %>%
            select(name, CAS, mz, rt, adduct, experiment_id, compound_id) %>%
            dplyr::collect() %>%
            distinct()

          if (!is.null(comp_names)) {
            if (!all(comp_names %in% suspects$name)) {
              stop("Chosen compounds not in DB")
            }
            suspects <- filter(suspects, name %in% comp_names)
          }

          susToCheck <- setdiff(unique(suspects$name), unique(results$comp_name))
          suspects <- suspects[suspects$name %in% susToCheck, ]
          # there must be an RT available, i.e. search based purely on mz not possible
          suspects <- suspects[!is.na(suspects$rt), ]
          eval_comp_mzrt <- function(compound) {
            # get spectrum of scan for RT of compound (if no RT given, return NULL)
            compMz <- compound$mz[1]
            compRt <- (compound$rt[1] + settings$rtoffset) * 60
            msScans <- which(abs(rawLink@scantime - compRt) <= settings$rttolm * 60)
            getSpecScan <- function(scanNr) {
              r <- xcms::getScan(rawLink, scanNr, wind(compMz, settings$mztolu_fine))
              suppressWarnings(cbind(r, scanNr))
            }
            specs <- lapply(msScans, getSpecScan)
            allSig <- Reduce(rbind, specs)
            # check to see if mass is present above treshold
            allSig <- allSig[allSig[, "intensity"] >= settings$baseline_noise_MS1, , drop = FALSE]
            if (nrow(allSig) == 0) {
              return(NULL)
            }
            allSig <- cbind(allSig, rawLink@scantime[allSig[, "scanNr"]])
            colnames(allSig)[4] <- "rt"
            allSig <- allSig[order(allSig[, "intensity"], decreasing = TRUE), , drop = FALSE]
            new_result <- data.frame(
              mz = allSig[, "mz"][1],
              rt = allSig[, "rt"][1],
              index = NA,
              score = NA,
              expID = NA,
              mz_ok = TRUE,
              real_mz = allSig[, "mz"][1],
              int_h = allSig[, "intensity"][1],
              peak = "A",
              rt_min = allSig[, "rt"][1] / 60,
              comp_id = compound$compound_id[1],
              samp = basename(datFile),
              comp_name = compound$name[1],
              comp_CAS = compound$CAS[1],
              mz_error_mDa = (allSig[, "mz"][1] - compMz) * 1000,
              rt_error_min = (allSig[, "rt"][1] - compRt) / 60,
              stringsAsFactors = FALSE
            )

            new_result
          }
          further_res <- by(suspects, suspects$name, eval_comp_mzrt, simplify = FALSE)
          further_res <- do.call("rbind", further_res)
          # need to clear class names so that rbind works
          if (nrow(further_res) >= 1) {
            class(results) <- "data.frame"
            results <- rbind(results, further_res)
          }
        }
        
        # Clear features which are under the MS1 baseline
        if (!is.null(results)) {
          results <- results[results$int_h >= settings$baseline_noise_MS1, , drop = FALSE]
        }
        
        # Continue with this rest only if peaks are found
        if (!is.null(results) && nrow(results) != 0) {
          # check for duplicate detections
          du <- results[, c("index", "samp")]
          du <- du[duplicated(du), ]
          if (nrow(du) == 0) {
            results$duplicate <- NA
          } else {
            # get current max duplicate id number from peakList
            if (nrow(peakList) > 0) {
              maxId <- max(peakList$duplicate, na.rm = TRUE)
            } else {
              maxId <- 0
            }
            for (i in seq_len(nrow(du))) {
              results[results$index == du$index[i] &
                results$samp == du$samp[i], "duplicate"] <- maxId + i
            }
          }
          # add peak IDs
          results$peakID <- NA
          for (i in seq_len(nrow(results))) {
            results$peakID[i] <- currentPeakID
            currentPeakID <<- currentPeakID + 1L
          }
          # need to clear class names so that rbind works
          class(results) <- "data.frame"
          results$int_a <- NA
          results$peak_start <- NA
          results$peak_end <- NA
          results$s_to_n <- NA
          results$real_rt_min <- NA
          results$eic_extraction_width <- NA
          
          # Add adduct to peak list
          # slib <- DBI::dbConnect(RSQLite::SQLite(), settings$db_path)
          # exptbl <- tbl(slib, "experiment") %>% 
          #   select(experiment_id, adduct) %>% 
          #   collect()
          # DBI::dbDisconnect(slib)
          # get_adduct <- function(expid) {
          #   exptbl %>% filter(experiment_id == !!expid) %>% 
          #     select(adduct) %>% unlist()
          # }
          # results$adduct <- vapply(results$expID, get_adduct, character(1))
          
          # bind results to existing results
          # export peaklist ####
          peakList <<- rbind(peakList, results)

          # get only results for this datafile
          subres <- peakList[peakList$samp == basename(datFile), ]

          if (nrow(subres) == 0) {
            message("No compounds found")
            next
          }
          getMS2 <- function(thisID, thisIndex) {
            if (is.na(thisIndex)) {
              return(NULL)
            }
            ms2Spec <- xcms::getMsnScan(rawLink, thisIndex)
            ms2Spec <- as.data.frame(ms2Spec)
            colnames(ms2Spec) <- c("mz", "int")

            ms2Spec$peakID <- thisID
            ms2Spec
          }

          MS2List <- Map(getMS2, subres$peakID, subres$index)
          MS2df <- do.call("rbind", MS2List)

          # export MS2 ####
          MS2 <<- rbind(MS2, MS2df)

          getMS1 <- function(thisID, thisRt, thisMz, thisInt) {
            ind <- which.min(abs(rawLink@scantime - thisRt))
            failvar <- FALSE
            tryCatch(
              ms1Spec <- xcms::getSpec(rawLink, scanrange = c(ind - 1, ind, ind + 1)),
              error = function(cnd) {
                log_warn("In getSpec() for peakID = {thisID}: {conditionMessage(cnd)}")
                failvar <<- TRUE
              }
            )
            if (failvar) {
              for(i in c(ind, ind-1, ind+1)) {  # check if any of them returns results with multiple rows
                ms1Spec <- xcms::getScan(rawLink, i)
                if (nrow(ms1Spec) > 1) {  # first one to return nrows > 1
                  # Test if m/z of peak is found in spectrum (within tolerance)
                  if (any(abs(ms1Spec[,1] - thisMz) <= settings$mztolu_fine)) {
                    log_info("Used single scan instead and found peak for peakID = {thisID}")
                    break
                  }
                }
              }
            }
            ms1Spec <- as.data.frame(ms1Spec)
            colnames(ms1Spec) <- c("mz", "int")
            ms1Spec <- ms1Spec[ms1Spec$int >= settings$baseline_noise_MS1, ]
            ms1Spec <- ms1Spec[ms1Spec$mz > thisMz - 2, ]
            ms1Spec <- ms1Spec[ms1Spec$mz < thisMz + 8, ]
            ms1Spec <- ms1Spec[ms1Spec$int >= thisInt * 0.05, ] # remove anything less than 5% of int
            ms1Spec$peakID <- thisID
            ms1Spec <- ms1Spec[!is.na(ms1Spec$mz) & !is.na(ms1Spec$int), ]
            ms1Spec
          }
          
          MS1List <- Map(getMS1, subres$peakID, subres$rt, subres$real_mz, subres$int_h)
          MS1df <- do.call("rbind", MS1List)
          # export MS1 ####
          MS1 <<- rbind(MS1, MS1df)

          getEic <- function(thisID, thisRt, thisMz) {
            ext <- settings$EIC_extraction
            thisEic <- xcms::rawEIC(rawLink, mzrange = c(thisMz - ext / 2, thisMz + ext / 2))
            thisEic <- as.data.frame(thisEic)
            colnames(thisEic) <- c("scan", "int")
            thisEic$time <- rawLink@scantime[thisEic$scan]
            thisEic <- thisEic[thisEic$int >= settings$baseline_noise_MS1, ]
            thisEic <- thisEic[abs(thisEic$time - thisRt) <= 200, ]

            if (nrow(thisEic) == 0) {
              thisEic <- data.frame(scan = NA, time = NA, int = NA)
            }

            thisEic$peakID <- thisID
            thisEic
          }
          EicList <- Map(getEic, subres$peakID, subres$rt, subres$real_mz)
          Eicdf <- do.call("rbind", EicList)
          # export eic ####
          EIC <<- rbind(EIC, Eicdf)

          # loop through all peaks and collect data ####
          for (j in seq_len(nrow(subres))) {
            thisID <- subres$peakID[j]

            minInd <- which.min(abs(rawLink@scantime - (subres$rt[j] - (settings$rtTolReinteg * 60) / 2)))
            maxInd <- which.min(abs(rawLink@scantime - (subres$rt[j] + (settings$rtTolReinteg * 60) / 2)))

            # calculate area for this peak

            comp_mz <- subres$real_mz[j]
            comp_rt <- subres$rt_min[j]

            # Peak integration ####
            # using the peak picking algorithm for this mass
            res <- .self$getPeak(
              rawLink, comp_mz, comp_rt, minInd, maxInd,
              settings$EIC_extraction,
              settings$mzTolReinteg,
              settings$rtTolReinteg
            )
            # Add data to peak list
            if (nrow(res) != 0) {
              # choose largest peak, normally there should only be one remaining peak
              res <- res[which.max(res$peakArea), ]
              peakList[peakList$peakID == thisID, "real_rt_min"] <<- round(res$scantime / 60, 2)
              peakList[peakList$peakID == thisID, "int_h"] <<- as.integer(round(res$peak_intens))
              peakList[peakList$peakID == thisID, "int_a"] <<- as.integer(round(res$peakArea))
              peakList[peakList$peakID == thisID, "peak_start"] <<- round(res$scantimeleft_end / 60, 2)
              peakList[peakList$peakID == thisID, "peak_end"] <<- round(res$scantimeright_end / 60, 2)
              peakList[peakList$peakID == thisID, "s_to_n"] <<- round(res$peak_intens / res$noisedeviation, 1)
              peakList[peakList$peakID == thisID, "eic_extraction_width"] <<- res$e_width
            }
          }
          # Remove any false positives
          if (alsoDeleteFP && nrow(falsePos) >= 1) {
            for (row in seq_len(nrow(falsePos))) {
              if (falsePos[row, "sampNum"] == 0) {
                .self$deleteFP(falsePos[row, "name"], 0)
              }
            }
          }
        }

        # Integrate IS ####
        for (k in seq_len(nrow(IS))) {
          thisCharge <- switch(settings$pol,
            pos = 1,
            neg = -1
          )
          isRtSec <- (IS[k, "rt"] + settings$IS_rtoffset) * 60
          minISInd <- which.min(abs(rawLink@scantime - (isRtSec - (settings$ISrttolm * 60) / 2)))
          maxISInd <- which.min(abs(rawLink@scantime - (isRtSec + (settings$ISrttolm * 60) / 2)))

          IS_mz <- ntsworkflow::get_mass(IS[k, "formula"],
            charge = thisCharge,
            adduct = IS[k, "adduct"]
          )
          IS_rt <- IS[k, "rt"] + settings$IS_rtoffset
          IS_res <- .self$getPeak(
            rawLink, IS_mz, IS_rt, minISInd, maxISInd,
            settings$EIC_extraction,
            settings$ISmztol,
            settings$ISrttolm
          )

          if (nrow(IS_res) == 0) {
            message(paste(IS[k, "name"], "not found."), appendLF = FALSE)
            next
          }
          # get highest peak
          IS_res <- IS_res[which.max(IS_res$peak_intens), ]

          IS_res2 <- data.frame(
            samp = basename(datFile), IS = IS[k, "name"],
            mz = IS_res$exactmass, rt = IS_res$scantime_min,
            int_h = as.integer(round(IS_res$peak_intens)),
            int_a = as.integer(round(IS_res$peakArea)),
            peak_start = round(IS_res$scantimeleft_end / 60, 2),
            peak_end = round(IS_res$scantimeright_end / 60, 2),
            eic_extraction_width = IS_res$e_width,
            ISpeakID = as.integer(currentISpeakID),
            stringsAsFactors = FALSE
          )
          # export IS results
          ISresults <<- rbind(ISresults, IS_res2)
          currentISpeakID <<- currentISpeakID + 1L
        }

        # update rawfiles completed list
        processed <- data.frame(path = datFile, date = date(), stringsAsFactors = FALSE)
        # update raw file complete
        rawFilesCompl <<- rbind(rawFilesCompl, processed)
        rm(rawLink)
        message("Complete")
      }
      if (settings$use_int_threshold %in% c("area", "height")) {
        .self$deleteBelowIntThresh()
      }

      if (settings$mustFindChromPeak) {
        .self$delPeakID(peakList[is.na(peakList$int_a), "peakID"])
      }

      .self$cleanPeakLetterCol()
      .self$reduceSizeEic()

      message(sprintf("Processing completed %i files", length(to_process)))
    },
    reprocess = function(indices = NULL, comp_names = NULL, alsoDeleteFP = TRUE) {
      "Reprocess files or specific compounds, will delete all previous information. As with
      process_all, previous FP will be deleted by default"
      # first delete results corresponding to results that are to be reprocessed
      # get peak IDs of all peaks to be deleted
      
      if (is.null(indices) && is.null(comp_names)) {
        indices <- seq_along(rawFiles)
        delIDs <- peakList[, "peakID"]
      } else if (!is.null(indices) && is.null(comp_names)) {
        stopifnot(all(indices %in% seq_along(rawFiles)))
        delIDs <- peakList[peakList$samp %in% basename(rawFiles[indices]), "peakID"]
      } else if (is.null(indices) && !is.null(comp_names)) {
        stopifnot(is.character(comp_names))
        delIDs <- peakList[peakList$comp_name %in% comp_names, "peakID"]
      } else if (!is.null(indices) && !is.null(comp_names)) {
        stopifnot(is.character(comp_names))
        stopifnot(all(indices %in% seq_along(rawFiles)))
        delIDs <- peakList[peakList$comp_name %in% comp_names &
          peakList$samp %in% basename(rawFiles[indices]), "peakID"]
      } else {
        stop("Could not find peaks to reprocess")
      }

      .self$delPeakID(delIDs)

      # if all files are to be reprocessed, need to clear some data
      if (is.null(indices)) {
        indices <- seq_along(rawFiles)
      }
      # in any case redo the IS, for now, it is just easier that way
      ISresults <<- ISresults[!(ISresults$samp %in% basename(rawFiles[indices])), ]

      # delete rows from rawFilesCompl table
      toKeep <- !(basename(rawFilesCompl$path) %in% basename(rawFiles[indices]))
      rawFilesCompl <<- rawFilesCompl[toKeep, ]

      # re-process files
      .self$loadData()
      .self$process_all(comp_names = comp_names, alsoDeleteFP = alsoDeleteFP)
    },

    # produce plot of EIC for a particular peak
    plotEIC = function(ID) {
      "Plot a particular EIC using the ID from the peak list table"
      stopifnot(length(ID) == 1)
      stopifnot(ID %in% peakList$peakID)
      toPlot <- EIC[EIC$peakID == ID, ]
      toPlot$time_min <- toPlot$time / 60
      rt_min <- peakList[peakList$peakID == ID, "rt_min"]
      comp_name <- peakList[peakList$peakID == ID, "comp_name"]
      comp_CAS <- peakList[peakList$peakID == ID, "comp_CAS"]
      int_i <- peakList[peakList$peakID == ID, "int_h"]
      mz <- peakList[peakList$peakID == ID, "real_mz"]
      samp <- peakList[peakList$peakID == ID, "samp"]
      e_width <- peakList[peakList$peakID == ID, "eic_extraction_width"]
      stopifnot(length(comp_CAS) == 1)
      start_t <- peakList[peakList$peakID == ID, "peak_start"]
      end_t <- peakList[peakList$peakID == ID, "peak_end"]
      thisPlot <- ggplot2::ggplot(toPlot, aes(x = time_min, y = int)) +
        geom_line() +
        theme_bw() +
        xlab("Time (min)") +
        ylab("Intensity (counts)") +
        coord_cartesian(
          xlim = c(rt_min - 5, rt_min + 5),
          ylim = c(0, max(toPlot[abs(toPlot$time_min - rt_min) <= 5, "int"]) * 1.1)
        ) +
        geom_vline(xintercept = rt_min, color = "red", alpha = 0.3) +
        annotate("text",
          x = rt_min - 1, y = Inf, vjust = 1, hjust = 1,
          label = sprintf(
            "\nName = %s\nCAS = %s\nRT = %.1f min\nm/z = %.4f\nwidth = %.3f Da",
            comp_name, comp_CAS, rt_min, mz, e_width
          ),
          color = "blue", alpha = 0.5, size = 3
        ) +
        annotate("text",
          x = rt_min + 1, y = Inf, label = paste0("\n", samp), color = "blue",
          alpha = 0.5, size = 3, vjust = 1, hjust = 0
        )

      if (!is.na(start_t) && !is.na(end_t)) {
        thisPlot <- thisPlot +
          annotate("segment",
            x = start_t, xend = start_t, y = 0, yend = int_i * 0.05, colour = "blue",
            alpha = 0.5, size = 0.5
          ) +
          annotate("segment",
            x = end_t, xend = end_t, y = 0, yend = int_i * 0.05, colour = "blue",
            alpha = 0.5, size = 0.5
          )
      }
      return(thisPlot)
    },
    plotMS1 = function(ID) {
      "Plot a MS1 spectrum using the ID from the peak list table. Includes comparison to calculated
      isotopic pattern."
      stopifnot(length(ID) == 1)
      stopifnot(ID %in% peakList$peakID)
      spec <- MS1[MS1$peakID == ID, ]
      expID <- peakList[peakList$peakID == ID, "expID"]
      compName <- peakList[peakList$peakID == ID, "comp_name"]
      newDbPath <- settings$db_path
      
      db <- DBI::dbConnect(RSQLite::SQLite(), newDbPath)
      expTbl <- tbl(db, "experiment")
      comTbl <- tbl(db, "compound")
      r <- if (is.na(expID)) {
        filter(comTbl, name == compName) %>%
          left_join(expTbl, by = "compound_id") %>%
          dplyr::collect() %>%
          .[1, ]
      } else {
        filter(expTbl, experiment_id == expID) %>%
          left_join(comTbl, by = "compound_id") %>%
          dplyr::collect()
      }

      # calculate theoretical MS1 from isotope composition for mol-formular mf
      formula <- r$formula
      adduct <- r$adduct

      DBI::dbDisconnect(db)

      # find the actual molecular formula
      form_charge <- ntsworkflow::correct_formula(r$formula, r$adduct)

      comp_mz <- peakList[peakList$peakID == ID, "real_mz"]
      comp_int <- peakList[peakList$peakID == ID, "int_h"]
      if (comp_int < settings$baseline_noise_MS1) {
        comp_int <- settings$baseline_noise_MS1
      }

      comp_mz_error <- peakList[peakList$peakID == ID, "mz_error_mDa"]
      # extract data for main peak
      mz_lab_data <- spec[abs(spec$mz - comp_mz) < settings$EIC_extraction, ]

      if (!is.na(formula)) {
        mf <- rcdk::get.formula(form_charge$form, charge = form_charge$charge)
        isotope_spec <- rcdk::get.isotopes.pattern(mf, minAbund = 0.01)

        isotope_spec <- as.data.frame(isotope_spec)
        colnames(isotope_spec) <- c("mz", "int")
        # mirror intensity of MS1 spec
        isotope_spec$int <- isotope_spec$int * -max(comp_int)

        ms1Plot <- ggplot2::ggplot(spec, aes(mz, int, label = round(mz, 4))) +
          geom_segment(aes(x = mz, xend = mz, y = 0, yend = int),
            stat = "identity", size = .5, alpha = 0.5
          ) +
          geom_segment(
            data = isotope_spec, aes(x = mz, xend = mz, y = 0, yend = int),
            stat = "identity", size = .5, alpha = 0.5
          ) +
          theme_bw() +
          # geom_text(data = mz_lab_data, vjust = -0.5, alpha = 0.4) +
          geom_text(data = spec[spec$int > comp_int * 0.2, ], check_overlap = TRUE, vjust = -0.5) +
          geom_text(
            data = isotope_spec[isotope_spec$int < -comp_int * 0.2, ],
            check_overlap = TRUE, vjust = 1.5
          ) +
          geom_vline(xintercept = comp_mz, color = "red", alpha = 0.2) +
          scale_y_continuous(expand = c(0, 0)) +
          coord_cartesian(
            xlim = c(comp_mz - 2, comp_mz + 8),
            ylim = c(-1.2 * max(comp_int), 1.2 * max(comp_int))
          ) +
          ylab("Intensity (abs.)") +
          xlab("m/z (u)") +
          geom_hline(yintercept = 0, color = "blue") +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = .1, fill = "orange") +
          annotate("text",
            x = -Inf, y = -Inf, label = "Calculated-Spec", vjust = -1, hjust = -.1,
            fontface = "bold.italic", alpha = .5
          ) +
          annotate("text",
            x = Inf, y = Inf, label = sprintf("\nmz-error (mDa): %.2f", comp_mz_error),
            vjust = .7, hjust = 1.1, fontface = "italic", alpha = .5
          ) +
          annotate("text",
            x = -Inf, y = Inf, label = "Data-Spec", vjust = 1.5, hjust = -.1,
            fontface = "bold.italic", alpha = .5
          )
      } else {
        ms1Plot <- ggplot2::ggplot(spec, aes(mz, int, label = round(mz, 4))) +
          geom_segment(aes(x = mz, xend = mz, y = 0, yend = int),
            stat = "identity", size = .5, alpha = 0.5
          ) +
          theme_bw() +
          geom_text(data = spec[spec$int > comp_int * 0.2, ], check_overlap = TRUE, vjust = -0.5) +
          geom_text(data = mz_lab_data, vjust = -0.5) +
          geom_vline(xintercept = comp_mz, color = "red", alpha = 0.2) +
          scale_y_continuous(expand = c(0, 0)) +
          coord_cartesian(
            xlim = c(comp_mz - 2, comp_mz + 8),
            ylim = c(0, 1.2 * max(comp_int))
          ) +
          ylab("Intensity (abs.)") +
          xlab("m/z (u)") +
          geom_hline(yintercept = 0, color = "blue") +
          annotate("text",
            x = Inf, y = Inf, label = sprintf("mz-error (mDa): %.2f", comp_mz_error),
            vjust = -1, hjust = 1.1, fontface = "italic", alpha = .5
          ) +
          annotate("text",
            x = -Inf, y = Inf, label = "Data-Spec", vjust = 1.5, hjust = -.1,
            fontface = "bold.italic", alpha = .5
          )
      }
      return(ms1Plot)
    },
    plotMS2 = function(ID) {
      "Plot an MS2 spectrum using the ID from the peak list table. Includes comparison to database
      spectrum used to make the identification."
      stopifnot(length(ID) == 1)
      stopifnot(ID %in% peakList$peakID)
      # get ms2 spec from database
      expID <- peakList[peakList$peakID == ID, "expID"]
      newDbPath <- settings$db_path
     
      db <- DBI::dbConnect(RSQLite::SQLite(), newDbPath)
      db_spec <- ntsworkflow::dbGetSpectrum(db, expID)
      DBI::dbDisconnect(db)
      # get MS2 from data
      data_spec <- MS2[MS2$peakID == ID, ]

      score <- peakList[peakList$peakID == ID, "score"]
      comp_mz <- peakList[peakList$peakID == ID, "real_mz"]
      name <- peakList[peakList$peakID == ID, "comp_name"]

      # get relative intensities
      data_spec <- data_spec[data_spec$mz < comp_mz + 0.015, ]
      data_spec$int <- data_spec$int / max(data_spec$int)
      db_spec$int <- db_spec$int / max(db_spec$int)
      db_spec$int <- db_spec$int * -1

      ms2Plot <- ggplot2::ggplot(data_spec, aes(mz, int, label = round(mz, 4))) +
        geom_segment(aes(x = mz, xend = mz, y = 0, yend = int),
          stat = "identity", size = .5, alpha = .5
        ) +
        geom_segment(
          data = db_spec, aes(x = mz, xend = mz, y = 0, yend = int),
          stat = "identity", size = .5, alpha = .5
        ) +
        theme_bw() +
        geom_text(data = data_spec[data_spec$int > 0.01, ], check_overlap = TRUE, vjust = -0.5) +
        geom_text(data = db_spec[db_spec$int < -0.01, ], check_overlap = TRUE, vjust = 1.5) +
        geom_vline(xintercept = comp_mz, color = "red", alpha = 0.2) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(-1.2, 1.2), xlim = c(0, comp_mz + 5)) +
        ylab("Intensity (relative)") +
        xlab("m/z (u)") +
        geom_hline(yintercept = 0, color = "blue") +
        annotate("rect",
          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, alpha = .1,
          fill = "orange"
        ) +
        annotate("text",
          x = -Inf, y = -1.1, label = "DB-Spec", vjust = 0.5, hjust = -.1,
          fontface = "bold.italic", alpha = .5
        ) +
        annotate("text",
          x = Inf, y = -1.1, label = name, vjust = 0.5, hjust = 1.1,
          fontface = "bold.italic", alpha = .5, colour = "darkorange4"
        ) +
        annotate("text",
          x = -Inf, y = 1.05, label = "Data-Spec", vjust = 0, hjust = -.1,
          fontface = "bold.italic", alpha = .5
        ) +
        annotate("text",
          x = -Inf, y = .95, label = paste("Score:", round(score)), vjust = 0.5,
          hjust = -.1, fontface = "bold.italic", alpha = .5, colour = "darkorange4"
        )
      return(ms2Plot)
    },
    delPeakID = function(delIDs) {
      peakList <<- peakList[!(peakList$peakID %in% delIDs), ]
      MS1 <<- MS1[!(MS1$peakID %in% delIDs), ]
      MS2 <<- MS2[!(MS2$peakID %in% delIDs), ]
      EIC <<- EIC[!(EIC$peakID %in% delIDs), ]

      .self$cleanPeakLetterCol()
    },
    deleteFP = function(substanceName, indices = 0) {
      "Delete all occurances of a compound given by name. Only one compound can be given. Name must
      match exactly. Compound is added to list of falsePos (see fields). If indices = 0, then
      compound is always deleted, in all samples and also in future processed samples.
      Otherwise it is only deleted in the specified samples (by index)."
      stopifnot(length(substanceName) == 1, is.numeric(indices))
      if (length(indices) > 1 && 0 %in% indices)
        stop("indices can be either 0 (all files) or the file index or indices as given by $rawFiles")

      if (substanceName %in% peakList$comp_name && length(indices) == 1 && indices[1] == 0) {
        # get all peakIDs corresponding to this substance
        delIDs <- peakList[peakList$comp_name == substanceName, "peakID"]

        .self$delPeakID(delIDs)
        integRes <<- integRes[integRes$comp_name != substanceName, ]
        # add fp to list
        if (substanceName %notin% falsePos$name |
          all(falsePos[falsePos$name == substanceName, "sampNum"] != 0)) {
          tempdf <- data.frame(name = substanceName, sampNum = 0, stringsAsFactors = FALSE)
          falsePos <<- rbind(falsePos, tempdf)
        }
      } else if (substanceName %in% peakList$comp_name && indices[1] != 0) {
        samps <- basename(rawFiles[indices])
        delIDs <- peakList[peakList$comp_name == substanceName &
          peakList$samp %in% samps, "peakID"]
        .self$delPeakID(delIDs)
        integRes <<- integRes[!(integRes$comp_name == substanceName &
          integRes$samp %in% samps), ]
        # add fp to list
        if (substanceName %notin% falsePos$name) {
          tempdf <- data.frame(name = substanceName, sampNum = indices, stringsAsFactors = FALSE)
          falsePos <<- rbind(falsePos, tempdf)
        } else {
          # indices which are missing from the current list
          toAdd <- indices[indices %notin% falsePos[falsePos$name == substanceName, "sampNum"]]
          tempdf <- data.frame(name = substanceName, sampNum = toAdd, stringsAsFactors = FALSE)
          falsePos <<- rbind(falsePos, tempdf)
        }
      } else {
        # Even if not found, still need to add this compound to the list
        if (substanceName %notin% falsePos$name) {
          tempdf <- data.frame(name = substanceName, sampNum = indices, stringsAsFactors = FALSE)
          falsePos <<- rbind(falsePos, tempdf)
        } else {
          missingIndices <- setdiff(indices, falsePos[falsePos$name == substanceName, "sampNum"])
          if (length(missingIndices) != 0) {
            tempdf <- data.frame(name = substanceName, sampNum = missingIndices, stringsAsFactors = FALSE)
            falsePos <<- rbind(falsePos, tempdf)
          }
        }
      }

      .self$cleanDuplicateCol()
    },
    cleanDuplicateCol = function() {
      "Check the peaklist for duplicates, if peaks previously indicated as
      duplicates are no longer so then set duplicate column back to NA"
      test <- peakList$duplicate
      ucount <- table(test)
      if (any(!is.na(test)) && any(ucount == 1)) {
        unary <- as.numeric(names(ucount[ucount == 1]))
        peakList[peakList$duplicate %in% unary, "duplicate"] <<- NA
      }
    },
    cleanPeakLetterCol = function() {
      "Check the peaklist for multiple peaks from one compound, if any B, C etc.
      peaks no longer have an A peak, these are renamed to start with A"
      if (nrow(peakList) == 0) {
        return(NULL)
      }
      # go through each sample in peakList
      pll <- split(peakList, peakList$samp)
      relabel <- function(pl) {
        # find which compounds have peaks with B or C... label
        notA <- unique(pl[pl$peak != "A", "comp_name"])
        for (checkComp in notA) { # checkComp <- notA[1]
          #  sort by intensity and rename peak col
          plComp <- pl[pl$comp_name == checkComp, ]
          idsToChange <- plComp[order(plComp$int_h, decreasing = TRUE), "peakID"]
          for (i in seq_along(idsToChange)) {
            pl[pl$peakID == idsToChange[i], "peak"] <- LETTERS[i]
          }
        }
        pl
      }
      pll <- lapply(pll, relabel)
      peakList <<- do.call("rbind", pll)
      rownames(peakList) <<- NULL
    },
    deleteBelowIntThresh = function() {
      switch(settings$use_int_threshold,
        "area" = {
          toThrow <- ifelse(is.na(peakList$int_a),
            peakList$int_h < settings$height_threshold,
            peakList$int_a < settings$area_threshold
          )
          delIDs <- peakList[toThrow, "peakID"]
          .self$delPeakID(delIDs)
          message(paste(
            length(delIDs),
            "peaks were below area threshold",
            settings$area_threshold, "and deleted."
          ))
        },
        "height" = {
          toThrow <- ifelse(is.na(peakList$int_h),
            TRUE,
            peakList$int_h < settings$height_threshold
          )
          delIDs <- peakList[toThrow, "peakID"]
          .self$delPeakID(delIDs)
          message(paste(
            length(delIDs),
            "peaks were below height threshold", settings$height_threshold, "and deleted."
          ))
        },
        "none" = {
          NULL
        },
        stop(
          "settings$use_int_threshold must be one of area, height or none"
        )
      )
    },
    reduceSizeEic = function() {
      "Will average out data points outside of the peak to reduce size of the EIC table"

      shrinkEic <- function(x) {
        thisId <- x$peakID[1]

        oo <- zoo::zoo(cbind(x$scan, x$int), x$time)
        if (nrow(oo) > 20) {
          oor <- zoo::rollapply(oo, 3, mean, by = 3)
        } else {
          oor <- oo
        }

        y <- data.frame(
          scan = round(oor[, 1]), int = oor[, 2], time = attr(oor, "index"),
          peakID = thisId, stringsAsFactors = FALSE
        )
        # if a peak was found, remove averaged points within peak and replace with originals
        plrow <- peakList[peakList$peakID == thisId, , drop = FALSE]
        if (!is.na(plrow[1, "int_a"])) {
          startTime <- plrow[1, "peak_start"] * 60
          endTime <- plrow[1, "peak_end"] * 60
          y <- rbind(
            y[y$time <= startTime | y$time >= endTime, ],
            x[x$time >= startTime & x$time <= endTime, ]
          )
        }
        y
      }

      if (nrow(EIC) > 0) {
        eicl <- by(EIC, EIC$peakID, shrinkEic, simplify = FALSE)
        EIC <<- do.call("rbind", eicl)
      }
    },
    deleteBackground = function(indicesData, indicesBlank, includeIS = FALSE) {
      "Delete background signals based on intensity. Use indices to indicate blanks and samples. If
      there is more than one blank then these will be combined. If IS should also be deleted use
      includeIS. Removal is based on intensity ratio, see settings."
      stopifnot(indicesData %in% seq_along(rawFiles))
      stopifnot(indicesBlank %in% seq_along(rawFiles))

      # using the indices, select the datafiles used for selection
      blankSamp <- basename(rawFiles[indicesBlank])
      dataSamp <- basename(rawFiles[indicesData])

      # find peakIDs and IS names of compounds found in blank files, but only those which are in the
      # data files specified
      # For a match, the intensity in sample must be within factor 3

      dPeaks <- peakList[peakList$samp %in% dataSamp, c("peakID", "comp_name", "int_a", "int_h")]
      bPeaks <- peakList[peakList$samp %in% blankSamp, c("peakID", "comp_name", "int_a", "int_h")]
      delIDs <- numeric()
      for (i in seq_len(nrow(dPeaks))) {
        cname <- dPeaks[i, "comp_name"]
        if (cname %in% bPeaks$comp_name) {
          # get average intensity in blank
          aveIntBlank <- mean(bPeaks[bPeaks$comp_name == cname, "int_a"], na.rm = TRUE)
          # use int_h for comparison if int_a not available
          if (is.na(dPeaks[i, "int_a"]) || is.na(aveIntBlank)) {
            aveIntBlank <- mean(bPeaks[bPeaks$comp_name == cname, "int_h"], na.rm = TRUE)

            if (dPeaks[i, "int_h"] <= (aveIntBlank * settings$blank_int_factor)) {
              delIDs <- append(delIDs, dPeaks[i, "peakID"])
            }
            next
          }

          if (dPeaks[i, "int_a"] <= (aveIntBlank * settings$blank_int_factor)) {
            delIDs <- append(delIDs, dPeaks[i, "peakID"])
          }
        }
      }
      if (length(delIDs) >= 1) {
        peakList <<- peakList[peakList$peakID %notin% delIDs, ]
        MS1 <<- MS1[MS1$peakID %notin% delIDs, ]
        MS2 <<- MS2[MS2$peakID %notin% delIDs, ]
        EIC <<- EIC[EIC$peakID %notin% delIDs, ]
      } else {
        message("No compounds found in blank")
      }

      if (includeIS) {
        # delete IS in the data files specified
        dISPeaks <- ISresults[ISresults$samp %in% dataSamp, c("ISpeakID", "IS", "int_a", "int_h")]
        bISPeaks <- ISresults[ISresults$samp %in% blankSamp, c("ISpeakID", "IS", "int_a", "int_h")]
        delISIDs <- numeric()
        for (j in seq_len(nrow(dISPeaks))) {
          ISname <- dISPeaks[j, "IS"]
          if (ISname %in% bISPeaks$IS) {
            aveIntISblk <- mean(bISPeaks[bISPeaks$IS == ISname, "int_a"], na.rm = FALSE)
            if (is.na(dISPeaks[j, "int_a"]) || is.na(aveIntISblk)) {
              aveIntISblk <- mean(bISPeaks[bISPeaks$IS == ISname, "int_h"], na.rm = FALSE)

              if (dISPeaks[j, "int_h"] <= (aveIntISblk * settings$blank_int_factor)) {
                delISIDs <- append(delISIDs, dISPeaks[j, "ISpeakID"])
              }
              next
            }

            if (dISPeaks[j, "int_a"] <= (aveIntISblk * settings$blank_int_factor)) {
              delISIDs <- append(delISIDs, dISPeaks[j, "ISpeakID"])
            }
          }
        }
        if (length(delISIDs) >= 1) {
          ISresults <<- ISresults[ISresults$ISpeakID %notin% delISIDs, ]
        } else {
          message("No IS found in blank")
        }
      }

      .self$cleanPeakLetterCol()
    },
    reIntegrate = function(indices) {
      "Reintegrate peak areas for all found compounds based only on m/z and rt. No MS2 comparison
      will be performed. The m/z and rt are taken as an average of those found in the peak list."
      # delete any previous information in integRes df for the indices
      # no argument implies all should be integrated, delete all previous reintegration results
      if (missing(indices)) {
        indices <- seq_along(rawFiles)
        integRes <<- data.frame(
          samp = character(), comp_name = character(),
          adduct = character(),
          isotopologue = character(),
          int_h = integer(),
          int_a = integer(), s_to_n = numeric(),
          rt_error_min = numeric(),
          eic_extraction_width = numeric(),
          real_mz = numeric(), real_rt_min = numeric(),
          stringsAsFactors = FALSE
        )
      }

      if (nrow(integRes) > 0) {
        integRes <<- integRes[integRes$samp %notin% basename(rawFiles[indices]), ]
      }

      # get all unique compound+adduct+isotop from peakList with their mz and rt
      
      tempcai <- paste(peakList$comp_name, peakList$adduct, peakList$isotopologue, sep = "|")
      allComps <- unique(tempcai)
      allSamps <- basename(rawFiles[indices])

      cc <- rep(allComps, times = length(allSamps))
      ss <- rep(allSamps, each = length(allComps))
      # check if results are present in peakList, if yes, copy these to integRes, if not
      # leave as NA, need to be processed later
      getPreDat <- function(nm, sa) {
        x <- strsplit(nm, "\\|")[[1]]
        cnm <- x[1]
        ad <- x[2]
        isot <- x[3]
        r <- peakList[
          peakList$comp_name == cnm & 
            peakList$adduct == ad &
            peakList$isotopologue == isot &
            peakList$samp == basename(sa), ]
        if (nrow(r) >= 1) {
          return(data.frame(
            samp = sa, comp_name = cnm, 
            adduct = ad, isotopologue = isot,
            int_h = r$int_h,
            int_a = r$int_a, s_to_n = r$s_to_n,
            rt_error_min = r$rt_error_min,
            eic_extraction_width = r$eic_extraction_width,
            real_mz = r$real_mz, real_rt_min = ifelse(is.na(r$real_rt_min), r$rt_min, r$real_rt_min),
            stringsAsFactors = FALSE
          ))
        } else {
          return(data.frame(
            samp = sa, comp_name = cnm, adduct = ad, isotopologue = isot, int_h = NA,
            int_a = NA, s_to_n = NA, rt_error_min = NA,
            eic_extraction_width = NA,
            real_mz = NA, real_rt_min = NA,
            stringsAsFactors = FALSE
          ))
        }
      }
      preDatL <- mcmapply(getPreDat, cc, ss,
        SIMPLIFY = FALSE,
        mc.cores = settings$numcores
      )
      preDat <- do.call("rbind", preDatL)

      # which compounds and samples have no area
      toProc <- preDat[is.na(preDat$int_a), c("comp_name", "samp", "adduct", "isotopologue")]

      if (nrow(toProc) == 0) {
        message("no additional peaks to integrate")
        integRes <<- rbind(integRes, preDat)
        return(NULL)
      }

      # Get the remaining results
      perDatFile <- function(proc) {
        thisSamp <- proc$samp[1]

        # if the file is not in memory, load this file, and remember to remove it again
        if (thisSamp %in% names(rawData)) {
          clearData <- FALSE
        } else {
          .self$loadData(indices = which(rawFiles == thisSamp))
          clearData <- TRUE
        }

        rawLink <- rawData[[thisSamp]]

        getCompDat <- function(thisCai) {
          x <- strsplit(thisCai, "\\|")[[1]]
          thisComp <- x[1]
          thisAd <- x[2]
          thisIsot <- x[3]
          # find best mz
          thisMz <- mean(
            peakList[
              peakList$comp_name == thisComp & 
              peakList$adduct == thisAd & 
              peakList$isotopologue == thisIsot
              , "real_mz"
            ], na.rm = TRUE)
          # find best rt
          allRt <- peakList[
            peakList$comp_name == thisComp & 
            peakList$adduct == thisAd & 
            peakList$isotopologue == thisIsot, "real_rt_min"]
          thisRt <- mean(allRt, na.rm = TRUE)
          # if no rt found, could be because no peak integrated, use rt from
          if (is.na(thisRt)) {
            allRt <- peakList[peakList$comp_name == thisComp, "rt_min"]
            thisRt <- mean(allRt, na.rm = TRUE)
          }
          newRtTol <- max(sd(allRt, na.rm = TRUE) * 3, settings$rtTolReinteg, na.rm = TRUE)

          minInd <- which.min(abs(rawLink@scantime - (thisRt * 60 - (settings$rtTolReinteg * 60) / 2)))
          maxInd <- which.min(abs(rawLink@scantime - (thisRt * 60 + (settings$rtTolReinteg * 60) / 2)))
          # integrate peak
          compRes <- .self$getPeak(
            rawLink, thisMz, thisRt, minInd, maxInd,
            settings$EIC_extraction,
            settings$mzTolReinteg, newRtTol
          )
          if (nrow(compRes) == 0) {
            return(NULL)
          }
          # get highest peak
          compRes <- compRes[which.max(compRes$peak_intens), ]
          data.frame(
            samp = basename(thisSamp), 
            comp_name = thisComp,
            adduct = thisAd,
            isotopologue = thisIsot,
            int_h = as.integer(round(compRes$peak_intens)),
            int_a = as.integer(round(compRes$peakArea)),
            s_to_n = round(compRes$peak_intens / compRes$noisedeviation, 1),
            rt_error_min = round(compRes$scantime / 60 - thisRt, 2),
            eic_extraction_width = compRes$e_width,
            real_mz = thisMz, real_rt_min = thisRt,
            stringsAsFactors = FALSE
          )
        }
        proc$cai <- paste(proc$comp_name, proc$adduct, proc$isotopologue, sep = "|")
        datL <- lapply(proc$cai, getCompDat)
        if (all(is.null(datL))) {
          return(NULL)
        }
        datL <- compact(datL)
        dat <- do.call("rbind", datL)
        rm(rawLink)
        if (clearData) {
          .self$clearData(indices = which(rawFiles == thisSamp))
        }
        dat
      }

      splitBySamp <- split(toProc, toProc$samp)
      newResL <- mclapply(splitBySamp, perDatFile, mc.cores = settings$numcores)
      if (all(is.null(newResL))) {
        message("no additional peaks found")
        return(NULL)
      }
      newResL <- compact(newResL)
      newRes <- do.call("rbind", newResL)

      # combine results with those already made
      preDat$samp <- basename(preDat$samp)
      allDat <- rbind(preDat, newRes)
      # remove NAs
      allDat <- allDat[!is.na(allDat$int_a), ]
      # remove peaks below threshold
      allDat <- switch(settings$use_int_threshold,
        height = allDat[allDat$int_h >= settings$height_threshold, ],
        area = allDat[allDat$int_a >= settings$area_threshold, ],
        stop("use_int_threshold can only be area or height")
      )

      integRes <<- rbind(integRes, allDat)
      rownames(integRes) <<- NULL
    },
    clearAndSave = function(dialog = TRUE, nameReport = NULL, clearData = TRUE) {
      "Clear data (optionally) from RAM and save report as .report file in the
      current working directory.
      Use nameReport to give a different location and different name as in
      *.clearAndSave(F, 'D:/exampleFolder/example'). Note: The folder must exist beforehand.
      To read the file again, use the function ntsworkflow::loadReport"

      if (dialog) {
        nameReport <- rstudioapi::selectFile("Save File as...",
          label = "Save", existing = FALSE,
          filter = "DBscreening report file (*.report)"
        )
        if (!grepl("\\.report$", nameReport)) {
          nameReport <- paste0(nameReport, ".report")
        }
      } else if (is.null(nameReport)) {
        nameReport <- stringr::str_match(deparse(sys.call()), "^(.*)\\$clearAndSave")[, 2]
        nameReport <- paste0(nameReport, ".report")
      }
      if (clearData)
        .self$clearData()
      saveRDS(.self, file = nameReport)
    },
    view = function() {
      "View results using shiny."

      require(shiny)
      require(shinyFiles)
      require(ggplot2)

      nameReport <- stringr::str_match(deparse(sys.call()), "^(.*)\\$view")[, 2]

      app <- shinyApp(
        ui = navbarPage(
          paste("dbscreening -", nameReport),
          # Overview ####
          tabPanel(
            "Overview",

            # Show a table of info

            verticalLayout(
              tableOutput("susSFileInfo"),
              hr(),
              DT::dataTableOutput("dataFileInfo")
            )
          ),
          tabPanel(
            "Peak list",
            verticalLayout(DT::dataTableOutput("peakList"))
          ),
          tabPanel(
            "View",
            verticalLayout(
              numericInput("peakIDfield", "Peak ID", 1, 1, NA, 1),
              plotOutput(
                "eicPlot",
                dblclick = "eicPlot_dblclick",
                brush = brushOpts(
                  id = "eicPlot_brush",
                  resetOnNew = TRUE
                )
              ),
              plotOutput(
                "ms1Plot",
                dblclick = "ms1Plot_dblclick",
                brush = brushOpts(
                  id = "ms1Plot_brush",
                  resetOnNew = TRUE
                )
              ),
              plotOutput(
                "ms2Plot",
                dblclick = "ms2Plot_dblclick",
                brush = brushOpts(
                  id = "ms2Plot_brush",
                  resetOnNew = TRUE
                )
              ),
              absolutePanel(plotOutput(outputId = "structure"),
                top = 0, height = 100,
                width = 350, left = "25%"
              )
            )
          ),
          tabPanel(
            "IS",
            verticalLayout(
              DT::dataTableOutput("ISlist"),
              DT::dataTableOutput("ISresults")
            )
          ),
          tabPanel(
            "View IS",
            verticalLayout(
              selectInput("ISnameChosen",
                label = "IS name",
                choices = unique(.self$IS$name)
              ),
              splitLayout(
                plotOutput("ISarea",
                  dblclick = "ISarea_dblclick",
                  brush = brushOpts(
                    id = "ISarea_brush",
                    resetOnNew = TRUE
                  )
                ),
                plotOutput("ISheight",
                  dblclick = "ISheight_dblclick",
                  brush = brushOpts(
                    id = "ISheight_brush",
                    resetOnNew = TRUE
                  )
                )
              ),
              splitLayout(
                plotOutput("ISwidth",
                  dblclick = "ISwidth_dblclick",
                  brush = brushOpts(
                    id = "ISwidth_brush",
                    resetOnNew = TRUE
                  )
                ),
                plotOutput("ISrt",
                  dblclick = "ISrt_dblclick",
                  brush = brushOpts(
                    id = "ISrt_brush",
                    resetOnNew = TRUE
                  )
                )
              )
            )
          ),
          tabPanel(
            "Trend",
            fluidRow(
              column(
                4,
                selectInput("chooseCompTrend", "Choose compound", choices = unique(.self$peakList$comp_name))
              ),
              column(
                4,
                selectInput("chooseDisplayTrend", "Metric", c(
                  "Rel._height",
                  "Rel._area",
                  "Height",
                  "Area"
                ))
              ),
              column(
                4,
                selectInput("chooseISTrend", "Choose IS", choices = unique(.self$IS$name))
              )
            ),
            fluidRow(
              column(
                12,
                plotOutput("trendPlot",
                  dblclick = "trend_dblclick",
                  brush = brushOpts(
                    id = "trend_brush",
                    resetOnNew = TRUE
                  )
                )
              )
            )
          ),
          tabPanel(
            "Settings",
            tags$pre(textOutput("settingsOut"))
          )
        ),
        server = function(input, output, session) {
          home <- c("Home" = "~")

          shinyFileChoose(input, "dbFile", roots = home, session = session, filetypes = "db")
          dbLocation <- reactive({
            newPath <- as.character(parseFilePaths(home, input$dbFile)$datapath)
            if (length(newPath) != 0) {
              .self$changeSettings("db_path", newPath)
            }
            newPath
          })

          output$dbLoc <- renderText(dbLocation())


          output$susSFileInfo <- renderTable(
            {
              if (nrow(.self$rawFilesCompl) >= 1) {
                lastDate <- max(.self$rawFilesCompl$date)
              } else {
                lastDate <- NA
              }
              # make dataframe with info
              data.frame(
                Last_modification = lastDate,
                No._datafiles = length(.self$rawFiles),
                No._processed_datafiles = nrow(.self$rawFilesCompl),
                No._peaks = nrow(.self$peakList), stringsAsFactors = FALSE
              )
            },
            align = "c",
            spacing = "s"
          )

          output$dataFileInfo <- DT::renderDataTable(
            {
              isComplete <- .self$rawFiles %in% .self$rawFilesCompl$path
              df <- data.frame(
                Files = .self$rawFiles,
                Complete = isComplete,
                stringsAsFactors = FALSE
              )
              df$Process_date <- ifelse(df$Complete,
                .self$rawFilesCompl[.self$rawFilesCompl$path == df$Files, "date"],
                NA
              )
              df$Files <- stringr::str_wrap(df$Files, 20)
              df
            },
            options = list(pageLength = 25)
          )

          output$peakList <- DT::renderDataTable(
            {
              data.frame(
                mz = round(.self$peakList$real_mz, 4),
                rt = round(.self$peakList$real_rt_min, 2),
                area = .self$peakList$int_a,
                hght = round(.self$peakList$int_h),
                name = .self$peakList$comp_name,
                CAS = .self$peakList$comp_CAS,
                MS2 = round(.self$peakList$score),
                dvmDa = round(.self$peakList$mz_error_mDa, 2),
                dvRt = round(.self$peakList$rt_error_min, 2),
                pkWdth = round(.self$peakList$peak_end - .self$peakList$peak_start, 2),
                sn = .self$peakList$s_to_n,
                sample = stringr::str_match(.self$peakList$samp, "^(.*)\\.mzX?ML$")[, 2],
                dup = .self$peakList$duplicate,
                stringsAsFactors = FALSE, row.names = .self$peakList$peakID
              )
            },
            selection = "single",
            options = list(pageLength = 50)
          )

          peakIdPlot <- reactive({
            sel <- input$peakList_rows_selected
            if (is.null(sel)) {
              sel <- 1
            } else {
              sel <- .self$peakList[sel, "peakID"]
            }
            sel
          })

          # get compound table as reactive value once
          compTable <- reactive({
            db <- DBI::dbConnect(RSQLite::SQLite(), .self$settings$db_path)
            result <- dplyr::collect(dplyr::tbl(db, "compound"))
            DBI::dbDisconnect(db)
            result
          })

          output$structure <- renderPlot({
            newId <- input$peakIDfield
            cpTbl <- compTable()
            compId <- .self$peakList[.self$peakList$peakID == newId, "comp_id"]
            newSmiles <- cpTbl[cpTbl$compound_id == compId, "SMILES"]
            mol <- rcdk::parse.smiles(as.character(newSmiles), kekulise = TRUE)[[1]]
            img <- rcdk::view.image.2d(molecule = mol)
            graphics::plot.window(c(0, 1), c(0, 1))
            graphics::rasterImage(img, 0, 0, 1, 1)
          })

          observeEvent(input$peakList_rows_selected, {
            updateNumericInput(session, "peakIDfield",
              value = peakIdPlot(),
              min = min(.self$peakList$peakID),
              max = max(.self$peakList$peakID)
            )
          })



          rangesEIC <- reactiveValues(x = NULL, y = NULL)
          rangesMS1 <- reactiveValues(x = NULL, y = NULL)
          rangesMS2 <- reactiveValues(x = NULL, y = NULL)

          output$eicPlot <- renderPlot({
            newPlot <- .self$plotEIC(input$peakIDfield)
            if (!is.null(rangesEIC$x)) {
              newPlot <- newPlot +
                coord_cartesian(xlim = rangesEIC$x, ylim = rangesEIC$y, expand = FALSE)
            }
            newPlot
          })
          output$ms1Plot <- renderPlot({
            newPlot <- .self$plotMS1(input$peakIDfield)
            if (!is.null(rangesMS1$x)) {
              newPlot <- newPlot +
                coord_cartesian(xlim = rangesMS1$x, ylim = rangesMS1$y, expand = FALSE)
            }
            newPlot
          })
          output$ms2Plot <- renderPlot({
            newPlot <- .self$plotMS2(input$peakIDfield)
            if (!is.null(rangesMS2$x)) {
              newPlot <- newPlot +
                coord_cartesian(xlim = rangesMS2$x, ylim = rangesMS2$y, expand = FALSE)
            }
            newPlot
          })

          observeEvent(input$eicPlot_dblclick, {
            brush <- input$eicPlot_brush
            if (!is.null(brush)) {
              rangesEIC$x <- c(brush$xmin, brush$xmax)
              rangesEIC$y <- c(brush$ymin, brush$ymax)
            } else {
              rangesEIC$x <- NULL
              rangesEIC$y <- NULL
            }
          })

          observeEvent(input$ms1Plot_dblclick, {
            brush <- input$ms1Plot_brush
            if (!is.null(brush)) {
              rangesMS1$x <- c(brush$xmin, brush$xmax)
              rangesMS1$y <- c(brush$ymin, brush$ymax)
            } else {
              rangesMS1$x <- NULL
              rangesMS1$y <- NULL
            }
          })

          observeEvent(input$ms2Plot_dblclick, {
            brush <- input$ms2Plot_brush
            if (!is.null(brush)) {
              rangesMS2$x <- c(brush$xmin, brush$xmax)
              rangesMS2$y <- c(brush$ymin, brush$ymax)
            } else {
              rangesMS2$x <- NULL
              rangesMS2$y <- NULL
            }
          })

          output$ISlist <- DT::renderDataTable(
            {
              data.frame(
                IS_name = .self$IS$name,
                Formula = .self$IS$formula,
                RT = .self$IS$rt,
                stringsAsFactors = FALSE
              )
            },
            selection = "single"
          )

          output$ISresults <- DT::renderDataTable(
            {
              data.frame(
                IS = .self$ISresults$IS,
                mz = round(.self$ISresults$mz, 4),
                rt = round(.self$ISresults$rt, 2),
                height = round(.self$ISresults$int_h),
                area = round(.self$ISresults$int_a),
                width_min = round(.self$ISresults$peak_end - .self$ISresults$peak_start, 2),
                sample = stringr::str_wrap(.self$ISresults$samp, 10),
                stringsAsFactors = FALSE
              )
            },
            options = list(pageLength = 25)
          )

          # get some general plots about the behavior of IS, barplots over all samples, similar to check_IS
          ISnameSelected <- reactive({
            selIS <- input$ISlist_rows_selected
            if (is.null(selIS)) {
              selIS <- .self$IS[1, "name"]
            } else {
              selIS <- .self$IS[selIS, "name"]
            }
            selIS
          })


          observeEvent(
            {
              input$ISlist_rows_selected
              input$dataFile
            },
            {
              updateSelectInput(session, "ISnameChosen", selected = ISnameSelected())
            }
          )

          trendBreaks <- function(limits) {
            if (length(limits) <= 40) {
              limits
            } else {
              limits[seq(1, length(limits), length.out = 40)]
            }
          }
          trendLabels <- function(breaks) {
            stringr::str_match(breaks, "^(.*)\\.mzX?ML$")[, 2]
          }
          rangesIS <- reactiveValues(x = NULL, y = NULL)
          output$ISarea <- renderPlot({
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            newPlot <- ggplot(toPlot, aes(x = samp, y = int_a)) +
              geom_bar(stat = "identity") +
              scale_x_discrete(
                limits = basename(.self$rawFiles), breaks = trendBreaks,
                labels = trendLabels
              ) +
              ylab("Intensity (area, cps)") +
              theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5),
                axis.title.x = element_blank()
              )
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(
                  limits = basename(.self$rawFiles)[from:to],
                  breaks = trendBreaks, labels = trendLabels
                )
            }
            newPlot
          })

          observeEvent(input$ISarea_dblclick, {
            brush <- input$ISarea_brush
            if (!is.null(brush)) {
              rangesIS$x <- c(brush$xmin, brush$xmax)
              rangesIS$y <- c(brush$ymin, brush$ymax)
            } else {
              rangesIS$x <- NULL
              rangesIS$y <- NULL
            }
          })
          rangesISheight <- reactiveValues(x = NULL, y = NULL)
          output$ISheight <- renderPlot({
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            newPlot <- ggplot(toPlot, aes(x = samp, y = int_h)) +
              geom_bar(stat = "identity") +
              scale_x_discrete(
                limits = basename(.self$rawFiles), breaks = trendBreaks,
                labels = trendLabels
              ) +
              ylab("Intensity (height, cps)") +
              theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5),
                axis.title.x = element_blank()
              )
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(
                  limits = basename(.self$rawFiles)[from:to],
                  breaks = trendBreaks, labels = trendLabels
                )
            }
            newPlot
          })
          observeEvent(input$ISheight_dblclick, {
            brush <- input$ISheight_brush
            if (!is.null(brush)) {
              rangesIS$x <- c(brush$xmin, brush$xmax)
              rangesIS$y <- c(brush$ymin, brush$ymax)
            } else {
              rangesIS$x <- NULL
              rangesIS$y <- NULL
            }
          })

          output$ISwidth <- renderPlot({
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            toPlot$width <- (toPlot$peak_end - toPlot$peak_start)

            newPlot <- ggplot(toPlot, aes(x = samp, y = width, group = 1)) +
              geom_line() +
              geom_point() +
              scale_x_discrete(
                limits = basename(.self$rawFiles), breaks = trendBreaks,
                labels = trendLabels
              ) +
              ylim(0, NA) +
              ylab("Peak width (min.)") +
              theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5),
                axis.title.x = element_blank()
              )
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(
                  limits = basename(.self$rawFiles)[from:to],
                  breaks = trendBreaks, labels = trendLabels
                )
            }
            newPlot
          })

          output$ISrt <- renderPlot({
            toPlot <- .self$ISresults[.self$ISresults$IS == input$ISnameChosen, ]
            newPlot <- ggplot(toPlot, aes(x = samp, y = rt, group = 1)) +
              geom_line() +
              geom_point() +
              scale_x_discrete(
                limits = basename(.self$rawFiles), breaks = trendBreaks,
                labels = trendLabels
              ) +
              ylim(0, NA) +
              ylab("RT (min., cps)") +
              theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5),
                axis.title.x = element_blank()
              )
            if (!is.null(rangesIS$x)) {
              from <- round(rangesIS$x[1])
              to <- round(rangesIS$x[2])
              newPlot <- newPlot +
                scale_x_discrete(
                  limits = basename(.self$rawFiles)[from:to],
                  breaks = trendBreaks, labels = trendLabels
                )
            }
            newPlot
          })
          observeEvent(input$ISrt_dblclick, {
            brush <- input$ISrt_brush
            if (!is.null(brush)) {
              rangesIS$x <- c(brush$xmin, brush$xmax)
              rangesIS$y <- c(brush$ymin, brush$ymax)
            } else {
              rangesIS$x <- NULL
              rangesIS$y <- NULL
            }
          })


          trendBreaks80 <- function(limits) {
            if (length(limits) <= 80) {
              limits
            } else {
              limits[seq(1, length(limits), length.out = 80)]
            }
          }
          rangesTrend <- reactiveValues(x = NULL, y = NULL)
          locShade <- reactiveValues(x = NULL)


          observeEvent(input$peakList_rows_selected, {
            updateSelectInput(session, "chooseCompTrend",
              selected = .self$peakList[.self$peakList$peakID == peakIdPlot(), "comp_name"]
            )
          })

          observeEvent(input$ISlist_rows_selected, {
            updateSelectInput(
              session, "chooseISTrend",
              selected = .self$IS[.self$IS$name == ISnameSelected(), "name"]
            )
          })

          prepareTrendPlot <- reactive({
            comp <- input$chooseCompTrend
            intStd <- input$chooseISTrend

            trendPlotRel <- function(overview, reInt, overviewAppr, reIntAppr, ycol, yname) {
              gg <- ggplot(overview, aes_(x = quote(samp), y = as.name(ycol))) +
                geom_col(width = 0.5) +
                geom_col(
                  data = reInt, aes_(quote(samp), as.name(ycol)), width = 0.5,
                  fill = "blue", alpha = .5
                ) +
                geom_col(
                  data = overviewAppr,
                  aes_(quote(samp), as.name(ycol)), width = 0.5, fill = "darkred"
                ) +
                geom_col(
                  data = reIntAppr, aes_(quote(samp), as.name(ycol)), width = 0.5,
                  fill = "darkred", alpha = .5
                ) +
                ylab(yname) +
                scale_x_discrete(
                  limits = basename(.self$rawFiles), breaks = trendBreaks80,
                  labels = trendLabels
                ) +
                theme(axis.title.x = element_blank())
              if (any(nchar(trendLabels(.self$rawFiles)) >= 3)) {
                gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
              }
              gg
            }

            trendPlotAbs <- function(df, reInt, ycol, yname) {
              gg <- ggplot(df, aes_(x = quote(samp), y = as.name(ycol))) +
                geom_col(width = 0.5) +
                geom_col(data = reInt, aes_(quote(samp), as.name(ycol)), width = 0.5, fill = "blue", alpha = .5) +
                ylab(yname) +
                scale_x_discrete(
                  limits = basename(.self$rawFiles), breaks = trendBreaks80,
                  labels = trendLabels
                ) +
                theme(axis.title.x = element_blank())
              if (any(nchar(trendLabels(.self$rawFiles)) >= 3)) {
                gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
              }
              gg
            }
            switch(input$chooseDisplayTrend,
              Rel._height = {
                comp_selected <- .self$peakList[.self$peakList$comp_name == comp, ]
                IS_selected <- .self$ISresults[.self$ISresults$IS == intStd, ]
                IS_all <- merge(
                  data.frame(samp = basename(.self$rawFiles)), IS_selected,
                  by = "samp", all.x = TRUE
                )
                # approximate IS concentration if some are missing
                IS_all$approx_int_h <- is.na(IS_all$int_h)
                if (any(IS_all$approx_int_h)) {
                  IS_all[is.na(IS_all$int_h), "int_h"] <- approx(
                    seq_along(IS_all$samp), IS_all$int_h,
                    xout = which(is.na(IS_all$int_h))
                  )$y
                }

                overview <- merge(comp_selected, IS_all, by = "samp", all.x = TRUE)
                overview$norm_int_h <- overview$int_h.x / overview$int_h.y

                # data for no MS2
                comp_reInt <- .self$integRes[.self$integRes$comp_name == comp, ]
                reInt <- merge(comp_reInt, IS_all, by = "samp", all.x = TRUE)
                reInt$norm_int_h <- reInt$int_h.x / reInt$int_h.y

                # from overview and reInt, cut out approximated values
                overviewAppr <- overview[overview$approx_int_h, ]
                reIntAppr <- reInt[reInt$approx_int_h, ]
                overview <- overview[!overview$approx_int_h, ]
                reInt <- reInt[!reInt$approx_int_h, ]

                # for samples where more than one peak is present, choose the one with the best MS2
                o <- by(overview, overview$samp, function(x) x[order(x$score, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                overview <- Reduce(rbind, o)
                # since there is no MS2 for reintegrated peaks, choose the most intense...
                r <- by(reInt, reInt$samp, function(x) x[order(x$int_h.x, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                reInt <- Reduce(rbind, r)

                trendPlotRel(overview, reInt, overviewAppr, reIntAppr,
                  ycol = "norm_int_h", yname = "Rel. intensity (height)"
                )
              },
              Rel._area = {
                comp_selected <- .self$peakList[.self$peakList$comp_name == comp, ]
                IS_selected <- .self$ISresults[.self$ISresults$IS == intStd, ]
                # approximate IS concentration if some are missing
                IS_all <- merge(
                  data.frame(samp = basename(.self$rawFiles)), IS_selected,
                  by = "samp", all.x = TRUE
                )
                IS_all$approx_int_a <- is.na(IS_all$int_a)
                if (any(IS_all$approx_int_a)) {
                  IS_all[is.na(IS_all$int_a), "int_a"] <- approx(
                    seq_along(IS_all$samp), IS_all$int_a,
                    xout = which(is.na(IS_all$int_a))
                  )$y
                }
                overview <- merge(comp_selected, IS_all, by = "samp", all.x = TRUE)
                overview$norm_int_a <- overview$int_a.x / overview$int_a.y

                # data for no MS2
                comp_reInt <- .self$integRes[.self$integRes$comp_name == comp, ]
                reInt <- merge(comp_reInt, IS_all, by = "samp", all.x = TRUE)
                reInt$norm_int_a <- reInt$int_a.x / reInt$int_a.y

                # from overview and reInt, cut out approximated values
                overviewAppr <- overview[overview$approx_int_a, ]
                reIntAppr <- reInt[reInt$approx_int_a, ]
                overview <- overview[!overview$approx_int_a, ]
                reInt <- reInt[!reInt$approx_int_a, ]

                # for samples where more than one peak is present, choose the one with the best MS2
                o <- by(overview, overview$samp, function(x) x[order(x$score, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                overview <- Reduce(rbind, o)
                # since there is no MS2 for reintegrated peaks, choose the most intense
                r <- by(reInt, reInt$samp, function(x) x[order(x$int_a.x, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                reInt <- Reduce(rbind, r)


                trendPlotRel(overview, reInt, overviewAppr, reIntAppr,
                  ycol = "norm_int_a", yname = "Rel. intensity (area)"
                )
              },
              Height = {
                overview <- .self$peakList
                df <- overview[overview$comp_name == comp, ]
                reInt <- .self$integRes
                reInt <- reInt[reInt$comp_name == comp, ]

                # for samples where more than one peak is present, choose the one with the best MS2
                o <- by(df, df$samp, function(x) x[order(x$score, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                df <- Reduce(rbind, o)
                # since there is no MS2 for reintegrated peaks, choose the most intense
                r <- by(reInt, reInt$samp, function(x) x[order(x$int_h, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                reInt <- Reduce(rbind, r)

                trendPlotAbs(df, reInt, "int_h", "Intensity (height)")
              },
              Area = {
                overview <- .self$peakList
                df <- overview[overview$comp_name == comp, ]
                reInt <- .self$integRes
                reInt <- reInt[reInt$comp_name == comp, ]

                # for samples where more than one peak is present, choose the one with the best MS2
                o <- by(df, df$samp, function(x) x[order(x$score, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                df <- Reduce(rbind, o)
                # since there is no MS2 for reintegrated peaks, choose the most intense
                r <- by(reInt, reInt$samp, function(x) x[order(x$int_a, decreasing = TRUE), ][1, ],
                  simplify = FALSE
                )
                reInt <- Reduce(rbind, r)

                trendPlotAbs(df, reInt, "int_a", "Intensity (area)")
              }
            )
          })

          output$trendPlot <- renderPlot({
            p <- prepareTrendPlot()

            if (!is.null(rangesTrend$x)) {
              from <- round(rangesTrend$x[1])
              to <- round(rangesTrend$x[2])
              p <- p + scale_x_discrete(
                limits = basename(.self$rawFiles)[from:to],
                breaks = trendBreaks80, labels = trendLabels
              )
              if (!is.null(locShade$x)) {
                pos <- round(length(from:to) / 2)
                p <- p + annotate("rect",
                  xmin = pos - 0.5, xmax = pos + 0.5,
                  ymin = -Inf, ymax = Inf, alpha = .1, fill = "orange"
                )
                locShade$x <- which(.self$rawFiles == .self$rawFiles[from:to][pos])
              }
            }
            p
          })
          observeEvent(input$trend_dblclick, {
            brush <- input$trend_brush
            if (!is.null(brush)) {
              rangesTrend$x <- c(brush$xmin, brush$xmax)
              rangesTrend$y <- c(brush$ymin, brush$ymax)
            } else {
              rangesTrend$x <- NULL
              rangesTrend$y <- NULL
            }

            locShade$x <- round(mean(c(brush$xmin, brush$xmax)))
          })

          peakIdTrend <- reactive({
            # find sample corresponding to x

            fl <- basename(.self$rawFiles[locShade$x])
            # find ID of peak
            pl <- .self$peakList
            # if more than one ID matches, the one with the highest MS2 score is shown.
            idRes <- pl[pl$samp == fl & pl$comp_name == input$chooseCompTrend, "peakID"]
            if (length(idRes) > 1) {
              plNew <- pl[pl$peakID == idRes, ]
              plNew <- plNew[order(plNew$score, decreasing = TRUE), ]
              plNew$peakID[1]
            } else {
              idRes
            }
          })

          observeEvent(input$trend_dblclick, {
            if (length(peakIdTrend()) != 0L) {
              updateNumericInput(session, "peakIDfield",
                value = peakIdTrend(),
                min = min(.self$peakList$peakID),
                max = max(.self$peakList$peakID)
              )
            }
          })


          output$trendSampleTable <- renderTable({
            data.frame(
              No. = seq_along(.self$rawFiles),
              Sample_Name = stringr::str_match(basename(.self$rawFiles), "^(.*)\\.mzX?ML$")[, 2],
              stringsAsFactors = FALSE
            )
          })

          output$settingsOut <- renderPrint({
            .self$settings
          })
        },
        onStart = function() {
          message(paste("Currently viewing", nameReport))
          onStop(function() {
            message("Visualization closed")
            detach("package:shinyFiles", unload = TRUE)
          })
        }
      )
      runApp(app)
    }
  )
)

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow
