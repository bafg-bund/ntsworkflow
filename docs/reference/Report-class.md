# A reference class for database screening.

This is a class for conducting non-target screening data evaluation
based on database searching by m/z, retention time and MS2 spectral
comparison. View the vignette "Database screening with Report in
ntsworkflow" for a step-by-step guide on how to use this function.

## Details

Initialize a Report by calling Report\$new() this starts a new Report
with the default settings (see settings field, below). Settings can be
changed by calling the `changeSettings` method. Load files by calling
the `addRawFiles` method and process files with the `process_all`
method. After processing the Report can be saved for viewing in the
visuallization tool (susS app).

## Fields

- `rawFiles`:

  A character vector with file paths in the order in which they are to
  be displayed.

- `settings`:

  A list with settings for all processing steps. The settings list
  contains the following fields and defaults:
  `list(db_path = NULL, rttolm = 1, mztolu = 0.05, mztolu_fine = 0.005, chromatography = "bfg_nts_rp1", pol = "pos", ce_s = 30:40, ces_s = 0:15, instr = "LC-ESI-QTOF TripleTOF 5600 SCIEX", ceunit = c("V", "eV"), comparison = "dot_product", threshold = 400, rt_res = 1.5, EIC_extraction = 0.05, baseline_noise_MS1 = 0.5, sn = 3, rtoffset = 0, IS_rtoffset = 0, ISrttolm = 1, blank_int_factor = 3, rtTolReinteg = 1, mzTolReinteg = 0.005, ndp_m = 0.4, ndp_n = 1, mztolu_ms2 = 0.015, area_threshold = 1, height_threshold = 1, use_int_threshold = "area"`.
  `db_path`: path to sqlite spectral-database (SDB). `rttolm`: retention
  time tolerance for the suspect search. `mztolu`: m/z tolerance for the
  spectrum extraction (MS2 precursor mass). `mztolu_fine`: m/z tolerance
  for precursor mass in MS1. `chromatography` Is the chromatographic
  method allowed from the SDB. `pol`: polarity. `ce_s`: allowed
  collision energies from SDB. `ces_s`: allowed collision energy spreads
  from SDB. `instr`: allowed instruments from SDB. `ceunit`: allowed
  collision energy unit from SDB. `comparison`: algorithm to use for
  spectral comparison, currently only the default is allowed.
  `threshold`: score of spectral comparison under which results are
  rejected automatically, 1000=perfect match. `rt_res`: Resolution of
  chrom. peaks (see
  [`ms2_search`](http://r.bafg.de/reference/ms2_search.md)) in min.
  Peaks with RT difference \> `rt_res` are considered to be from
  different substances. `EIC_extraction`: Extraction width to produce
  EIC, affects integration. `baseline_noise_MS1`: Signals under this
  intensity are ignored and also not included in the saved spectra.
  `sn`: signal-to-noise ratio limit for peak integration. `rtoffset`:
  Retention time offset between samples and database. `ISrttolm`: RT
  tolerance for IS peak finding and integration. `blank_int_factor`:
  Intensity factor to remove compounds found in blank. Peak in sample
  must be within `blank_int_factor`: \* intensity in blank.
  `rtTolReinteg`: is the retention time tolerance for reintegration.
  `mzTolReinteg`: is the m/z tolerance for reintegration. `ndp_m`: Peak
  intensity weighting factor for dot-product, `ndp_n`: m/z weighting
  factor for dot-product, `mztolu_ms2`: m/z tolerance for dot-product
  fragment mass binning, `area_threshold`: Peak area intensity
  threshold. `height_threshold`: Peak height intensity threshold.
  `use_int_threshold`: can be either "area", "height", or "none"
  `peaksPerPeak`: `mustFindChromPeak`: default FALSE, if TRUE, peaks
  without area are deleted

- `rawFilesCompl`:

  `data.frame` with processed files and date of processing.

- `peakList`:

  `data.frame` with all suspect search results.

- `currentPeakID`:

  numeric to keep track of peakIDs, do not change.

- `MS1`:

  `data.frame` holding all MS1 spectra.

- `EIC`:

  `data.frame` holding all chromatograms.

- `MS2`:

  `data.frame` holding all MS2 spectra.

- `IS`:

  `data.frame` with list of IS to process, this has to be imported from
  a csv file (comma sep) using this method `addIS()`. The csv file
  should have the columns name, formula, rt, adduct. Where formula is in
  the form e.g. for CBZ: C14 13CH12N 15NO, and adduct is in the form
  \[M+H\]+ or \[M\]+.

- `ISresults`:

  `data.frame` with results of IS processings.

- `falsePos`:

  `data.frame` recording which false positives should be deleted in
  which files

## Methods

- `addIS(dialog = TRUE, fileName = NULL)`:

  Include a list of internal standards in the report. Dialog indicates
  use of interactive file choosing dialog. The file should be a csv with
  4 columns: 'name', 'formula', 'rt', 'adduct'.

- `addRawFiles(dialog = TRUE, file_list = NULL)`:

  Add raw files, (mzML or mzXML) Dialog indicates use of interactive
  file choosing dialog. Does not load files to RAM, only marks file
  path.

- `addRawFilesDir(dialog = TRUE, dir_path = NULL)`:

  Add all raw files in the chosen directory

- `changeSettings(parameter, value)`:

  Change any setting by name and then value. See settings field for more
  details.

- `cleanDuplicateCol()`:

  Check the peaklist for duplicates, if peaks previously indicated as
  duplicates are no longer so then set duplicate column back to NA

- `cleanPeakLetterCol()`:

  Check the peaklist for multiple peaks from one compound, if any B, C
  etc. peaks no longer have an A peak, these are renamed to start with A

- `clearAndSave(dialog = TRUE, nameReport = NULL, clearData = TRUE)`:

  Clear data (optionally) from RAM and save report as .report file in
  the current working directory. Use nameReport to give a different
  location and different name as in \*.clearAndSave(F,
  'D:/exampleFolder/example'). Note: The folder must exist beforehand.
  To read the file again, use the function ntsworkflow::loadReport

- `clearData(indices = NULL)`:

  Remove data from RAM to clear memory, use indices of raw files

- `deleteBackground(indicesData, indicesBlank, includeIS = FALSE)`:

  Delete background signals based on intensity. Use indices to indicate
  blanks and samples. If there is more than one blank then these will be
  combined. If IS should also be deleted use includeIS. Removal is based
  on intensity ratio, see settings.

- `deleteFP(substanceName, indices = 0)`:

  Delete all occurances of a compound given by name. Only one compound
  can be given. Name must match exactly. Compound is added to list of
  falsePos (see fields). If indices = 0, then compound is always
  deleted, in all samples and also in future processed samples.
  Otherwise it is only deleted in the specified samples (by index).

- `getPeak(rawLinki, comp_mzi, comp_rti, minIndi, maxIndi, width, mztoli, rttoli)`:

  Internal function to integrate peaks during processing

- `loadData(all = FALSE, indices = NULL)`:

  Load files that still need to be processed into RAM for fast access,
  if all = TRUE then all files are loaded regardless, if indices is an
  integer vector, only these samps will be loaded

- `loadSettings()`:

  Load previously settings file from current working directory. This
  will fail if there is more than one file present in the directory.

- `moveRawFile(index, direction)`:

  Move a raw file in direction 'up' or 'down'

- `plotEIC(ID)`:

  Plot a particular EIC using the ID from the peak list table

- `plotMS1(ID)`:

  Plot a MS1 spectrum using the ID from the peak list table. Includes
  comparison to calculated isotopic pattern.

- `plotMS2(ID)`:

  Plot an MS2 spectrum using the ID from the peak list table. Includes
  comparison to database spectrum used to make the identification.

- `process_all(comp_names = NULL, alsoDeleteFP = TRUE)`:

  Process all currently unprocessed files in the report object.
  Previously recorded false positives (without specified sample) are
  deleted by default.

- `reIntegrate(indices)`:

  Reintegrate peak areas for all found compounds based only on m/z and
  rt. No MS2 comparison will be performed. The m/z and rt are taken as
  an average of those found in the peak list.

- `reduceSizeEic()`:

  Will average out data points outside of the peak to reduce size of the
  EIC table

- `remRawFiles(indices)`:

  Delete files based on their indices

- `reprocess(indices = NULL, comp_names = NULL, alsoDeleteFP = TRUE)`:

  Reprocess files or specific compounds, will delete all previous
  information. As with process_all, previous FP will be deleted by
  default

- `saveSettings(path = getwd())`:

  Save settings to a json file. Name of file is generated from current
  time.

- `view()`:

  View results using shiny.
