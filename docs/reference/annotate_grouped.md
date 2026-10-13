# Annotate alignment table with compound names by spectral library screening

Using the SQLite spectral library (CSL), this function will pass over
the alignment table, open the measurement files stored in datenList and
compare m/z, retention time and MS² spectra found in the data and the
library.

## Usage

``` r
annotate_grouped(
  sampleListLocal,
  peakListList,
  alignmentTable,
  db_path,
  threshold_score = 500,
  mztolu = 0.005,
  rttol = 1,
  polarity = "pos",
  CE = c(30, 40),
  CES = c(0, 15),
  instrument = c("LC-ESI-QTOF TripleTOF 5600 SCIEX", "LC-ESI-QTOF TripleTOF 6600 SCIEX",
    "LC-ESI-Orbitrap QExactive Thermo", "LC-ESI-QTOF Agilent 6500 Series Q-TOF",
    "LC-ESI-QTOF TripleTOF X500R SCIEX", "LC-ESI-QTOF ZenoTOF 7600 SCIEX",
    "LC-ESI-Orbitrap Exploris 240 Thermo"),
  chrom_method = "bfg_nts_rp1",
  ndp_m = 2,
  ndp_n = 1,
  mztolu_ms2 = 0.015,
  rtoffset = 0,
  intCutData = 0,
  numcores = 1,
  datenListLocal,
  expGroups = "bfg"
)
```

## Arguments

- peakListList:

  list of data.frames of the peak-lists

- alignmentTable:

  Knonw as "grouped" in the app

- db_path:

  Path to the spectral library (CSL in SQLite format)

- threshold_score:

  dot product threshold score

- mztolu:

  m/z tolerance in Da

- rttol:

  RT tolerance in min.

- polarity:

  Polarity of the measurement

- CE:

  Collision energy

- CES:

  Collision energy spread

- instrument:

  Instruments allowed, must match the instruments listed in the library

- chrom_method:

  Chromatography name to choose the correct retention times from the
  library

- ndp_m:

  Peak intensity weighting factor for dot-product

- ndp_n:

  m/z weighting factor for dot-product

- mztolu_ms2:

  m/z window for dot-product

- rtoffset:

  Offset the database to make it match your values (in min)

- intCutData:

  Cut off intensity under which fragments are ignored (in data spectra)

- numcores:

  Number of cores for parallelization (currently not used)

- sampleList:

  data.frame of filenames and paths

## Value

\`data.frame\` with the annotations for the alignment table referenced
by the column "alignementID" (compound identifiers and results of the
spectral comparison)

## Details

This function operates either in conjunction with the non-target app but
can also work in stand-alone scripts, see examples in ntsportal. If data
files are not in memory (xcmsRaw objects), they will be temporary
loaded. \`datenListLocal\` is a shallow copy. Therefore, changes to the
underlying xcmsRaw objects affect the calling environment.

\## Spectral library The spectral library is an SQLite file (see Jewell
et al. 2019).

\## References Jewell, K. S., Kunkel, U., Ehlig, B., Thron, F.,
Schlüsener, M., Dietrich, C., . . . Ternes, T. A. (2019). Comparing
mass, retention time and MS2 spectra as criteria for the automated
screening of small molecules in aqueous environmental samples analyzed
by LC-QToF-MS/MS. Rapid Communications in Mass Spectrometry, 34, e8541.
doi:10.1002/rcm.8541
