# Library screening without prior peak-picking

A direct library screening approach using measurement files and the CSL
spectral library in SQLite format.

## Usage

``` r
ms2_search(
  data_path,
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
)
```

## Arguments

- data_path:

  Path to rawfiles, can be a vector of file locations or a list of
  xcmsRaw objects

- csl_path:

  Path to the Collective Spectral Library (CSL)

- chromatography:

  Chromatography name to choose the correct retention times from the
  library

- pol:

  Polarity

- ce_s:

  Collision energy

- ces_s:

  Collision energy spread

- instr:

  Instruments allowed, must match the instruments listed in the library

- rtoffset:

  Retention time offset (will be added to library RTs)

- rt_res:

  Resolution of two peaks in chromatography (min)

- rttolm:

  Retention time tolerance in minutes

- mztolu:

  m/z tolerance in Da (MS2 precursor mass)

- mztolu_fine:

  m/z tolerance in Da (second stage)

- mztolu_ms2:

  m/z window for dot-product

- threshold:

  Threshold score for spectral comparison (for dot_product 1000 is max)

- ndp_m:

  Peak intensity weighting factor for dot-product

- ndp_n:

  m/z weighting factor for dot-product

- compounds:

  character vector of compounds which should be processed

## Value

a \`data.frame\` with hits from the CSL

## Details

The function is used in the Report RC to obtain the peak list. No prior
peak picking is needed for this function. The function starts from the
list of masses provided in the spectral library and searches for these
in the measurement data. The measurement data is in mz(X)ML format and
contains DDA MS² spectra. **Multiple peaks for one compound** Will be
sorted by intensity and labelled A-Z \# TODO: Check if this is already
the case

\## References Jewell, K. S., Kunkel, U., Ehlig, B., Thron, F.,
Schlüsener, M., Dietrich, C., . . . Ternes, T. A. (2019). Comparing
mass, retention time and MS2 spectra as criteria for the automated
screening of small molecules in aqueous environmental samples analyzed
by LC-QToF-MS/MS. Rapid Communications in Mass Spectrometry, 34, e8541.
doi:10.1002/rcm.8541

## See also

\[annotate_grouped()\]
