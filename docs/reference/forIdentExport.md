# Export data to For-Ident.

For batch processing in For-Ident, txt file is created for import on
For-Ident platform.

## Usage

``` r
forIdentExport(
  align_matrix,
  sampleList,
  rawDataList,
  intThreshRel = 0.01,
  report_nm
)
```

## Arguments

- align_matrix:

  Alignment `matrix` from Christian's App

- sampleList:

  sampleList `matrix` from App

- rawDataList:

  datenList `list` from App

- intThreshRel:

  Intensity threshold for including peaks in MS2 spectra

- report_nm:

  Name of output file

## Value

A txt file with each found mass, retention time and MS2 spectrum (if
available).

## Details

MS2 spectra are included for Met-frag or Massbank search. The spectrum
for each row is taken from the sample with the highest intensity in this
row. The function takes approx. 1 min to process 100 rows, so limit the
size of the alignment matrix passed to the app. If no MS2 is available,
row is skipped.
