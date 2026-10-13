# Peak-picking for an m/z range

Function will pick chromatographic peaks in the raw data file by binning
and analyzing the extracted ion chromatograms.

## Usage

``` r
pickPeaksMzRange(
  daten,
  mz_min,
  mz_max,
  mz_step,
  rt_min,
  rt_max,
  sn,
  int_threshold,
  peak_NoiseScans,
  precursormzTol,
  peakwidth_min,
  peakwidth_max,
  maxPeaksPerSignal
)
```

## Arguments

- daten:

  Measurement data of class \`xcms::xcmsRAW\`

- mz_min:

  m/z range to peak-picking (lower) (Da)

- mz_max:

  m/z range to peak-picking (upper) (Da)

- mz_step:

  m/z width for extracted ion chromatogram (Da)

- rt_min:

  Retention time range minimum in which to look for peaks (s)

- rt_max:

  Retention time range maximum in which to look for peaks (s)

- sn:

  Minimum signal-to-noise ratio (apex peak height over noise spread
  before and after peak)

- int_threshold:

  Minimum peak intensity (at peak apex)

- peak_NoiseScans:

  Number of scans before and after peak to measure noise

- precursormzTol:

  m/z tolerance for linking MS2 fragment spectra by the precursor m/z
  (ppm)

- peakwidth_min:

  Minimum peak width given (s)

- peakwidth_max:

  Maximum peak width given (s)

- maxPeaksPerSignal:

  Maximum number of sub-peaks within a peak (direction changes) within a
  peak

## Value

\`data.frame\` of the peak inventory (peak table)
