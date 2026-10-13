# Peakpeaking algorithm using C++

This function wraps the c++ peakpicking algorithm. The function picks
peaks in a specific mz range (defined by i and mz_step). Used by the
function FindPeaks_BfG.

## Usage

``` r
pickPeaksOneEic(
  i,
  rawData,
  mz_step,
  rt_min_scan,
  rt_max_scan,
  sn,
  int_threshold,
  NoiseScans,
  peakwidth_min,
  peakwidth_max,
  precursormzTol,
  maxPeaksPerSignal
)
```

## Arguments

- i:

  Lower m/z for extracted ion chromatogram in Da

- rawData:

  Measurement data of class \`xcms::xcmsRAW\`

- mz_step:

  m/z width for extracted ion chromatogram in Da

- rt_min_scan:

  Retention time range minimum in which to look for peaks, in scans

- rt_max_scan:

  Retention time range maximum in which to look for peaks, in scans

- sn:

  Minimum signal-to-noise ratio (apex peak height over noise spread
  before and after peak)

- int_threshold:

  Minimum peak intensity (at peak apex)

- NoiseScans:

  Number of scans before and after peak to measure noise

- peakwidth_min:

  Minimum peak width given in seconds

- peakwidth_max:

  Maximum peak width given in seconds

- precursormzTol:

  m/z tolerance for linking MS2 fragment spectra by the precursor m/z,
  in ppm

- maxPeaksPerSignal:

  Maximum number of sub-peaks within a peak (direction changes) within a
  peak

## Value

'matrix' of peaks detected in the extracted ion chromatogram, each peak
is one row with the columns representing various parameters for the
peak.
