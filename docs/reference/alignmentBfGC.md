# Align peaks in separate peak lists by m/z and retention time tolerance

This function will align the peaks in multiple peak lists by comparing
the m/z and RT.

## Usage

``` r
alignmentBfGC(peaklistR, mz_dev, DeltaRT, mz_dev_unit)
```

## Arguments

- peaklistR:

  list of peak tables

- mz_dev:

  m/z tolerance

- DeltaRT:

  RT tolerance (s)

- mz_dev_unit:

  Units for the m/z tolerance 1 = ppm, 2 = mDa

## Value

A matrix of class Rcpp::integerMatrix giving the .
