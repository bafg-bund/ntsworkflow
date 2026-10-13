# Calculate Normalized Dot Product Similarity Score Between MS2 Spectra

Calculate Normalized Dot Product Similarity Score Between MS2 Spectra

## Usage

``` r
calc_ndp(data_spec, csl_spec, ndp_m = 2, ndp_n = 1, mztolu_ms2 = 0.015)
```

## Arguments

- data_spec:

  Spectrum from data

- csl_spec:

  Spectrum from CSL database

- ndp_m:

  Peak intensity weighting factor

- ndp_n:

  m/z weighting factor

- mztolu:

  m/z tolerance for binning peaks (u)

## Value

Similarity score from 0 to 1000, where 1000 indicates a perfect match

## Details

See MetCirc::NDP for more information
