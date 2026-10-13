# Group features into components

Componentization/ Grouping of adducts, isotopologues, in-source
fragments, etc. Either by: fixed thresholds for RT, RT_FWHM_left and
RT_FWHM_right or automatically adjusted thresholds based on the 13C
isotopologue.

## Usage

``` r
componentization_BfG(
  Liste,
  daten,
  ppm = 10,
  Grenzwert_RT = 10,
  Grenzwert_FWHM_left = 10,
  Grenzwert_FWHM_right = 10,
  Summe_all = 10,
  adjust_tolerance = TRUE
)
```

## Arguments

- Liste:

  peaklist from the peak-picking step

- daten:

  (datenList)

- ppm:

  (mass deviation in ppm)

- Grenzwert_RT:

  (Threshold value for the retention time)

- Grenzwert_FWHM_left:

  (Limit value for the retention time of the left FWHM)

- Grenzwert_FWHM_right:

  (Limit value for the retention time of the right FWHM)

- Summe_all:

  (Threshold value for the sum of the three retention times (apex, FWHM
  left, FWHM right)

- adjust_tolerance:

  (adjust thresholds with 13C isotopologue)

## Value

peaklist with added cols groupleader and group for each component
