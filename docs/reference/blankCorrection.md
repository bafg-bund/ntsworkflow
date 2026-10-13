# Remove features also found in the blanks

Which samples are blanks is given in the sample table. The features
found both in samples and blanks are removed from the alignment table,
unless the intensity in the samples is higher than the blank by the
intensityFactor.

## Usage

``` r
blankCorrection(
  alignment,
  samplels,
  intensityFactor = 10,
  deleteGrouped = TRUE
)
```

## Arguments

- alignment:

  Alignment table (matrix)

- samplels:

  sample table

- intensityFactor:

  Intensity difference between samples and blanks beyond which features
  are not deleted

- deleteGrouped:

  If true (default) then when a blank feature is found, it will be
  removed along with the whole component group (isotopologues, adducts,
  in-source fragments)

## Value

An alignment table with rows removed

## Details

If more than one blank sample exists then the highest intensity will be
used for the comparison.
