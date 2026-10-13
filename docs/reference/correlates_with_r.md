# Find Pearson's correlation between intensity trends

Wrapper for the cpp function correlates_with. Takes a matrix of
intensities (each row representing one feature) and a row, returns which
rows have correlation coefficients within threshold

## Usage

``` r
correlates_with_r(aligned_intensities, zeile, koeffizient)
```

## Arguments

- aligned_intensities:

  Matrix of intensities

- zeile:

  Reference row of matrix to which all other are compared to

- koeffizient:

  Pearson correlation coefficient (minimum)

## Value

Boolean. True for those rows of the matrix which have a correlation
coefficient higher than the threshold.
