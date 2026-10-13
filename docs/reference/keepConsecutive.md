# Keep only rows found in at least X consecutive samples

Keep only rows found in at least X consecutive samples

## Usage

``` r
keepConsecutive(alignment, samples, consecutive)
```

## Arguments

- alignment:

  Alignment Table

- samples:

  Samples to consider (everything else will be ignored), by ID

- consecutive:

  Number of consecutive samples a peak should be present in

## Value

Alignment table with possibly some rows removed
