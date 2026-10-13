# Only retain peaks found in all replicates

This function will check, for a given number of replicate injections,
that a feature was found in all replicates. Features will be set to zero
if they are not found in replicates. Aligned features not found in any
samples (everything zero) will be removed.

## Usage

``` r
keepReps(alignment, samples, reps, least)
```

## Arguments

- alignment:

  alignment table

- samples:

  numeric vector sample IDs to merge, length must be divisible by reps

- reps:

  number of replicates

- least:

  peak found in at least this many replicates, must be \<= reps

## Value

An alignment table with the same number of columns but a reduced number
of rows.
