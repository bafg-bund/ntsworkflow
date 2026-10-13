# Get mean intensities for replicate injections

This function will take an alignment table and average the intensities
of replicate injections. The first sample in each replicate will remain
in the table with average int. Zeros will be ignored (unless all
intensities are zero).

## Usage

``` r
averageReplicates(alignment, samples, reps)
```

## Arguments

- alignment:

  alignment table

- samples:

  numeric sample IDs to merge, length must be divisable by repNum

- reps:

  number of replicates

## Value

An alignment table with a reduced number of samples but the same number
of rows.
