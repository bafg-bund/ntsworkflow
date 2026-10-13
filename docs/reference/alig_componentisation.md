# Second stage componentisation using the alignment table

Based on four criteria: rt, peak shape, intensity correlation and common
mz differences. Will replace group column in the alignment table.
Distance matrices for all four criteria are combined and then DBSCAN
clusters the features.

## Usage

``` r
alig_componentisation(
  altable,
  rttols = 3,
  fracComponMatch = 0.5,
  mztol = 0.005,
  pearsonCorr = 0.5,
  pol = "pos",
  numcores = 6
)
```

## Arguments

- altable:

  alignment table matrix ("grouped" in the app)

- rttols:

  rt tolerance in s

- fracComponMatch:

  fraction of samples with matching peak shape

- mztol:

  mz tolerance

- pearsonCorr:

  minimum pearson's r for comparing intensity trends

- pol:

  polarity, must be either "pos" or "neg"

- numcores:

  number of cores to use for parallel distance matrix computations

## Value

Alignment table with the group column "Gruppe" replaced with the new
values
