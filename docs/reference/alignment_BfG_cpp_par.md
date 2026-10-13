# Align peaks with parallel processing

To align very large datasets it is necessary to split the alignment
process into parallel jobs. This function slips the peak lists and runs
the alignment_BfG_cpp on each chunk in parallel. The results are
combined. Uses forking therefore not available on Windows.

## Usage

``` r
alignment_BfG_cpp_par(
  peaklists,
  ppm_dev,
  DeltaRT,
  mz_dev_unit,
  mDa_split = 100,
  numSplits = 16
)
```

## Arguments

- peaklists:

  list of peaklists

- ppm_dev:

  mz tolerance

- DeltaRT:

  rt tolerance in s

- mz_dev_unit:

  units for mz tolerance (ppm or mDa)

- mDa_split:

  mz gap between chunks

- numSplits:

  Number parallel chunks to create

## Value

matrix of aligned peaks, one aligned feature per row
