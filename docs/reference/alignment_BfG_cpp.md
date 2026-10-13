# Alignment of peaks from difference samples

This a wrapper for the C++ function alignmentBfGC. It takes a list of
peaklists and produces an alignment table.

## Usage

``` r
alignment_BfG_cpp(peaklists, mz_dev, DeltaRT, mz_dev_unit)
```

## Arguments

- peaklists:

  list of peaklists

- mz_dev:

  m/z tolerance

- DeltaRT:

  rt tolerance in s

- mz_dev_unit:

  Units for m/z tolerance (ppm or mDa)

## Value

Matrix of aligned peaks, one aligned feature per row
