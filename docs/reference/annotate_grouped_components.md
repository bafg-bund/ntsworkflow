# Annotate alignment table with adduct and isotopologue information

This function is for the annotation section of the app. If no db is
chosen, it will annotate the alignmentTable with peaklist component info
(Cl and Br and some adducts)

## Usage

``` r
annotate_grouped_components(
  sampleListLocal,
  peakListList,
  alignmentTable,
  numcores
)
```

## Arguments

- sampleListLocal:

  (\`data.frame\` of processed samples)

- peakListList:

  (from peakpicking)

- alignmentTable:

  (from alignment)

- numcores:

  (number of CPU cores to be used)

## Value

alignmentTable with annotated isotopologues based on componentization
