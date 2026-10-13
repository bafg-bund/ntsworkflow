# Keep only features which are the group leaders

This will go through the alignment table and check every feature if it
is a group leader. Features which are marked as group leader in ANY
sample will kept. In other words, if the componentization is
inconsistent and a feature is sometimes marked as groupleader and
sometimes not, the feature will be kept.

## Usage

``` r
onlyGroupLeaders(alignment, pll, samples)
```

## Arguments

- alignment:

  alignment table

- pll:

  collection of peaklists as a \`list\`

- samples:

  sample IDs of the samples which will be looked at, all other samples
  will be ignored

## Value

An alignment table with reduced number of rows
