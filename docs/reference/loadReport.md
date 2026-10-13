# Load a previously saved Report

A dialog box is used or the path can be given as an argument.

## Usage

``` r
loadReport(dialog = TRUE, path = NULL)
```

## Arguments

- dialog:

  logical should a dialog window open to select the file (only works in
  rstudio)

- path:

  if dialog is false a path to the file. Must be given if rstudioapi is
  not available

## Value

an \`ntsworkflow::Report\` object

## Details

rstudioapi is used for the dialog. readRDS can be used directly to open
saved Report objects as well.

## See also

readRDS
