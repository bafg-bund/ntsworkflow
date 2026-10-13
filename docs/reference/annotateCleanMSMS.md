# Annotate Cleaned MSMS

This takes the output from
[`RunGenForm`](http://r.bafg.de/reference/RunGenForm.md) and creates a
peaklist with the annotations for further use.

## Usage

``` r
annotateCleanMSMS(cleanMSMS, GenForm_an_out, GenForm_loss_out, lowestPpm=TRUE, round=TRUE,
removeAnnoSpace=FALSE)
```

## Arguments

- cleanMSMS:

  File name for Cleaned MSMS text file output from
  [`RunGenForm`](http://r.bafg.de/reference/RunGenForm.md)

- GenForm_an_out:

  File name for annotated MSMS text file output from
  [`RunGenForm`](http://r.bafg.de/reference/RunGenForm.md)

- GenForm_loss_out:

  File name for annotated "loss" MSMS text file output from
  [`RunGenForm`](http://r.bafg.de/reference/RunGenForm.md)

- lowestPpm:

  If `TRUE`, choses the lowest ppm subformula for assignment of formula
  to fragments. If `FALSE`, takes the first formula.

- round:

  If `TRUE`, rounds the output to 4 (mz) and 1 (Int) decimal places for
  aesthetics. If `FALSE`, no rounding is done.

- removeAnnoSpace:

  If `TRUE`, removes the space in the `"no loss"` annotation.

## Value

Returns annotated peaks for further use

## Author

Emma Schymanski (`emma.schymanski@uni.lu`)
