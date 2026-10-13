# Fix Whitespace Error in GenForm Loss Annotation

This takes a filename with e.g. the output from
[`annotateCleanMSMS`](http://r.bafg.de/reference/annotateCleanMSMS.md)
and removes a space in the `"no loss"` annotation to allow for easier
handling.

## Usage

``` r
fixNoLossSpace(anno_MS_file, replace=TRUE)
```

## Arguments

- anno_MS_file:

  File name for annotated cleaned MSMS text file output from
  [`annotateCleanMSMS`](http://r.bafg.de/reference/annotateCleanMSMS.md)

- replace:

  If `TRUE`, replace the input file. If `FALSE`, saved under the new
  name `"*_ed.txt"`.

## Value

File name of the "fixed" file.

## See also

[`annotateCleanMSMS`](http://r.bafg.de/reference/annotateCleanMSMS.md)

## Author

Emma Schymanski (`emma.schymanski@uni.lu`)
