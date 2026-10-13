# Annotate alignment table based only on m/z and RT comparison

Enables also annotation based only on m/z. Uses a csv file substance
library.

## Usage

``` r
annotate_grouped_mz_rt(alig, compLibPath, mztol, rttol)
```

## Arguments

- alig:

  alignment table from shiny app

- compLibPath:

  path to csv file

- mztol:

  tolerance in u

- rttol:

  tolerance in minutes

## Value

Annotation table (\`data.frame\`) which gives the annotations of the
alignment table linked by \`alignmentID\`
