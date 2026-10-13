# Extraction of a MS2-Spectrum for a given ExperimentID from spectral library

Extraction of a MS2-Spectrum for a given ExperimentID from spectral
library

## Usage

``` r
dbGetSpectrum(db, Exp.ID)
```

## Arguments

- db:

  Either a DBI connection object OR a csl_data list as returned by
  get_csl_data() (containing \$exptbl, \$comptbl, \$fragtbl)

- Exp.ID:

  integer, a single experiment ID

## Value

MS2 spectrum as a data.frame (mz, int)
