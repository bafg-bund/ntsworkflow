# Get a recommendation for the m/z-step parameter (binning width)

Calculates the optimum mz step based on mass differences. Takes a sample
of 100 spectra to make the calculation

## Usage

``` r
optimumMzStep(daten, probs)
```

## Arguments

- daten:

  raw data as list of xcmsRAW objects

- probs:

  quantile to use (see quantile function)

## Value

Recommended m/z width (Da)
