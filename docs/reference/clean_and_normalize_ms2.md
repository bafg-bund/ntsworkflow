# Filter and normalize MS2 spectra by removing precursor peaks and noise.

This function filters out peaks close to the precursor m/z and
low-intensity noise, normalizes intensities to a maximum of 1, and
retains relevant metadata. Returns NULL if no valid peaks remain.

## Usage

``` r
clean_and_normalize_ms2(spec)
```
