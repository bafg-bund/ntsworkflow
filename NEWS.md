
# ntsworkflow 


# 0.2.11

- When performing background correction within a Report, (`deleteBackground` method), this now also 
removes all rows of the reintegration table (IntegRes) of all the compounds (name, adduct, isotopologue combination) 
which are deleted from the peakList.
- Bugfixes in mergeReport()

## Dependency version updates

- R=4.6.1
- C++=C++17

# 0.2.10

Moved to CSL version 25.7. Older CSL version will no longer work.

- Refactor variables and update fixtures for CSL version 25.7
- Refactor library screening functions
- Add and update integration and unit tests for library screening functions
- Cleanup and update documentation

## Other dependency version updates

- R=4.5.2
- Bioconductor=3.22
- xcms=4.8.0
- Rcpp=1.1.1
- tidyverse=2.0.0
- g++=13.3.1


# 0.2.9

Moved to CSL version 25.5. Older CSL versions will no longer work.

# 0.2.8

Removed dependency on `shinyBS`, as it no longer works with `DT` and is no longer maintained.

# 0.2.7 

In `Report$integRes`, the column `samp` contains only the sample name (previously sample path)
