
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ntsworkflow

The goal of ntsworkflow is to process non-target LC/GC-HRMS data and
carry out simple data analysis. It is designed as a basic peak-picking
and annotation algorithm which feeds into more specialized workflows for
different tasks.

If you use the package or any derivative thereof, please cite the work.
To cite this work please use the following citation or run
`citation("ntsworkflow")`.

> Kevin S. Jewell, Christian Dietrich, Toni Köppe, Franziska Thron, Arne
> Wick and Thomas A. Ternes (2025). ntsworkflow: A Non-Target Screening
> Data Evaluation Tool. R package version 0.2.9.

### Peer-reviewed articles reporting the development of ntsworkflow

- K. S. Jewell, U. Kunkel, B. Ehlig, F. Thron, M. Schlüsener, C.
  Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in Mass
  Spectrometry*, **2019**, 34, e8541.
- T. Köppe, K. S. Jewell, C. Dietrich, A. Wick, and T. A. Ternes, *Water
  Research*, **2020**, 178, 115703.
- C. Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in Mass
  Spectrometry*, **2021**, e9206.

## Installation

R 4.5.1 and Java JRE and JDK (v1.8) (for `rcdk`) are needed. Compilation
of C++ code requires GCC v13.

#### Install the Bioconductor package `xcms`

It is recommended to install `xcms` first. Follow the instructions on
[bioconductor - xcms](https://doi.org/doi:10.18129/B9.bioc.xcms).

#### Install ntsworkflow

``` r
remotes::install_github("bafg-bund/ntsworkflow")
```

#### Install Genform (optional)

The Windows .exe can be downloaded directly from Sourceforge. It must be
put on the system path.

For Linux, the program needs to be compiled from source. Download the
source files and in the folder with the .cpp files run the following:
`g++ main.cpp ms*.cpp -o genform`. This will compile the program and
create the `genform` executable, which must be moved to the `~/bin`
folder.

#### Installation on Ubuntu

We recommend using the PPA `r2u` and installing dependencies via `apt`.
The bash command for installing dependencies is:

    sudo apt install r-cran-tidyverse r-cran-rsqlite r-cran-shiny r-cran-shinyfiles r-cran-shinybs r-cran-dt r-cran-tcltk2 r-cran-foreach r-cran-doparallel r-cran-iterators r-cran-rcpp r-cran-rcpparmadillo r-cran-yaml r-cran-jsonlite r-cran-zoo r-cran-rcdk r-cran-future r-cran-furrr r-cran-paralleldist r-bioc-xcms

## Help files

After loading the package with `library(ntsworkflow)`, You can find
documentation by typing `browseVignettes("ntsworkflow")`.

## License

Copyright 2025 Bundesanstalt für Gewässerkunde (Federal Institute of
Hydrology)

ntsworkflow is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

ntsworkflow is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with ntsworkflow. If not, see <https://www.gnu.org/licenses/>.
