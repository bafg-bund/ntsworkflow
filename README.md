
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
> Data Evaluation Tool. R package version 0.2.7.

### Peer-reviewed articles reporting the development of ntsworkflow

- K. S. Jewell, U. Kunkel, B. Ehlig, F. Thron, M. Schlüsener, C.
  Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in Mass
  Spectrometry*, **2019**, 34, e8541.
- T. Köppe, K. S. Jewell, C. Dietrich, A. Wick, and T. A. Ternes, *Water
  Research*, **2020**, 178, 115703.
- C. Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in Mass
  Spectrometry*, **2021**, e9206.

## Installation

R 4.5.0 and Java JRE and JDK (v1.8) (for rcdk package) are needed.
Compilation of C++ code requires GCC v13.

#### Install the Bioconductor package xcms

Follow the instructions on [bioconductor -
xcms](https://doi.org/doi:10.18129/B9.bioc.xcms).

#### Install Genform (optional)

The Windows .exe can be downloaded directly from Sourceforge. It must be
put on the system path.

For Linux, the program needs to be compiled from source. Download the
source files and in the folder with the .cpp files run the following:
`g++ main.cpp ms*.cpp -o genform`. This will compile the program and
create the `genform` executable, which must be moved to the `~/bin`
folder.

#### Install ntsworkflow

``` r
remotes::install_github("bafg-bund/ntsworkflow")
```

#### Installation on Ubuntu

Using `apt`, the packages `r-base`, `default-jre` and `default-jdk` are
needed. Installing Rstudio server for WSL2 was done through the website
instructions (deb package).

xcms requires *mzR* which in turn requires *ncdf4* which needs the
`nc-config` script to be installed. Prior to installing xcms you
therefore need to install *ncdf4*. Use `sudo apt install r-cran-ncdf4`.
This installs *ncdf4* into `/usr/lib/R/site-library`.

Installation of `gert` for `devtools` requires `libgit2-dev` on Ubuntu.

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

## Change log

### 0.2.7

In `Report$integRes`, the column `samp` contains only the sample name
(previously sample path)
