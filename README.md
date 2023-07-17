
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ntsworkflow

The goal of ntsworkflow is to process non-target LC/GC-HRMS data and
carry out simple data analysis. It is designed as a basic peak-picking
and annotation algorithm which feeds into more specialized workflows for
different tasks.

## License

Copyright 2016-2023 Bundesanstalt für Gewässerkunde (Federal Institute
of Hydrology)

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

If you use the package or any derivative thereof, please cite the work.
To cite this work please use the following citation or run
`citation("ntsworkflow")`.

> Kevin S. Jewell, Christian Dietrich, Toni Köppe, Franziska Thron, Arne
> Wick and Thomas A. Ternes (2023). ntsworkflow: A Non-Target Screening
> Data Evaluation Tool. R package version 0.2.2.

### Peer-reviewed articles reporting the development of ntsworkflow

- K. S. Jewell, U. Kunkel, B. Ehlig, F. Thron, M. Schlüsener, C.
  Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in Mass
  Spectrometry*, **2019**, 34, e8541.
- T. Köppe, K. S. Jewell, C. Dietrich, A. Wick, and T. A. Ternes, *Water
  Research*, **2020**, 178, 115703.
- C. Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in Mass
  Spectrometry*, **2021**, e9206.

## Installation

ntsworkflow is designed for use on servers and currently only works on
Linux-based systems. It is usually run on CentOS or Ubuntu.

### Step 1: Install dependencies

To use the package, you must first install some dependencies. R (v4.3.0)
and Java JRE and JDK (v1.8) (for rcdk package).

#### Installation on Ubuntu

Using `apt`, the packages r-base, default-jre and default-jdk are
needed. Installing Rstudio server for WSL2 was done through the website
instructions (deb package).

#### Install CRAN packages

``` r
 install.packages(
    c(
      "tidyverse", 
      
      # shiny and co:
      "shiny",  
      "shinyBS",  
      "shinyFiles", 
      "tcltk2",  
      "DT",  
      
      # Database
      "RSQLite",  
      "DBI",  
      
      # base extensions
      "foreach",  
      "iterators", 
      "doParallel", 
      "Rcpp", 
      "RcppArmadillo",  
      
      # chemistry
      "rcdk",  
     
      # clustering
      "dtw", 
      "parallelDist",
      
      # utils
      "htmltools" 
      "yaml"  
      )
    )
```

#### Install Bioconductor package xcms

Installation of xcms takes a long time due to all the dependencies. If
you can, use a pre-installed library folder.

The bioconductor version for R 4.3.0 is Bioconductor 3.17.

``` r
install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install("xcms")  
```

xcms requires mzR which in turn requires ncdf4 which needs the nc-config
script to be installed. Prior to installing xcms you therefore need to
install ncdf4. If, for some reason, using `BiocManager::install` or
`install.packages` does not work, you can use `apt`
(`sudo apt install r-cran-ncdf4`) on Ubuntu. This installs ncdf4 into
`/usr/lib/R/site-library`.

#### Install Genform

The Windows .exe can be downloaded directly from Sourceforge. For Linux,
the program needs to be compiled from source. Download the source files
and in the folder with the .cpp files run the following:
`g++ main.cpp ms*.cpp -o genform`

This will compile the program and create the `genform` executable, which
must be moved to the `~/bin` folder.

### Step 2: Install the ntsworkflow package

Install the ntsworkflow package by using RStudio menu Tools –\> Install
Packages… Under “Install from:” select “Package Archive File” and choose
the `.tar.gz` file containing the package.

#### Installing from source

For this you need devtools.

Installation of gert for devtools requires libgit2-dev on Ubuntu.

## Help files

After loading the package with `library(ntsworkflow)`, You can find
documentation by typing `browseVignettes("ntsworkflow")`.

## Known bugs

When peak-picking multiple samples, if the first sample is marked as a
blank the app will hang with an out-of-bounds index error. To avoid this
do not mark samples as blanks until after the peak-picking is complete.
