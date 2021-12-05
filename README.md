# ntsworkflow 

=======

<!-- README.md is generated from README.Rmd. Please edit that file -->

# ntsworkflow

The goal of ntsworkflow is to process non-target LC/GC-HRMS data and
carry out simple data analysis. It is designed as a basic peak-picking
and annotation algorithm which feeds into more specialized workflows for
different tasks.

## Usage rights

The ntsworkflow package is shared under the GNU GPL-3.0 license. If you
use the package or any derivative thereof, please cite the work. To cite
this work please use the following citation or run
`citation("ntsworkflow")`.

> Christian Dietrich, Kevin S. Jewell, Toni Köppe, Arne Wick and Thomas
> A. Ternes (2021). ntsworkflow: A Non-Target Screening Data Evaluation
> Tool. R package version 0.2.2.

### Peer-reviewed articles reporting the development of ntsworkflow

-   K. S. Jewell, U. Kunkel, B. Ehlig, F. Thron, M. Schlüsener, C.
    Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in Mass
    Spectrometry*, **2019**, 34, e8541.
-   T. Köppe, K. S. Jewell, C. Dietrich, A. Wick, and T. A. Ternes,
    *Water Research*, **2020**, 178, 115703.
-   C. Dietrich, A. Wick, and T. A. Ternes, *Rapid Communications in
    Mass Spectrometry*, **2021**, e9206.

## Installation

ntsworkflow is designed for use on servers and is only available for
Linux-based systems. It is usually run on centOS or Ubuntu.

### Step 1: Install dependencies

To use the package, you must first install some dependencies. R (v4.0.4)
and RStudio (v1.1) are needed firstly (www.r-project.org,
www.rstudio.com). Java JRE and JDK (v1.8) will also be needed
(www.java.com/download). Once you have these installed, please execute
the following code in the R console.

#### Installation on Ubuntu

Using `apt`, the packages r-base, default-jre and default-jdk are
needed. Installing Rstudio server for WSL2 was done through the website
instructions (deb package).

#### Install packages

This may take some time on Linux OSs since everything needs to be
recompiled.

**Important:**

-   Do not allow updating of packages, if R asks to update packages,
    answer with no.
-   If R asks “Do you want to install from sources packages which need
    compilation?” answer with no.
-   This should work off the bat but make sure the right package
    versions have been installed. I try to keep this updated with the
    newest package versions for R 4.0.4.

#### installing rcdk

This package is particularly tricky since it depends on java. Make sure
you have a JRE and a JDK available. If it does not work try installing
rJava first.

``` r
 install.packages(
    c(
      "tidyverse",  # 1.3.0
      
      # shiny and co:
      "shiny",  # 1.5.0
      "shinyBS",  # 0.61
      "shinyFiles",  # 0.9
      "tcltk2",  # 1.2-11
      "DT",  # 0.16
      
      # Database
      "RSQLite",  # 2.2.1
      "DBI",  # 1.1.0
      
      # base extensions
      "foreach",  # 1.4.8
      "iterators", # 1.0.12
      "doParallel", # 1.0.15
      "Rcpp", # 1.0.5
      "RcppArmadillo",  # 0.10.1.0.0
      
      # chemistry
      "rcdk",  # 3.5.0
     
      # clustering
      "dtw",  # 1.21-3
      "parallelDist", # 0.2.4
      
      # utils
      "htmltools" # 0.5.0
      "yaml"  # 2.2.1
      )
    )
```

#### Install Bioconductor packages

I like to keep these packages separate from the cran packages. So I keep
them in a separate directory.

First update the .Rprofile and move the bioclib directory to the front.
Restart R.

##### Installation of xcms

xcms requires mzR which in turn requires ncdf4 which needs the nc-config
script to be installed. Prior to installing xcms you therefore need to
install ncdf4. For some reason, using `BiocManager::install` or
`install.packages` does not work, you must use `apt`
(`sudo apt install r-cran-ncdf4`). This installs ncdf4 into
`/usr/lib/R/site-library`.

Installation of xcms takes a long time due to all the dependencies. If
you can use a pre-installed library folder.

The bioconductor version for R 4.0.4 is Bioconductor 3.12.

``` r
install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("xcms")  
```

#### Install Genform

The Windows .exe can be downloaded directly from Sourceforge. For Linux,
the program needs to be compiled from source. Download the source files
and in the folder with the cpp files run the following:
`g++ main.cpp ms*.cpp -o genform`

This will compile the program and create the `genform` executable, which
must be moved to the `~/bin` folder.

### Step 2: Install the ntsworkflow package

Install the ntsworkflow package by using RStudio menu Tools –> Install
Packages… Under “Install from:” select “Package Archive File” and choose
the `.zip` (Windows) or `.tar.gz` (Linux) file containing the package.

#### Installing from source

For this you need devtools and on Windows also Rtools.

Installation of gert for devtools requires libgit2-dev on Ubuntu (use
`apt`)

## Help files

After loading the package with `library(ntsworkflow)`, You can find
documentation by typing `browseVignettes("ntsworkflow")`.

## Known bugs

When picking multiple samples, if the first sample is marked as a blank
the app will hang with an out-of-bounds index error. To avoid this do
not mark samples as blanks until after the peak-picking is complete.
