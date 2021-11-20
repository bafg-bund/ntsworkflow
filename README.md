
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
> A. Ternes (2019). ntsworkflow: A Non-Target Screening Data Evaluation
> Tool. R package version 0.1.6.

## Installation

### Step 1: Install dependencies

To use the package, you must first install some dependencies. R (v3.6.0)
and RStudio (v1.1) are needed firstly (www.r-project.org,
www.rstudio.com). Java (v1.8) will also be needed
(www.java.com/download). Once you have these installed, please execute
the following code in the R console.

**Important:**

  - Do not allow updating of packages, if R asks to update packages,
    answer with no.
  - If R asks “Do you want to install from sources packages which need
    compilation?” answer with no.
  - This should work off the bat but make sure the right package
    versions have been installed. I try to keep this updated with the
    newest package versions for R 3.6.0.

<!-- end list -->

``` r
 install.packages(
    c(
      "tidyverse",  # 1.3.0
      
      # shiny and co:
      "shiny",  # 1.4.0.2
      "shinyBS",  # 0.61
      "shinyFiles",  # 0.8
      "tcltk2",  # 1.2-11
      "DT",  # 0.13
      
      # Database
      "RSQLite",  # 2.2.0
      "DBI",  # 1.1.0
      
      # base extensions
      "data.table", # 1.12.8
      "foreach",  # 1.4.8
      "iterators", # 1.0.12
      "doParallel", # 1.0.15
      "Rcpp", # 1.0.4.6
      "RcppArmadillo",  # 0.9.900.1.0
      
      # chemistry
      "rcdk",  # 3.5.0
     
      # clustering
      "dtw",  # 1.21-3
      "parallelDist", # 0.2.4
      
      "htmltools" # 0.5.0
      )
    )
```

#### Install Bioconductor packages

The bioconductor version for R 3.6.0 is Bioconductor 3.9.

``` r
install.packages("BiocManager")
BiocManager::install(version = "3.9", lib = "~/R/myBiocLib")
BiocManager::install("MetCirc")  # version 1.14.0
BiocManager::install("xcms")  # version 3.6.2
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("xcms", "MetCirc"), suppressUpdates = TRUE) # xcms version 3.2.0, MetCirc version 1.10.0
```

#### Install packages which need special versions

For Windows users: This works only if you have R development tools
installed. If you do not, you can install the `.zip` files of the
additional packages if you have them using the same proceedure used in
step 2 below.

``` r
devtools::install_version("amap", "0.8-16")
```

### Step 2: Install the ntsworkflow package

Install the ntsworkflow package by using RStudio menu Tools –\> Install
Packages… Under “Install from:” select “Package Archive File” and choose
the `.zip` (Windows) or `.tar.gz` (Linux) file containing the package.

## Help files

After loading the package with `library(ntsworkflow)`, You can find
documentation by typing `browseVignettes("ntsworkflow")`.

## Known bugs

When picking multiple samples, if the first sample is marked as a blank
the app will hang with an out-of-bounds index error. To avoid this do
not mark samples as blanks until after the peak-picking is complete.
