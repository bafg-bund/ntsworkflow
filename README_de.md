
<!-- README.md is generated from README.Rmd. Please edit that file -->
ntsworkflow
===========

Das Ziel von "ntsworkflow" ist es non-target LC/GC-HRMS Daten zu verarbeiten und eine einfache Datenanalyse durchzuführen. Das Paket ist als ein grundlegender Peak-Picking- und Erkennungs-Algorithmus erstellt, welcher zu weiteren differenzierten Workflows für verschiedene Aufgaben führt.

Nutzungsrechte
--------------

Das "ntsworkflow" Paket ist unter der GNU GPL-3.0 Lizenz verfügbar. Bitte zitieren Sie das Werk, wenn das Paket oder Varianten davon genutzt werden. Um das Werk zu zitieren, nutzen Sie bitte das folgende Zitat oder führen Sie `citation("ntsworkflow")` aus.

> Christian Dietrich, Kevin S. Jewell, Ute Thorenz, Toni Köppe, Arne Wick and Thomas A. Ternes (2019). ntsworkflow: A Non-Target Screening Data Evaluation Tool. R package version 0.1.4.

Installation
------------

### Schritt 1: Grundvoraussetzungen

Um das Paket nutzen zu können, ist es notwendig einige Grundvorrausetzungen zu schaffen. Zunächst werden R (Version 3.5.x) und RStudio benötigt (www.r-project.org, www.rstudio.com). Weiterhin wird Java gebraucht (www.java.com/download). Sobald Sie diese Programme installiert haben, führen Sie bitte den folgenden Code in der R Konsole aus. Lassen Sie NICHT zu, dass die Pakete aktualisiert werden. Wenn R nachfragt, ob die Pakete aktualisiert werden sollen, antworten Sie mit "Nein".

``` r
 install.packages(
    c(
      "flux",  # 0.3-0
      "ggplot2",  # 3.1.0
      "plyr",  # 1.8.4
      "dplyr",  # 0.7.7 
      "dbplyr",  # 1.2.2
      "tidyr",  # 0.8.2
      "RSQLite",  # 2.1.1
      "rcdk",  # 3.4.7.1
      "webchem",
      "data.table",
      "doParallel",
      "shiny",
      "RSQLite",
      "DBI",
      "rcdk",
      "webchem",
      "OrgMassSpecR",
      "DT",
      "roxygen2",
      "Rcpp",
      "RcppArmadillo",
      "shinyBS",
      "pryr",
      "foreach",
      "lubridate",
      "iterators",
      "tcltk2",
      "knitr",
      "rmarkdown",
      "dtw",
      "shinyFiles",
      "parallelDist",
      "shinyjs"
      )
    )
```

Installieren von bioconductor packages

``` r
source("https://bioconductor.org/biocLite.R")
biocLite(c("xcms", "MetCirc"), suppressUpdates = TRUE)

# Für R 3.6.x

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MetCirc")

BiocManager::install("xcms")
```

Installieren von Paketen, die besondere Versionen benötigen

(Für MS-Windows: Dieser Schritt funktioniert nur mit Hilfe von "R Development Tools". Wenn diese nicht verfügbar sind, kann man alternativ die .zip Dateien der zusätzliche Pakete installieren anhand der Anleitung in Schritt 2, unten)

``` r
devtools::install_version("amap", "0.8-16")
```

### Step 2: Installieren des Pakets

Installieren Sie das "ntsworklflow" Paket, indem Sie den RStudio Menüpunkt Tools --&gt; Install Packages... anwählen. Unter "Install from:" wählen Sie "Package Archieve File" aus und dann nutzen Sie die `.zip` (Windows) oder `.tar.gz` (Linux) Datei um das Paket auszuführen.

Hilfe Dateien
-------------

Nachdem das Paket mit `library(ntsworkflow)` gestartet wurde, lässt sich die Dokumentation durch Eingabe von `browseVignettes("ntsworkflow")` finden.
