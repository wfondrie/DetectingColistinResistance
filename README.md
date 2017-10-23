Machine learning enables pathogen identification and detection of antimicrobial resistance from mass spectrometry analysis of membrane glycolipids
================
William E Fondrie,
October 23, 2017

-   [Introduction](#introduction)
-   [Directory structure](#directory-structure)
-   [Code to run the analysis](#code-to-run-the-analysis)
    -   [Prepare workspace](#prepare-workspace)
    -   [Run Analysis Scripts](#run-analysis-scripts)
-   [Session Info](#session-info)

Introduction
------------

This repository contains the code and data required to reproduce the analysis presented in [*Machine learning enables pathogen identification and detection of antimicrobial resistance from mass spectrometry analysis of membrane glycolipids*]().

**WARNING:** When attempting to reproduce this analysis, it should be noted that `Rmd/simulateCompleSpectra.Rmd` requires approximately 50 Gb of memory to execute. Additionally, the entire analysis time may take a number of hours depending on the hardware.

Order of analysis (see the Rmd subdirectory for details):
1. [modelTraining.Rmd](%22Rmd/modelTraining.md%22)
2. [simulateComplexSpectra.Rmd](%22Rmd/simulateComplexSpectra.md%22)
3. [mixtureAnalysis.Rmd](%22Rmd/mixtureAnalysis.md%22)
4. [evaluateModels.Rmd](%22Rmd/evaluateModels.md%22)
5. [makeMiscFigures.Rmd](%22Rmd/makeMiscFigures.md%22)
6. [illustrations.Rmd](%22Rmd/illustrations.md%22)
7. [makeTables.Rmd](%22Rmd/makeTables.md%22)

Directory structure
-------------------

To reproduce this analysis, the project directory must be structured as follows:

    |- README.Rmd
    |- data
    |  `- data.zip
    |- R
    |  |- createNewFeatureTbl.R
    |  |- ggplotTheme.R
    |  |- utilityFunctions.R
    |  |- prepareData.R
    |  |- extract.R
    |  `- preProcessSpec.R
    `- Rmd
       |- evaluateModels.Rmd
       |- illustrations.Rmd
       |- makeMiscFigures.Rmd
       |- mixtureAnalysis.Rmd
       |- modelTraining.Rmd
       |- simulateComplexSpectra.Rmd
       `- writeTables.Rmd

Code to run the analysis
------------------------

With a correctly prepared directory, the entire analysis can be run with a single line of code. Also, take note of the `instPckgs` variable in the "Prepare workspace" section below. If you system does not have all of these packages already installed, change this variable to `TRUE`.

To run the analysis:

``` r
render("path/to/README.Rmd", envir = new.env())
```

When this rmarkdown file is rendered, the following code is executed:

### Prepare workspace

``` r
instPckgs <- FALSE # change to "TRUE" to install all packages used in the analyses
unzipData <- FALSE # change to "TRUE" to unzip "data.zip" in the data directory

if(instPckgs) {
    install.packages(c("tidyverse",
                       "stringr",
                       "forcats",
                       "MALDIquant",
                       "MALDIquantForeign",
                       "caret",
                       "PRROC",
                       "xgboost",
                       "tictoc",
                       "devtools",
                       "openxlsx",
                       "knitr",
                       "rmarkdown"))
}

if(unzipData) {
    unzip("data/data.zip", overwrite = T, exdir = "../data")
}

library(rmarkdown)
library(tictoc)

tic()

# Make results and temp directories:
dir.create("temp")
dir.create("results")

# Set global parameters for analysis
mzTol <- 1.5 # m/z tolerance for feature extraction in Da
saveRDS(mzTol, "temp/mzTol.rds")
```

### Run Analysis Scripts

``` r
render("Rmd/modelTraining.Rmd", envir = new.env()) # Trains xgboost models
render("Rmd/simulateComplexSpectra.Rmd", envir = new.env()) # simulates polymicrobial spectra
render("Rmd/mixtureAnalysis.Rmd", envir = new.env()) # Imports spectra from experimental two-species mixtures
render("Rmd/evaluateModels.Rmd", envir = new.env()) # Calculate performance metrics and makes figures
render("Rmd/makeMiscFigures.Rmd", envir = new.env()) # Creates additional figures
render("Rmd/illustrations.Rmd", envir = new.env()) # Creates additional figures for illustrations
render("Rmd/writeTables.Rmd", envir = new.env()) # Creates supplementarly Tables
```

Session Info
------------

``` r
times <- toc()

cat(c("Execution time:", round((times$toc - times$tic)/60, 0), "min\n\n"))

session_info()
```
