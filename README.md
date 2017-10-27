Machine learning enables pathogen identification and detection of antimicrobial resistance from mass spectrometry analysis of membrane glycolipids
================
William E Fondrie,
October 26, 2017

-   [Introduction](#introduction)
-   [Directory structure](#directory-structure)
-   [Install necessary packages](#install-necessary-packages)
-   [Code to run the analysis](#code-to-run-the-analysis)
    -   [Prepare workspace](#prepare-workspace)
    -   [Run Analysis Scripts](#run-analysis-scripts)
-   [Session Info](#session-info)

Introduction
------------

This repository contains the code and data required to reproduce the analysis presented in [*Machine learning enables pathogen identification and detection of antimicrobial resistance from mass spectrometry analysis of membrane glycolipids*]().

**WARNING:** When attempting to reproduce this analysis, it should be noted that `Rmd/simulateCompleSpectra.Rmd` requires approximately 50 Gb of memory to execute. Additionally, the entire analysis time may take a number of hours depending on the hardware.

Order of analysis (see the Rmd subdirectory for details):
1. [modelTraining.Rmd](Rmd/mdOutput/modelTraining.md)
2. [simulateComplexSpectra.Rmd](Rmd/mdOutput/simulateComplexSpectra.md)
3. [mixtureAnalysis.Rmd](Rmd/mdOutput/mixtureAnalysis.md)
4. [evaluateModels.Rmd](Rmd/mdOutput/evaluateModels.md)
5. [makeMiscFigures.Rmd](Rmd/mdOutput/makeMiscFigures.md)
6. [illustrations.Rmd](Rmd/mdOutput/illustrations.md)
7. [makeTables.Rmd](Rmd/mdOutput/makeTables.md)

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

Install necessary packages
--------------------------

This analysis uses a number of R packages. To ensure that all of them are installed on your machine, you can execute the following:

``` r
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
                   "e1071",
                   "knitr",
                   "rmarkdown"))
```

Code to run the analysis
------------------------

With a correctly prepared directory, the entire analysis can be run with a single line of code. If it is your fist time running the analysis, you can unzip the `data/data.zip` into the `data` directory manually, or change `unzipData` to `TRUE` in the next section.

To run the analysis:

``` r
render("path/to/README.Rmd", envir = new.env())
```

When this rmarkdown file is rendered, the following code is executed:

### Prepare workspace

``` r
unzipData <- FALSE # change to "TRUE" to unzip "data.zip" in the data directory

if(unzipData) {
    unzip("data/data.zip", overwrite = T, exdir = "../data")
}

library(rmarkdown)

# Make results, temp and mdOutput directories:
dir.create("temp")
dir.create("results")
dir.create("Rmd/mdOutput")

# Set global parameters for analysis
mzTol <- 1.5 # m/z tolerance for feature extraction in Da
saveRDS(mzTol, "temp/mzTol.rds")
```

### Run Analysis Scripts

The `runAnalysis()` function takes a list of `.Rmd` files and knits them into GitHub compatible markdown files.

``` r
runAnalysis <- function(rmd) {
    lapply(rmd, 
           render,
           envir = new.env(),
           output_dir = "Rmd/mdOutput")
}

files <- c("modelTraining", # Trains xgboost models
           "simulateComplexSpectra", # simulates polymicrobial spectra
           "mixtureAnalysis", # Imports spectra from experimental two-species
           "evaluateModels", # Calculate performance metrics and makes figures
           "makeMiscFigures", # Creates additional figures
           "illustrations", # Creates additional figures for illustrations
           "writeTables") # Creates supplementarly Tables

files <- paste0("Rmd/", files, ".Rmd")
```

Now run the analysis:

``` r
runAnalysis(files)
```

Session Info
------------

``` r
session_info()
```

    ##  setting  value                       
    ##  version  R version 3.4.2 (2017-09-28)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  tz       America/New_York            
    ##  date     2017-10-26                  
    ## 
    ##  package            * version  date       source        
    ##  assertthat           0.2.0    2017-04-11 CRAN (R 3.4.2)
    ##  backports            1.1.1    2017-09-25 CRAN (R 3.4.1)
    ##  base               * 3.4.2    2017-09-28 local         
    ##  base64enc            0.1-3    2015-07-28 CRAN (R 3.4.1)
    ##  bindr                0.1      2016-11-13 CRAN (R 3.4.2)
    ##  bindrcpp           * 0.2      2017-06-17 CRAN (R 3.4.2)
    ##  broom                0.4.2    2017-02-13 CRAN (R 3.4.0)
    ##  caret              * 6.0-77   2017-09-07 CRAN (R 3.4.2)
    ##  cellranger           1.1.0    2016-07-27 CRAN (R 3.4.2)
    ##  class                7.3-14   2015-08-30 CRAN (R 3.4.2)
    ##  codetools            0.2-15   2016-10-05 CRAN (R 3.4.2)
    ##  colorspace           1.3-2    2016-12-14 CRAN (R 3.4.2)
    ##  compiler             3.4.2    2017-09-28 local         
    ##  CVST                 0.2-1    2013-12-10 CRAN (R 3.4.2)
    ##  data.table           1.10.4-2 2017-10-12 CRAN (R 3.4.2)
    ##  datasets           * 3.4.2    2017-09-28 local         
    ##  ddalpha              1.3.1    2017-09-27 CRAN (R 3.4.2)
    ##  DEoptimR             1.0-8    2016-11-19 CRAN (R 3.4.1)
    ##  devtools           * 1.13.3   2017-08-02 CRAN (R 3.4.2)
    ##  digest               0.6.12   2017-01-27 CRAN (R 3.4.2)
    ##  dimRed               0.1.0    2017-05-04 CRAN (R 3.4.2)
    ##  dplyr              * 0.7.4    2017-09-28 CRAN (R 3.4.2)
    ##  DRR                  0.0.2    2016-09-15 CRAN (R 3.4.2)
    ##  e1071                1.6-8    2017-02-02 CRAN (R 3.4.2)
    ##  evaluate             0.10.1   2017-06-24 CRAN (R 3.4.2)
    ##  forcats            * 0.2.0    2017-01-23 CRAN (R 3.4.2)
    ##  foreach              1.4.3    2015-10-13 CRAN (R 3.4.2)
    ##  foreign              0.8-69   2017-06-22 CRAN (R 3.4.2)
    ##  ggplot2            * 2.2.1    2016-12-30 CRAN (R 3.4.2)
    ##  glue                 1.1.1    2017-06-21 CRAN (R 3.4.2)
    ##  gower                0.1.2    2017-02-23 CRAN (R 3.4.2)
    ##  graphics           * 3.4.2    2017-09-28 local         
    ##  grDevices          * 3.4.2    2017-09-28 local         
    ##  grid                 3.4.2    2017-09-28 local         
    ##  gtable               0.2.0    2016-02-26 CRAN (R 3.4.2)
    ##  haven                1.1.0    2017-07-09 CRAN (R 3.4.2)
    ##  hms                  0.3      2016-11-22 CRAN (R 3.4.2)
    ##  htmltools            0.3.6    2017-04-28 CRAN (R 3.4.2)
    ##  httr                 1.3.1    2017-08-20 CRAN (R 3.4.2)
    ##  ipred                0.9-6    2017-03-01 CRAN (R 3.4.2)
    ##  iterators            1.0.8    2015-10-13 CRAN (R 3.4.1)
    ##  jsonlite             1.5      2017-06-01 CRAN (R 3.4.2)
    ##  kernlab              0.9-25   2016-10-03 CRAN (R 3.4.1)
    ##  knitr                1.17     2017-08-10 CRAN (R 3.4.2)
    ##  labeling             0.3      2014-08-23 CRAN (R 3.4.1)
    ##  lattice            * 0.20-35  2017-03-25 CRAN (R 3.4.2)
    ##  lava                 1.5.1    2017-09-27 CRAN (R 3.4.2)
    ##  lazyeval             0.2.0    2016-06-12 CRAN (R 3.4.2)
    ##  lubridate            1.6.0    2016-09-13 CRAN (R 3.4.2)
    ##  magrittr             1.5      2014-11-22 CRAN (R 3.4.2)
    ##  MALDIquant         * 1.16.4   2017-09-01 CRAN (R 3.4.2)
    ##  MALDIquantForeign  * 0.11     2017-09-06 CRAN (R 3.4.2)
    ##  MASS                 7.3-47   2017-02-26 CRAN (R 3.4.2)
    ##  Matrix               1.2-11   2017-08-21 CRAN (R 3.4.2)
    ##  memoise              1.1.0    2017-04-21 CRAN (R 3.4.2)
    ##  methods            * 3.4.2    2017-09-28 local         
    ##  mnormt               1.5-5    2016-10-15 CRAN (R 3.4.1)
    ##  ModelMetrics         1.1.0    2016-08-26 CRAN (R 3.4.2)
    ##  modelr               0.1.1    2017-07-24 CRAN (R 3.4.2)
    ##  munsell              0.4.3    2016-02-13 CRAN (R 3.4.2)
    ##  nlme                 3.1-131  2017-02-06 CRAN (R 3.4.2)
    ##  nnet                 7.3-12   2016-02-02 CRAN (R 3.4.2)
    ##  openxlsx           * 4.0.17   2017-03-23 CRAN (R 3.4.2)
    ##  parallel             3.4.2    2017-09-28 local         
    ##  pkgconfig            2.0.1    2017-03-21 CRAN (R 3.4.2)
    ##  plyr                 1.8.4    2016-06-08 CRAN (R 3.4.2)
    ##  prodlim              1.6.1    2017-03-06 CRAN (R 3.4.2)
    ##  PRROC              * 1.3      2017-04-21 CRAN (R 3.4.2)
    ##  psych                1.7.8    2017-09-09 CRAN (R 3.4.2)
    ##  purrr              * 0.2.4    2017-10-18 CRAN (R 3.4.2)
    ##  R6                   2.2.2    2017-06-17 CRAN (R 3.4.2)
    ##  Rcpp                 0.12.13  2017-09-28 CRAN (R 3.4.2)
    ##  RcppRoll             0.2.2    2015-04-05 CRAN (R 3.4.2)
    ##  readBrukerFlexData   1.8.5    2017-04-22 CRAN (R 3.4.1)
    ##  readMzXmlData        2.8.1    2015-09-16 CRAN (R 3.4.2)
    ##  readr              * 1.1.1    2017-05-16 CRAN (R 3.4.2)
    ##  readxl               1.0.0    2017-04-18 CRAN (R 3.4.2)
    ##  recipes              0.1.0    2017-07-27 CRAN (R 3.4.2)
    ##  reshape2             1.4.2    2016-10-22 CRAN (R 3.4.2)
    ##  rlang                0.1.2    2017-08-09 CRAN (R 3.4.2)
    ##  rmarkdown          * 1.6      2017-06-15 CRAN (R 3.4.2)
    ##  robustbase           0.92-7   2016-12-09 CRAN (R 3.4.2)
    ##  rpart                4.1-11   2017-03-13 CRAN (R 3.4.2)
    ##  rprojroot            1.2      2017-01-16 CRAN (R 3.4.2)
    ##  rstudioapi           0.7      2017-09-07 CRAN (R 3.4.2)
    ##  rvest                0.3.2    2016-06-17 CRAN (R 3.4.2)
    ##  scales               0.5.0    2017-08-24 CRAN (R 3.4.2)
    ##  sfsmisc              1.1-1    2017-06-08 CRAN (R 3.4.2)
    ##  splines              3.4.2    2017-09-28 local         
    ##  stats              * 3.4.2    2017-09-28 local         
    ##  stats4               3.4.2    2017-09-28 local         
    ##  stringi              1.1.5    2017-04-07 CRAN (R 3.4.1)
    ##  stringr            * 1.2.0    2017-02-18 CRAN (R 3.4.2)
    ##  survival             2.41-3   2017-04-04 CRAN (R 3.4.2)
    ##  tibble             * 1.3.4    2017-08-22 CRAN (R 3.4.2)
    ##  tictoc             * 1.0      2014-06-17 CRAN (R 3.4.1)
    ##  tidyr              * 0.7.2    2017-10-16 CRAN (R 3.4.2)
    ##  tidyselect           0.2.2    2017-10-10 CRAN (R 3.4.2)
    ##  tidyverse          * 1.1.1    2017-01-27 CRAN (R 3.4.2)
    ##  timeDate             3012.100 2015-01-23 CRAN (R 3.4.1)
    ##  tools                3.4.2    2017-09-28 local         
    ##  utils              * 3.4.2    2017-09-28 local         
    ##  withr                2.0.0    2017-07-28 CRAN (R 3.4.2)
    ##  xgboost            * 0.6-4    2017-01-05 CRAN (R 3.4.2)
    ##  XML                  3.98-1.9 2017-06-19 CRAN (R 3.4.1)
    ##  xml2                 1.1.1    2017-01-24 CRAN (R 3.4.2)
    ##  yaml                 2.1.14   2016-11-12 CRAN (R 3.4.2)
