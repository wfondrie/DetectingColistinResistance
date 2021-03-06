Pathogen Identification Direct from Polymicrobial Specimens Using Membrane Glycolipids
================
William E Fondrie,
January 18, 2018

-   [Introduction](#introduction)
-   [Directory structure](#directory-structure)
-   [Install necessary packages](#install-necessary-packages)
-   [Code to run the analysis](#code-to-run-the-analysis)
    -   [Prepare workspace](#prepare-workspace)
    -   [Run Analysis Scripts](#run-analysis-scripts)
-   [Session Info](#session-info)

Introduction
------------

This repository contains the code required to reproduce the analysis presented in [*Pathogen Identification Direct from Polymicrobial Specimens Using Membrane Glycolipids*](https://doi.org/10.1038/s41598-018-33681-8). The data required for this analysis (`data.zip` below) can be downloaded from [here](). *We are currently in the process of working with the UMB Office of Technology Transfer to provide a download link. In the interim, please contact dgoodlett@rx.umaryland.edu and the data will be provided to you promptly.*

**WARNING:** When attempting to reproduce this analysis, it should be noted that `Rmd/simulateCompleSpectra.Rmd` requires approximately 36 Gb of memory to execute. Additionally, the entire analysis time may take a number of hours depending on the hardware.

Order of analysis (see the Rmd subdirectory for details):  

1. [modelTraining.Rmd](Rmd/mdOutput/modelTraining.md)
2. [simulateComplexSpectra.Rmd](Rmd/mdOutput/simulateComplexSpectra.md)  
3. [mixtureAnalysis.Rmd](Rmd/mdOutput/mixtureAnalysis.md)  
4. [evaluateModels.Rmd](Rmd/mdOutput/evaluateModels.md)  
5. [baselineClassifiers.Rmd](Rmd/mdOutput/baselineClassifiers.md)  
6. [makeMiscFigures.Rmd](Rmd/mdOutput/makeMiscFigures.md)  
7. [illustrations.Rmd](Rmd/mdOutput/illustrations.md)  
8. [makeTables.Rmd](Rmd/mdOutput/makeTables.md)  

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
       |- baselineClassifiers.Rmd
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
rmarkdown::render("path/to/README.Rmd", envir = new.env())
```

When this rmarkdown file is rendered, the following code is executed:

### Prepare workspace

``` r
tictoc::tic()
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
           "baselineClassifiers", # Implements simple baseline classifiers
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
times <- tictoc::toc()
cat(c("Execution time:", round((times$toc - times$tic)/60, 0), "min\n\n"))
session_info()
```

    ## 3490.38 sec elapsed
    ## Execution time: 58 min
    ## 
    ##  setting  value                       
    ##  version  R version 3.4.0 (2017-04-21)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  tz       America/New_York            
    ##  date     2018-01-18                  
    ## 
    ##  package            * version  date       source        
    ##  assertthat           0.2.0    2017-04-11 CRAN (R 3.3.3)
    ##  backports            1.1.0    2017-05-22 CRAN (R 3.4.0)
    ##  base               * 3.4.0    2017-04-21 local         
    ##  base64enc            0.1-3    2015-07-28 CRAN (R 3.2.1)
    ##  broom              * 0.4.2    2017-02-13 CRAN (R 3.4.0)
    ##  car                  2.1-4    2016-12-02 CRAN (R 3.4.0)
    ##  caret              * 6.0-76   2017-04-18 CRAN (R 3.4.0)
    ##  cellranger           1.1.0    2016-07-27 CRAN (R 3.4.0)
    ##  class                7.3-14   2015-08-30 CRAN (R 3.4.0)
    ##  codetools            0.2-15   2016-10-05 CRAN (R 3.4.0)
    ##  colorspace           1.3-2    2016-12-14 CRAN (R 3.4.0)
    ##  compiler             3.4.0    2017-04-21 local         
    ##  data.table           1.10.4   2017-02-01 CRAN (R 3.4.0)
    ##  datasets           * 3.4.0    2017-04-21 local         
    ##  DBI                  0.6-1    2017-04-01 CRAN (R 3.4.0)
    ##  devtools           * 1.13.2   2017-06-02 CRAN (R 3.4.0)
    ##  digest               0.6.12   2017-01-27 CRAN (R 3.4.0)
    ##  dplyr              * 0.5.0    2016-06-24 CRAN (R 3.4.0)
    ##  e1071                1.6-8    2017-02-02 CRAN (R 3.4.0)
    ##  evaluate             0.10     2016-10-11 CRAN (R 3.4.0)
    ##  forcats            * 0.2.0    2017-01-23 CRAN (R 3.4.0)
    ##  foreach              1.4.3    2015-10-13 CRAN (R 3.4.0)
    ##  foreign              0.8-68   2017-04-24 CRAN (R 3.4.0)
    ##  ggplot2            * 2.2.1    2016-12-30 CRAN (R 3.4.0)
    ##  graphics           * 3.4.0    2017-04-21 local         
    ##  grDevices          * 3.4.0    2017-04-21 local         
    ##  grid                 3.4.0    2017-04-21 local         
    ##  gtable               0.2.0    2016-02-26 CRAN (R 3.4.0)
    ##  haven                1.0.0    2016-09-23 CRAN (R 3.4.0)
    ##  hms                  0.3      2016-11-22 CRAN (R 3.4.0)
    ##  htmltools            0.3.6    2017-04-28 CRAN (R 3.4.0)
    ##  httr                 1.2.1    2016-07-03 CRAN (R 3.4.0)
    ##  iterators            1.0.8    2015-10-13 CRAN (R 3.4.0)
    ##  jsonlite             1.4      2017-04-08 CRAN (R 3.4.0)
    ##  knitr                1.16     2017-05-18 CRAN (R 3.4.0)
    ##  labeling             0.3      2014-08-23 CRAN (R 3.4.0)
    ##  lattice            * 0.20-35  2017-03-25 CRAN (R 3.4.0)
    ##  lazyeval             0.2.0    2016-06-12 CRAN (R 3.4.0)
    ##  lme4                 1.1-13   2017-04-19 CRAN (R 3.4.0)
    ##  lubridate            1.6.0    2016-09-13 CRAN (R 3.4.0)
    ##  magrittr             1.5      2014-11-22 CRAN (R 3.4.0)
    ##  MALDIquant         * 1.16.2   2017-04-04 CRAN (R 3.4.0)
    ##  MALDIquantForeign  * 0.10     2015-11-01 CRAN (R 3.4.0)
    ##  MASS                 7.3-47   2017-02-26 CRAN (R 3.4.0)
    ##  Matrix               1.2-10   2017-04-28 CRAN (R 3.4.0)
    ##  MatrixModels         0.4-1    2015-08-22 CRAN (R 3.4.0)
    ##  memoise              1.1.0    2017-04-21 CRAN (R 3.4.0)
    ##  methods            * 3.4.0    2017-04-21 local         
    ##  mgcv                 1.8-17   2017-02-08 CRAN (R 3.4.0)
    ##  minqa                1.2.4    2014-10-09 CRAN (R 3.4.0)
    ##  mnormt               1.5-5    2016-10-15 CRAN (R 3.4.0)
    ##  ModelMetrics         1.1.0    2016-08-26 CRAN (R 3.4.0)
    ##  modelr               0.1.0    2016-08-31 CRAN (R 3.4.0)
    ##  munsell              0.4.3    2016-02-13 CRAN (R 3.4.0)
    ##  nlme                 3.1-131  2017-02-06 CRAN (R 3.4.0)
    ##  nloptr               1.0.4    2014-08-04 CRAN (R 3.4.0)
    ##  nnet                 7.3-12   2016-02-02 CRAN (R 3.4.0)
    ##  openxlsx           * 4.0.17   2017-03-23 CRAN (R 3.4.2)
    ##  parallel             3.4.0    2017-04-21 local         
    ##  pbkrtest             0.4-7    2017-03-15 CRAN (R 3.4.0)
    ##  plyr                 1.8.4    2016-06-08 CRAN (R 3.4.0)
    ##  PRROC              * 1.3      2017-04-21 CRAN (R 3.4.0)
    ##  psych                1.7.5    2017-05-03 CRAN (R 3.4.0)
    ##  purrr              * 0.2.2.2  2017-05-11 CRAN (R 3.4.0)
    ##  quantreg             5.33     2017-04-18 CRAN (R 3.4.0)
    ##  R6                   2.2.1    2017-05-10 CRAN (R 3.4.0)
    ##  Rcpp                 0.12.11  2017-05-22 CRAN (R 3.4.0)
    ##  readBrukerFlexData   1.8.5    2017-04-22 CRAN (R 3.4.0)
    ##  readMzXmlData        2.8.1    2015-09-16 CRAN (R 3.4.0)
    ##  readr              * 1.1.1    2017-05-16 CRAN (R 3.4.0)
    ##  readxl               1.0.0    2017-04-18 CRAN (R 3.4.0)
    ##  reshape2             1.4.2    2016-10-22 CRAN (R 3.4.0)
    ##  rlang                0.1.1    2017-05-18 CRAN (R 3.4.0)
    ##  rmarkdown          * 1.6      2017-06-15 CRAN (R 3.4.2)
    ##  rprojroot            1.2      2017-01-16 CRAN (R 3.4.0)
    ##  rstudioapi           0.6      2016-06-27 CRAN (R 3.4.0)
    ##  rvest                0.3.2    2016-06-17 CRAN (R 3.4.0)
    ##  scales               0.4.1    2016-11-09 CRAN (R 3.4.0)
    ##  SparseM              1.77     2017-04-23 CRAN (R 3.4.0)
    ##  splines              3.4.0    2017-04-21 local         
    ##  stats              * 3.4.0    2017-04-21 local         
    ##  stats4               3.4.0    2017-04-21 local         
    ##  stringi              1.1.5    2017-04-07 CRAN (R 3.4.0)
    ##  stringr            * 1.2.0    2017-02-18 CRAN (R 3.4.0)
    ##  tibble             * 1.3.1    2017-05-17 CRAN (R 3.4.0)
    ##  tictoc               1.0      2014-06-17 CRAN (R 3.4.0)
    ##  tidyr              * 0.6.3    2017-05-15 CRAN (R 3.4.0)
    ##  tidyverse          * 1.1.1    2017-01-27 CRAN (R 3.4.0)
    ##  tools                3.4.0    2017-04-21 local         
    ##  utils              * 3.4.0    2017-04-21 local         
    ##  withr                1.0.2    2016-06-20 CRAN (R 3.4.0)
    ##  xgboost            * 0.6-4    2017-01-05 CRAN (R 3.4.1)
    ##  XML                  3.98-1.7 2017-05-03 CRAN (R 3.4.0)
    ##  xml2                 1.1.1    2017-01-24 CRAN (R 3.4.0)
    ##  yaml                 2.1.14   2016-11-12 CRAN (R 3.4.0)
