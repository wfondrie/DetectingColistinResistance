Insert Paper Title Here
================
William E Fondrie

Order of analysis (see the Rmd subdirectory for details):
1. [modelTraining.Rmd](%22Rmd/modelTraining.md%22)
2. [simulateComplexSpectra.Rmd](%22Rmd/simulateComplexSpectra.md%22)
3. [mixtureAnalysis.Rmd](%22Rmd/mixtureAnalysis.md%22)
4. [evaluateModels.Rmd](%22Rmd/evaluateModels.md%22)
5. [makeMiscFigures.Rmd](%22Rmd/makeMiscFigures.md%22)
6. [illustrations.Rmd](%22Rmd/illustrations.md%22)
7. [makeTables.Rmd](%22Rmd/makeTables.md%22)

Code to run the analysis
------------------------

### Prepare workspace

``` r
instPckgs <- F # change to "T" to install all packages used in the analyses

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
                       "knitr"))
}

library(knitr)
library(tictoc)

tic()

# Make results and temp directories:
dir.create("../temp")
dir.create("../results")

# Set global parameters for analysis
mzTol <- 1.5 # m/z tolerance for feature extraction in Da
saveRDS(mzTol, "../temp/mzTol.rds")

knit("modelTraining.Rmd") # Trains xgboost models
knit("simulateComplexSpectra.Rmd") # simulates polymicrobial spectra
knit("mixtureAnalysis.Rmd") # Imports spectra from experimental two-species mixtures
knit("evaluateModels.Rmd") # Calculate performance metrics and makes figures
knit("makeMiscFigures.Rmd") # Creates additional figures
knit("illustrations.Rmd") # Creates additional figures for illustrations
knit("writeTables.Rmd") # Creates supplementarly Tables
```

### Copy .md files to top directory

``` r
file.copy("README.md", "..", overwrite=T)
```

    ## [1] FALSE
