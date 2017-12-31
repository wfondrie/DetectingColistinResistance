Write Supplementary Tables
================
William E Fondrie

-   [Load Libraries](#load-libraries)
-   [Supplementary Table 1](#supplementary-table-1)
-   [Supplementary Table 2](#supplementary-table-2)
-   [Supplementary Table 4](#supplementary-table-4)
-   [Session Info](#session-info)

Load Libraries
--------------

``` r
# data manipulation
suppressMessages(library(tidyverse, quietly = T))
library(stringr, quietly = T)
library(forcats, quietly = T)
library(devtools)

# writing to excel files
library(openxlsx, quietly = T)

# other functions
source("../R/utilityFunctions.R")

# for reproducibility
set.seed(937426)
```

Supplementary Table 1
---------------------

``` r
testList <- readRDS("../temp/testDatList.rds")
abRes <- readRDS("../temp/abRes.rds")
kpRes <- readRDS("../temp/kpRes.rds")

createTbl <- function(summDf, resDf) {
    
    summDf %>%
    cbind(resDf) %>%
    mutate(Class = ifelse(truth == "pos", "Colistin-Resistant", "Other Species"),
           Class = ifelse(truth == "neg", "Colistin-Susceptible", Class),
           type = str_replace(type, " - .+$", ""),
           type = str_replace(type, " [^ ]*[0-9].*$", "")) %>%
    rename(Species = type, 
           `Colistin-Resistant Score` = pos,
           `Colistin-Susceptible Score` = neg,
           `Other Species Score` = other,
           `Species Score` = speciesVsOther,
           `Susceptible or Other Species Score` = posVsAll) %>%
    select(Species,
           Class,
           `Colistin-Resistant Score`,
           `Colistin-Susceptible Score`,
           `Other Species Score`,
           `Species Score`,
           `Susceptible or Other Species Score`,
           starts_with("mz"))
    
}

abTest <- createTbl(testList$Ab, abRes)
kpTest <- createTbl(testList$Kp, kpRes)

testSetRes <- list(`A. baumannii Test Set` = abTest,
                   `K. pneumoniae Test Set` = kpTest)

write.xlsx(testSetRes, "../results/testSetResults.xlsx")
```

Supplementary Table 2
---------------------

``` r
mixedSumm <- readRDS("../temp/complexSpectraSummary.rds") %>%
    mutate(Ab2 = ifelse(Ab == "pos", "Colistin-Resistant", "Other Species"),
           Ab2 = ifelse(Ab == "neg", "Colistin-Susceptible", Ab2),
           Kp2 = ifelse(Kp == "pos", "Colistin-Resistant", "Other Species"),
           Kp2 = ifelse(Kp == "neg", "Colistin-Susceptible", Kp2)) %>%
    select(-Ab, -Kp) %>%
    rename(`Spectrum ID` = spec_id,
           `Number of Species` = n,
           `A. baumannii Class` = Ab2,
           `A. baumannii Weight` = AbCoeff,
           `K. pneumoniae Class` = Kp2,
           `K. pneumoniae Weight` = KpCoeff)


mixedComponents <- readRDS("../temp/complexComponents.rds") %>%
    select(`Spectrum ID` = spec_id,
           `Species` = species,
           `Weight` = coeff,
           `mzXML File` = fname)
    

mixedList <- readRDS("../temp/mixtureDatList.rds")
mixedAbRes <- readRDS("../temp/mixedAbRes.rds")
mixedKpRes <- readRDS("../temp/mixedKpRes.rds")

createMixedTbl <- function(summDf, resDf) {
    
    summDf %>%
    left_join(resDf) %>%
    mutate(Class = ifelse(truth == "pos", "Colistin-Resistant", "Other Species"),
           Class = ifelse(truth == "neg", "Colistin-Susceptible", Class)) %>%
    rename(`Spectrum ID` = id, 
           `Colistin-Resistant Score` = pos,
           `Colistin-Susceptible Score` = neg,
           `Other Species Score` = other,
           `Species Score` = speciesVsOther,
           `Susceptible or Other Species Score` = posVsAll) %>%
    select(`Spectrum ID`,
           Class,
           `Colistin-Resistant Score`,
           `Colistin-Susceptible Score`,
           `Other Species Score`,
           `Species Score`,
           `Susceptible or Other Species Score`,
           starts_with("mz"))
    
}

abMixed <- createMixedTbl(mixedList$Ab, mixedAbRes)
kpMixed <- createMixedTbl(mixedList$Kp, mixedKpRes)

simMixTbl <- list(`Sim Spectra Summary` = mixedSumm,
                  `Sim Spectra Weights` = mixedComponents,
                  `A. baummannii Results` = abMixed,
                  `K. pneumoniae Results` = kpMixed)

write.xlsx(simMixTbl, "../results/simMixResults.xlsx")
```

Supplementary Table 4
---------------------

``` r
twoSumm <- readRDS("../temp/twoSpeciesSpecInfo.rds")
twoList <- readRDS("../temp/twoSpeciesDatList.rds")
twoAbRes <- readRDS("../temp/twoAbRes.rds")
twoKpRes <- readRDS("../temp/twoKpRes.rds")

twoAb <- twoList$Ab %>%
    left_join(twoAbRes) %>%
    filter(!is.na(percentEc),
           Kp == "other") %>%
    mutate(Class = ifelse(truth == "pos", "Colistin-Resistant", "Other Species"),
           Class = ifelse(truth == "neg", "Colistin-Susceptible", Class),
           Class = ifelse(percentEc == 100, "Other Species", Class)) %>%
    rename(`Percent A. baumannii` = percentTarget,
           `Percent E. coli` = percentEc,
           `mzXML file` = fname,
           `Colistin-Resistant Score` = pos,
           `Colistin-Susceptible Score` = neg,
           `Other Species Score` = other,
           `Species Score` = speciesVsOther,
           `Susceptible or Other Species Score` = posVsAll) %>%
    select(`Percent A. baumannii`,
           `Percent E. coli`,
           Class,
           `Colistin-Resistant Score`,
           `Colistin-Susceptible Score`,
           `Other Species Score`,
           `Species Score`,
           `Susceptible or Other Species Score`,
           `mzXML file`,
           starts_with("mz"))

twoKp <- twoList$Kp %>%
    left_join(twoKpRes) %>%
    filter(!is.na(percentEc),
           Ab == "other") %>%
    mutate(Class = ifelse(truth == "pos", "Colistin-Resistant", "Other Species"),
           Class = ifelse(truth == "neg", "Colistin-Susceptible", Class),
           Class = ifelse(percentEc == 100, "Other Species", Class)) %>%
    rename(`Percent K. pneumoniae` = percentTarget,
           `Percent E. coli` = percentEc,
           `mzXML file` = fname,
           `Colistin-Resistant Score` = pos,
           `Colistin-Susceptible Score` = neg,
           `Other Species Score` = other,
           `Species Score` = speciesVsOther,
           `Susceptible or Other Species Score` = posVsAll) %>%
    select(`Percent K. pneumoniae`,
           `Percent E. coli`,
           Class,
           `Colistin-Resistant Score`,
           `Colistin-Susceptible Score`,
           `Other Species Score`,
           `Species Score`,
           `Susceptible or Other Species Score`,
           `mzXML file`,
           starts_with("mz"))

twoMixTbl <- list(`A. baummannii Two-Species Mixtures` = twoAb,
               `K. pneumoniae Two-Species Mixtures` = twoKp)

write.xlsx(twoMixTbl, "../results/twoMixResults.xlsx")
```

Session Info
------------

``` r
session_info()
```

    ##  setting  value                       
    ##  version  R version 3.4.0 (2017-04-21)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  tz       America/New_York            
    ##  date     2017-12-30                  
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
    ##  tictoc             * 1.0      2014-06-17 CRAN (R 3.4.0)
    ##  tidyr              * 0.6.3    2017-05-15 CRAN (R 3.4.0)
    ##  tidyverse          * 1.1.1    2017-01-27 CRAN (R 3.4.0)
    ##  tools                3.4.0    2017-04-21 local         
    ##  utils              * 3.4.0    2017-04-21 local         
    ##  withr                1.0.2    2016-06-20 CRAN (R 3.4.0)
    ##  xgboost            * 0.6-4    2017-01-05 CRAN (R 3.4.1)
    ##  XML                  3.98-1.7 2017-05-03 CRAN (R 3.4.0)
    ##  xml2                 1.1.1    2017-01-24 CRAN (R 3.4.0)
    ##  yaml                 2.1.14   2016-11-12 CRAN (R 3.4.0)
