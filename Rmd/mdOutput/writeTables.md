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
    ##  version  R version 3.4.2 (2017-09-28)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  tz       America/New_York            
    ##  date     2017-10-25                  
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
