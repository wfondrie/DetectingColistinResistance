Simulating Glycolipid MALDI-TOF-MS of Complex Samples
================

-   [Introduction](#introduction)
-   [Load Libraries & Prepare Workspace](#load-libraries-prepare-workspace)
-   [Retrieve File Information](#retrieve-file-information)
    -   [Select files](#select-files)
    -   [Import Spectra](#import-spectra)
    -   [Summarize Metadata](#summarize-metadata)
-   [Select Linear Combinations](#select-linear-combinations)
-   [Create Complex Spectra from Selections And Extract Features](#create-complex-spectra-from-selections-and-extract-features)
-   [Session Information](#session-information)

Introduction
============

This file contains the code to simulate glycolipid MALDI-TOF mass spectra of organism mixtures using linear combinations of the real glycolipid mass spectra in our library. **Warning**: Running this code requires a minimum of ~36 Gb of RAM to execute.

Load Libraries & Prepare Workspace
==================================

``` r
suppressMessages(library(tidyverse, quietly = T))
library(MALDIquant, quietly = T)
library(MALDIquantForeign, quietly = T)
library(stringr, quietly = T)
library(forcats, quietly = T)
library(tictoc, quietly = T)
library(devtools, quietly = T)

source("../R/preProcessSpec.R")
source("../R/extract.R")
source("../R/ggplotTheme.R")
source("../R/createNewFeatureTbl.R")

set.seed(657854)
tic()

featTol <- readRDS("../temp/mzTol.rds") # load feature extraction tolerance
```

Retrieve File Information
=========================

Select files
------------

``` r
files <- list.files("../data/fullLib", 
                    full.names = T, 
                    pattern = "mzXML$",
                    recursive = T)


specList <- unlist(map(files, preProcessSpec, hws = 80)) # Same preprocessing

spec <- map_df(specList, extractSpectra)

spec <- spec %>%
    group_by(id) %>%
    summarize(keep = min(mz) <= 1100 & max(mz) >= 2100) %>%
    filter(keep) %>%
    left_join(spec) %>%
    select(-keep)

fileDf <- tibble(fname = files) %>%
    mutate(id = str_match(fname, "([^\\\\^/]+).mzXML$")[ , 2]) %>%
    filter(id %in% unique(spec$id))

rm(specList)
rm(spec)
```

Import Spectra
--------------

``` r
# list of features, created in rmd/modelTraining.rmd
features <- readRDS("../temp/features.rds")

specList <- unlist(map(fileDf$fname, preProcessSpec, hws = 80)) # Same preprocessing

spec <- map_df(specList, extractSpectra)

# Free some memory
rm(specList)
```

Summarize Metadata
------------------

``` r
# fileDf <- tibble(fname = files) %>%
#     mutate(id = str_match(fname, "([^\\\\^/]+).mzXML$")[ , 2])

specInfo <- fileDf %>%
    mutate(type = str_match(fname, "([^\\^/]+)[\\/][^\\^/]+mzXML$")[ , 2],
           species = str_match(type, "^[^ ]+ [^ ]+")[ , 1],
           Ab = ifelse(str_detect(type, "Acinetobacter baumannii - res"), 
                       "pos", "other"),
           Kp = ifelse(str_detect(type, "Klebsiella pneumoniae - res"), 
                       "pos", "other"),
           Ab = as.factor(ifelse(str_detect(type, "Acinetobacter baumannii - sen"), 
                                 "neg", Ab)),
           Kp = as.factor(ifelse(str_detect(type, "Klebsiella pneumoniae - sen"), 
                                 "neg", Kp)))

summary(specInfo)
```

Select Linear Combinations
==========================

``` r
selectLinearCombos <- function(n, info) {
    
    regex <- "(Acinetobacter baumannii|Klebsiella pneumoniae)"
    pAbKp <- 0.3/n
    res <- sample(c(" - res", " - sen"), 1) # 50:50 shot of res or sen
    
    speciesSelect <- info %>%
        group_by(species) %>%
        summarize(num = length(id)) %>%
        mutate(weights = ifelse(str_detect(species, regex), 
                                pAbKp, 
                                (1-2*pAbKp)/(length(species)-2))) %>%
        sample_n(size = n, weight = weights)
    
    selected <- speciesSelect %>%
        left_join(info, by = "species") %>%
        filter(type != paste0(species, res)) %>%
        group_by(species) %>%
        sample_n(size = 1) %>%
        ungroup() %>%
        sample_frac(size = 1) %>%
        mutate(coeff = c(1, sample(c(1, 0.5, 0.25, 0.1), n-1, replace = T)))
    
    return(selected)
}

numSpecies <- rep(c(1:5), 1000) # number of spectra to simulate
#numSpecies <- rep(c(1:5), 1) # test line

combos <- tibble(spec_id = as.factor(paste0("spec_", 1:length(numSpecies))),
                 comp = map(numSpecies, selectLinearCombos, info = specInfo)) %>%
    unnest()

fctOrder <- c("neg", "other", "pos")

comboSummary <- combos %>%
    group_by(spec_id) %>%
    summarize(n = length(id),
              AbCoeff = ifelse(any(Ab != "other"), coeff[Ab != "other"] / max(coeff), 0),
              KpCoeff = ifelse(any(Kp != "other"), coeff[Kp != "other"] / max(coeff), 0),
              AbPres = ifelse(any(Ab == "pos"), "pos", "other"),
              KpPres = ifelse(any(Kp == "pos"), "pos", "other"),
              AbPres = as.factor(ifelse(any(Ab == "neg"), "neg", AbPres)),
              KpPres = as.factor(ifelse(any(Kp == "neg"), "neg", KpPres))) %>%
    rename(Ab = AbPres, Kp = KpPres) %>%
    mutate(Ab = fct_relevel(Ab, fctOrder),
           Kp = fct_relevel(Kp, fctOrder))

saveRDS(comboSummary, file = "../temp/complexSpectraSummary.rds")
saveRDS(combos, file = "../temp/complexComponents.rds")

summary(comboSummary)

comboSummary %>%
    group_by(n) %>%
    summarize(numAb = sum(Ab != "other"),
              numKp = sum(Kp != "other"),
              numOther = sum(Ab == "other" & Kp == "other"),
              numPosAb = sum(Ab == "pos"),
              numPosKp = sum(Kp == "neg"))
```

    ##       spec_id           n        AbCoeff         KpCoeff      
    ##  spec_1   :   1   Min.   :1   Min.   :0.000   Min.   :0.0000  
    ##  spec_10  :   1   1st Qu.:2   1st Qu.:0.000   1st Qu.:0.0000  
    ##  spec_100 :   1   Median :3   Median :0.000   Median :0.0000  
    ##  spec_1000:   1   Mean   :3   Mean   :0.209   Mean   :0.1932  
    ##  spec_1001:   1   3rd Qu.:4   3rd Qu.:0.250   3rd Qu.:0.1000  
    ##  spec_1002:   1   Max.   :5   Max.   :1.000   Max.   :1.0000  
    ##  (Other)  :4994                                               
    ##      Ab           Kp      
    ##  neg  : 750   neg  : 685  
    ##  other:3503   other:3658  
    ##  pos  : 747   pos  : 657  
    ##                           
    ##                           
    ##                           
    ##                           
    ## # A tibble: 5 x 6
    ##       n numAb numKp numOther numPosAb numPosKp
    ##   <int> <int> <int>    <int>    <int>    <int>
    ## 1     1   303   292      405      162      146
    ## 2     2   318   248      480      156      128
    ## 3     3   298   267      500      150      126
    ## 4     4   291   278      497      139      155
    ## 5     5   287   257      522      140      130

Create Complex Spectra from Selections And Extract Features
===========================================================

``` r
simSpectra <- function(df, specDat) {
    gc()
    
    df %>%
        left_join(specDat, by = c("type", "id")) %>%
        group_by(spec_id, id, type) %>%
        do(massSpecObj = createMassSpectrum(mass = .$mz, intensity = .$relInt,
                                            metaData = list(file = paste0("/cmb/",
                                                                          .$spec_id[1],
                                                                          ".mzXML")))) %>%
        group_by(spec_id) %>%
        do(massSpecObj = averageMassSpectra(.$massSpecObj, method = "sum")) %>%
        do(extractSpectra(.$massSpecObj)) %>%
        select(-id)
    
}

gc()

# Takes ~30 min
comboSpec <- combos %>% 
    group_by(spec_id) %>%
    do(simSpectra(df = ., specDat = filter(spec, id %in% .$id))) %>%
    rename(id = spec_id)

#saveRDS(comboSpec, "../temp/comboSpec.rds")

trainIdx <- readRDS("../temp/trainIdx.rds")


createMixtureFeatureTbl <- function(trainList, specDf, summaryDat, mzTol) {
    suffix <- if(str_detect(trainList$regex, "Acineto")) "_Ab" else "_Kp"
    
    # extract specified features
    cat(paste0("\nCreating features for ", suffix, "...\n"))
    
    feattbl <- specDf %>%
      group_by(id) %>%
      do(extractFeatures(., featureVec = trainList$features, tol = mzTol))
    
    saveRDS(feattbl, file = paste0("../temp/mixtureFeatureInfo", suffix, ".rds"))
    
    mlFeat <- feattbl %>%
        select(id, feat, relInt) %>%
        group_by(id) %>%
        mutate(relInt = relInt / max(relInt)) %>%
        ungroup() %>%
        spread(key = feat, value = relInt, fill = 0) %>%
        left_join(summaryDat, by = c("id" = "spec_id")) 
    
    return(mlFeat)
}

# extract features for prediction
mixtureDatList <- map(trainIdx,createMixtureFeatureTbl, 
                      specDf = comboSpec,
                      summaryDat = comboSummary,
                      mzTol = featTol)


saveRDS(mixtureDatList, file = "../temp/mixtureDatList.rds")
```

Session Information
===================

``` r
times <- toc()

cat(c("Execution time:", round((times$toc - times$tic)/60, 0), "min\n\n"))

session_info()
```

    ## 10653.45 sec elapsed
    ## Execution time: 178 min
    ## 
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
    ##  broom                0.4.2    2017-02-13 CRAN (R 3.4.0)
    ##  car                  2.1-4    2016-12-02 CRAN (R 3.4.0)
    ##  caret              * 6.0-76   2017-04-18 CRAN (R 3.4.0)
    ##  cellranger           1.1.0    2016-07-27 CRAN (R 3.4.0)
    ##  codetools            0.2-15   2016-10-05 CRAN (R 3.4.0)
    ##  colorspace           1.3-2    2016-12-14 CRAN (R 3.4.0)
    ##  compiler             3.4.0    2017-04-21 local         
    ##  data.table           1.10.4   2017-02-01 CRAN (R 3.4.0)
    ##  datasets           * 3.4.0    2017-04-21 local         
    ##  DBI                  0.6-1    2017-04-01 CRAN (R 3.4.0)
    ##  devtools           * 1.13.2   2017-06-02 CRAN (R 3.4.0)
    ##  digest               0.6.12   2017-01-27 CRAN (R 3.4.0)
    ##  dplyr              * 0.5.0    2016-06-24 CRAN (R 3.4.0)
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
