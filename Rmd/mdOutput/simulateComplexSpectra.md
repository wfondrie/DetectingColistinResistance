Simulating Glycolipid MALDI-TOF-MS of Complex Samples
================

-   [Introduction](#introduction)
-   [Load Libraries & Prepare Workspace](#load-libraries-prepare-workspace)
-   [Retrieve File Information](#retrieve-file-information)
-   [Select Linear Combinations](#select-linear-combinations)
-   [Create Complex Spectra from Selections And Extract Features](#create-complex-spectra-from-selections-and-extract-features)
-   [Session Information](#session-information)

Introduction
============

This file contains the code to simulate glycolipid MALDI-TOF mass spectra of organism mixtures using linear combinations of the real glycolipid mass spectra in our library. **Warning**: Running this code requires a minimum of ~25 Gb of RAM to execute. The machine this was run on had 64 Gb of RAM and 2x AMD Opteron 6272 processors (2.10 GHz).

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

set.seed(6578548)
tic()

featTol <- readRDS("../temp/mzTol.rds") # load feature extraction tolerance
```

Retrieve File Information
=========================

``` r
files <- list.files("../data/fullLib", 
                    full.names = T, 
                    pattern = "mzXML$",
                    recursive = T)
#files <- sample(files, 50)

specInfo <- tibble(fname = files) %>%
    mutate(type = str_match(fname, "([^\\^/]+)[\\/][^\\^/]+mzXML$")[ , 2],
           id = str_match(fname, "([^\\^/]+).mzXML$")[ , 2],
           species = str_match(type, "^[^ ]+ [^ ]+")[ , 1],
           Ab = ifelse(str_detect(type, "Acinetobacter baumannii - res"), "pos", "other"),
           Kp = ifelse(str_detect(type, "Klebsiella pneumoniae - res"), "pos", "other"),
           Ab = as.factor(ifelse(str_detect(type, "Acinetobacter baumannii - sen"), "neg", Ab)),
           Kp = as.factor(ifelse(str_detect(type, "Klebsiella pneumoniae - sen"), "neg", Kp)))

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

    ##       spec_id           n        AbCoeff          KpCoeff      
    ##  spec_1   :   1   Min.   :1   Min.   :0.0000   Min.   :0.0000  
    ##  spec_10  :   1   1st Qu.:2   1st Qu.:0.0000   1st Qu.:0.0000  
    ##  spec_100 :   1   Median :3   Median :0.0000   Median :0.0000  
    ##  spec_1000:   1   Mean   :3   Mean   :0.1966   Mean   :0.2015  
    ##  spec_1001:   1   3rd Qu.:4   3rd Qu.:0.1000   3rd Qu.:0.1000  
    ##  spec_1002:   1   Max.   :5   Max.   :1.0000   Max.   :1.0000  
    ##  (Other)  :4994                                                
    ##      Ab           Kp      
    ##  neg  : 674   neg  : 692  
    ##  other:3641   other:3581  
    ##  pos  : 685   pos  : 727  
    ##                           
    ##                           
    ##                           
    ##                           
    ## # A tibble: 5 x 6
    ##       n numAb numKp numOther numPosAb numPosKp
    ##   <int> <int> <int>    <int>    <int>    <int>
    ## 1     1   277   330      393      140      158
    ## 2     2   268   266      518      129      143
    ## 3     3   302   274      496      148      138
    ## 4     4   270   274      528      145      120
    ## 5     5   242   275      537      123      133

Create Complex Spectra from Selections And Extract Features
===========================================================

``` r
# list of features, created in rmd/modelTraining.rmd
features <- readRDS("../temp/features.rds")

specList <- preProcessSpec(files, hws = 80) # Same preprocessing

# rm files that contain multiple spectra
multiSpecIdx <- map_lgl(specList, ~ metaData(.)$num > 1)
specList <- specList[!multiSpecIdx]

spec <- map_df(specList, extractSpectra)

gc() # garbage collection. Free as much memory as possible.

# Takes ~30 min
comboSpec <- combos %>% 
    left_join(spec, by = c("type", "id")) %>%
    group_by(spec_id) %>%
    group_by(spec_id, id, type) %>%
    do(massSpecObj = createMassSpectrum(mass = .$mz, intensity = .$relInt, 
                                        metaData = list(file = paste0("/cmb/", .$spec_id[1], ".mzXML")))) %>%
    group_by(spec_id) %>%
    do(massSpecObj = averageMassSpectra(.$massSpecObj, method = "sum")) %>%
    do(extractSpectra(.$massSpecObj))

saveRDS(comboSpec, "../temp/comboSpec.rds")

trainIdx <- readRDS("../temp/trainIdx.rds")


createMixtureFeatureTbl <- function(trainList, specDf, summaryDat, mzTol) {
    suffix <- if(str_detect(trainList$regex, "Acineto")) "_Ab" else "_Kp"
    
    # extract specified features
    cat(paste0("Creating features for ", suffix, "...\n"))
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

    ## 15352.9 sec elapsed
    ## Execution time: 256 min
    ## 
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
