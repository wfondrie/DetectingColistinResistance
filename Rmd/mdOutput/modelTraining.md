Training XGBoost Models
================

-   [Introduction](#introduction)
-   [Load Libraries & Prepare Workspace](#load-libraries-prepare-workspace)
-   [Partition Data](#partition-data)
-   [Select Features From Training Data](#select-features-from-training-data)
-   [Extract features from Spectra](#extract-features-from-spectra)
-   [Train the XGBoost Tree Models](#train-the-xgboost-tree-models)
    -   [Prepare the Data for Model Training](#prepare-the-data-for-model-training)
    -   [Tune Hyperparameter by Grid Search](#tune-hyperparameter-by-grid-search)
    -   [Refine Models](#refine-models)
    -   [Final Model](#final-model)
-   [Prepare Test Set](#prepare-test-set)
-   [Session Information](#session-information)

Introduction
------------

This file contains the code used to build the XGBoost machine learning models.

Load Libraries & Prepare Workspace
----------------------------------

``` r
# data manipulation
suppressMessages(library(tidyverse, quietly = T))
library(stringr, quietly = T)
library(forcats, quietly = T)

# handling MALDI spectra
library(MALDIquant, quietly = T)
library(MALDIquantForeign, quietly = T)

# machine learning
library(caret, quietly = T)
library(PRROC, quietly = T)
library(xgboost, quietly = T)

# timing and session information
library(tictoc)
library(devtools)

# import helper functions
source("../R/preProcessSpec.R")
source("../R/extract.R")
source("../R/prepareData.R")

# for reproducibility
set.seed(1832939)
tic()

featTol <- readRDS("../temp/mzTol.rds") # load feature extraction tolerance
```

Partition Data
--------------

``` r
files <- list.files("../data/fullLib", 
                    full.names = T, 
                    pattern = "mzXML$",
                    recursive = T)

specInfo <- tibble(fname = files) %>%
    mutate(type = str_match(fname, "([^\\^/]+)[\\/][^\\^/]+mzXML$")[ , 2],
           id = str_match(fname, "([^\\^/]+).mzXML$")[ , 2],
           Ab = ifelse(str_detect(type, "Acinetobacter baumannii - res"), "pos", "other"),
           Kp = ifelse(str_detect(type, "Klebsiella pneumoniae - res"), "pos", "other"),
           Ab = as.factor(ifelse(str_detect(type, "Acinetobacter baumannii - sen"), "neg", Ab)),
           Kp = as.factor(ifelse(str_detect(type, "Klebsiella pneumoniae - sen"), "neg", Kp)))

trainIdx <- list(Ab = list(idx = createDataPartition(specInfo$Ab, p = 0.6, list = F),
                           regex = "Acinetobacter baumannii - res"),
                 Kp = list(idx = createDataPartition(specInfo$Kp, p = 0.6, list = F),
                           regex = "Klebsiella pneumoniae - res"))
```

Select Features From Training Data
----------------------------------

``` r
determineFeatures <- function(trainList, n = 50) {
    
    lab <- gsub(" - res", "", trainList$regex)
    
    # Preprocess spectra
    cat(paste0("Importing spectra for ",  lab, "...\n"))
    spec <- preProcessSpec(files[trainList$idx], hws = 80)
    multiSpecIdx <- map_lgl(spec, ~ metaData(.)$num > 1)
    spec <- spec[!multiSpecIdx]

    # Detect peaks
    cat(paste0("Peak picking for ",  lab, "...\n"))
    peaks <- detectPeaks(spec, halfWindowSize = 80, SNR = 5)
    peaks <- binPeaks(peaks, tolerance = 0.5)
    
    # extract peaks
    peakDat <- extractPeaks(peaks, spec)
    
    # extract features
    cat(paste0("Determining top ", n, " features for ", lab, "...\n"))
    features <- peakDat %>%
        mutate(sel = str_detect(type, trainList$regex)) %>%
        filter(sel) %>%
        group_by(sel, mz) %>%
        summarize(relInt = sum(relInt)) %>%
        group_by(sel) %>%
        mutate(relInt = relInt / max(relInt)) %>%
        arrange(desc(relInt)) %>%
        filter(relInt >= relInt[n]) %>%
        arrange(desc(mz))
    
    return(features)
}

featureList <- map(trainIdx, determineFeatures)
trainIdx$Ab$features <- unique(featureList$Ab$mz)
trainIdx$Kp$features <- unique(featureList$Kp$mz)

saveRDS(trainIdx, "../temp/trainIdx.rds")

features <- tibble(response = names(featureList), featList = featureList) %>%
    unnest() %>%
    select(response, mz)

features <- unique(features)
saveRDS(features, file = "../temp/features.rds")
```

Extract features from Spectra
-----------------------------

``` r
createFeatureTbl <- function(trainList, mzTol, testSet = F) {
    idx <- if(testSet) -trainList$idx else trainList$idx
    suffix <- if(testSet) "_test" else "_train"
    suffix2 <- if(str_detect(trainList$regex, "Acineto")) "_Ab" else "_Kp"
    
    # Preprocess spectra
    cat(paste0("\nImporting spectra for ", suffix2, "...\n"))
    spec <- preProcessSpec(files[idx], hws = 80)
    multiSpecIdx <- map_lgl(spec, ~ metaData(.)$num > 1)
    spec <- spec[!multiSpecIdx]
    
    # extract spectra 
    cat(paste0("Extracting spectra for ", suffix2, "...\n"))
    spectbl <- map_df(spec, extractSpectra) %>% 
      mutate(type = as.factor(type), id = as.factor(id))
    
    # extract specified features
    cat(paste0("Creating features for ", suffix2, "...\n"))
    feattbl <- spectbl %>%
      group_by(id) %>%
      do(extractFeatures(., featureVec = trainList$features, tol = mzTol))
    
    saveRDS(feattbl, file = paste0("../temp/featureInfo", suffix2, suffix))
    
    mlFeat <- feattbl %>%
      select(id, type, feat, relInt) %>%
      group_by(id) %>%
      mutate(relInt = relInt / max(relInt)) %>%
      ungroup() %>%
      spread(key = feat, value = relInt, fill = 0) %>%
      mutate(Ab = ifelse(str_detect(type, "i - r"), "pos", "other"),
             Kp = ifelse(str_detect(type, "e - r"), "pos", "other"),
             Ab = as.factor(ifelse(str_detect(type, "i - s"), "neg", Ab)),
             Kp = as.factor(ifelse(str_detect(type, "e - s"), "neg", Kp)))
    
    return(mlFeat)
}

trainDatList <- map(trainIdx, createFeatureTbl, 
                    mzTol = featTol,
                    testSet = F)

saveRDS(trainDatList, "../temp/trainDatList.rds")
```

Train the XGBoost Tree Models
-----------------------------

### Prepare the Data for Model Training

``` r
dmatList <- prepareData(trainDatList, orgLabs = c("Ab", "Kp"))
saveRDS(dmatList, "../temp/dmatList.rds")
```

### Tune Hyperparameter by Grid Search

``` r
paramGrid <- expand.grid(max_depth = seq(1, 10, by = 2),
                              min_child_weight = c(0, 1, 2, 3),
                              gamma = c(0, 0.01, 0.05, 0.1))

paramGrid$id <- 1:nrow(paramGrid)

xgbParams <- list(objective = "multi:softprob",
                  booster = "gbtree",
                  num_class = 3)

optimizeXgbParams <- function(params, data, eta, verbose = F) {
    tune <- xgb.cv(params = xgbParams,
                   data = data,
                   nrounds = 20000,
                   eta = eta,
                   nfold = 10,
                   early_stopping_rounds = 50,
                   max_depth = params$max_depth,
                   min_child_weight = params$min_child_weight,
                   gamma = params$gamma,
                   metrics = "mlogloss",
                   verbose = verbose)
    
    bst <- as.data.frame(tune$evaluation_log) %>%
        filter(test_mlogloss_mean == min(test_mlogloss_mean))
    
    return(cbind(params, bst))
}

paramRes <- map_df(dmatList, function(org) {
    paramGrid %>%
        group_by(id) %>%
        do(optimizeXgbParams(., org$dat, 0.3))
}, .id = "org")

bstParams <- paramRes %>%
    group_by(org) %>%
    filter(test_mlogloss_mean == min(test_mlogloss_mean))
```

### Refine Models

``` r
finalParams <- bstParams %>%
                group_by(org) %>%
                do(optimizeXgbParams(., dmatList[[.$org]]$dat, 0.01, verbose = T))
```

### Final Model

``` r
orgs <- c("Ab", "Kp")

mods <- map(orgs, function(org){
    cat(paste0("Training ", org, "\n"))
    params <- finalParams[finalParams$org == org, ]
    
    mod <- xgb.train(params = xgbParams,
              data = dmatList[[org]]$dat,
              nrounds = params$iter,
              max_depth = params$max_depth,
              min_child_weight = params$min_child_weight,
              gamma = params$gamma,
              eta = 0.01)
    xgb.save.raw(mod)
})
```

    ## Training Ab
    ## Training Kp

``` r
names(mods) <- orgs
saveRDS(mods, "../temp/finalModels.rds")
```

Prepare Test Set
----------------

``` r
testDatList <- map(trainIdx, createFeatureTbl, 
                   mzTol = featTol,
                   testSet = T)

saveRDS(testDatList, "../temp/testDatList.rds")
```

Session Information
-------------------

``` r
times <- toc()
cat(c("Executions time:", round((times$toc - times$tic)/60, 0), "min"))
session_info()
```

    ## 23074.7 sec elapsed
    ## Executions time: 385 min setting  value                       
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
