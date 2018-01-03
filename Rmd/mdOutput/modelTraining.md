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
source("../R/utilityFunctions.R")

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
           Ab = ifelse(str_detect(type, "Acinetobacter baumannii - res"), 
                       "pos", "other"),
           Kp = ifelse(str_detect(type, "Klebsiella pneumoniae - res"), 
                       "pos", "other"),
           Ab = as.factor(ifelse(str_detect(type, "Acinetobacter baumannii - sen"), 
                                 "neg", Ab)),
           Kp = as.factor(ifelse(str_detect(type, "Klebsiella pneumoniae - sen"), 
                                 "neg", Kp)),
           strain = str_replace(id, " *[0-9]{8}", ""),
           strain = str_replace(strain, "K_p_ ", "K pneumoniae -"),
           strain = str_replace(strain, "mut_", "mut"),
           strain = str_replace_all(strain, "[ _][0-9]([ _]|$)", ""),
           strain = str_replace(strain, "[0-9]{4}-[0-9]{4}", ""),
           strain = str_replace(strain, "\\([0-9]\\)", ""),
           strain = str_replace(strain, "[0-9]{2}_[0-9]{2}", ""),
           strain = str_replace(strain, "\\+.*", ""),
           strain = str_replace(strain, " (col med|lyo pellet).*$", ""),
           strain = str_replace(strain, "[_ ]+$", ""),
           strain = ifelse(Kp == "other" & Ab == "other", "Other Species", strain))

# specInfo %>%
#     filter(Ab != "other" | Kp != "other") %>%
#     group_by(type) %>%
#     summarize(n = length(unique(strain)))

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
    spec <- unlist(map(files[trainList$idx], preProcessSpec, hws = 80))

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
        #summarize(relInt = median(relInt)) %>%
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
    spec <- unlist(map(files[idx], preProcessSpec, hws = 80))
    #spec <- preProcessSpec(files[idx], hws = 80)
    #multiSpecIdx <- map_lgl(spec, ~ metaData(.)$num > 1)
    #spec <- spec[!multiSpecIdx]
    
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
    
    # Mergin to Table
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
    filter(test_mlogloss_mean == min(test_mlogloss_mean)) %>%
    filter(iter == min(iter))
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

    ## 3686.18 sec elapsed
    ## Executions time: 61 min setting  value                       
    ##  version  R version 3.4.0 (2017-04-21)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  tz       America/New_York            
    ##  date     2018-01-02                  
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
