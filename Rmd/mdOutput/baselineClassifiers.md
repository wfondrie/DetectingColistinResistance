Baseline Classifiers for Comparison
================
William E Fondrie

-   [Load libraries and prepare workspace](#load-libraries-and-prepare-workspace)
-   [Load Data](#load-data)
-   [Baseline Classifiers on Isolate Spectra](#baseline-classifiers-on-isolate-spectra)
    -   [Single Feature](#single-feature)
    -   [Nearest Neighbor](#nearest-neighbor)
-   [Baseline Classifiers on Simulated Spectra](#baseline-classifiers-on-simulated-spectra)
    -   [Single Feature](#single-feature-1)
    -   [Nearest Neighbor](#nearest-neighbor-1)
-   [Baseline Classifiers on in vitro Mixtures](#baseline-classifiers-on-in-vitro-mixtures)
    -   [Single feature](#single-feature-2)
    -   [Nearest Neighbor](#nearest-neighbor-2)

Load libraries and prepare workspace
------------------------------------

Load Data
---------

``` r
trainList <- readRDS("../temp/trainDatList.rds")
testList <- readRDS("../temp/testDatList.rds")
```

Baseline Classifiers on Isolate Spectra
---------------------------------------

### Single Feature

``` r
# select and format data -------------------------------------------------------
featAb <- tibble(col = names(testList$Ab)) %>%
    filter(str_detect(col, "mz")) %>%
    mutate(mz = as.numeric(str_match(col, "mz(.*)")[ , 2])) %>%
    filter(nom2feat(1912, mz) == mz |
               nom2feat(2036, mz) == mz)

singleAb <- testList$Ab %>%
    select_("type", "Ab", featAb$col[1], featAb$col[2]) #%>%

names(singleAb) <- c("type", "truth", "speciesVsOther", "pos")


featKp <- tibble(col = names(testList$Kp)) %>%
    filter(str_detect(col, "mz")) %>%
    mutate(mz = as.numeric(str_match(col, "mz(.*)")[ , 2])) %>%
    filter(nom2feat(1841, mz) == mz |
               nom2feat(1973, mz) == mz)

singleKp <- testList$Kp %>%
    select_("type", "Kp", featKp$col[1], featKp$col[2])

names(singleKp) <- c("type", "truth", "speciesVsOther", "pos")

# Calculate curves -------------------------------------------------------------
sfAbRoc <- makeRocCurves(singleAb)
sfAbPr <- makePrCurves(singleAb)

sfKpRoc <- makeRocCurves(singleKp)
sfKpPr <- makePrCurves(singleKp)

sfRes <- list(sfAbRoc = sfAbRoc,
              sfAbPr = sfAbPr,
              sfKpRoc = sfKpRoc,
              sfKpPr = sfKpPr)

saveRDS(sfRes, "../temp/sfRes.rds")


# Calulate curve statistics ----------------------------------------------------
sfAbStat <- singleAb %>%
    bootstrap(nboot) %>%
    do(calcAUC(.)) %>%
    mutate(org = "Ab")

sfKpStat <- singleKp %>%
    bootstrap(nboot) %>%
    do(calcAUC(.)) %>%
    mutate(org = "Kp")

curveStat <- sfAbStat %>%
    full_join(sfKpStat) %>% 
    group_by(org, type, level) %>%
    summarize(avg_AUC = mean(AUC),
              ci_low = quantile(AUC, 0.05 / 2),
              ci_high = quantile(AUC, 1 - 0.05 / 2))

saveRDS(curveStat, "../temp/sfCurveStat.rds")
```

### Nearest Neighbor

``` r
# Scales the rows of a matrix to unit vectors
scaleVectors <- function(m) {
    t(apply(m, 1, function(x) x/as.vector(sqrt(x %*% x))))
}

# picks the maximum column of a matrix for each row
pickMaxCol <- function(m, colname = "test") {
    dat <- apply(m, 1, function(x) max(x))
    df <- data.frame(id = names(dat),
                     x = dat)
    
    names(df)[2] <- colname
    
    return(df)
}

# The actual classifier

# for testing:
trainDat <- trainList
testDat <- testList
species <- "Kp"

nnClassifier <- function(trainDat, testDat, species = "Ab") {
    
    # format data
    trainDat <- trainDat[[species]]
    testDat <- testDat[[species]]
    
    names(trainDat)[names(trainDat) == species] <- "truth"
    names(testDat)[names(testDat) == species] <- "truth"
    
    # create matrices
    pos <- trainDat %>%
        filter(truth == "pos") %>%
        select(starts_with("mz")) %>%
        as.matrix(.) %>%
        scaleVectors(.)
    
    pos <- pos[!is.na(rowSums(pos)) != 0, ]
    
    neg <- trainDat %>%
        filter(truth == "neg") %>%
        select(starts_with("mz")) %>%
        as.matrix(.) %>%
        scaleVectors(.)
    
    neg <- neg[!is.na(rowSums(neg)), ]
    
    other <- trainDat %>%
        filter(truth == "other") %>%
        select(starts_with("mz")) %>%
        as.matrix(.) %>%
        scaleVectors(.)
    
    other <- other[!is.na(rowSums(other)), ]
    
    testMat <- testDat %>%
        select(starts_with("mz")) %>%
        as.matrix(.) %>%
        scaleVectors(.)
    
    row.names(testMat) <- testDat$id
    
    # calculate distances
    pos <- pickMaxCol(testMat %*% t(pos), "pos")
    pos$pos[is.na(pos$pos)] <- 0
    
    neg <- pickMaxCol(testMat %*% t(neg), "neg")
    neg$neg[is.na(neg$neg)] <- 0
    
    other <- pickMaxCol(testMat %*% t(other), "other")
    other$other[is.na(other$other)] <- 1
    
    # merge scores
    scores <- pos %>%
        full_join(neg) %>%
        full_join(other) %>%
        mutate(summed = neg + pos + other,
               pos = pos / summed,
               neg = neg / summed,
               other = other / summed,
               speciesVsOther = pos + neg)
    
    ret <- testDat %>%
        select(type, id, truth) %>%
        full_join(scores) %>%
        select(-summed)
    
   return(ret) 
}

# Classify test set
nnAb <- nnClassifier(trainList, testList, species = "Ab")
nnKp <- nnClassifier(trainList, testList, species = "Kp")

# ROC and PR curves
nnAbRoc <- makeRocCurves(nnAb)
nnAbPr <- makePrCurves(nnAb)

nnKpRoc <- makeRocCurves(nnKp)
nnKpPr <- makePrCurves(nnKp)

nnRes <- list(nnAbRoc = nnAbRoc,
              nnAbPr = nnAbPr,
              nnKpRoc = nnKpRoc,
              nnKpPr = nnKpPr)

saveRDS(nnRes, "../temp/nnRes.rds")
```

Baseline Classifiers on Simulated Spectra
-----------------------------------------

### Single Feature

### Nearest Neighbor

Baseline Classifiers on in vitro Mixtures
-----------------------------------------

### Single feature

### Nearest Neighbor
