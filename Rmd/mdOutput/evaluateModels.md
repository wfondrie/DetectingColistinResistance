Evaluate XGBoost Models
================
William E Fondrie

-   [Load Libraries & Prepare Workspace](#load-libraries-prepare-workspace)
-   [Evaluation of Model Performance on Isolate Spectra in Library](#evaluation-of-model-performance-on-isolate-spectra-in-library)
    -   [Import Acinetobacter baumannii](#import-acinetobacter-baumannii)
    -   [Import Klebsiella pneumoniae](#import-klebsiella-pneumoniae)
    -   [Plot PR and ROC Curves](#plot-pr-and-roc-curves)
    -   [Plot Curve Stats](#plot-curve-stats)
-   [Select 97% Sensitivity Threshold](#select-97-sensitivity-threshold)
    -   [Print FPRs](#print-fprs)
-   [Performance at 97% Sensitivity](#performance-at-97-sensitivity)
    -   [Performance on Isolate Spectra](#performance-on-isolate-spectra)
-   [Variable Importance](#variable-importance)
-   [Performance on Simulated Mixture Spectra](#performance-on-simulated-mixture-spectra)
    -   [Acinetobacter baumannii](#acinetobacter-baumannii)
    -   [Klebsiella pneumoniae](#klebsiella-pneumoniae)
    -   [Plot PR and ROC curves for Simulated Mixtures by Number of Spectra in Mixture](#plot-pr-and-roc-curves-for-simulated-mixtures-by-number-of-spectra-in-mixture)
    -   [Plot the AUC of PR and ROC Curves With Increasing Matrix Complexity](#plot-the-auc-of-pr-and-roc-curves-with-increasing-matrix-complexity)
    -   [Plot statistics from simulated mixtures](#plot-statistics-from-simulated-mixtures)
-   [Two-Species UTI Mixtures](#two-species-uti-mixtures)
    -   [Plot scores vs cutoff](#plot-scores-vs-cutoff)
-   [Session Info](#session-info)

Load Libraries & Prepare Workspace
----------------------------------

``` r
# data manipulation
suppressMessages(library(tidyverse, quietly = T))
library(stringr, quietly = T)
library(forcats, quietly = T)
library(broom, quietly = T)

# machine learning tools
library(PRROC, quietly = T)
library(caret, quietly = T)
library(xgboost, quietly = T)

# for session info
library(devtools, quietly = T)

# ggplot2 theme
source("../R/ggplotTheme.R")
theme_set(coolTheme)

# other functions
source("../R/prepareData.R")
source("../R/utilityFunctions.R")

# for reproducibility
set.seed(0847532)

# Bootstrap
nboot <- 2000 #2000

# for barplots
dodge <- position_dodge(width = 0.9)
```

Evaluation of Model Performance on Isolate Spectra in Library
-------------------------------------------------------------

``` r
testList <- readRDS("../temp/testDatList.rds")
testdmatList <- prepareData(testList, orgLabs = c("Ab", "Kp"))

mods <- readRDS("../temp/finalModels.rds")

formatResults <- function(model, dataList) {
    results <- as.data.frame(predict(model, dataList$dat, reshape = T))
    names(results) <- levels(unique(dataList$encoding)[[1]])
    
    results <- results %>%
        mutate(truth = dataList$encoding[[1]],
               negVsAll = pos + other,
               posVsAll = neg + other,
               speciesVsOther = pos + neg,
               pred = ifelse(neg > pos & neg > other, "neg", NA),
               pred = ifelse(pos > neg & pos > other, "pos", pred),
               pred = ifelse(other > neg & other > pos, "other", pred))
    
    return(results)
}
```

### Import Acinetobacter baumannii

``` r
abMod <- xgb.load(mods$Ab)
abRes <- formatResults(abMod, testdmatList$Ab)

# Ab PR curves
abPr <- makePrCurves(abRes)

# Ab ROC curves
abRoc <- makeRocCurves(abRes)

saveRDS(abRes, file = "../temp/abRes.rds")

# Ab curve statistics
abCurveStat <- abRes %>%
    bootstrap(nboot) %>%
    do(calcAUC(.)) %>%
    mutate(org = "Ab")
```

### Import Klebsiella pneumoniae

``` r
kpMod <- xgb.load(mods$Kp)
kpRes <- formatResults(kpMod, testdmatList$Kp)

# Kp PR curves
kpPr <- makePrCurves(kpRes)

# Kp ROC curves
kpRoc <- makeRocCurves(kpRes)

saveRDS(kpRes, file = "../temp/kpRes.rds")

# Kp curve statistics
kpCurveStat <- kpRes %>%
    bootstrap(nboot) %>%
    do(calcAUC(.)) %>%
    mutate(org = "Kp")


# Aggregate stats 
curveStat <- abCurveStat %>%
    full_join(kpCurveStat) %>% 
    group_by(org, type, level) %>%
    summarize(avg_AUC = mean(AUC),
              ci_low = quantile(AUC, 0.05 / 2),
              ci_high = quantile(AUC, 1 - 0.05 / 2)) %>%
    mutate(classifier = "XGBoost")

curveStat    
```

    ## Source: local data frame [8 x 7]
    ## Groups: org, type [4]
    ## 
    ## # A tibble: 8 x 7
    ##     org   type   level   avg_AUC    ci_low   ci_high classifier
    ##   <chr> <fctr>  <fctr>     <dbl>     <dbl>     <dbl>      <chr>
    ## 1    Ab     PR     Res 0.8040917 0.6831564 0.9069533    XGBoost
    ## 2    Ab     PR Species 0.9986373 0.9965608 0.9998725    XGBoost
    ## 3    Ab    ROC     Res 0.9742474 0.9536570 0.9905575    XGBoost
    ## 4    Ab    ROC Species 0.9994182 0.9985615 0.9999420    XGBoost
    ## 5    Kp     PR     Res 0.9827733 0.9631255 0.9967303    XGBoost
    ## 6    Kp     PR Species 0.9888387 0.9690220 0.9997858    XGBoost
    ## 7    Kp    ROC     Res 0.9979750 0.9958071 0.9995842    XGBoost
    ## 8    Kp    ROC Species 0.9985387 0.9962585 0.9999620    XGBoost

### Plot PR and ROC Curves

#### Plot function

``` r
# dat should be the ouput of PRROC::pr.curve(..., curve = T) or roc.curve(..., curve = T)
plotCurve <- function(dat, type = "ROC", title = NULL, mar = 0.5) {
    dat <- dat %>%
        mutate(classifier = fct_rev(as.factor(classifier)),
               Level = fct_rev(as.factor(Level)))
    
    # Main Plotting
    p <- ggplot(dat, aes(x = V1, y = V2, linetype = classifier, color = Level))
    
    if(type == "ROC") {
        p <-  p + geom_abline(intercept = 0, slope = 1) +
            xlab("1 - Specificity (FPR)") +
            ylab("Sensitivity (TPR)")
    } else {
        p <- p + xlab("Recall") +
            ylab("Precision")
    }
    
    # Common theming
    p <- p + geom_path(size = 0.5) +
        ylim(c(0, 1)) +
        xlim(c(0, 1)) +
        scale_color_discrete(name = "Level") +
        scale_linetype_discrete(name = "Classifier") +
        coord_equal()+
        labs(title = title) +
        theme(legend.key.height = unit(0.6, "lines"),
              legend.key.width = unit(0.6, "lines"),
              legend.margin = margin(l = -0.5, 
                                     r = -0.5,
                                     unit = "lines"),
              plot.title = element_text(margin = margin(b = mar, unit = "lines")))
    
    return(p)
}
```

#### Import baseline classifier results

``` r
sfRes <- readRDS("../temp/sfRes.rds")
sfCurveStat <- readRDS("../temp/sfCurveStat.rds") %>%
    mutate(classifier = "Single Feature")

sfCurveStat
```

    ## Source: local data frame [8 x 7]
    ## Groups: org, type [4]
    ## 
    ## # A tibble: 8 x 7
    ##     org   type   level   avg_AUC     ci_low   ci_high     classifier
    ##   <chr> <fctr>  <fctr>     <dbl>      <dbl>     <dbl>          <chr>
    ## 1    Ab     PR     Res 0.1377097 0.07224094 0.2263167 Single Feature
    ## 2    Ab     PR Species 0.8885745 0.85123769 0.9199739 Single Feature
    ## 3    Ab    ROC     Res 0.7287866 0.62848079 0.8205407 Single Feature
    ## 4    Ab    ROC Species 0.9421489 0.92337676 0.9585021 Single Feature
    ## 5    Kp     PR     Res 0.6093802 0.49565210 0.7254965 Single Feature
    ## 6    Kp     PR Species 0.9373398 0.89200400 0.9740877 Single Feature
    ## 7    Kp    ROC     Res 0.9534337 0.93627215 0.9687551 Single Feature
    ## 8    Kp    ROC Species 0.9901142 0.98401603 0.9953598 Single Feature

#### Function to Format ROC and PR results

``` r
combineResults <- function(mlResults, 
                           sfResults, 
                           type = "species") {
    ml <- as.tibble(mlResults[[type]]$curve) %>% mutate(classifier = "XGBoost")
    sf <- as.tibble(sfResults[[type]]$curve) %>% mutate(classifier = "Single Feature")
    
    r <-  ml %>% rbind(sf) %>% 
        mutate(type = type,
               Level = ifelse(type == "pos", "Resistance", "Species"))
    
    return(r)
}
```

#### Acinetobacter ROC & PR Curves

``` r
abTitle <- expression(italic("A. baumannii"))

abFullRoc <- combineResults(abRoc, sfRes$sfAbRoc, type = "species") %>%
    full_join(combineResults(abRoc, sfRes$sfAbRoc, type = "pos"))
plotCurve(abFullRoc, type = "ROC", title = abTitle)
```

``` r
ggsave("../results/AbRocCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

abFullPr <- combineResults(abPr, sfRes$sfAbPr, type = "species") %>%
    full_join(combineResults(abPr, sfRes$sfAbPr, type = "pos"))
plotCurve(abFullPr, type = "PR", title = abTitle)
```

``` r
ggsave("../results/AbPrCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)
```

#### Klebsiella ROC & PR Curves

``` r
kpTitle <- expression(italic("K. pneumoniae"))

kpFullRoc <- combineResults(kpRoc, sfRes$sfKpRoc, type = "species") %>%
    full_join(combineResults(kpRoc, sfRes$sfKpRoc, type = "pos"))
plotCurve(kpFullRoc, type = "ROC", title = kpTitle, mar = 0)
```

``` r
ggsave("../results/KpRocCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)


# Species PR
kpFullPr <- combineResults(kpPr, sfRes$sfKpPr, type = "species") %>%
    full_join(combineResults(kpPr, sfRes$sfKpPr, type = "pos"))
```

``` r
plotCurve(kpFullPr, type = "PR", title = kpTitle, mar = 0)
```

``` r
ggsave("../results/KpPrCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)
```

### Plot Curve Stats

``` r
curveStatPlot <- curveStat %>%
    full_join(sfCurveStat) %>%
    ungroup() %>%
    mutate(org = ifelse(org == "Ab", "A. baumannii", "K. pneumoniae"),
           level = ifelse(level == "Res", "Resistance", "Species"),
           classifier = str_replace(classifier, " ", "\n"))

curveStatPlot %>%
    filter(org == "A. baumannii")  %>%
    ggplot(aes(x = fct_rev(classifier), y = avg_AUC, fill = fct_rev(level),
               ymin = ci_low, ymax = ci_high)) +
    geom_col(color = "black", position = dodge) +
    geom_errorbar(position = dodge, width = 0.25) +
    facet_wrap(~fct_rev(type), nrow = 1) +
    ylab("Area Under the Curve") +
    xlab("Classifier") +
    scale_fill_discrete(name = "Level") +
    theme(legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(l = -0.5, 
                                 r = -0.2,
                                 unit = "lines"),
          plot.title = element_text(margin = margin(b = 0.5, unit = "lines"))) +
    ggtitle(expression(italic("A. baumannii")))
```

``` r
ggsave("../results/abCurveSummary.pdf", width = 105, height = 50, units = "mm", useDingbats = F)

curveStatPlot %>%
    filter(org == "K. pneumoniae")  %>%
    ggplot(aes(x = fct_rev(classifier), y = avg_AUC, fill = fct_rev(level),
               ymin = ci_low, ymax = ci_high)) +
    geom_col(color = "black", position = dodge) +
    geom_errorbar(position = dodge, width = 0.25) +
    facet_wrap(~fct_rev(type), nrow = 1) +
    ylab("Area Under the Curve") +
    xlab("Classifier") +
    scale_fill_discrete(name = "Level") +
    theme(legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(l = -0.5, 
                                 r = -0.2,
                                 unit = "lines"),
          plot.title = element_text(margin = margin(b = 0, unit = "lines"))) +
    ggtitle(expression(italic("K. pneumoniae")))
```

``` r
ggsave("../results/kpCurveSummary.pdf", width = 105, height = 50, units = "mm", useDingbats = F)
```

Select 97% Sensitivity Threshold
--------------------------------

``` r
pickThold <- function(curve, minSens) {
    dat <- as.data.frame(curve$curve)
    dat <- dat[dat[ , 2] >= minSens, ]
    dat <- dat[dat[ , 3] == max(dat[ , 3]), ]
    names(dat) <- c("FPR", "TPR", "threshold")
    
    return(dat)
} 

# For XGBoost Classifiers
thold <- tibble(type = c("Ab_species", "Ab_resistant","Kp_species", "Kp_resistant"),
                curves = list(abRoc$species, 
                              abRoc$pos,
                              kpRoc$species, 
                              kpRoc$pos)) %>%
    group_by(type) %>%
    do(pickThold(curve = .$curves[[1]], 0.97)) %>%
    mutate(classifier = "XGBoost")

abRes2 <- abRes %>% 
    mutate(posPred = ifelse(pos >= thold$threshold[thold$type == "Ab_resistant"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= thold$threshold[thold$type == "Ab_species"],
                              "Ab", "other"),
           specTruth = ifelse(truth != "other", "Ab", "other"))


kpRes2 <- kpRes %>%
    mutate(posPred = ifelse(pos >= thold$threshold[thold$type == "Kp_resistant"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= thold$threshold[thold$type == "Kp_species"],
                             "Kp", "other"),
           specTruth = ifelse(truth != "other", "Kp", "other"))


# For sf Classifiers
sfthold <- tibble(type = c("Ab_species", "Ab_resistant","Kp_species", "Kp_resistant"),
                  curves = list(sfRes$sfAbRoc$species, 
                                sfRes$sfAbRoc$pos,
                                sfRes$sfKpRoc$species, 
                                sfRes$sfKpRoc$pos)) %>%
    group_by(type) %>%
    do(pickThold(curve = .$curves[[1]], 0.97)) %>%
    mutate(classifier = "Single\nFeature")


# plot for 97% Sensitivity
ratePlot <- thold %>%
    full_join(sfthold) %>%
    mutate(FNR = 1 - TPR,
           TNR = 1 - FPR) %>%
    gather(metric, value, FNR, TNR, FPR, TPR) %>%
    mutate(metricType = ifelse(metric %in% c("TPR", "FNR"), "Postives", "Negatives"),
           fills = ifelse(str_detect(metric, "T"), "good", "bad"),
           org = ifelse(str_detect(type, "Ab_"), "A. baumannii", "K. pneumoniae"),
           Level = ifelse(str_detect(type, "_resistant"), "Resistance", "Species"),
           metric = parse_factor(metric, c("TPR", "FNR", "TNR", "FPR")))


ratePlot %>%
    filter(metric == "FPR") %>%
    ggplot(aes(x = fct_rev(classifier), y = value, fill = Level)) +
    geom_col(color = "black", position = "dodge") +
    facet_wrap(~org, nrow = 1) +
    xlab("Classifier") +
    ylab("False Positive Rate (FPR)") +
    theme(legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(l = -0.5, 
                                 r = -0.2,
                                 unit = "lines"),
          panel.spacing = unit(0, "lines"))
```

``` r
ggsave("../results/FPR.pdf", width = 70, height = 50, units = "mm", useDingbats = F)


ratePlot %>%
    filter(org == "A. baumannii",
           metric %in% c("TNR", "FPR")) %>%
    ggplot(aes(x = fct_rev(classifier), y = value, fill = fct_rev(metric))) + 
    geom_col(color = "black") +
    facet_wrap(~fct_rev(Level)) +
    scale_fill_discrete(name = "Metric") +
    theme(legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(l = -0.5, 
                                 r = -0.2,
                                 unit = "lines"),
          axis.title.y = element_blank(),
          plot.title = element_text(margin = margin(b = 0.5, unit = "lines"))) +
    ggtitle(expression(italic("A. baumannii"))) +
    xlab("Classifier")
```

``` r
ggsave("../results/abFPR.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

ratePlot %>%
    filter(org == "K. pneumoniae",
           metric %in% c("TNR", "FPR")) %>%
    ggplot(aes(x = fct_rev(classifier), y = value, fill = fct_rev(metric))) + 
    geom_col(color = "black") +
    facet_wrap(~fct_rev(Level)) +
    scale_fill_discrete(name = "Metric") +
    theme(legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(l = -0.5, 
                                 r = -0.2,
                                 unit = "lines"),
          axis.title.y = element_blank(),
          plot.title = element_text(margin = margin(b = 0, unit = "lines"))) +
    ggtitle(expression(italic("K. pneumoniae"))) +
    xlab("Classifier")
```

``` r
ggsave("../results/kpFPR.pdf", width = 70, height = 50, units = "mm", useDingbats = F)
```

### Print FPRs

``` r
ratePlot %>%
    filter(metric == "FPR")
```

    ## Source: local data frame [8 x 9]
    ## Groups: type [4]
    ## 
    ## # A tibble: 8 x 9
    ##           type   threshold        classifier metric       value
    ##          <chr>       <dbl>             <chr> <fctr>       <dbl>
    ## 1 Ab_resistant 0.001703996           XGBoost    FPR 0.163291139
    ## 2   Ab_species 0.201784899           XGBoost    FPR 0.007042254
    ## 3 Kp_resistant 0.150716931           XGBoost    FPR 0.012178620
    ## 4   Kp_species 0.414723821           XGBoost    FPR 0.002857143
    ## 5 Ab_resistant 0.001792086 "Single\nFeature"    FPR 0.930379747
    ## 6   Ab_species 0.037191286 "Single\nFeature"    FPR 0.385563380
    ## 7 Kp_resistant 0.048104094 "Single\nFeature"    FPR 0.263870095
    ## 8   Kp_species 0.281802944 "Single\nFeature"    FPR 0.100000000
    ## # ... with 4 more variables: metricType <chr>, fills <chr>, org <chr>,
    ## #   Level <chr>

Performance at 97% Sensitivity
------------------------------

``` r
calcIsoStats <- function(resDf, species, level) {
    
    if(level == "species") {
        cm <- confusionMatrix(data = resDf$specPred, 
                              reference = resDf$specTruth, 
                              positive = species, mode = "everything")
        
        
    } else if(level == "res") {
        cm <- confusionMatrix(data = resDf$posPred, 
                              reference = resDf$posTruth, 
                              positive = "pos", mode = "everything")
    }
    
    ret <- as.tibble(as.list(c(cm$overall, cm$byClass)))
    ret$orgResults <- species
    ret$level <- level
    
    ret <- select(ret, level, orgResults, Accuracy,
                  Sensitivity, Specificity, Precision)
    
    return(ret)
}

abSpecStats <- abRes2 %>%
    bootstrap(nboot) %>%
    do(calcIsoStats(., "Ab", "species"))

abResStats <- abRes2 %>%
    bootstrap(nboot) %>%
    do(calcIsoStats(., "Ab", "res"))

kpSpecStats <- kpRes2 %>%
    bootstrap(nboot) %>%
    do(calcIsoStats(., "Kp", "species"))

kpResStats <- kpRes2 %>%
    bootstrap(nboot) %>%
    do(calcIsoStats(., "Kp", "res"))

statOrder <- c("Accuracy", "Sensitivity", "Specificity", "Precision")

isoStats <- abSpecStats %>%
    full_join(abResStats) %>%
    full_join(kpSpecStats) %>%
    full_join(kpResStats) %>%
    gather(statistic, value, -level, -orgResults, -replicate) %>%
    group_by(level, orgResults, statistic) %>%
    summarize(avg = mean(value),
              ci_low = quantile(value, 0.05 / 2),
              ci_high = quantile(value, 1 - 0.05 / 2)) %>%
    ungroup() %>%
    mutate(statistic = fct_relevel(statistic, statOrder),
           species = str_replace(orgResults, "Ab", "A. baumannii"),
           species = str_replace(species, "Kp", "K. pneumoniae"),
           level = str_replace(level, "res", "Resistance"),
           level = str_replace(level, "species", "Species"),
           level = fct_rev(as.factor(level)))

isoStats %>%
    ggplot(aes(x = statistic, y = avg, 
               ymax = ci_high, ymin = ci_low, fill = level)) +
    geom_col(color = "black", position = dodge) +
    geom_errorbar(position = dodge, width = 0.25) +
    scale_fill_discrete(name = "Level") +
    facet_wrap(~species, ncol = 1) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(l = -0.5, 
                                 r = -0.2,
                                 unit = "lines"))
```

``` r
ggsave("../results/isoStats.pdf", width = 140, height = 100, units = "mm", useDingbats = F)
```

### Performance on Isolate Spectra

``` r
isoStats
```

    ## # A tibble: 16 x 7
    ##         level orgResults   statistic       avg    ci_low   ci_high
    ##        <fctr>      <chr>      <fctr>     <dbl>     <dbl>     <dbl>
    ##  1 Resistance         Ab    Accuracy 0.8412028 0.8171913 0.8656174
    ##  2 Resistance         Ab   Precision 0.2105517 0.1491622 0.2743902
    ##  3 Resistance         Ab Sensitivity 0.9721873 0.9090909 1.0000000
    ##  4 Resistance         Ab Specificity 0.8352831 0.8098147 0.8602052
    ##  5 Resistance         Kp    Accuracy 0.9878372 0.9794189 0.9951574
    ##  6 Resistance         Kp   Precision 0.9049463 0.8378294 0.9592041
    ##  7 Resistance         Kp Sensitivity 0.9880953 0.9611619 1.0000000
    ##  8 Resistance         Kp Specificity 0.9878093 0.9793360 0.9946452
    ##  9    Species         Ab    Accuracy 0.9877924 0.9794189 0.9939467
    ## 10    Species         Ab   Precision 0.9843691 0.9665388 0.9962121
    ## 11    Species         Ab Sensitivity 0.9764080 0.9563449 0.9925387
    ## 12    Species         Ab Specificity 0.9929499 0.9852546 0.9982729
    ## 13    Species         Kp    Accuracy 0.9927579 0.9866828 0.9975787
    ## 14    Species         Kp   Precision 0.9757616 0.9435484 1.0000000
    ## 15    Species         Kp Sensitivity 0.9767331 0.9469027 1.0000000
    ## 16    Species         Kp Specificity 0.9956365 0.9900143 1.0000000
    ## # ... with 1 more variables: species <chr>

Variable Importance
-------------------

``` r
vImpAb <- xgb.importance(colnames(testdmatList$Ab$dat), model = abMod) %>%
    mutate(org = "Ab")

vImpKp <- xgb.importance(colnames(testdmatList$Kp$dat), model = kpMod) %>%
    mutate(org = "Kp")

abFeats <- str_replace(vImpAb$Feature, "mz", "")
kpFeats <- str_replace(vImpKp$Feature, "mz", "")

basePeaks <- c(nom2feat(1910, abFeats), # Ab
               nom2feat(1894, abFeats),
               nom2feat(1882, abFeats),
               nom2feat(1840, kpFeats), # Kp
               nom2feat(1824, kpFeats),
               nom2feat(2063, kpFeats),
               nom2feat(2079, kpFeats))

resPeaks <- c(nom2feat(2033, abFeats), # Ab
              nom2feat(2005, abFeats),
              nom2feat(1955, kpFeats), # Kp
              nom2feat(1971, kpFeats))#,
              #nom2feat(1895, kpFeats))

vImp <- vImpAb %>%
    full_join(vImpKp) %>%
    mutate(featureMz = as.factor(str_replace(Feature, "mz", "")),
           featureMz = fct_reorder(featureMz, Gain, .desc = T),
           color = ifelse(featureMz %in% basePeaks, "Published Base Lipid A Structure", NA),
           color = ifelse(featureMz %in% resPeaks, 
                          "Published Resistance Lipid A Structure", color)) %>%
    group_by(org) %>%
    arrange(desc(Gain)) %>%
    filter(Gain >= Gain[20])

# Acinetobacter
vImp %>%
    filter(org == "Ab") %>%
    ggplot(aes(x = featureMz, y = Gain, ymax = Gain, ymin = 0)) +
    geom_linerange(size = 0.75) +
    geom_point(size = 2, aes(color = color)) +
    xlab(expression(paste("Feature (", italic("m/z"), ")"))) +
    scale_color_discrete(breaks = c("Published Base Lipid A Structure", 
                                    "Published Resistance Lipid A Structure")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = c(0.996, 0.992),
          legend.justification = c(1, 1),
          legend.title = element_blank(),
          legend.background = element_rect(color = "black"),
          legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(t = 0, b = 0.2, l = 0.2, r = 0.2,  unit = "lines"),
          plot.title = element_text(margin = margin(b = 0.5, unit = "lines"))) +
    ggtitle(expression(italic("A. baumannii")))
```

``` r
ggsave("../results/vImpAbPlot.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

saveRDS(vImpAb, "../temp/vImpAb.rds")
```

``` r
# Klebsiella
vImp %>%
    filter(org == "Kp") %>%
    ggplot(aes(x = featureMz, y = Gain, ymax = Gain, ymin = 0)) +
    geom_linerange(size = 0.75) +
    geom_point(size = 2, aes(color = color)) +
    xlab(expression(paste("Feature (", italic("m/z"), ")"))) +
    scale_color_discrete(breaks = c("Published Base Lipid A Structure", 
                                    "Published Resistance Lipid A Structure")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = c(0.996, 0.992),
          legend.justification = c(1, 1),
          legend.title = element_blank(),
          legend.background = element_rect(color = "black"),
          legend.key.size = unit(0.5, "lines"),
          legend.margin = margin(t = 0, b = 0.2, l = 0.2, r = 0.2,  unit = "lines"),
          plot.title = element_text(margin = margin(b = 0, unit = "lines"))) +
    ggtitle(expression(italic("K. pneumoniae")))
```

``` r
ggsave("../results/vImpKpPlot.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

saveRDS(vImpKp, "../temp/vImpKp.rds")
```

Performance on Simulated Mixture Spectra
----------------------------------------

``` r
# Import simulated mixture spectra features and information
mixedSumm <- readRDS("../temp/complexSpectraSummary.rds")
mixedComponents <- readRDS("../temp/complexComponents.rds")

mixedList <- readRDS("../temp/mixtureDatList.rds")

mixedDmatList <- prepareData(mixedList, orgLabs = c("Ab", "Kp"))
```

### Acinetobacter baumannii

``` r
mixedAbRes <- formatResults(abMod, mixedDmatList$Ab) %>%
    cbind(., mixedList$Ab) %>%
    select(-starts_with("mz")) %>%
    mutate(orgResults = "A. baumannii") %>%
    as_tibble(.)

mixedAbPRall <- mixedAbRes %>%
    group_by(orgResults) %>%
    do(prPos = pr.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       prSpecies = pr.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T),
       rocPos = roc.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       rocSpecies = roc.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T)) %>%
    mutate(prPosAUC = prPos$auc.integral,
           prSpeciesAUC = prSpecies$auc.integral,
           rocPosAUC = rocPos$auc,
           rocSpeciesAUC = rocSpecies$auc,
           n = "Overall")

mixedAbPRByn <- mixedAbRes %>%
    group_by(orgResults, n) %>%
    do(prPos = pr.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       prSpecies = pr.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T),
       rocPos = roc.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       rocSpecies = roc.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T)) %>%
    mutate(prPosAUC = prPos$auc.integral,
           prSpeciesAUC = prSpecies$auc.integral,
           rocPosAUC = rocPos$auc,
           rocSpeciesAUC = rocSpecies$auc)

mixedAbCurves <- mixedAbPRall %>% rbind(mixedAbPRByn)

saveRDS(mixedAbRes, file = "../temp/mixedAbRes.rds")
```

### Klebsiella pneumoniae

``` r
mixedKpRes <- formatResults(kpMod, mixedDmatList$Kp) %>%
    cbind(., mixedList$Kp) %>%
    select(-starts_with("mz"))%>%
    mutate(orgResults = "K. pneumoniae") %>%
    as_tibble(.)

mixedKpAll <- mixedKpRes %>%
    group_by(orgResults) %>%
    do(prPos = pr.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       prSpecies = pr.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T),
       rocPos = roc.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       rocSpecies = roc.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T)) %>%
    mutate(prPosAUC = prPos$auc.integral,
           prSpeciesAUC = prSpecies$auc.integral,
           rocPosAUC = rocPos$auc,
           rocSpeciesAUC = rocSpecies$auc,
           n = "Overall")

mixedKpByn <- mixedKpRes %>%
    group_by(orgResults, n) %>%
    do(prPos = pr.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       prSpecies = pr.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T),
       rocPos = roc.curve(.$pos[.$truth == "pos"], .$pos[.$truth != "pos"], curve = T),
       rocSpecies = roc.curve(.$speciesVsOther[.$truth != "other"], 
                            .$speciesVsOther[.$truth ==  "other"], curve = T)) %>%
    mutate(prPosAUC = prPos$auc.integral,
           prSpeciesAUC = prSpecies$auc.integral,
           rocPosAUC = rocPos$auc,
           rocSpeciesAUC = rocSpecies$auc)

mixedKpCurves <- mixedKpAll %>% rbind(mixedKpByn)

# All PR and ROC results
mixedCurves <- mixedKpCurves %>% rbind(mixedAbCurves)

saveRDS(mixedKpRes, file = "../temp/mixedKpRes.rds")
```

### Plot PR and ROC curves for Simulated Mixtures by Number of Spectra in Mixture

``` r
mixedCurveDat <- mixedCurves %>%
    select(-ends_with("AUC")) %>%
    gather(curveType, curveDat, -n, -orgResults) %>%
    group_by(orgResults, n, curveType) %>% 
    do(curveDf = as_tibble(.$curveDat[[1]]$curve)) %>%
    ungroup() %>%
    unnest(curveDf) %>%
    rename("Rec_Spec" = V1, "Prec_Sens" = V2, "Threshold" = V3) %>%
    mutate(PrOrRoc = ifelse(str_detect(curveType, "pr"), "PR", "ROC"),
           type = ifelse(str_detect(curveType, "Pos"), "Resistance", "Species"), #type),
           type = as.factor(type))

mixedCurveDat %>%
    filter(PrOrRoc == "PR") %>%
    ggplot(aes(x = Rec_Spec, y = Prec_Sens, color = n)) +
    geom_path() +
    coord_equal() +
    facet_grid(orgResults ~ type) +
    xlab("Recall") +
    ylab("Precision") +
    ylim(0, 1) +
    xlim(0, 1) +
    scale_color_discrete(name = "Number of Species") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

``` r
ggsave("../results/mixturePRCurves.pdf", width = 70, height = 80, units = "mm", useDingbats = F)


mixedCurveDat %>%
    filter(PrOrRoc == "ROC") %>%
    ggplot(aes(x = Rec_Spec, y = Prec_Sens, color = n)) +
    geom_abline(xintercept = 0, slope = 1) +
    geom_path() +
    coord_equal() +
    facet_grid(orgResults ~ type) +
    xlab("1 - Specificity (FPR)") +
    ylab("Sensitivity (TPR)") +
    ylim(0, 1) +
    xlim(0, 1) +
    scale_color_discrete(name = "Number of Species") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

``` r
ggsave("../results/mixtureROCCurves.pdf", width = 70, height = 80, units = "mm", useDingbats = F)
```

### Plot the AUC of PR and ROC Curves With Increasing Matrix Complexity

``` r
mixedCurveStats <- mixedKpRes %>%
    full_join(mixedAbRes) %>%
    nest(-n, -orgResults) %>%
    mutate(AUC = map(data, ~ bootstrap(., nboot) %>% do(calcAUC(.)))) %>%
    unnest(AUC) %>%
    group_by(n, orgResults, type, level) %>%
    summarize(avg_AUC = mean(AUC),
              ci_low = quantile(AUC, 0.05 / 2),
              ci_high = quantile(AUC, 1 - 0.05 / 2)) %>%
    mutate(level = str_replace(level, "Res", "Resistance"))
    


labOrder <- c("Species", 
              "Resistance")

mixedCurveStats %>%
    ggplot(aes(x = as.numeric(n), y = avg_AUC, color = fct_rev(level))) +
    geom_line() +
    geom_point() + 
    geom_errorbar(aes(ymax = ci_high, ymin = ci_low), width = 0.25) +
    xlab("Number of Species in Mixture") +
    ylab("Area Under the Curve") +
    ylim(c(0, 1)) +
    coord_fixed(ratio = 4) +
    scale_color_discrete(name = "Level") +
    facet_grid(fct_rev(type) ~ orgResults) +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom")
```

``` r
ggsave("../results/mixedCurvesAUCPlot.pdf", width = 70, height = 80, units = "mm", useDingbats = F)
```

``` r
mixedThold <- mixedCurveDat %>%
    filter(PrOrRoc == "ROC",
           n == "Overall") %>%
    group_by(orgResults, type) %>%
    summarise(threshold = max(Threshold[Prec_Sens >= 0.97])) %>%
    mutate(type2 = tolower(str_replace(type, "Colistin-", "")),
           species = str_replace(orgResults, "\\. ", ""),
           species = str_match(species, "^.."),
           type = paste0(species, "_", type2)) %>%
    select(-species, -type2)

mixedThold # 97% Sensitivity Thresholds for Simulated Mixtures
```

    ## Source: local data frame [4 x 3]
    ## Groups: orgResults [2]
    ## 
    ## # A tibble: 4 x 3
    ##      orgResults          type    threshold
    ##           <chr>         <chr>        <dbl>
    ## 1  A. baumannii Ab_resistance 0.0002114323
    ## 2  A. baumannii    Ab_species 0.0007294051
    ## 3 K. pneumoniae Kp_resistance 0.0033867187
    ## 4 K. pneumoniae    Kp_species 0.0061715187

``` r
mixedAbRes <- mixedAbRes %>%
    mutate(posPred = ifelse(pos >= mixedThold$threshold[mixedThold$type == "Ab_resistance"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= mixedThold$threshold[mixedThold$type == "Ab_species"],
                             "Ab", "other"),
           specTruth = ifelse(truth != "other", "Ab", "other"))


mixedKpRes <- mixedKpRes %>%
    mutate(posPred = ifelse(pos >= mixedThold$threshold[mixedThold$type == "Kp_resistance"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= mixedThold$threshold[mixedThold$type == "Kp_species"],
                             "Kp", "other"),
           specTruth = ifelse(truth != "other", "Kp", "other"))

# Function to calculate Species Stats
specStats <- function(resDf, posClass) {
    d1 <- map_df(1:5, function(n) {
        cm <- confusionMatrix(data = resDf$specPred[resDf$n == n], 
                              reference = resDf$specTruth[resDf$n == n], 
                              positive = posClass, mode = "everything")
        
        data.frame(as.list(c(cm$overall, cm$byClass)))
    })
    
    cm2 <- confusionMatrix(data = resDf$specPred, 
                           reference = resDf$specTruth, 
                           positive = posClass, mode = "everything")
    
    ret <- rbind(data.frame(as.list(c(cm2$overall, cm2$byClass))), d1)
    ret$n <- c("Overall", 1, 2, 3, 4, 5)
    ret$orgResults <- resDf$orgResults[1]
    
    ret <- select(ret, orgResults, n, Accuracy,
                  Sensitivity, Specificity, Precision)

    
    return(ret)
}

# Species Stats
mixedAbSpecStats <- mixedAbRes %>%
    bootstrap(nboot) %>%
    do(specStats(., "Ab"))

mixedKpSpecStats <- mixedKpRes %>%
    bootstrap(nboot) %>%
    do(specStats(., "Kp"))

mixedSpecStats <- mixedAbSpecStats %>%
    full_join(mixedKpSpecStats) %>%
    gather(statistic, value, -orgResults, -n, -replicate) %>%
    group_by(orgResults, n, statistic) %>%
    summarize(avg = mean(value),
              ci_low = quantile(value, 0.05 / 2),
              ci_high = quantile(value, 1 - 0.05 / 2))


mixedSpecStats
```

    ## Source: local data frame [48 x 6]
    ## Groups: orgResults, n [?]
    ## 
    ## # A tibble: 48 x 6
    ##      orgResults     n   statistic       avg    ci_low   ci_high
    ##           <chr> <chr>       <chr>     <dbl>     <dbl>     <dbl>
    ##  1 A. baumannii     1    Accuracy 0.6387172 0.6093237 0.6686752
    ##  2 A. baumannii     1   Precision 0.4564523 0.4182856 0.4932799
    ##  3 A. baumannii     1 Sensitivity 1.0000000 1.0000000 1.0000000
    ##  4 A. baumannii     1 Specificity 0.4813829 0.4450950 0.5195947
    ##  5 A. baumannii     2    Accuracy 0.6595055 0.6289980 0.6885661
    ##  6 A. baumannii     2   Precision 0.4826549 0.4444401 0.5214545
    ##  7 A. baumannii     2 Sensitivity 1.0000000 1.0000000 1.0000000
    ##  8 A. baumannii     2 Specificity 0.5009819 0.4630402 0.5383503
    ##  9 A. baumannii     3    Accuracy 0.7112633 0.6830169 0.7384189
    ## 10 A. baumannii     3   Precision 0.5079206 0.4659656 0.5485759
    ## # ... with 38 more rows

``` r
# Function to calculate Resistance Stats
resStats <- function(resDf, posClass = "pos") {
    d1 <- map_df(1:5, function(n) {
        cm <- confusionMatrix(data = resDf$posPred[resDf$n == n], 
                              reference = resDf$posTruth[resDf$n == n], 
                              positive = posClass, mode = "everything")
        
        data.frame(as.list(c(cm$overall, cm$byClass)))
    })
    
    cm2 <- confusionMatrix(data = resDf$posPred, 
                           reference = resDf$posTruth, 
                           positive = posClass, mode = "everything")
    
    ret <- rbind(data.frame(as.list(c(cm2$overall, cm2$byClass))), d1)
    ret$n <- c("Overall", 1, 2, 3, 4, 5)
    ret$orgResults <- resDf$orgResults[1]
    
    ret <- select(ret, orgResults, n, Accuracy,
                  Sensitivity, Specificity, Precision)
    
    return(ret)
}

# Resistance Stats
mixedAbResStats <- mixedAbRes %>%
    bootstrap(nboot) %>%
    do(resStats(.))

mixedKpResStats <- mixedKpRes %>%
    bootstrap(nboot) %>%
    do(resStats(.))

mixedResStats <- mixedAbResStats %>%
    full_join(mixedKpResStats) %>%
    gather(statistic, value, -orgResults, -n, -replicate) %>%
    group_by(orgResults, n, statistic) %>%
    summarize(avg = mean(value),
              ci_low = quantile(value, 0.05 / 2),
              ci_high = quantile(value, 1 - 0.05 / 2))

mixedResStats
```

    ## Source: local data frame [48 x 6]
    ## Groups: orgResults, n [?]
    ## 
    ## # A tibble: 48 x 6
    ##      orgResults     n   statistic       avg    ci_low   ci_high
    ##           <chr> <chr>       <chr>     <dbl>     <dbl>     <dbl>
    ##  1 A. baumannii     1    Accuracy 0.4623663 0.4332521 0.4921818
    ##  2 A. baumannii     1   Precision 0.2313406 0.1994379 0.2621399
    ##  3 A. baumannii     1 Sensitivity 1.0000000 1.0000000 1.0000000
    ##  4 A. baumannii     1 Specificity 0.3585691 0.3271125 0.3911548
    ##  5 A. baumannii     2    Accuracy 0.4460136 0.4153962 0.4756292
    ##  6 A. baumannii     2   Precision 0.2196744 0.1902632 0.2500083
    ##  7 A. baumannii     2 Sensitivity 1.0000000 1.0000000 1.0000000
    ##  8 A. baumannii     2 Specificity 0.3436398 0.3124200 0.3747160
    ##  9 A. baumannii     3    Accuracy 0.4810933 0.4499987 0.5132591
    ## 10 A. baumannii     3   Precision 0.2214491 0.1899262 0.2527825
    ## # ... with 38 more rows

### Plot statistics from simulated mixtures

``` r
mixedSpecStats %>%
    ggplot(aes(x = n, y = avg, ymax = ci_high, ymin = ci_low)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(width = 0.25) +
    facet_grid(orgResults ~ statistic) +
    ylab("Value") +
    xlab("Number of Species in Mixture") +
    ggtitle("Species-Level Detection")
```

``` r
ggsave("../results/simSpecStats.pdf", width = 210, height = 80, units = "mm", useDingbats = F)

mixedResStats %>%
    ggplot(aes(x = n, y = avg, ymax = ci_high, ymin = ci_low)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(width = 0.25) +
    facet_grid(orgResults ~ statistic) +
    ylab("Value") +
    xlab("Number of Species in Mixture") +
    ggtitle("Colistin Resistance Detection")
```

``` r
ggsave("../results/simResStats.pdf", width = 210, height = 80, units = "mm", useDingbats = F)
```

Two-Species UTI Mixtures
------------------------

``` r
# Import simulated mixture spectra features and information
twoSpeciesSumm <- readRDS("../temp/twoSpeciesSpecInfo.rds")

twoSpeciesList <- readRDS("../temp/twoSpeciesDatList.rds")

twoSpeciesDmatList <- prepareData(twoSpeciesList, orgLabs = c("Ab", "Kp"))

# predict 
twoAbRes <- formatResults(abMod, twoSpeciesDmatList$Ab) %>%
    cbind(., twoSpeciesList$Ab) %>%
    select(-starts_with("mz"))%>%
    mutate(orgResults = "A. baumannii",
           percentTarget = ifelse(truth == "other", 0, 100 - percentEc),
           Colistin = ifelse(percentTarget > 0 & truth != "other", as.character(truth), NA),
           Colistin = str_replace(Colistin, "pos", "resistant"),
           Colistin = str_replace(Colistin, "neg", "susceptible")) %>%
    as_tibble(.)

twoKpRes <- formatResults(kpMod, twoSpeciesDmatList$Kp) %>%
    cbind(., twoSpeciesList$Kp) %>%
    select(-starts_with("mz"))%>%
    mutate(orgResults = "K. pneumoniae",
           percentTarget = ifelse(truth == "other", 0, 100 - percentEc),
           Colistin = ifelse(percentTarget > 0 & truth != "other", as.character(truth), NA),
           Colistin = str_replace(Colistin, "pos", "resistant"),
           Colistin = str_replace(Colistin, "neg", "susceptible")) %>%
    as_tibble(.)


twoRes <- twoAbRes %>% 
    full_join(twoKpRes) %>%
    mutate(Colistin = fct_rev(as.factor(Colistin))) %>%
    filter(!(truth == "other" & percentEc != 100)) %>%
    group_by(Colistin, percentTarget, orgResults) %>%
    mutate(avg_speciesVsOther = mean(log10(speciesVsOther)),
           avg_pos = mean(log10(pos)),
           ci_lo_speciesVsOther = t.test(log10(speciesVsOther))$conf.int[1],
           ci_hi_speciesVsOther = t.test(log10(speciesVsOther))$conf.int[2],
           ci_lo_pos = t.test(log10(pos))$conf.int[1],
           ci_hi_pos = t.test(log10(pos))$conf.int[2])

twoResSumm <- twoRes %>%
    group_by(Colistin, percentTarget, orgResults) %>%
    summarize(avg_speciesVsOther = mean(log10(speciesVsOther)),
              avg_pos = mean(log10(pos)),
              ci_lo_speciesVsOther = t.test(log10(speciesVsOther))$conf.int[1],
              ci_hi_speciesVsOther = t.test(log10(speciesVsOther))$conf.int[2],
              ci_lo_pos = t.test(log10(pos))$conf.int[1],
              ci_hi_pos = t.test(log10(pos))$conf.int[2])

saveRDS(twoAbRes, file = "../temp/twoAbRes.rds")
saveRDS(twoKpRes, file = "../temp/twoKpRes.rds")
```

### Plot scores vs cutoff

``` r
tholdSpecies <- mixedThold %>%
    filter(str_detect(type, "species"))

twoRes %>%
    ggplot(aes(x = as.factor(percentTarget), y = log10(speciesVsOther),
               fill = Colistin, group = Colistin)) +
    geom_errorbar(aes(ymax = ci_hi_speciesVsOther, 
                      ymin = ci_lo_speciesVsOther, 
                      color = Colistin), 
                  width = 0.25, position = dodge) +
    geom_crossbar(aes(ymax = avg_speciesVsOther, ymin = avg_speciesVsOther, 
                      y = avg_speciesVsOther, color = Colistin),
                  width = 0.5, fatten = 1, 
                  position = dodge, show.legend = F) +
    geom_point(width = 0.25, shape = 21, position = position_jitterdodge()) +
    facet_wrap(~orgResults, ncol = 2) +
    scale_fill_discrete(breaks = c("resistant", "susceptible")) +
    scale_color_discrete(guide = "none") +
    xlab("Percent Target Species by Volume") +
    ylab("Species log(Score)") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom")
```

``` r
ggsave("../results/twoSpeciesSpecRes.pdf", width = 105, height = 60, units = "mm", useDingbats = F)


tholdRes <- mixedThold %>%
    filter(str_detect(type, "res"))
```

``` r
twoRes %>%
    ggplot(aes(x = as.factor(percentTarget), y = log10(pos), 
               fill = Colistin, group = Colistin)) +
    geom_errorbar(aes(ymax = ci_hi_pos, ymin = ci_lo_pos, color = Colistin), 
                  width = 0.25, position = dodge) +
    geom_crossbar(aes(ymax = avg_pos, ymin = avg_pos, 
                      y = avg_pos, color = Colistin),
                  width = 0.5, fatten = 1, 
                  position = dodge, show.legend = F) +
    geom_point(width = 0.25, shape = 21, position = position_jitterdodge()) +
    facet_wrap(~orgResults, ncol = 2) +
    scale_fill_discrete(breaks = c("resistant", "susceptible")) +
    scale_color_discrete(guide = "none") +
    xlab("Percent Target Species by Volume") +
    ylab("Colistin Resistance log(Score)") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom")
```

``` r
    #ylim(c(0, 1))

ggsave("../results/twoSpeciesResRes.pdf", width = 105, height = 60, units = "mm", useDingbats = F)
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
    ##  date     2018-01-18                  
    ## 
    ##  package      * version date       source        
    ##  assertthat     0.2.0   2017-04-11 CRAN (R 3.3.3)
    ##  backports      1.1.0   2017-05-22 CRAN (R 3.4.0)
    ##  base         * 3.4.0   2017-04-21 local         
    ##  broom        * 0.4.2   2017-02-13 CRAN (R 3.4.0)
    ##  car            2.1-4   2016-12-02 CRAN (R 3.4.0)
    ##  caret        * 6.0-76  2017-04-18 CRAN (R 3.4.0)
    ##  cellranger     1.1.0   2016-07-27 CRAN (R 3.4.0)
    ##  class          7.3-14  2015-08-30 CRAN (R 3.4.0)
    ##  codetools      0.2-15  2016-10-05 CRAN (R 3.4.0)
    ##  colorspace     1.3-2   2016-12-14 CRAN (R 3.4.0)
    ##  compiler       3.4.0   2017-04-21 local         
    ##  data.table     1.10.4  2017-02-01 CRAN (R 3.4.0)
    ##  datasets     * 3.4.0   2017-04-21 local         
    ##  DBI            0.6-1   2017-04-01 CRAN (R 3.4.0)
    ##  devtools     * 1.13.2  2017-06-02 CRAN (R 3.4.0)
    ##  digest         0.6.12  2017-01-27 CRAN (R 3.4.0)
    ##  dplyr        * 0.5.0   2016-06-24 CRAN (R 3.4.0)
    ##  e1071          1.6-8   2017-02-02 CRAN (R 3.4.0)
    ##  evaluate       0.10    2016-10-11 CRAN (R 3.4.0)
    ##  forcats      * 0.2.0   2017-01-23 CRAN (R 3.4.0)
    ##  foreach        1.4.3   2015-10-13 CRAN (R 3.4.0)
    ##  foreign        0.8-68  2017-04-24 CRAN (R 3.4.0)
    ##  ggplot2      * 2.2.1   2016-12-30 CRAN (R 3.4.0)
    ##  graphics     * 3.4.0   2017-04-21 local         
    ##  grDevices    * 3.4.0   2017-04-21 local         
    ##  grid           3.4.0   2017-04-21 local         
    ##  gtable         0.2.0   2016-02-26 CRAN (R 3.4.0)
    ##  haven          1.0.0   2016-09-23 CRAN (R 3.4.0)
    ##  hms            0.3     2016-11-22 CRAN (R 3.4.0)
    ##  htmltools      0.3.6   2017-04-28 CRAN (R 3.4.0)
    ##  httr           1.2.1   2016-07-03 CRAN (R 3.4.0)
    ##  iterators      1.0.8   2015-10-13 CRAN (R 3.4.0)
    ##  jsonlite       1.4     2017-04-08 CRAN (R 3.4.0)
    ##  knitr          1.16    2017-05-18 CRAN (R 3.4.0)
    ##  labeling       0.3     2014-08-23 CRAN (R 3.4.0)
    ##  lattice      * 0.20-35 2017-03-25 CRAN (R 3.4.0)
    ##  lazyeval       0.2.0   2016-06-12 CRAN (R 3.4.0)
    ##  lme4           1.1-13  2017-04-19 CRAN (R 3.4.0)
    ##  lubridate      1.6.0   2016-09-13 CRAN (R 3.4.0)
    ##  magrittr       1.5     2014-11-22 CRAN (R 3.4.0)
    ##  MASS           7.3-47  2017-02-26 CRAN (R 3.4.0)
    ##  Matrix         1.2-10  2017-04-28 CRAN (R 3.4.0)
    ##  MatrixModels   0.4-1   2015-08-22 CRAN (R 3.4.0)
    ##  memoise        1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods      * 3.4.0   2017-04-21 local         
    ##  mgcv           1.8-17  2017-02-08 CRAN (R 3.4.0)
    ##  minqa          1.2.4   2014-10-09 CRAN (R 3.4.0)
    ##  mnormt         1.5-5   2016-10-15 CRAN (R 3.4.0)
    ##  ModelMetrics   1.1.0   2016-08-26 CRAN (R 3.4.0)
    ##  modelr         0.1.0   2016-08-31 CRAN (R 3.4.0)
    ##  munsell        0.4.3   2016-02-13 CRAN (R 3.4.0)
    ##  nlme           3.1-131 2017-02-06 CRAN (R 3.4.0)
    ##  nloptr         1.0.4   2014-08-04 CRAN (R 3.4.0)
    ##  nnet           7.3-12  2016-02-02 CRAN (R 3.4.0)
    ##  parallel       3.4.0   2017-04-21 local         
    ##  pbkrtest       0.4-7   2017-03-15 CRAN (R 3.4.0)
    ##  plyr           1.8.4   2016-06-08 CRAN (R 3.4.0)
    ##  PRROC        * 1.3     2017-04-21 CRAN (R 3.4.0)
    ##  psych          1.7.5   2017-05-03 CRAN (R 3.4.0)
    ##  purrr        * 0.2.2.2 2017-05-11 CRAN (R 3.4.0)
    ##  quantreg       5.33    2017-04-18 CRAN (R 3.4.0)
    ##  R6             2.2.1   2017-05-10 CRAN (R 3.4.0)
    ##  Rcpp           0.12.11 2017-05-22 CRAN (R 3.4.0)
    ##  readr        * 1.1.1   2017-05-16 CRAN (R 3.4.0)
    ##  readxl         1.0.0   2017-04-18 CRAN (R 3.4.0)
    ##  reshape2       1.4.2   2016-10-22 CRAN (R 3.4.0)
    ##  rlang          0.1.1   2017-05-18 CRAN (R 3.4.0)
    ##  rmarkdown    * 1.6     2017-06-15 CRAN (R 3.4.2)
    ##  rprojroot      1.2     2017-01-16 CRAN (R 3.4.0)
    ##  rvest          0.3.2   2016-06-17 CRAN (R 3.4.0)
    ##  scales         0.4.1   2016-11-09 CRAN (R 3.4.0)
    ##  SparseM        1.77    2017-04-23 CRAN (R 3.4.0)
    ##  splines        3.4.0   2017-04-21 local         
    ##  stats        * 3.4.0   2017-04-21 local         
    ##  stats4         3.4.0   2017-04-21 local         
    ##  stringi        1.1.5   2017-04-07 CRAN (R 3.4.0)
    ##  stringr      * 1.2.0   2017-02-18 CRAN (R 3.4.0)
    ##  tibble       * 1.3.1   2017-05-17 CRAN (R 3.4.0)
    ##  tictoc         1.0     2014-06-17 CRAN (R 3.4.0)
    ##  tidyr        * 0.6.3   2017-05-15 CRAN (R 3.4.0)
    ##  tidyverse    * 1.1.1   2017-01-27 CRAN (R 3.4.0)
    ##  tools          3.4.0   2017-04-21 local         
    ##  utils        * 3.4.0   2017-04-21 local         
    ##  withr          1.0.2   2016-06-20 CRAN (R 3.4.0)
    ##  xgboost      * 0.6-4   2017-01-05 CRAN (R 3.4.1)
    ##  xml2           1.1.1   2017-01-24 CRAN (R 3.4.0)
    ##  yaml           2.1.14  2016-11-12 CRAN (R 3.4.0)
