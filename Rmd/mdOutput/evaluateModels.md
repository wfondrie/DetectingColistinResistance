Evaluate XGBoost Models
================
William E Fondrie

-   [Load Libraries & Prepare Workspace](#load-libraries-prepare-workspace)
-   [Evaluation of Model Performance on Isolate Spectra in Library](#evaluation-of-model-performance-on-isolate-spectra-in-library)
    -   [Acinetobacter baumannii](#acinetobacter-baumannii)
    -   [Klebsiella pneumoniae](#klebsiella-pneumoniae)
    -   [Plot PR and ROC Curves](#plot-pr-and-roc-curves)
-   [New Confusion Matrices](#new-confusion-matrices)
-   [Variable Importance](#variable-importance)
-   [Performance on Simulated Mixture Spectra](#performance-on-simulated-mixture-spectra)
    -   [Acinetobacter baumannii](#acinetobacter-baumannii-1)
    -   [Klebsiella pneumoniae](#klebsiella-pneumoniae-1)
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

### Acinetobacter baumannii

``` r
abMod <- xgb.load(mods$Ab)
abRes <- formatResults(abMod, testdmatList$Ab)

# Ab PR curves
abPosPR <- pr.curve(abRes$pos[abRes$truth == "pos"], 
                    abRes$pos[abRes$truth != "pos"], curve = T)

abSpeciesPR <- pr.curve(abRes$speciesVsOther[abRes$truth != "other"], 
                        abRes$speciesVsOther[abRes$truth ==  "other"], curve = T)

# Ab ROC curves
abPosROC <- roc.curve(abRes$pos[abRes$truth == "pos"], 
                    abRes$pos[abRes$truth != "pos"], curve = T)

abSpeciesROC <- roc.curve(abRes$speciesVsOther[abRes$truth != "other"], 
                        abRes$speciesVsOther[abRes$truth ==  "other"], curve = T)

saveRDS(abRes, file = "../temp/abRes.rds")
```

### Klebsiella pneumoniae

``` r
kpMod <- xgb.load(mods$Kp)
kpRes <- formatResults(kpMod, testdmatList$Kp)

# Kp PR curves
kpPosPR <- pr.curve(kpRes$pos[kpRes$truth == "pos"], 
                    kpRes$pos[kpRes$truth != "pos"], curve = T)

kpSpeciesPR <- pr.curve(kpRes$speciesVsOther[kpRes$truth != "other"], 
                        kpRes$speciesVsOther[kpRes$truth ==  "other"], curve = T)

# Kp ROC curves
kpPosROC <- roc.curve(kpRes$pos[kpRes$truth == "pos"], 
                    kpRes$pos[kpRes$truth != "pos"], curve = T)

kpSpeciesROC <- roc.curve(kpRes$speciesVsOther[kpRes$truth != "other"], 
                        kpRes$speciesVsOther[kpRes$truth ==  "other"], curve = T)

saveRDS(kpRes, file = "../temp/kpRes.rds")
```

### Plot PR and ROC Curves

#### Acinetobacter PR & ROC Curves

``` r
# PR
prnames <- c("Recall", "Precision", "threshold")
prAb <- as.tibble(rbind(abPosPR$curve, abSpeciesPR$curve))
names(prAb) <- prnames

labsAb <- c(paste0("Resistant\nAUC = ", sigRound(abPosPR$auc.integral, 3)),
            paste0("Species\nAUC = ", sigRound(abSpeciesPR$auc.integral, 3)))

AbPrCurve <- prAb %>%
    mutate(cond = c(rep(labsAb[1], nrow(abPosPR$curve)),
                    rep(labsAb[2], nrow(abSpeciesPR$curve))),
           cond = fct_rev(fct_relevel(as.factor(cond), labsAb))) %>%
    ggplot(aes(x = Recall, y = Precision, color = cond)) +
    geom_path(size = 1)  +
    coord_equal() +
    ylim(c(0, 1)) +
    xlim(c(0, 1)) +
    theme(legend.title = element_blank(),
          legend.key.height = unit(1.2, "lines"),
          legend.key.width = unit(0.6, "lines"),
          legend.margin = margin(l = -0.5, unit = "lines"),
          plot.title = element_text(margin = margin(b = 0.5, unit = "lines"))) +
    ggtitle(expression(italic("A. baumannii")))

AbPrCurve
```

``` r
ggsave("../results/AbPrCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

# ROC
rocnames <- c("1 - Specificity (FPR)", "Sensitivity (TPR)", "threshold")
```

``` r
rocAb <- as.tibble(rbind(abPosROC$curve, 
                         abSpeciesROC$curve))
names(rocAb) <- rocnames

roclabsAb <- c(paste0("Resistant\nAUC = ", sigRound(abPosROC$auc, 3)),
               paste0("Species\nAUC = ", sigRound(abSpeciesROC$auc, 3)))


AbRocCurve <- rocAb %>%
    mutate(cond = c(rep(roclabsAb[1], nrow(abPosROC$curve)),
                    rep(roclabsAb[2], nrow(abSpeciesROC$curve))),
           cond = fct_rev(fct_relevel(as.factor(cond), roclabsAb))) %>%
    ggplot(aes(x = `1 - Specificity (FPR)`, y = `Sensitivity (TPR)`, color = cond)) +
    geom_abline(xintercept = 0, slope = 1) +
    geom_path(size = 1) +
    ylim(c(0, 1)) +
    xlim(c(0, 1)) +
    coord_equal() +
    theme(legend.title = element_blank(),
          legend.key.height = unit(1.2, "lines"),
          legend.key.width = unit(0.6, "lines"),
          legend.margin = margin(l = -0.5, unit = "lines"),
          plot.title = element_text(margin = margin(b = 0.5, unit = "lines"))) +
    ggtitle(expression(italic("A. baumannii")))

AbRocCurve
```

``` r
ggsave("../results/AbRocCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)
```

#### Klebsiella pneumoniae PR & ROC Curves

``` r
# PR
prnames <- c("Recall", "Precision", "threshold")
prKp <- as.tibble(rbind(kpPosPR$curve, 
                        kpSpeciesPR$curve))
names(prKp) <- prnames

labsKp <- c(paste0("Resistant\nAUC = ", sigRound(kpPosPR$auc.integral, 3)),
            paste0("Species\nAUC = ", sigRound(kpSpeciesPR$auc.integral, 3)))

prKp %>%
    mutate(cond = c(rep(labsKp[1], nrow(kpPosPR$curve)),
                    rep(labsKp[2], nrow(kpSpeciesPR$curve))),
           cond = fct_rev(fct_relevel(as.factor(cond), labsKp))) %>%
    ggplot(aes(x = Recall, y = Precision, color = cond)) +
    geom_path(size = 1)  +
    ylim(c(0, 1)) +
    xlim(c(0, 1)) +
    coord_equal() +
    theme(legend.title = element_blank(),
          legend.key.height = unit(1.2, "lines"),
          legend.key.width = unit(0.6, "lines"),
          legend.margin = margin(l = -0.5, unit = "lines"),
          plot.title = element_text(margin = margin(b = -0, unit = "lines"))) +
    ggtitle(expression(italic("K. pneumoniae")))
```

``` r
ggsave("../results/KpPrCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

# ROC
rocnames <- c("1 - Specificity (FPR)", "Sensitivity (TPR)", "threshold")
```

``` r
rocKp <- as.tibble(rbind(kpPosROC$curve, 
                         kpSpeciesROC$curve))
names(rocKp) <- rocnames

roclabsKp <- c(paste0("Resistant\nAUC = ", sigRound(kpPosROC$auc, 3)),
               paste0("Species\nAUC = ", sigRound(kpSpeciesROC$auc, 3)))

rocKp %>%
    mutate(cond = c(rep(roclabsKp[1], nrow(kpPosROC$curve)),
                    rep(roclabsKp[2], nrow(kpSpeciesROC$curve))),
           cond = fct_rev(fct_relevel(as.factor(cond), roclabsKp))) %>%
    ggplot(aes(x = `1 - Specificity (FPR)`, y = `Sensitivity (TPR)`, color = cond)) +
    geom_abline(xintercept = 0, slope = 1) +
    geom_path(size = 1) +
    ylim(c(0, 1)) +
    xlim(c(0, 1)) +
    coord_equal() + 
    theme(legend.title = element_blank(),
          legend.key.height = unit(1.2, "lines"),
          legend.key.width = unit(0.6, "lines"),
          legend.margin = margin(l = -0.5, unit = "lines"),
          plot.title = element_text(margin = margin(b = 0, unit = "lines"))) +
    ggtitle(expression(italic("K. pneumoniae")))
```

``` r
ggsave("../results/KpRocCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)
```

New Confusion Matrices
----------------------

``` r
pickThold <- function(curve, minSens) {
    dat <- as.data.frame(curve$curve)
    dat <- dat[dat[ , 2] >= minSens, ]
    dat <- dat[dat[ , 3] == max(dat[ , 3]), ]
    names(dat) <- c("FPR", "TPR", "threshold")
    
    return(dat)
} 

thold <- tibble(type = c("Ab_species", "Ab_resistant","Kp_species", "Kp_resistant"),
                curves = list(abSpeciesROC, abPosROC, kpSpeciesROC, kpPosROC)) %>%
    group_by(type) %>%
    do(pickThold(curve = .$curves[[1]], 0.97))

abRes2 <- abRes %>% 
    mutate(posPred = ifelse(pos >= thold$threshold[thold$type == "Ab_resistant"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= thold$threshold[thold$type == "Ab_species"],
                              "Ab", "other"),
           specTruth = ifelse(truth != "other", "Ab", "other"))

# Ab Species
confusionMatrix(data = abRes2$specPred, reference = abRes2$specTruth, 
                positive = "Ab", mode = "everything")
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  Ab other
    ##      Ab    252     0
    ##      other   6   580
    ##                                           
    ##                Accuracy : 0.9928          
    ##                  95% CI : (0.9845, 0.9974)
    ##     No Information Rate : 0.6921          
    ##     P-Value [Acc > NIR] : < 2e-16         
    ##                                           
    ##                   Kappa : 0.9831          
    ##  Mcnemar's Test P-Value : 0.04123         
    ##                                           
    ##             Sensitivity : 0.9767          
    ##             Specificity : 1.0000          
    ##          Pos Pred Value : 1.0000          
    ##          Neg Pred Value : 0.9898          
    ##               Precision : 1.0000          
    ##                  Recall : 0.9767          
    ##                      F1 : 0.9882          
    ##              Prevalence : 0.3079          
    ##          Detection Rate : 0.3007          
    ##    Detection Prevalence : 0.3007          
    ##       Balanced Accuracy : 0.9884          
    ##                                           
    ##        'Positive' Class : Ab              
    ## 

``` r
# Colistin-Resistant Ab
confusionMatrix(data = abRes2$posPred, reference = abRes2$posTruth, 
                positive = "pos", mode = "everything")
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction neg pos
    ##        neg 778   1
    ##        pos  24  35
    ##                                           
    ##                Accuracy : 0.9702          
    ##                  95% CI : (0.9563, 0.9806)
    ##     No Information Rate : 0.957           
    ##     P-Value [Acc > NIR] : 0.03169         
    ##                                           
    ##                   Kappa : 0.722           
    ##  Mcnemar's Test P-Value : 1.083e-05       
    ##                                           
    ##             Sensitivity : 0.97222         
    ##             Specificity : 0.97007         
    ##          Pos Pred Value : 0.59322         
    ##          Neg Pred Value : 0.99872         
    ##               Precision : 0.59322         
    ##                  Recall : 0.97222         
    ##                      F1 : 0.73684         
    ##              Prevalence : 0.04296         
    ##          Detection Rate : 0.04177         
    ##    Detection Prevalence : 0.07041         
    ##       Balanced Accuracy : 0.97115         
    ##                                           
    ##        'Positive' Class : pos             
    ## 

``` r
kpRes2 <- kpRes %>%
    mutate(posPred = ifelse(pos >= thold$threshold[thold$type == "Kp_resistant"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= thold$threshold[thold$type == "Kp_species"],
                             "Kp", "other"),
           specTruth = ifelse(truth != "other", "Kp", "other"))

# Kp Species
confusionMatrix(data = kpRes2$specPred, reference = kpRes2$specTruth, 
                positive = "Kp", mode = "everything")
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  Kp other
    ##      Kp    123     7
    ##      other   3   704
    ##                                           
    ##                Accuracy : 0.9881          
    ##                  95% CI : (0.9781, 0.9943)
    ##     No Information Rate : 0.8495          
    ##     P-Value [Acc > NIR] : <2e-16          
    ##                                           
    ##                   Kappa : 0.9539          
    ##  Mcnemar's Test P-Value : 0.3428          
    ##                                           
    ##             Sensitivity : 0.9762          
    ##             Specificity : 0.9902          
    ##          Pos Pred Value : 0.9462          
    ##          Neg Pred Value : 0.9958          
    ##               Precision : 0.9462          
    ##                  Recall : 0.9762          
    ##                      F1 : 0.9609          
    ##              Prevalence : 0.1505          
    ##          Detection Rate : 0.1470          
    ##    Detection Prevalence : 0.1553          
    ##       Balanced Accuracy : 0.9832          
    ##                                           
    ##        'Positive' Class : Kp              
    ## 

``` r
# Colistin-Resistant Kp
confusionMatrix(data = kpRes2$posPred, reference = kpRes2$posTruth, 
                positive = "pos", mode = "everything")
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction neg pos
    ##        neg 651   2
    ##        pos  99  85
    ##                                           
    ##                Accuracy : 0.8793          
    ##                  95% CI : (0.8553, 0.9006)
    ##     No Information Rate : 0.8961          
    ##     P-Value [Acc > NIR] : 0.9473          
    ##                                           
    ##                   Kappa : 0.5661          
    ##  Mcnemar's Test P-Value : <2e-16          
    ##                                           
    ##             Sensitivity : 0.9770          
    ##             Specificity : 0.8680          
    ##          Pos Pred Value : 0.4620          
    ##          Neg Pred Value : 0.9969          
    ##               Precision : 0.4620          
    ##                  Recall : 0.9770          
    ##                      F1 : 0.6273          
    ##              Prevalence : 0.1039          
    ##          Detection Rate : 0.1016          
    ##    Detection Prevalence : 0.2198          
    ##       Balanced Accuracy : 0.9225          
    ##                                           
    ##        'Positive' Class : pos             
    ## 

Variable Importance
-------------------

``` r
vImpAb <- xgb.importance(colnames(testdmatList$Ab$dat), model = abMod) %>%
    mutate(org = "Ab")

vImpKp <- xgb.importance(colnames(testdmatList$Kp$dat), model = kpMod) %>%
    mutate(org = "Kp")

basePeaks <- c("1912.5448", # Ab
               "1841.4497", # Kp
               "1826.3531",
               "2080.3158")
resPeaks <- c("2036.4927", # Ab
              "1957.0536",
              "1957.1088", # Kp
              "1973.1961",
              "1895.7716")

vImp <- vImpAb %>%
    full_join(vImpKp) %>%
    mutate(featureMz = as.factor(str_replace(Feature, "mz", "")),
           featureMz = fct_reorder(featureMz, Gain, .desc = T),
           color = ifelse(featureMz %in% basePeaks, "Base Lipid A Structure", NA),
           color = ifelse(featureMz %in% resPeaks, "Resistance Lipid A Structure", color)) %>%
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
    scale_color_discrete(breaks = c("Base Lipid A Structure", "Resistance Lipid A Structure")) +
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
    scale_color_discrete(breaks = c("Base Lipid A Structure", "Resistance Lipid A Structure")) +
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
           type = ifelse(str_detect(curveType, "Pos"), "Colistin-Resistant", "Species"), #type),
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
labOrder <- c("Species", 
              "Colistin-Resistant")

mixedCurvesAUC <- mixedCurves %>%
    select(ends_with("AUC"), n, orgResults) %>%
    gather(curveType, AUC, -n, -orgResults) %>%
    mutate(PrOrRoc = ifelse(str_detect(curveType, "pr"), "PR", "ROC"),
           type = ifelse(str_detect(curveType, "Neg"), "Colistin-Susceptible", "Species"),
           type = ifelse(str_detect(curveType, "Pos"), "Colistin-Resistant", type),
           type = fct_relevel(as.factor(type), labOrder))

mixedCurvesAUC %>%
    filter(n != "Otherwise") %>%
    ggplot(aes(x = as.numeric(n), y = AUC, color = type)) +
    geom_line() +
    geom_point() + 
    xlab("Number of Species in Mixture") +
    ylab("Area Under the Curve") +
    ylim(c(0, 1)) +
    coord_fixed(ratio = 4) +
    scale_color_discrete(name = NULL) +
    facet_grid(fct_rev(PrOrRoc) ~ orgResults) +
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

    ## # A tibble: 4 x 3
    ## # Groups:   orgResults [2]
    ##      orgResults         type    threshold
    ##           <chr>        <chr>        <dbl>
    ## 1  A. baumannii Ab_resistant 0.0004897423
    ## 2  A. baumannii   Ab_species 0.0016128461
    ## 3 K. pneumoniae Kp_resistant 0.0098284408
    ## 4 K. pneumoniae   Kp_species 0.0143815272

``` r
mixedAbRes <- mixedAbRes %>%
    mutate(posPred = ifelse(pos >= mixedThold$threshold[mixedThold$type == "Ab_resistant"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= mixedThold$threshold[mixedThold$type == "Ab_species"],
                             "Ab", "other"),
           specTruth = ifelse(truth != "other", "Ab", "other"))


mixedKpRes <- mixedKpRes %>%
    mutate(posPred = ifelse(pos >= mixedThold$threshold[mixedThold$type == "Kp_resistant"],
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
mixedSpecStats <- specStats(mixedAbRes, "Ab") %>%
    full_join(specStats(mixedKpRes, "Kp"))
mixedSpecStats
```

    ##       orgResults       n  Accuracy Sensitivity Specificity Precision
    ## 1   A. baumannii Overall 0.7198000   0.9705666   0.6262016 0.4921642
    ## 2   A. baumannii       1 0.6800000   1.0000000   0.5573997 0.4639866
    ## 3   A. baumannii       2 0.6580000   0.9925373   0.5355191 0.4389439
    ## 4   A. baumannii       3 0.7350000   0.9701987   0.6332378 0.5336976
    ## 5   A. baumannii       4 0.7480000   0.9444444   0.6753425 0.5182927
    ## 6   A. baumannii       5 0.7780000   0.9421488   0.7255937 0.5229358
    ## 7  K. pneumoniae Overall 0.8956133   0.9700000   0.8660613 0.7420765
    ## 8  K. pneumoniae       1 0.9498495   0.9909091   0.9295352 0.8743316
    ## 9  K. pneumoniae       2 0.9040816   0.9884615   0.8736111 0.7385057
    ## 10 K. pneumoniae       3 0.8886619   0.9776952   0.8549296 0.7185792
    ## 11 K. pneumoniae       4 0.8676171   0.9630996   0.8312236 0.6850394
    ## 12 K. pneumoniae       5 0.8671400   0.9259259   0.8449721 0.6925208

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
mixedResStats <- resStats(mixedAbRes) %>%
    full_join(resStats(mixedKpRes))

mixedResStats
```

    ##       orgResults       n  Accuracy Sensitivity Specificity Precision
    ## 1   A. baumannii Overall 0.5142000   0.9708029   0.4417149 0.2163305
    ## 2   A. baumannii       1 0.4350000   1.0000000   0.3430233 0.1985816
    ## 3   A. baumannii       2 0.4700000   1.0000000   0.3915040 0.1957511
    ## 4   A. baumannii       3 0.5400000   0.9662162   0.4659624 0.2391304
    ## 5   A. baumannii       4 0.5390000   0.9310345   0.4725146 0.2303754
    ## 6   A. baumannii       5 0.5870000   0.9593496   0.5347777 0.2243346
    ## 7  K. pneumoniae Overall 0.7916328   0.9707521   0.7610556 0.4095182
    ## 8  K. pneumoniae       1 0.8726179   0.9651163   0.8533333 0.5783972
    ## 9  K. pneumoniae       2 0.7938776   0.9833333   0.7674419 0.3710692
    ## 10 K. pneumoniae       3 0.7701736   0.9849624   0.7364066 0.3700565
    ## 11 K. pneumoniae       4 0.7515275   0.9803922   0.7092883 0.3836317
    ## 12 K. pneumoniae       5 0.7687627   0.9428571   0.7399527 0.3750000

### Plot statistics from simulated mixtures

``` r
mixedSpecStats %>%
    gather(statistic, value, -orgResults, -n) %>%
    ggplot(aes(x = n, y = value)) +
    geom_col(position = "dodge") +
    facet_grid(orgResults ~ statistic) +
    ylab("Value") +
    xlab("Number of Species in Mixture") +
    ggtitle("Species-Level Detection")
```

``` r
ggsave("../results/simSpecStats.pdf", width = 210, height = 80, units = "mm", useDingbats = F)

mixedResStats %>%
    gather(statistic, value, -orgResults, -n) %>%
    ggplot(aes(x = n, y = value)) +
    geom_col(position = "dodge") +
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
    filter(!(truth == "other" & percentEc != 100))

saveRDS(twoAbRes, file = "../temp/twoAbRes.rds")
saveRDS(twoKpRes, file = "../temp/twoKpRes.rds")
```

### Plot scores vs cutoff

``` r
tholdSpecies <- mixedThold %>%
    filter(str_detect(type, "species"))

twoRes %>%
    ggplot(aes(x = as.factor(percentTarget), y = log10(speciesVsOther), fill = Colistin)) +
    #geom_hline(data = tholdSpecies, aes(yintercept = threshold), 
    #           size = 1, linetype = "dashed") +
    geom_jitter(width = 0.25, shape = 21) +
    facet_wrap(~orgResults, ncol = 2) +
    scale_fill_discrete(breaks = c("resistant", "susceptible")) +
    xlab("Percent Target Species by Volume") +
    ylab("Species log(Score)") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom") #+
```

``` r
    #ylim(c(0, 1))

ggsave("../results/twoSpeciesSpecRes.pdf", width = 105, height = 60, units = "mm", useDingbats = F)


tholdRes <- mixedThold %>%
    filter(str_detect(type, "res"))
```

``` r
twoRes %>%
    ggplot(aes(x = as.factor(percentTarget), y = log10(pos), fill = Colistin)) +
    #geom_hline(data = tholdRes, aes(yintercept = log10(threshold)), 
    #           size = 1, linetype = "dashed") +
    geom_jitter(width = 0.25, shape = 21) +
    facet_wrap(~orgResults, ncol = 2) +
    scale_fill_discrete(breaks = c("resistant", "susceptible")) +
    xlab("Percent Target Species by Volume") +
    ylab("Colistin Resistance log[Score]") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom") #+
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
    ##  version  R version 3.4.2 (2017-09-28)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  tz       America/New_York            
    ##  date     2017-10-27                  
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
    ##  e1071                1.6-8    2017-02-02 CRAN (R 3.4.2)
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
