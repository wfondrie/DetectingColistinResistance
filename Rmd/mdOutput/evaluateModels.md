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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/AbPR-1.png)

``` r
ggsave("../results/AbPrCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

# ROC
rocnames <- c("1 - Specificity (FPR)", "Sensitivity (TPR)", "threshold")
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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/AbPR-2.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/KpPR-1.png)

``` r
ggsave("../results/KpPrCurve.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

# ROC
rocnames <- c("1 - Specificity (FPR)", "Sensitivity (TPR)", "threshold")
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/KpPR-2.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/KpPR-3.png)

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

# Colistin-Resistant Ab
confusionMatrix(data = abRes2$posPred, reference = abRes2$posTruth, 
                positive = "pos", mode = "everything")


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

# Colistin-Resistant Kp
confusionMatrix(data = kpRes2$posPred, reference = kpRes2$posTruth, 
                positive = "pos", mode = "everything")
```

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/varImp-1.png)

``` r
ggsave("../results/vImpAbPlot.pdf", width = 70, height = 50, units = "mm", useDingbats = F)

saveRDS(vImpAb, "../temp/vImpAb.rds")
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/varImp-2.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/varImp-3.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/simMixPlot-1.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/simMixPlot-2.png)![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/simMixPlot-3.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/simAuc-1.png)

``` r
ggsave("../results/mixedCurvesAUCPlot.pdf", width = 70, height = 80, units = "mm", useDingbats = F)
```

``` r
mixedAbRes <- mixedAbRes %>%
    mutate(posPred = ifelse(pos >= thold$threshold[thold$type == "Ab_resistant"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= thold$threshold[thold$type == "Ab_species"],
                             "Ab", "other"),
           specTruth = ifelse(truth != "other", "Ab", "other"))


mixedKpRes <- mixedKpRes %>%
    mutate(posPred = ifelse(pos >= thold$threshold[thold$type == "Kp_resistant"],
                            "pos", "neg"),
           posTruth = ifelse(truth == "pos", "pos", "neg"),
           specPred = ifelse(speciesVsOther >= thold$threshold[thold$type == "Kp_species"],
                             "Kp", "other"),
           specTruth = ifelse(truth != "other", "Kp", "other"))


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
    
    return(ret)
}

# Ab Species
specStats(mixedAbRes, "Ab")
```

    ##   Accuracy     Kappa AccuracyLower AccuracyUpper AccuracyNull
    ## 1   0.8832 0.6593830     0.8739730     0.8919769       0.7282
    ## 2   0.9980 0.9949956     0.9927942     0.9997577       0.7230
    ## 3   0.9140 0.7566332     0.8948797     0.9306375       0.7320
    ## 4   0.8340 0.5335192     0.8094632     0.8565513       0.6980
    ## 5   0.8270 0.4501303     0.8021034     0.8499536       0.7300
    ## 6   0.8430 0.4527557     0.8189474     0.8650123       0.7580
    ##   AccuracyPValue McnemarPValue Sensitivity Specificity Pos.Pred.Value
    ## 1  3.560493e-158 7.411441e-127   0.5717439   0.9994507      0.9974326
    ## 2  1.013410e-136  4.795001e-01   0.9927798   1.0000000      1.0000000
    ## 3   1.235919e-47  3.550956e-19   0.6828358   0.9986339      0.9945652
    ## 4   2.891881e-23  1.508154e-37   0.4503311   1.0000000      1.0000000
    ## 5   3.268622e-13  4.461673e-39   0.3592593   1.0000000      1.0000000
    ## 6   3.195514e-11  1.018261e-34   0.3553719   0.9986807      0.9885057
    ##   Neg.Pred.Value Precision    Recall        F1 Prevalence Detection.Rate
    ## 1      0.8621180 0.9974326 0.5717439 0.7268475     0.2718         0.1554
    ## 2      0.9972414 1.0000000 0.9927798 0.9963768     0.2770         0.2750
    ## 3      0.8958333 0.9945652 0.6828358 0.8097345     0.2680         0.1830
    ## 4      0.8078704 1.0000000 0.4503311 0.6210046     0.3020         0.1360
    ## 5      0.8084164 1.0000000 0.3592593 0.5286104     0.2700         0.0970
    ## 6      0.8291347 0.9885057 0.3553719 0.5227964     0.2420         0.0860
    ##   Detection.Prevalence Balanced.Accuracy       n
    ## 1               0.1558         0.7855973 Overall
    ## 2               0.2750         0.9963899       1
    ## 3               0.1840         0.8407349       2
    ## 4               0.1360         0.7251656       3
    ## 5               0.0970         0.6796296       4
    ## 6               0.0870         0.6770263       5

``` r
# Kp Species
specStats(mixedKpRes, "Kp")
```

    ##    Accuracy     Kappa AccuracyLower AccuracyUpper AccuracyNull
    ## 1 0.9451665 0.8598316     0.9384377     0.9513612    0.7156783
    ## 2 0.9869609 0.9705804     0.9778059     0.9930394    0.6690070
    ## 3 0.9755102 0.9376392     0.9637792     0.9842473    0.7346939
    ## 4 0.9560776 0.8859750     0.9412902     0.9680338    0.7252298
    ## 5 0.9205703 0.7850371     0.9018575     0.9367104    0.7240326
    ## 6 0.8864097 0.6804704     0.8649326     0.9055459    0.7261663
    ##   AccuracyPValue McnemarPValue Sensitivity Specificity Pos.Pred.Value
    ## 1   0.000000e+00  1.236630e-27   0.8392857   0.9872304      0.9631148
    ## 2  1.405346e-149  1.000000e+00   0.9818182   0.9895052      0.9788520
    ## 3   1.180997e-94  3.074342e-01   0.9653846   0.9791667      0.9436090
    ## 4   5.659585e-80  1.955081e-05   0.8661710   0.9901408      0.9708333
    ## 5   4.844728e-54  1.841850e-13   0.7343173   0.9915612      0.9707317
    ## 6   7.841891e-35  8.062992e-18   0.6222222   0.9860335      0.9438202
    ##   Neg.Pred.Value Precision    Recall        F1 Prevalence Detection.Rate
    ## 1      0.9392549 0.9631148 0.8392857 0.8969466  0.2843217      0.2386271
    ## 2      0.9909910 0.9788520 0.9818182 0.9803328  0.3309930      0.3249749
    ## 3      0.9873950 0.9436090 0.9653846 0.9543726  0.2653061      0.2561224
    ## 4      0.9512855 0.9708333 0.8661710 0.9155206  0.2747702      0.2379980
    ## 5      0.9073359 0.9707317 0.7343173 0.8361345  0.2759674      0.2026477
    ## 6      0.8737624 0.9438202 0.6222222 0.7500000  0.2738337      0.1703854
    ##   Detection.Prevalence Balanced.Accuracy       n
    ## 1            0.2477660         0.9132581 Overall
    ## 2            0.3319960         0.9856617       1
    ## 3            0.2714286         0.9722756       2
    ## 4            0.2451481         0.9281559       3
    ## 5            0.2087576         0.8629393       4
    ## 6            0.1805274         0.8041279       5

``` r
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
    
    return(ret)
}

# Ab Species
resStats(mixedAbRes)
```

    ##   Accuracy     Kappa AccuracyLower AccuracyUpper AccuracyNull
    ## 1    0.904 0.6065774     0.8954970     0.9120279        0.863
    ## 2    0.980 0.9198204     0.9692800     0.9877417        0.860
    ## 3    0.903 0.6198166     0.8829522     0.9206355        0.871
    ## 4    0.882 0.5637643     0.8603723     0.9013475        0.852
    ## 5    0.866 0.4403375     0.8433072     0.8865117        0.855
    ## 6    0.889 0.4570322     0.8678739     0.9078021        0.877
    ##   AccuracyPValue McnemarPValue Sensitivity Specificity Pos.Pred.Value
    ## 1   6.043112e-19  1.992162e-02   0.6875912   0.9383546      0.6390773
    ## 2   2.087008e-40  1.390630e-02   0.9714286   0.9813953      0.8947368
    ## 3   1.051712e-03  4.878252e-05   0.7829457   0.9207807      0.5941176
    ## 4   3.514773e-03  2.136697e-02   0.6891892   0.9154930      0.5862069
    ## 5   1.731055e-01  3.419826e-01   0.4965517   0.9286550      0.5413534
    ## 6   1.334245e-01  1.839070e-01   0.4878049   0.9452680      0.5555556
    ##   Neg.Pred.Value Precision    Recall        F1 Prevalence Detection.Rate
    ## 1      0.9498006 0.6390773 0.6875912 0.6624473      0.137         0.0942
    ## 2      0.9952830 0.8947368 0.9714286 0.9315068      0.140         0.1360
    ## 3      0.9662651 0.5941176 0.7829457 0.6755853      0.129         0.1010
    ## 4      0.9443099 0.5862069 0.6891892 0.6335404      0.148         0.1020
    ## 5      0.9158016 0.5413534 0.4965517 0.5179856      0.145         0.0720
    ## 6      0.9293722 0.5555556 0.4878049 0.5194805      0.123         0.0600
    ##   Detection.Prevalence Balanced.Accuracy       n
    ## 1               0.1474         0.8129729 Overall
    ## 2               0.1520         0.9764120       1
    ## 3               0.1700         0.8518632       2
    ## 4               0.1740         0.8023411       3
    ## 5               0.1330         0.7126033       4
    ## 6               0.1080         0.7165364       5

``` r
# Kp Species
resStats(mixedKpRes)
```

    ##    Accuracy     Kappa AccuracyLower AccuracyUpper AccuracyNull
    ## 1 0.7106011 0.3570065     0.6977110     0.7232411    0.8541836
    ## 2 0.7933801 0.5043871     0.7668973     0.8181176    0.8274824
    ## 3 0.7244898 0.3438151     0.6953613     0.7522625    0.8775510
    ## 4 0.6956078 0.3288010     0.6657157     0.7243180    0.8641471
    ## 5 0.6720978 0.3233238     0.6417440     0.7014155    0.8441955
    ## 6 0.6663286 0.3001299     0.6359266     0.6957335    0.8580122
    ##   AccuracyPValue McnemarPValue Sensitivity Specificity Pos.Pred.Value
    ## 1      1.0000000 2.318525e-299   0.9805014   0.6645269      0.3328605
    ## 2      0.9976757  1.032448e-43   0.9825581   0.7539394      0.4543011
    ## 3      1.0000000  1.638073e-58   0.9833333   0.6883721      0.3056995
    ## 4      1.0000000  1.300360e-64   0.9849624   0.6501182      0.3067916
    ## 5      1.0000000  5.517788e-69   0.9803922   0.6151990      0.3198294
    ## 6      1.0000000  1.168826e-69   0.9714286   0.6158392      0.2950108
    ##   Neg.Pred.Value Precision    Recall        F1 Prevalence Detection.Rate
    ## 1      0.9950160 0.3328605 0.9805014 0.4969996  0.1458164      0.1429732
    ## 2      0.9952000 0.4543011 0.9825581 0.6213235  0.1725176      0.1695085
    ## 3      0.9966330 0.3056995 0.9833333 0.4664032  0.1224490      0.1204082
    ## 4      0.9963768 0.3067916 0.9849624 0.4678571  0.1358529      0.1338100
    ## 5      0.9941520 0.3198294 0.9803922 0.4823151  0.1558045      0.1527495
    ## 6      0.9923810 0.2950108 0.9714286 0.4525790  0.1419878      0.1379310
    ##   Detection.Prevalence Balanced.Accuracy       n
    ## 1            0.4295288         0.8225141 Overall
    ## 2            0.3731194         0.8682488       1
    ## 3            0.3938776         0.8358527       2
    ## 4            0.4361593         0.8175403       3
    ## 5            0.4775967         0.7977956       4
    ## 6            0.4675456         0.7936339       5

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
tholdSpecies <- thold %>%
    filter(str_detect(type, "species")) %>%
    mutate(orgResults = ifelse(str_detect(type, "Ab"),
                               "A. baumannii",
                               "K. pneumoniae"))

twoRes %>%
    ggplot(aes(x = as.factor(percentTarget), y = speciesVsOther, fill = Colistin)) +
    geom_hline(data = tholdSpecies, aes(yintercept = threshold), 
               size = 1, linetype = "dashed") +
    geom_jitter(width = 0.25, shape = 21) +
    facet_wrap(~orgResults, ncol = 2) +
    scale_fill_discrete(breaks = c("resistant", "susceptible")) +
    xlab("Percent Target Species by Volume") +
    ylab("Species Score") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom") +
    ylim(c(0, 1))
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/twoSpecScores-1.png)

``` r
ggsave("../results/twoSpeciesSpecRes.pdf", width = 105, height = 60, units = "mm", useDingbats = F)


tholdRes <- thold %>%
    filter(str_detect(type, "res")) %>%
    mutate(orgResults = ifelse(str_detect(type, "Ab"),
                               "A. baumannii",
                               "K. pneumoniae"))
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/twoSpecScores-2.png)

``` r
twoRes %>%
    ggplot(aes(x = as.factor(percentTarget), y = pos, fill = Colistin)) +
    geom_hline(data = tholdRes, aes(yintercept = threshold), 
               size = 1, linetype = "dashed") +
    geom_jitter(width = 0.25, shape = 21) +
    facet_wrap(~orgResults, ncol = 2) +
    scale_fill_discrete(breaks = c("resistant", "susceptible")) +
    xlab("Percent Target Species by Volume") +
    ylab("Colistin-Resistance Score") +
    theme(legend.key.size = unit(0.75, "lines"),
          legend.margin = margin(t = -0.5, unit = "lines"),
          legend.position = "bottom") +
    ylim(c(0, 1))
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\evaluateModels_files/figure-markdown_github-ascii_identifiers/twoSpecScores-3.png)

``` r
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
    ##  date     2017-10-25                  
    ## 
    ##  package      * version  date       source        
    ##  assertthat     0.2.0    2017-04-11 CRAN (R 3.4.2)
    ##  backports      1.1.1    2017-09-25 CRAN (R 3.4.1)
    ##  base         * 3.4.2    2017-09-28 local         
    ##  bindr          0.1      2016-11-13 CRAN (R 3.4.2)
    ##  bindrcpp     * 0.2      2017-06-17 CRAN (R 3.4.2)
    ##  broom          0.4.2    2017-02-13 CRAN (R 3.4.0)
    ##  caret        * 6.0-77   2017-09-07 CRAN (R 3.4.2)
    ##  cellranger     1.1.0    2016-07-27 CRAN (R 3.4.2)
    ##  class          7.3-14   2015-08-30 CRAN (R 3.4.2)
    ##  codetools      0.2-15   2016-10-05 CRAN (R 3.4.2)
    ##  colorspace     1.3-2    2016-12-14 CRAN (R 3.4.2)
    ##  compiler       3.4.2    2017-09-28 local         
    ##  CVST           0.2-1    2013-12-10 CRAN (R 3.4.2)
    ##  data.table     1.10.4-2 2017-10-12 CRAN (R 3.4.2)
    ##  datasets     * 3.4.2    2017-09-28 local         
    ##  ddalpha        1.3.1    2017-09-27 CRAN (R 3.4.2)
    ##  DEoptimR       1.0-8    2016-11-19 CRAN (R 3.4.1)
    ##  devtools     * 1.13.3   2017-08-02 CRAN (R 3.4.2)
    ##  digest         0.6.12   2017-01-27 CRAN (R 3.4.2)
    ##  dimRed         0.1.0    2017-05-04 CRAN (R 3.4.2)
    ##  dplyr        * 0.7.4    2017-09-28 CRAN (R 3.4.2)
    ##  DRR            0.0.2    2016-09-15 CRAN (R 3.4.2)
    ##  e1071          1.6-8    2017-02-02 CRAN (R 3.4.2)
    ##  evaluate       0.10.1   2017-06-24 CRAN (R 3.4.2)
    ##  forcats      * 0.2.0    2017-01-23 CRAN (R 3.4.2)
    ##  foreach        1.4.3    2015-10-13 CRAN (R 3.4.2)
    ##  foreign        0.8-69   2017-06-22 CRAN (R 3.4.2)
    ##  ggplot2      * 2.2.1    2016-12-30 CRAN (R 3.4.2)
    ##  glue           1.1.1    2017-06-21 CRAN (R 3.4.2)
    ##  gower          0.1.2    2017-02-23 CRAN (R 3.4.2)
    ##  graphics     * 3.4.2    2017-09-28 local         
    ##  grDevices    * 3.4.2    2017-09-28 local         
    ##  grid           3.4.2    2017-09-28 local         
    ##  gtable         0.2.0    2016-02-26 CRAN (R 3.4.2)
    ##  haven          1.1.0    2017-07-09 CRAN (R 3.4.2)
    ##  hms            0.3      2016-11-22 CRAN (R 3.4.2)
    ##  htmltools      0.3.6    2017-04-28 CRAN (R 3.4.2)
    ##  httr           1.3.1    2017-08-20 CRAN (R 3.4.2)
    ##  ipred          0.9-6    2017-03-01 CRAN (R 3.4.2)
    ##  iterators      1.0.8    2015-10-13 CRAN (R 3.4.1)
    ##  jsonlite       1.5      2017-06-01 CRAN (R 3.4.2)
    ##  kernlab        0.9-25   2016-10-03 CRAN (R 3.4.1)
    ##  knitr          1.17     2017-08-10 CRAN (R 3.4.2)
    ##  labeling       0.3      2014-08-23 CRAN (R 3.4.1)
    ##  lattice      * 0.20-35  2017-03-25 CRAN (R 3.4.2)
    ##  lava           1.5.1    2017-09-27 CRAN (R 3.4.2)
    ##  lazyeval       0.2.0    2016-06-12 CRAN (R 3.4.2)
    ##  lubridate      1.6.0    2016-09-13 CRAN (R 3.4.2)
    ##  magrittr       1.5      2014-11-22 CRAN (R 3.4.2)
    ##  MASS           7.3-47   2017-02-26 CRAN (R 3.4.2)
    ##  Matrix         1.2-11   2017-08-21 CRAN (R 3.4.2)
    ##  memoise        1.1.0    2017-04-21 CRAN (R 3.4.2)
    ##  methods      * 3.4.2    2017-09-28 local         
    ##  mnormt         1.5-5    2016-10-15 CRAN (R 3.4.1)
    ##  ModelMetrics   1.1.0    2016-08-26 CRAN (R 3.4.2)
    ##  modelr         0.1.1    2017-07-24 CRAN (R 3.4.2)
    ##  munsell        0.4.3    2016-02-13 CRAN (R 3.4.2)
    ##  nlme           3.1-131  2017-02-06 CRAN (R 3.4.2)
    ##  nnet           7.3-12   2016-02-02 CRAN (R 3.4.2)
    ##  parallel       3.4.2    2017-09-28 local         
    ##  pkgconfig      2.0.1    2017-03-21 CRAN (R 3.4.2)
    ##  plyr           1.8.4    2016-06-08 CRAN (R 3.4.2)
    ##  prodlim        1.6.1    2017-03-06 CRAN (R 3.4.2)
    ##  PRROC        * 1.3      2017-04-21 CRAN (R 3.4.2)
    ##  psych          1.7.8    2017-09-09 CRAN (R 3.4.2)
    ##  purrr        * 0.2.4    2017-10-18 CRAN (R 3.4.2)
    ##  R6             2.2.2    2017-06-17 CRAN (R 3.4.2)
    ##  Rcpp           0.12.13  2017-09-28 CRAN (R 3.4.2)
    ##  RcppRoll       0.2.2    2015-04-05 CRAN (R 3.4.2)
    ##  readr        * 1.1.1    2017-05-16 CRAN (R 3.4.2)
    ##  readxl         1.0.0    2017-04-18 CRAN (R 3.4.2)
    ##  recipes        0.1.0    2017-07-27 CRAN (R 3.4.2)
    ##  reshape2       1.4.2    2016-10-22 CRAN (R 3.4.2)
    ##  rlang          0.1.2    2017-08-09 CRAN (R 3.4.2)
    ##  rmarkdown    * 1.6      2017-06-15 CRAN (R 3.4.2)
    ##  robustbase     0.92-7   2016-12-09 CRAN (R 3.4.2)
    ##  rpart          4.1-11   2017-03-13 CRAN (R 3.4.2)
    ##  rprojroot      1.2      2017-01-16 CRAN (R 3.4.2)
    ##  rvest          0.3.2    2016-06-17 CRAN (R 3.4.2)
    ##  scales         0.5.0    2017-08-24 CRAN (R 3.4.2)
    ##  sfsmisc        1.1-1    2017-06-08 CRAN (R 3.4.2)
    ##  splines        3.4.2    2017-09-28 local         
    ##  stats        * 3.4.2    2017-09-28 local         
    ##  stats4         3.4.2    2017-09-28 local         
    ##  stringi        1.1.5    2017-04-07 CRAN (R 3.4.1)
    ##  stringr      * 1.2.0    2017-02-18 CRAN (R 3.4.2)
    ##  survival       2.41-3   2017-04-04 CRAN (R 3.4.2)
    ##  tibble       * 1.3.4    2017-08-22 CRAN (R 3.4.2)
    ##  tidyr        * 0.7.2    2017-10-16 CRAN (R 3.4.2)
    ##  tidyselect     0.2.2    2017-10-10 CRAN (R 3.4.2)
    ##  tidyverse    * 1.1.1    2017-01-27 CRAN (R 3.4.2)
    ##  timeDate       3012.100 2015-01-23 CRAN (R 3.4.1)
    ##  tools          3.4.2    2017-09-28 local         
    ##  utils        * 3.4.2    2017-09-28 local         
    ##  withr          2.0.0    2017-07-28 CRAN (R 3.4.2)
    ##  xgboost      * 0.6-4    2017-01-05 CRAN (R 3.4.2)
    ##  xml2           1.1.1    2017-01-24 CRAN (R 3.4.2)
    ##  yaml           2.1.14   2016-11-12 CRAN (R 3.4.2)
