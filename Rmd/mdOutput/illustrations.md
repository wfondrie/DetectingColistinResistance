Strategy Illustrations
================
William E Fondrie

-   [Load Libraries and Prepare Workspace](#load-libraries-and-prepare-workspace)
-   [Feature Creation](#feature-creation)
-   [Real vs Simulated Mixture Spectra](#real-vs-simulated-mixture-spectra)
    -   [Feature extraction performance](#feature-extraction-performance)
-   [Session Info](#session-info)

Load Libraries and Prepare Workspace
------------------------------------

``` r
suppressMessages(library(tidyverse, quietly = T))
library(MALDIquant, quietly = T)
library(MALDIquantForeign, quietly = T)
library(caret, quietly = T)
library(PRROC, quietly = T)
library(stringr, quietly = T)
library(forcats, quietly = T)
library(devtools, quietly = T)

# ggplot2 theme
source("../R/ggplotTheme.R")
theme_set(coolTheme)

source("../R/preProcessSpec.R")
source("../R/extract.R")

set.seed(6554)
```

Feature Creation
----------------

``` r
features <- readRDS("../temp/features.rds")
features <- filter(features, response == "Ab")

featFiles <- c(list.files("../data/fullLib/Acinetobacter baumannii - res", 
                    full.names = T, 
                    pattern = "mzXML$",
                    recursive = T)[1],
                list.files("../data/fullLib/Acinetobacter baumannii - sen", 
                    full.names = T, 
                    pattern = "mzXML$",
                    recursive = T)[1])

spec <- preProcessSpec(featFiles, hws = 80)

specDat <- map_df(spec, extractSpectra) %>%
    mutate(class = str_replace(type, "Acinetobacter", "A."),
           class = str_replace(class, " - res", " EAS004, Colistin-Resistant"),
           class = str_replace(class, " - sen", " EAS001, Colistin-Susceptible"),
           relInt = relInt * 100)


featurePlot <- specDat %>%
    ggplot(aes(x = mz, y = relInt)) +
    geom_segment(data = features, aes(x = mz, xend = mz, y = -1, yend = -10), 
                 color = ggColors(3)[3]) +
    geom_line() +
    facet_wrap(~ class, nrow = 2) +
    xlab(expression(italic("m/z"))) +
    ylab("Relative Intensity")

featurePlot 
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/createFeatures-1.png)

``` r
ggsave("../results/featureSelect_1.pdf", width = 105, height = 100, unit = "mm", useDingbats = F)

singlePlot <- specDat %>%
    filter(str_detect(type, "res")) %>%
    ggplot(aes(x = mz, y = relInt)) +
    geom_line() +
    xlab(expression(italic("m/z"))) +
    ylab("Relative Intensity")

singlePlot
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/createFeatures-2.png)

``` r
ggsave("../results/singlePlot.pdf", width = 105, height = 50, unit = "mm", useDingbats = F)

selectedMz <- features$mz[which.min(abs(features$mz - 2033))] 

maxDf <- specDat %>%
    filter(mz > selectedMz - 1.5,
           mz < selectedMz + 1.5) %>%
    group_by(class) %>%
    summarize(mz = mz[relInt == max(relInt)],
              relInt = max(relInt))

featurePlotZoom1 <- specDat %>%
    filter(str_detect(class, "Resistant")) %>%
    ggplot(aes(x = mz, y = relInt)) +
    annotate("rect", xmin = selectedMz - 1.5, xmax = selectedMz + 1.5, 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = ggColors(3)[3], color = ggColors(3)[3]) +
    geom_vline(data = features,  aes(xintercept = mz), color = ggColors(3)[3], linetype = "dashed") +
    geom_line() +
    geom_point(data = filter(maxDf, str_detect(class, "Resistant")), 
               aes(x = mz, y = relInt + 2), color = "black", shape = 8, size = 0.5) + 
    xlim(selectedMz + c(-5, 5)) +
    ylim(c(0, 35)) +
    theme(axis.title = element_blank())

featurePlotZoom1
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/createFeatures-3.png)

``` r
ggsave("../results/featureSelect_2_1.pdf",  width = 30, height = 20, unit = "mm", useDingbats = F)


# Zoom 2
featurePlotZoom2 <- specDat %>%
    filter(str_detect(class, "Susceptible")) %>%
    ggplot(aes(x = mz, y = relInt)) +
    annotate("rect", xmin = selectedMz - 1.5, xmax = selectedMz + 1.5, 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = ggColors(3)[3], color = ggColors(3)[3]) +
    geom_vline(data = features,  aes(xintercept = mz), color = ggColors(3)[3], linetype = "dashed") +
    geom_line() +
    geom_point(data = filter(maxDf, str_detect(class, "Susceptible")), 
               aes(x = mz, y = relInt + 2), color = "black", shape = 8, size = 0.5) + 
    xlim(selectedMz + c(-5, 5)) +
    ylim(c(0, 35)) +
    theme(axis.title = element_blank())

featurePlotZoom2
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/createFeatures-4.png)

``` r
ggsave("../results/featureSelect_2_2.pdf",  width = 30, height = 20, unit = "mm", useDingbats = F)
```

Real vs Simulated Mixture Spectra
---------------------------------

``` r
colRanges <- tibble(minMz = c(1300, 1440, 1810),
                    maxMz = c(1430, 1640, 2100),
                    lab = c("S. aureus",
                            "P. aeruginosa",
                            "K. pneumoniae"))

shortenName <- function(origName) {
  paste0(str_match(origName, "^."),
         ". ",
         str_match(origName, " ([^ ]+)")[ , 2])
}


# import and plot individual spec
soi <- c("S aureus NRS384 1 10282016",
         "P aeruginosa BE399 4 10282016",
         "K pneumoniae TBE818 5 10132016")
files <- list.files("../data/fullLib", 
                    full.names = T, 
                    pattern = "mzXML$",
                    recursive = T)

files <- map_chr(soi, ~files[str_detect(files, .)])

spec <- preProcessSpec(files, hws = 80)

specDat <- map_df(spec, extractSpectra) 

spPlot <- specDat %>%
  mutate(spc = shortenName(type),
         relInt = relInt * 100,
         mzRange = ifelse(mz > colRanges$minMz[1] & mz < colRanges$maxMz[1] & spc == "S. aureus", 
                          colRanges$lab[1], NA),
         mzRange = ifelse(mz > colRanges$minMz[2] & mz < colRanges$maxMz[2] & spc == "P. aeruginosa", 
                          colRanges$lab[2], mzRange),
         mzRange = ifelse(mz > colRanges$minMz[3] & mz < colRanges$maxMz[3] & spc == "K. pneumoniae", 
                          colRanges$lab[3], mzRange)) 


spPlot %>%
  ggplot(aes(x = as.numeric(mz), y = relInt)) +
  geom_line(size = 0.25) +
  geom_line(data = spPlot[!is.na(spPlot$mzRange), ], aes(color = mzRange), size = 0.30) +
  coolTheme +
  theme(legend.position = "none") +
  labs(x = expression(italic("m/z")),
       y = "Relative Intensity") +
  facet_wrap(~ spc, ncol = 1)
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/realVsSim-1.png)

``` r
ggsave("../results/individualSpec.pdf", width = 105, height = 100, unit = "mm", useDingbats = F)


 # Make combined Spectrum ----------------------------------------------------------------------
cmb <- specDat %>% 
  mutate(spec_id = "spec_1") %>%
  group_by(spec_id ,id, type) %>%
  do(massSpecObj = createMassSpectrum(mass = .$mz, intensity = .$relInt, 
                                      metaData = list(file = "a"))) %>%
  group_by(spec_id) %>%
  do(massSpecObj = averageMassSpectra(.$massSpecObj, method = "sum"))
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/realVsSim-2.png)

``` r
cmbSpec <- extractSpectra(cmb$massSpecObj[[1]])

cmbSpecNorm <-  cmbSpec %>%
  mutate(relInt = relInt / max(relInt) * 100,
         type = "combined",
         id = "combined") %>%
  group_by(mz) %>%
  mutate(mzRange = ifelse(mz > colRanges$minMz[1] & mz < colRanges$maxMz[1], 
                          colRanges$lab[1], NA),
         mzRange = ifelse(mz > colRanges$minMz[2] & mz < colRanges$maxMz[2], 
                          colRanges$lab[2], mzRange),
         mzRange = ifelse(mz > colRanges$minMz[3] & mz < colRanges$maxMz[3], 
                          colRanges$lab[3], mzRange))

cmbSpecNorm %>%
  ggplot(aes(x = as.numeric(mz), y = relInt)) +
  geom_line(size = 0.25) +
  geom_line(data = cmbSpecNorm[!is.na(cmbSpecNorm$mzRange), ], aes(color = mzRange), size = 0.30) +
  theme(legend.position = "none") +
  labs(title = "Simulated",
       x = expression(italic("m/z")),
       y = "Relative Intensity")
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/realVsSim-3.png)

``` r
ggsave("../results/inSilicoSpec.pdf", width = 105, height = 50, unit = "mm", useDingbats = F)

# Import and plot experimental mixed spectrum --------------------------------------------------------
mixFile <- list.files("../data/mixSpec", pattern = ".mzXML", full.names = T)[1]
origSpec <- preProcessSpec(mixFile, hws = 80)
origSpecDat <- map_df(origSpec, extractSpectra)

origSpecNorm <- origSpecDat %>%
  group_by(mz) %>%
  summarize(relInt = sum(relInt)) %>%
  ungroup() %>%
  mutate(relInt = relInt / max(relInt) * 100,
         type = "combined",
         id = "combined") %>%
  mutate(mzRange = ifelse(mz > colRanges$minMz[1] & mz < colRanges$maxMz[1], 
                          colRanges$lab[1], NA),
         mzRange = ifelse(mz > colRanges$minMz[2] & mz < colRanges$maxMz[2], 
                          colRanges$lab[2], mzRange),
         mzRange = ifelse(mz > colRanges$minMz[3] & mz < colRanges$maxMz[3], 
                          colRanges$lab[3], mzRange))

origSpecNorm %>%
  ggplot(aes(x = as.numeric(mz), y = relInt)) +
  geom_line(size = 0.25) +
  geom_line(data = origSpecNorm[!is.na(origSpecNorm$mzRange), ], aes(color = mzRange), size = 0.30) +
  coolTheme +
  theme(legend.position = "none") +
  labs(title = "Experimental",
       x = expression(italic("m/z")),
       y = "Relative Intensity")
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/realVsSim-4.png)

``` r
ggsave("../results/experimentalSpec.pdf", width = 105, height = 50, unit = "mm", useDingbats = F)
```

### Feature extraction performance

``` r
vImpKp <- readRDS("../temp/vImpKp.rds") %>%
    arrange(desc(Gain)) %>%
    mutate(featureMz = as.factor(str_replace(Feature, "mz", "")),
           featureMz = fct_reorder(featureMz, Gain, .desc = T),
           cols = c(rep("top", 5), rep("bot", length(Feature) - 5))) %>%
    select(featureMz, Feature, Gain, cols) %>%
    rename(feat = Feature)

vImpShort <- vImpKp[1:20, ]

singleSpec <- specDat %>%
    filter(type == "Klebsiella pneumoniae - sen")

kpFeat <- readRDS("../temp/features.rds") %>%
    filter(response == "Kp")

singleFeat <- extractFeatures(singleSpec, kpFeat$mz, tol = 1.5) %>%
    mutate(type = "Isolate K. pneumoniae",
           id = "isolate")
mixFeat <- extractFeatures(cmbSpec, kpFeat$mz, tol = 1.5) %>%
    mutate(type = "Simulated Mixture with K. pneumoniae",
           id = "mixture")

feats <- singleFeat %>%
    full_join(mixFeat)


feats %>%
    right_join(vImpShort) %>%
    ggplot(aes(x = featureMz, y = relInt, ymax = relInt, ymin = 0)) +
    geom_linerange(size = 0.75) +
    geom_point(size = 2) +
    ylim(c(0, 1.05)) +
    ylab("Relative Intensity") +
    xlab(expression(paste("Feature (", italic("m/z"), ")"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    facet_wrap(~ type, ncol = 2)
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/mixtureFeatureExtraction-1.png)

``` r
ggsave("../results/featureComparison.pdf", width = 160, height = 50, units = "mm", useDingbats = F)

feats %>%
    right_join(vImpKp) %>%
    select(featureMz, id, relInt, Gain) %>%
    spread(id, relInt) %>%
    mutate(diff = abs(isolate - mixture),
           Interference = ifelse(diff <= 0.05, "No", "Yes")) %>%
    ggplot(aes(x = isolate, y = mixture, size = Gain, fill = Interference)) +
    geom_point(shape = 21) +
    scale_size_continuous(guide = "none") +
    ylim(c(0, 1.05)) +
    xlim(c(0, 1.05)) +
    theme(legend.position = c(0.99, 0.01),
          legend.justification = c(1, 0),
          legend.background = element_rect(color = "black"),
          legend.key.size = unit(0.5, "lines")) +
    coord_equal() +
    xlab("Feature Intensity in Isolate") +
    ylab("Feature Intensity in Mixture")
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/mixtureFeatureExtraction-2.png)![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\illustrations_files/figure-markdown_github-ascii_identifiers/mixtureFeatureExtraction-3.png)

``` r
ggsave("../results/interferences.pdf", width = 50, height = 50, units = "mm", useDingbats = F)
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
