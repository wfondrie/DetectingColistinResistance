Make Miscellaneous Figures
================
William E Fondrie

-   [Load Libraries and Prepare Workspace](#load-libraries-and-prepare-workspace)
-   [Load Data](#load-data)
-   [Plot of Resistance Peaks and Their Variation](#plot-of-resistance-peaks-and-their-variation)
-   [Mass Accuracy Plot](#mass-accuracy-plot)
-   [Plot of the Isolate Dataset](#plot-of-the-isolate-dataset)
-   [Session Info](#session-info)

Load Libraries and Prepare Workspace
------------------------------------

``` r
# data manipulation
suppressMessages(library(tidyverse, quietly = T))
library(stringr, quietly = T)
library(forcats, quietly = T)
library(devtools, quietly = T)

# handling MALDI spectra
library(MALDIquant, quietly = T)
library(MALDIquantForeign, quietly = T)

# ggplot2 theme
source("../R/ggplotTheme.R")
theme_set(coolTheme)

set.seed(09813734)

# import helper functions
source("../R/preProcessSpec.R")
source("../R/extract.R")
source("../R/utilityFunctions.R")
```

Load Data
---------

``` r
files <- list.files("../data/fullLib", 
                     full.names = T, 
                     pattern = "mzXML$",
                     recursive = T)

specInfo <- tibble(fname = files) %>%
    mutate(type = str_match(fname, "([^\\^/]+)[\\/][^\\^/]+mzXML$")[ , 2],
           id = str_match(fname, "([^\\^/]+).mzXML$")[ , 2],
           species = str_match(type, "^[^ ]+ [^ ]+"),
           Ab = ifelse(str_detect(type, "Acinetobacter baumannii - res"), "pos", "other"),
           Kp = ifelse(str_detect(type, "Klebsiella pneumoniae - res"), "pos", "other"),
           Ab = as.factor(ifelse(str_detect(type, "Acinetobacter baumannii - sen"), "neg", Ab)),
           Kp = as.factor(ifelse(str_detect(type, "Klebsiella pneumoniae - sen"), "neg", Kp)))

summary(specInfo)

features <- readRDS("../temp/features.RDS")

specList <- unlist(map(files, preProcessSpec, hws = 80)) # Same preprocessing

spec <- map_df(specList, extractSpectra)

feattbl <- spec %>%
  group_by(id) %>%
  do(extractFeatures(., featureVec = unique(features$mz), tol = 1.5))
```

Plot of Resistance Peaks and Their Variation
--------------------------------------------

``` r
# Ab
abFeat <- feattbl$feat[feattbl$mz == nom2feat(2036, feattbl$mz)]
# Kp
kpFeat <- c(feattbl$feat[feattbl$mz == nom2feat(1956, feattbl$mz)],
            feattbl$feat[feattbl$mz == nom2feat(1972, feattbl$mz)])

feattbl %>%
  filter((str_detect(type, "Klebsiella pneumoniae") & feat %in% kpFeat) |
           (str_detect(type, "Acinetobacter baumannii") & feat %in% abFeat)) %>%
  mutate(res = ifelse(str_detect(type, "res"), "Colistin\nResistant", "Colistin\nSusceptible"),
         species = ifelse(str_detect(type, "Acineto"), "A. baumannii", "K. pneumoniae"),
         feat = str_replace(feat, "mz", "Feature m/z ")) %>%
  ggplot(aes(y = relInt, x = res)) +
  geom_jitter(alpha = 0.5, size = 1) +
  stat_summary(fun.ymax = median, fun.ymin = median, 
               fun.y = median, color = ggColors(3)[3], geom="crossbar") +
  facet_wrap(~species + feat) +
  ylab("Relative Intensity") +
  theme(axis.title.x = element_blank())
```

``` r
ggsave("../results/resPeaks.pdf", width = 105, height = 50, unit = "mm", useDingbats = F)
```

Mass Accuracy Plot
------------------

``` r
AbSpec <- nom2feat(2033, feattbl$mz)[1]

feattblWide <- spec %>%
  group_by(id) %>%
  do(extractFeatures(., featureVec = AbSpec, tol = 1.5))


specSelect <- c("Acinetobacter baumannii - res",
                "Salmonella minnesota")

feattblWide %>%
    filter(type %in% specSelect) %>%
    mutate(species = ifelse(str_detect(type, "Acinetobacter baumannii - res"),
                             "Resistant\nA. baumannii", "S. minnesota"),
           species = as.factor(species)) %>%
    mutate(feat = str_replace(feat, "mz", "Feature m/z ")) %>%
    ggplot(aes(y = real_mz, x = species)) +
    geom_dotplot(aes(color = fct_rev(species), 
                     fill = fct_rev(species)), 
                 binaxis = "y", stackdir = "center", dotsize = 0.75) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    geom_hline(yintercept = AbSpec, 
               size = 0.5, linetype = "dashed") +
    ylim(c(AbSpec - 1.5, AbSpec + 1.5)) +
    ylab(expression("Extracted"~italic("m/z"))) +
    theme(legend.position = "none",
          axis.title.y = element_blank()) + 
    coord_flip()
```

``` r
ggsave("../results/peakExtraction.pdf", width = 70, height = 50, unit = "mm", useDingbats = F)


twoSpecSel <- c("Ab ACAT PM3850 2 05182016",
                "Salmonella minnesota R595 TBE1121 3 06242016")
```

``` r
twoSpec <- spec %>%
    filter(id %in% twoSpecSel) %>%
    mutate(type = str_replace(type, "Acinetobacter", "Resistant A."),
           type = str_replace(type, " - res", ""),
           type = str_replace(type, "Salmonella", "S."))


twoSpec %>%
    ggplot(aes(x = mz, y = relInt)) +
    geom_line() +
    facet_wrap(~ type, ncol = 1) +
    xlab(expression(italic("m/z"))) +
    ylab("Relative Intensity")
```

``` r
ggsave("../results/minnesotaSpec.pdf", width = 140, height = 100, unit = "mm", useDingbats = F)


twoSpec %>%
    filter(mz >= AbSpec - 5 & mz <= AbSpec + 5) %>%
    group_by(type) %>%
    mutate(relInt = relInt / max(relInt)) %>%
    ungroup() %>%
    mutate(type = str_replace(type, "Resistant ", "Resistant\n")) %>%
    ggplot(aes(x = mz, y = relInt, color = fct_rev(type))) +
    annotate("rect", xmin = AbSpec - 1.5, xmax = AbSpec + 1.5, 
             ymin = -Inf, ymax = Inf, alpha = 1, fill = "lightgrey", color = "lightgrey") +
    geom_vline(xintercept = AbSpec, size = 0.5, linetype = "dashed") +
    geom_vline(xintercept = c(AbSpec - 1.5, AbSpec + 1.5), size = 0.5) +
    geom_line() +
    xlab(expression(italic("m/z"))) +
    ylab("Relative Intensity") +
    theme(legend.title = element_blank(),
          legend.margin = margin(l = -0.5, r = -0.2, unit = "lines"),
          legend.key.size = unit(0.5, "lines"))
```

``` r
ggsave("../results/minnesotaZoom.pdf", width = 70, height = 50, unit = "mm", useDingbats = F)
```

Plot of the Isolate Dataset
---------------------------

``` r
specInfo %>% 
  group_by(species) %>%
  mutate(num = length(id)) %>%
  ungroup() %>%
  mutate(`Colistin\nResistant` = 
           ifelse(Ab == "pos" | Kp == "pos", "Colistin-Resistant", NA),
         `Colistin\nResistant` = 
           ifelse(Ab == "neg" | Kp == "neg", "Colistin-Susceptible", `Colistin\nResistant`),
         `Colistin\nResistant` = fct_rev(`Colistin\nResistant`),
         species = fct_lump(species, 20),
         species = fct_reorder(species, num)) %>%
  ggplot(aes(x = species, fill = `Colistin\nResistant`)) +
  geom_bar(color = "black", width = 0.75) +
  ylab("Number of Mass Spectra") +
  scale_fill_discrete(name = NULL, breaks = c("Colistin-Resistant","Colistin-Susceptible")) +
  coord_flip() +
  coolTheme +
  theme(legend.position = c(0.95, 0.25),
        legend.justification = c(1,0),
        axis.title.y = element_blank(),
        legend.background = element_rect(color = "black"),
        legend.key.size = unit(0.75, "lines"),
        legend.title = element_blank())
```

``` r
ggsave("../results/library.pdf", height = 125, width = 105, unit = "mm", useDingbats = F)


specInfo %>% 
  group_by(species) %>%
  mutate(num = length(id)) %>%
  ungroup() %>%
  mutate(`Colistin\nResistant` = 
           ifelse(Ab == "pos" | Kp == "pos", "Colistin-Resistant", NA),
         `Colistin\nResistant` = 
           ifelse(Ab == "neg" | Kp == "neg", "Colistin-Susceptible", `Colistin\nResistant`),
         `Colistin\nResistant` = fct_rev(`Colistin\nResistant`),
         species = fct_lump(species, 10),
         species = fct_reorder(species, num)) %>%
  ggplot(aes(x = species, fill = `Colistin\nResistant`)) +
  geom_bar(color = "black", width = 0.75) +
  ylab("Number of Mass Spectra") +
  scale_fill_discrete(name = NULL, breaks = c("Colistin-Resistant","Colistin-Susceptible")) +
  coord_flip() +
  theme(legend.position = c(0.95, 0.25),
        legend.justification = c(1,0),
        axis.title.y = element_blank(),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.75, "lines"),
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10))
```

``` r
ggsave("../results/librarySmall.pdf", height = 100, width = 150, unit = "mm", useDingbats = F)

summary(specInfo$species)
length(unique(specInfo$species))
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
    ##  date     2017-12-30                  
    ## 
    ##  package            * version  date       source        
    ##  assertthat           0.2.0    2017-04-11 CRAN (R 3.3.3)
    ##  backports            1.1.0    2017-05-22 CRAN (R 3.4.0)
    ##  base               * 3.4.0    2017-04-21 local         
    ##  base64enc            0.1-3    2015-07-28 CRAN (R 3.2.1)
    ##  broom              * 0.4.2    2017-02-13 CRAN (R 3.4.0)
    ##  car                  2.1-4    2016-12-02 CRAN (R 3.4.0)
    ##  caret              * 6.0-76   2017-04-18 CRAN (R 3.4.0)
    ##  cellranger           1.1.0    2016-07-27 CRAN (R 3.4.0)
    ##  class                7.3-14   2015-08-30 CRAN (R 3.4.0)
    ##  codetools            0.2-15   2016-10-05 CRAN (R 3.4.0)
    ##  colorspace           1.3-2    2016-12-14 CRAN (R 3.4.0)
    ##  compiler             3.4.0    2017-04-21 local         
    ##  data.table           1.10.4   2017-02-01 CRAN (R 3.4.0)
    ##  datasets           * 3.4.0    2017-04-21 local         
    ##  DBI                  0.6-1    2017-04-01 CRAN (R 3.4.0)
    ##  devtools           * 1.13.2   2017-06-02 CRAN (R 3.4.0)
    ##  digest               0.6.12   2017-01-27 CRAN (R 3.4.0)
    ##  dplyr              * 0.5.0    2016-06-24 CRAN (R 3.4.0)
    ##  e1071                1.6-8    2017-02-02 CRAN (R 3.4.0)
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
    ##  labeling             0.3      2014-08-23 CRAN (R 3.4.0)
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
