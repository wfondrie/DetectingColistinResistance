Make Miscellaneous Figures
================
William E Fondrie

-   [Load Libraries and Prepare Workspace](#load-libraries-and-prepare-workspace)
-   [Load Data](#load-data)
-   [Plot of Resistance Peaks and Their Variation](#plot-of-resistance-peaks-and-their-variation)
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

# import helper functions
source("../R/preProcessSpec.R")
source("../R/extract.R")
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

specList <- preProcessSpec(files, hws = 80) # Same preprocessing

# rm files that contain multiple spectra
multiSpecIdx <- map_lgl(specList, ~ metaData(.)$num > 1)
specList <- specList[!multiSpecIdx]

spec <- map_df(specList, extractSpectra)

feattbl <- spec %>%
  group_by(id) %>%
  do(extractFeatures(., featureVec = unique(features$mz), tol = 1.5))
```

Plot of Resistance Peaks and Their Variation
--------------------------------------------

``` r
# Ab
abFeat <- c("mz2036.4927")
# Kp
kpFeat <- c("mz1957.0536", "mz1973.1961")

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\makeMiscFigures_files/figure-markdown_github-ascii_identifiers/peaks-1.png)

``` r
ggsave("../results/resPeaks.pdf", width = 105, height = 50, unit = "mm", useDingbats = F)
```

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\makeMiscFigures_files/figure-markdown_github-ascii_identifiers/peaks-2.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\makeMiscFigures_files/figure-markdown_github-ascii_identifiers/libraryDistro-1.png)

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

![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\makeMiscFigures_files/figure-markdown_github-ascii_identifiers/libraryDistro-2.png)![](C:\Users\WEF\DetectingColistinResistance\Rmd\mdOutput\makeMiscFigures_files/figure-markdown_github-ascii_identifiers/libraryDistro-3.png)

``` r
ggsave("../results/librarySmall.pdf", height = 100, width = 150, unit = "mm", useDingbats = F)

summary(specInfo$species)
```

    ##                        V1     
    ##  Acinetobacter baumannii:647  
    ##  Klebsiella pneumoniae  :317  
    ##  Pseudomonas aeruginosa :149  
    ##  Serratia marcescens    :111  
    ##  Staphylococcus aureus  :110  
    ##  Enterobacter cloacae   :100  
    ##  (Other)                :663

``` r
length(unique(specInfo$species))
```

    ## [1] 50

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
