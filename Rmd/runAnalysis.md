R Notebook
================

Introduction
============

Prepare Workspace
=================

Importing Data
==============

``` r
files <- list.files("data/fullLib", 
                    full.names = T, 
                    pattern = "mzXML$",
                    recursive = T)
files <- sample(files, 50) # uncomment for testing

specInfo <- tibble(fname = files) %>%
    mutate(type = str_match(fname, "([^\\^/]+)[\\/][^\\^/]+mzXML$")[ , 2],
           id = str_match(fname, "([^\\^/]+).mzXML$")[ , 2],
           Ab = ifelse(str_detect(type, "Acinetobacter baumannii - res"), "pos", "other"),
           Kp = ifelse(str_detect(type, "Klebsiella pneumoniae - res"), "pos", "other"),
           Ab = as.factor(ifelse(str_detect(type, "Acinetobacter baumannii - sen"), "neg", Ab)),
           Kp = as.factor(ifelse(str_detect(type, "Klebsiella pneumoniae - sen"), "neg", Kp)))

trainIdx <- list(Ab = createDataPartition(specInfo$Ab, p = 0.6, list = F),
                 Kp = createDataPartition(specInfo$Kp, p = 0.6, list = F))
```
