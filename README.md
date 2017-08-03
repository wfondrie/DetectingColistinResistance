Insert Paper Title Here
================
William E Fondrie

Order of analysis:
1. [modelTraining.Rmd](%22./Rmd/modelTraining.md%22) 2. [modelTraining.Rmd](%22./Rmd/modelTraining.md%22) 3. [modelTraining.Rmd](%22./Rmd/modelTraining.md%22) 4. [evaluateModels.Rmd](%22./Rmd/evaluateModels.md%22) 5. [modelTraining.Rmd](%22./Rmd/modelTraining.md%22)

``` r
library(knitr)
# Set global parameters for analysis
mzTol <- 1.5 # m/z tolerance for feature extraction in Da
saveRDS(mzTol, "../temp/mzTol.rds")

# Run analysis scripts
knit("evaluateModels.Rmd")
```

    ## 
    ## 
    ## processing file: evaluateModels.Rmd

    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |.....                                                            |   8%
    ##   ordinary text without R code
    ## 
    ## 
      |                                                                       
      |...........                                                      |  17%
    ## label: loadLibraries (with options) 
    ## List of 1
    ##  $ results: chr "hide"
    ## 
    ## 
      |                                                                       
      |................                                                 |  25%
    ##   ordinary text without R code
    ## 
    ## 
      |                                                                       
      |......................                                           |  33%
    ## label: singles
    ## 
      |                                                                       
      |...........................                                      |  42%
    ##   ordinary text without R code
    ## 
    ## 
      |                                                                       
      |................................                                 |  50%
    ## label: AbSingleSpec
    ## 
      |                                                                       
      |......................................                           |  58%
    ##   ordinary text without R code
    ## 
    ## 
      |                                                                       
      |...........................................                      |  67%
    ## label: KpSingleSpec
    ## 
      |                                                                       
      |.................................................                |  75%
    ##   ordinary text without R code
    ## 
    ## 
      |                                                                       
      |......................................................           |  83%
    ## label: AbPR

    ## 
      |                                                                       
      |............................................................     |  92%
    ##   ordinary text without R code
    ## 
    ## 
      |                                                                       
      |.................................................................| 100%
    ## label: KpPR

    ## output file: evaluateModels.md

![](README_files/figure-markdown_github/runAnalysis-1.png)

    ## [1] "evaluateModels.md"

``` r
# Copy GitHub markdown file to top directory
file.copy("README.md", "..", overwrite=T)
```

    ## Warning in file.copy("README.md", "..", overwrite = T): problem copying .
    ## \README.md to ..\README.md: No such file or directory

    ## [1] FALSE
