prepareData <- function(trainDatList, orgLabs) {
    dmlist <- map(orgLabs, function(o){
        dat <- trainDatList[[o]]
        
        y <- dat[ , o] %>% mutate(encoding = as.numeric(.[[o]]) - 1)
        
        x <- as.matrix(select(dat, starts_with("mz")))
        
        return(list(encoding = y,
                    dat = xgb.DMatrix(x, label = y$encoding)))
    })
    
    names(dmlist) <- orgLabs
    return(dmlist)
    
} 