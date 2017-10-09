# Creates a feature list as input for the ML models
#
# trainList - Should be trainIdx.rds. Tells which features belong with each model
# specDf - extracted spectra. The output from extractSpectra()


createNewFeatureTbl <- function(trainList, specDf, summaryDat, mzTol, fileName) {
    suffix <- if(str_detect(trainList$regex, "Acineto")) "_Ab" else "_Kp"
    
    # extract specified features
    cat(paste0("Creating features for ", suffix, "...\n"))
    feattbl <- specDf %>%
        group_by(id) %>%
        do(extractFeatures(., featureVec = trainList$features, tol = mzTol))
    
    saveRDS(feattbl, file = paste0("../temp/", fileName, suffix, ".rds"))
    
    mlFeat <- feattbl %>%
        select(id, feat, relInt) %>%
        group_by(id) %>%
        mutate(relInt = relInt / max(relInt)) %>%
        ungroup() %>%
        spread(key = feat, value = relInt, fill = 0) %>%
        left_join(summaryDat, by = c("id")) 
    
    return(mlFeat)
}