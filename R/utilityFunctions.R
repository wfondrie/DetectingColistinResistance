# A simple rounding function
sigRound <- function(x, digits = 0) {
    format(round(x = x, digits = digits), nsmall = digits)
}


# nom2feat - selects the closest feature greater than or equal to the provided 
# nominal mass
# Input: 
#    mz - the nominal m/z of the ion you want to select
#    feat - a vector of features m/z values
#    maxDiff - The maximum difference between the feature and mz.
nom2feat <- function(mz, feat, maxDiff = 4) {
    feat <- as.numeric(feat)
    feat <- feat[which(feat - mz >= 0 & feat - mz <= maxDiff)]
    feat <- feat[which(abs(feat - mz) == min(abs(feat - mz)))]
    
    return(feat)
}

# Require PRROC package --------------------------------------------------------

# Make ROC curves
# Returns a list of resistance and species-level ROC curves
# Input:
#   dat - a dataframe containing a minimum of 3 named columns:
#        -> truth - c("pos", "neg", "other")
#        -> pos - Score for colistin resistance (truth == "pos")
#        -> speciesVsOther - Score for species (truth != "other)

makeRocCurves <- function(dat) {
    posRoc <- roc.curve(dat$pos[dat$truth == "pos"], 
                        dat$pos[dat$truth != "pos"], curve = T)
    
    speciesRoc <- roc.curve(dat$speciesVsOther[dat$truth != "other"], 
                            dat$speciesVsOther[dat$truth ==  "other"], curve = T)
    
    ret <- list(pos = posRoc,
                species = speciesRoc)
    
    return(ret)
}

# Make PR curves
# Returns a list of resistance and species-level PR curves
# Input:
#   dat - a dataframe containing a minimum of 3 named columns:
#        -> truth - c("pos", "neg", "other")
#        -> pos - Score for colistin resistance (truth == "pos")
#        -> speciesVsOther - Score for species (truth != "other)

makePrCurves <- function(dat) {
    posPr <- pr.curve(dat$pos[dat$truth == "pos"], 
                       dat$pos[dat$truth != "pos"], curve = T)
    
    speciesPr <- pr.curve(dat$speciesVsOther[dat$truth != "other"], 
                           dat$speciesVsOther[dat$truth ==  "other"], curve = T)
    
    ret <- list(pos = posPr,
                species = speciesPr)
    
    return(ret)
}

# calcAUC 
# Returns a list of AUC for PR and ROC curves. Useful for Boostrapping
calcAUC <- function(dat) {
        posPr <- pr.curve(dat$pos[dat$truth == "pos"], 
                          dat$pos[dat$truth != "pos"])
        
        speciesPr <- pr.curve(dat$speciesVsOther[dat$truth != "other"], 
                              dat$speciesVsOther[dat$truth ==  "other"])
        
        posRoc <- roc.curve(dat$pos[dat$truth == "pos"], 
                            dat$pos[dat$truth != "pos"])
        
        speciesRoc <- roc.curve(dat$speciesVsOther[dat$truth != "other"], 
                                dat$speciesVsOther[dat$truth ==  "other"])
        
        data.frame(type = c("ROC", "ROC", "PR", "PR"),
                   level = c("Species", "Res", "Species", "Res"),
                   AUC = c(speciesRoc$auc,
                           posRoc$auc,
                           speciesPr$auc.integral,
                           posPr$auc.integral)) 
}
