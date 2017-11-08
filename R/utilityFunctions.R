# A simple rounding function
sigRound <- function(x, digits = 0) {
    format(round(x = x, digits = digits), nsmall = digits)
}


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

