# preProcessSpec()

# fileList - A vector of file names to analyze
# hws - Half window size for smoothing

# Output: A list of MALDIquant MassSpectra objects

preProcessSpec <- function(fileList, hws = 100) {
  spec <- importMzXml(fileList, verbose = F)
  spec <- trim(spec, c(1000, 2400))
  #spec <- trim(spec, c(1100, 2100))
  spec <- transformIntensity(spec, method="sqrt")
  spec <- smoothIntensity(spec, method="SavitzkyGolay", halfWindowSize = hws)
  spec <- removeBaseline(spec, method="SNIP", iterations = 60)
  spec <- calibrateIntensity(spec, method = "TIC")
  
  spec <- list(spec[[1]])
  
  return(spec)
}