library(tidyverse)
library(stringr)

# extractPeaks()

# massPeaks - a MassPeaks object from MALDIquant
# massSpec  - a MassSpectra object from MALDIquant

# Output - A tibble with 4 columns:
#   type   - the species or resistance + species of the spectrum
#   mz     - the m/z
#   relInt - the relative intensity at the specified m/z, relative to the 
#            base peak
#   id     - a unique identifier for the spectrum

extractPeaks <- function(massPeaks, massSpec) {
  feat <- as_tibble(intensityMatrix(massPeaks, massSpec))
  fname <- map_chr(massPeaks, ~ metaData(.)$file)
  type <- str_match(fname, "([^\\\\^/]+)[\\\\/][^\\\\^/]+mzXML$")[ , 2]
  id <- str_match(fname, "([^\\\\^/]+).mzXML$")[ , 2]
  
  feat$id <- as.factor(id)
  feat$type <- as.factor(type)
  
  feat <- feat %>%
    gather(mz, relInt, -id, -type) %>%
    group_by(type, id) %>% 
    mutate(relInt = relInt / max(relInt),
           mz = as.numeric(mz)) %>%
    ungroup()
  
  return(feat)
}

#-------------------------------------------------------------------------------
# extractSpectra()

# massSpec - a MassSpectra object from MALDIquant

# Output - A long formated tibble with 4 columns:
#   type   - the species or resistance + species of the spectrum
#   mz     - the m/z
#   relInt - the relative intensity at the specified m/z, relative to the 
#            base peak
#   id     - a unique identifier for the spectrum

extractSpectra <- function(massSpec) {
  s <- tibble(mz = as.numeric(massSpec@mass),
              relInt = as.numeric(massSpec@intensity))
  
  fname <- metaData(massSpec)$file
  type <- str_match(fname, "([^\\\\^/]+)[\\\\/][^\\\\^/]+mzXML$")[ , 2]
  id <- str_match(fname, "([^\\\\^/]+).mzXML$")[ , 2]
  
  s$type <- type
  s$id <- id
  
  s <- s %>%
    group_by(type, id) %>% 
    mutate(relInt = relInt / max(relInt),
           mz = mz) %>%
    ungroup()
  
  return(s)
}

#-------------------------------------------------------------------------------
# extractFeatures()

# Input:
#   specDf     - a mass spectrum in dataframe format, such as a single spectrum 
#                from the output of "extractSpectra()"
#   featureVec - a numeric vector of mz values to extract
#   tol        - window around specified mz values. For example, 2.5 specifies a 
#                5 m/z window. The function extracts the maximum intensity in this
#                window

# Output:  A long formated tibble with 6 columns:
#   mz      - the selected feature m/z
#   relInt  - the relative intensity at the specified m/z, relative to the 
#             base peak
#   type    - the species or resistance + species of the spectrum
#   id      - a unique identifier for the spectrum
#   real_mz - The true m/z of the extracted intensity
#   feat    - A character vector version of "mz", suitable for column names

extractFeatures <- function(specDf, featureVec, tol) {
  featureVec <- unique(featureVec)
  map_df(featureVec, function(feat) {
    win <- filter(specDf, mz <= feat + tol, mz >= feat - tol)
    m <- filter(win, relInt == max(relInt)) %>% 
        filter(mz == median(mz)) %>%
        mutate(real_mz = mz,
               mz = feat,
               feat = paste0("mz", round(mz, 4)))
    return(m)
  })
}
