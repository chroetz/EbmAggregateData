source("stats_gswp3-w5w5_year.R") # Make sure the working dir is set correctly.

if (!interactive()) {
  args <- commandArgs(TRUE)
  varName <- args[1]
} else {
  varName <- allVarNames[1]
}

cat("Variable:", varName, "\n")
stopifnot(varName %in% allVarNames)

doseMaskPath <- file.path(rootDir, "projects/ebm/aggregated_impacts/DOSE/fromShapefile")
doseMaskFile <- "regionMasksDose30arcmin.nc"

maxLen <- 500
nParts <- 7 # 3000 < # DOSE regions < 3500; TODO: should be loaded from data.

for (partNr in seq_len(nParts)) {
  
  cat("Mask Part Nr", partNr, "\n")
  
  mask <- loadWeights(
    file.path(doseMaskPath, doseMaskFile), 
    normalize = TRUE,
    rescaleOcean = TRUE,
    varSubset = seq_len(maxLen) + (partNr-1)*maxLen)
  
  calculateYearlyAggregationAndWriteAsCsv(
    varName,
    varPath = varPath, 
    stats = stats, 
    aggregationMask = mask, 
    population = pop, 
    cellArea = area,
    outName = paste0("DoseYear_", partNr))
}

