source("stats_gswp3-w5w5_year.R") # Make sure the working dir is set correctly.

if (!interactive()) {
  args <- commandArgs(TRUE)
  varName <- args[1]
} else {
  varName <- allVarNames[1]
}

cat("Variable:", varName, "\n")
stopifnot(varName %in% allVarNames)

countryMaskPath <- file.path(isimipInputPath, "geo_conditions/countrymasks")
countryMaskFile <- "countrymasks_fractional.nc"

mask <- loadWeights(
  file.path(countryMaskPath, countryMaskFile), 
  rescaleOcean = TRUE,
  normalize = TRUE)

calculateYearlyAggregationAndWriteAsCsv(
  varName,
  varPath = varPath, 
  stats = stats, 
  aggregationMask = mask, 
  population = pop, 
  cellArea = area,
  outName = "CountryYear")
