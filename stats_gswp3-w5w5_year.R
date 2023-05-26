options(warnPartialMatchDollar = TRUE)
source("netCDFfunctions.R") # Make sure the working dir is set correctly.

rootDir <- 
  if (dir.exists("//clusterfs.pik-potsdam.de/projects")) {
    "//clusterfs.pik-potsdam.de" 
  } else if (dir.exists("/p/projects")) {
    "/p"
  }




# Masks -------------------------------------------------------------

isimipInputPath <- file.path(
  rootDir,
  "projects/isimip/isimip/ISIMIP3a/InputData")

areaPath <- "."
areaFile <- "cellArea.nc"
area <- loadWeights(file.path(areaPath, areaFile), normalize = FALSE, rescaleOcean = FALSE)

popPath <- file.path(isimipInputPath, "socioeconomic/pop/histsoc")
popFile <- "population_histsoc_30arcmin_annual_1901_2021.nc"
pop <- loadGlobalYearlyNc(file.path(popPath, popFile), "total-population")

varPath <- file.path(
  isimipInputPath, 
  "climate/atmosphere/obsclim/global/daily/historical/GSWP3-W5E5")
allNcFiles <- dir(
  varPath, 
  pattern = paste0("^gswp3-w5e5_obsclim_.+_global_daily_\\d{4}_\\d{4}\\.nc$"))
allVarNames <-
  substr(
    allNcFiles, 
    start = nchar("gswp3-w5e5_obsclim_") + 1, 
    stop = nchar(allNcFiles)-nchar("_global_daily_0000_0000.nc")) |> 
  unique()




# Statistics -------------------------------------------------------------

# Apply these functions for the daily time series of each year and grid cell.
statsBasics <- list(
  N = length,
  Na = \(x) sum(is.na(x)),
  Mean = \(x) mean(x, na.rm=TRUE),
  Min = \(x) min(x, na.rm=TRUE),
  Max = \(x) max(x, na.rm=TRUE),
  Sd = \(x) sqrt(mean((x-mean(x, na.rm=TRUE))^2, na.rm=TRUE)))
probs <- seq(0.05, 0.95, 0.05)
statsQuantiles <- lapply(
  probs,
  \(p) \(x) quantile(x, p, na.rm = TRUE, names = FALSE))
names(statsQuantiles) <- sprintf("Q%02d", as.integer(probs*100))
stats <- c(statsBasics, statsQuantiles)




# Calculations -------------------------------------------------------------

calculateYearlyAggregationAndWriteAsCsv <- function(
    varName, varPath, stats, aggregationMask, population, cellArea, outName
) {
  
  cat(paste0("start calculateYearlyAggregationAndWriteAsCsv: ", varName, "\n"))
  
  files <- dir(
    varPath, 
    pattern = paste0("gswp3-w5e5_obsclim_",varName,"_global_daily_\\d{4}_\\d{4}\\.nc"))
  
  varStats <- stats
  names(varStats) <- paste0(varName, names(stats))
  
  res <- lapply(
    files,
    \(file) aggregateYearly(
      file.path(varPath, file),
      aggregationMask = aggregationMask, 
      population = population,
      cellArea = cellArea,
      funs = varStats,
      isNaFun = \(x) !is.finite(x) | x > 1e19 # 1e20 is officially NA
    ))
  
  cat(paste0("write calculateYearlyAggregationAndWriteAsCsv: ", varName, "\n"))
  
  output <- 
    bind_rows(res) |> 
    unite(name, c(fun, weight)) |> 
    pivot_wider()
  readr::write_csv(output, paste0(varName, outName, ".csv"))
  
  cat(paste0("end calculateYearlyAggregationAndWriteAsCsv: ", varName, "\n"))
}
