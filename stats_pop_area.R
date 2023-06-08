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
area <- loadWeights(
  file.path(areaPath, areaFile), 
  normalize = FALSE,
  rescaleOcean = FALSE)

popPath <- file.path(isimipInputPath, "socioeconomic/pop/histsoc")
popFile <- "population_histsoc_30arcmin_annual_1901_2021.nc"
pop <- loadGlobalYearlyNc(file.path(popPath, popFile), "total-population")

countryMaskPath <- file.path(isimipInputPath, "geo_conditions/countrymasks")
countryMaskFile <- "countrymasks_fractional.nc"

mask <- loadWeights(
  file.path(countryMaskPath, countryMaskFile), 
  normalize = FALSE,
  rescaleOcean = TRUE)

countryAreaMat <- crossprod(mask$weightMatrix, area$weightMatrix)
rownames(countryAreaMat) <- stringr::str_remove(rownames(countryAreaMat), "^m_")
colnames(countryAreaMat) <- "area"
countryArea <- as_tibble(countryAreaMat, rownames = "iso")

readr::write_csv(countryArea, "countryArea.csv")

popValues <- pop$values
popValues[is.na(popValues)] <- 0
countryPopMat <- crossprod(mask$weightMatrix, popValues)
colnames(countryPopMat) <- as.character(pop$years$year)
rownames(countryPopMat) <- stringr::str_remove(rownames(countryPopMat), "^m_")
countryPop <- 
  as_tibble(countryPopMat, rownames = "iso") |> 
  pivot_longer(-iso, names_to="year", values_to="pop", names_transform=list(year = as.integer))

readr::write_csv(countryPop, "countryPop.csv")
