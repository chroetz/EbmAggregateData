library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(lubridate)


#' Load netCDF file of dimensions lon, lat, time
#' @param filePath path to a netCDF data file
#' @return A list of following objects.
#' \itemize{
#'   \item nc: object from RNetCDF package
#'   \item lon: longitude values,
#'   \item lat: latitude values,
#'   \item time: time (days) values,
#'   \item days: table of days in the time dimension,
#'   \item months: table of months in the time dimension,
#'   \item years: table of years in the time dimension
#' }
#' @examples
#' \dontrun{
#' data <- loadGlobalDailyNc("path/to/data.nc")
#' RNetCDF::close.nc(data$nc) # Need to do this after work is done.
#' }
loadGlobalDailyNc <- function(filePath) {
  
  stopifnot(is.character(filePath), length(filePath) == 1)
  
  nc <- RNetCDF::open.nc(filePath)
  ncInfo <- RNetCDF::file.inq.nc(nc)
  stopifnot(ncInfo$ndims == 3)
  dimList <- list(
    RNetCDF::dim.inq.nc(nc, 0),
    RNetCDF::dim.inq.nc(nc, 1),
    RNetCDF::dim.inq.nc(nc, 2))
  names(dimList) <- c(dimList[[1]]$name, dimList[[2]]$name, dimList[[3]]$name)
  stopifnot(sort(names(dimList)) == c("lat", "lon", "time"))
  dimLonValues <- RNetCDF::var.get.nc(nc, "lon")
  dimLatValues <- RNetCDF::var.get.nc(nc, "lat")
  dimTimeValues <- RNetCDF::var.get.nc(nc, "time")
  orderLon <- order(dimLonValues, decreasing=FALSE)
  orderLat <- order(dimLatValues, decreasing=TRUE)
  rightLonOrder <- all(orderLon == seq_len(dimList[["lon"]]$length))
  rightLatOrder <- all(orderLat == seq_len(dimList[["lat"]]$length))
  if (!rightLonOrder) dimLonValues <- dimLonValues[orderLon]
  if (!rightLatOrder) dimLatValues <- dimLatValues[orderLat]

  varInfo <- RNetCDF::var.inq.nc(nc, 3) # First variable of dims.
  
  # TODO: Would need to permute if either of the following assertions does not hold.
  stopifnot(varInfo$dimids == c(dimList[["lon"]]$id, dimList[["lat"]]$id, dimList[["time"]]$id)) 
  stopifnot(rightLonOrder)
  stopifnot(rightLatOrder)
  
  timeUnits <- RNetCDF::att.get.nc(nc, "time", "units")
  
  startDate <- str_match(
    timeUnits, 
    "^days since (\\d{4}-\\d{1,2}-\\d{1,2})")[2]
  stopifnot(!is.na(startDate))
  
  dayTable <- 
    tibble(
      index = seq_along(dimTimeValues),
      daysSinceStart = dimTimeValues,
      date = ymd(startDate) + days(daysSinceStart),
      month = month(date),
      year = year(date),
      day = day(date))
  monthTable <- 
    dayTable |> 
    group_by(year, month) |> 
    summarize(startIndex = min(index), nDays = n(), .groups = "drop")
  yearTable <- 
    dayTable |> 
    group_by(year) |> 
    summarize(startIndex = min(index), nDays = n())
  
  return(list(
    nc = nc,
    lon = dimLonValues,
    lat = dimLatValues,
    time = dimTimeValues,
    days = dayTable,
    months = monthTable,
    years = yearTable))

}


loadGlobalYearlyNc <- function(filePath, varName) {
  
  stopifnot(is.character(filePath), length(filePath) == 1)
  
  nc <- RNetCDF::open.nc(filePath)
  ncInfo <- RNetCDF::file.inq.nc(nc)
  stopifnot(ncInfo$ndims == 3)
  dimList <- list(
    RNetCDF::dim.inq.nc(nc, 0),
    RNetCDF::dim.inq.nc(nc, 1),
    RNetCDF::dim.inq.nc(nc, 2))
  names(dimList) <- c(dimList[[1]]$name, dimList[[2]]$name, dimList[[3]]$name)
  stopifnot(sort(names(dimList)) == c("lat", "lon", "time"))
  dimLonValues <- RNetCDF::var.get.nc(nc, "lon")
  dimLatValues <- RNetCDF::var.get.nc(nc, "lat")
  dimTimeValues <- RNetCDF::var.get.nc(nc, "time")
  orderLon <- order(dimLonValues, decreasing=FALSE)
  orderLat <- order(dimLatValues, decreasing=TRUE)
  rightLonOrder <- all(orderLon == seq_len(dimList[["lon"]]$length))
  rightLatOrder <- all(orderLat == seq_len(dimList[["lat"]]$length))
  if (!rightLonOrder) dimLonValues <- dimLonValues[orderLon]
  if (!rightLatOrder) dimLatValues <- dimLatValues[orderLat]

  varInfo <- RNetCDF::var.inq.nc(nc, varName)
  values <- RNetCDF::var.get.nc(nc, varName)
  
  # TODO: Would need to permute if either of the following assertions does not hold.
  stopifnot(varInfo$dimids == c(dimList[["lon"]]$id, dimList[["lat"]]$id, dimList[["time"]]$id)) 
  stopifnot(rightLonOrder)
  stopifnot(rightLatOrder)
  
  dim(values) <- c(prod(dim(values)[1:2]), dim(values)[3])
  
  timeUnits <- RNetCDF::att.get.nc(nc, "time", "units")
  
  RNetCDF::close.nc(nc)
 
  startDate <- str_match(
    timeUnits, 
    "^days since (\\d{4}-\\d{1,2}-\\d{1,2})")[2]
  stopifnot(!is.na(startDate))
  
  yearTable <- 
    tibble(
      index = seq_along(dimTimeValues),
      daysSinceStart = dimTimeValues,
      date = ymd(startDate) + days(daysSinceStart),
      day = day(date),
      month = month(date),
      year = year(date))
  
  stopifnot(
    yearTable$day == 1,
    yearTable$month == 1)
  
  return(list(
    nc = nc,
    lon = dimLonValues,
    lat = dimLatValues,
    time = dimTimeValues,
    years = yearTable |> select(-c(day, month)),
    values = values))
}

#' Load a weight matrix from a netCDF mask array file
#'
#' @param filePath Path to a netCDF file with dimension `lon` and `lat`.
#' @param normalize Should the weights be normalized to colSum to 1?
#' @param rescaleOcean Should the weights be normalized to rowSum to 1? Done
#'   before normalize.
#' @param varSubset `NULL` or an integer vector of indices of the variables to
#'   be read (starting at 1).
#' @return A list of following objects.
#' \itemize{
#'   \item lon: longitude values (increasing order),
#'   \item lat: latitude values (decreasing order),
#'   \item weightMatrix: A matrix of with one column for each variable in the netCDF file and number of rows equal to the product of number of lon and lat values.
#' }
loadWeights <- function(
    filePath, 
    normalize = FALSE, 
    varSubset = NULL, 
    rescaleOcean = FALSE, 
    removePattern = "world"
) {
  
  stopifnot(is.character(filePath), length(filePath) == 1)
  
  nc <- RNetCDF::open.nc(filePath)
  ncInfo <- RNetCDF::file.inq.nc(nc)
  stopifnot(ncInfo$ndims == 2)
  dim0 <- RNetCDF::dim.inq.nc(nc, 0)
  dim1 <- RNetCDF::dim.inq.nc(nc, 1)
  stopifnot(sort(c(dim0$name, dim1$name)) == c("lat", "lon"))
  if (dim0$name == "lon") {
    dimLon <- dim0
    dimLat <- dim1
  } else {
    dimLon <- dim1
    dimLat <- dim0
  }
  dimLonValues <- RNetCDF::var.get.nc(nc, dimLon$id)
  dimLatValues <- RNetCDF::var.get.nc(nc, dimLat$id)
  orderLon <- order(dimLonValues, decreasing=FALSE)
  orderLat <- order(dimLatValues, decreasing=TRUE)
  rightLonOrder <- all(orderLon == seq_len(dimLon$length))
  rightLatOrder <- all(orderLat == seq_len(dimLat$length))
  if (!rightLonOrder) {
    warning("Lon values seem to have the wrong order. Reordering. Please check that this is correct.")
    dimLonValues <- dimLonValues[orderLon]
  }
  if (!rightLatOrder) {
    dimLatValues <- dimLatValues[orderLat]
    warning("Lat values seem to have the wrong order. Reordering. Please check that this is correct.")
  }
  proto <- matrix(NA_real_, nrow = dimLon$length, dimLat$length)
  indices <- if (is.null(varSubset)) {
      seq_len(ncInfo$nvars - 2) 
    } else {
      intersect(seq_len(ncInfo$nvars - 2), varSubset)
    }
  weightMatrix <- 
    vapply(
      indices,
      \(i) {
        varInfo <- RNetCDF::var.inq.nc(nc, i + 1)
        values <- RNetCDF::var.get.nc(nc, i + 1)
        if (all(varInfo$dimids == c(dimLon$id, dimLat$id))) {
        } else if (all(varInfo$dimids == c(dimLon$id, dimLat$id))) {
          values <- t(values)
        } else {
          stop(paste0("Cannot deal with dimids ", paste0(dimids, collapse=",")))
        }
        stopifnot(dim(values) == c(dimLon$length, dimLat$length))
        if (!rightLonOrder) values <- values[orderLon, ]
        if (!rightLatOrder) values <- values[, orderLat]
        return(values)
      },
      FUN.VALUE = proto)
  varNames <- 
    vapply(
      indices, 
      \(i) RNetCDF::var.inq.nc(nc, i + 1)$name,
      FUN.VALUE = "")
  RNetCDF::close.nc(nc)
  
  dim(weightMatrix) <- c(prod(dim(weightMatrix)[1:2]), dim(weightMatrix)[3])
  colnames(weightMatrix) <- varNames
  
  for (regex in removePattern) {
    sel <- !grepl(regex, colnames(weightMatrix))
    weightMatrix <- weightMatrix[, sel]
  }
  
  if (rescaleOcean) {
    weightMatrix <- sweep(weightMatrix, 1, rowSums(weightMatrix), "/")
    weightMatrix[is.na(weightMatrix)] <- 0
  }
  
  if (normalize) {
    weightMatrix <- sweep(weightMatrix, 2, colSums(weightMatrix), "/")
    weightMatrix[is.na(weightMatrix)] <- 0
  }
  
  return(list(
    weightMatrix = weightMatrix,
    lon = dimLonValues,
    lat = dimLatValues))
}

#' Apply a function to a slice of an netCDF data with weights
#'
#' @param start Start index of the time period.
#' @param count Number of consecutive entries in the time period.
#' @param nc An netCDF data object from the RNetCDF package.
#' @param funs A list of functions objects. Will be applied to the vector of
#'   values of the variable `varName` form the time period specified by `start`
#'   and `count`.
#' @param weightList A list of matrices of row dimension equal to the product of
#'   first two dimensions of the array `nc`. Used to aggregate the function
#'   values by taking a weighted mean.
#' @return A vector of length `ncol(weightMatrix)` containing the aggregated
#'   function values. The vector is named respectively if `weightMatrix` has
#'   column names.
#' @param varName The variable name to be read form `nc`. If set to `NA` and
#'   `nc` contains only one variable, then this variable is used.
applyPeriod <- function(start, count, nc, funs, weightList, varName = NA, isNaFun = \(x) !is.finite(x)) {
  
  data <- RNetCDF::var.get.nc(
    nc, 
    if (is.na(varName)) 3 else varName, 
    start = c(1, 1, start),
    count = c(NA, NA, count))
  
  data[isNaFun(data)] <- NA
        
  nCells <- prod(dim(data)[1:2])
  dim(data) <- c(nCells, count)
  
  results <- 
    sapply(
      funs, 
      \(fun) {
        valuesGrid <- apply(data, 1, fun)
        sapply(weightList, \(weightMatrix) {
          valuesOnMask <- as.vector(crossprod(valuesGrid, weightMatrix))
          names(valuesOnMask) <- colnames(weightMatrix)
          return(valuesOnMask)
        })},
      simplify = "array")
  
  return(results)
}

getMaskedPopWeights <- function(population, year, mask) {
  yearIndex <- population$years |> filter(.data$year == .env$year) |> pull(index)
  popWeights <- population$values[,yearIndex]
  popWeights[is.na(popWeights)] <- 0 # `population$values` is NA where the used country mask says there is no country. But a different country mask may be used as `mask`. To avoid NA-results, these values are set to 0, which probably is not far from the truth. 
  popWeightsMasked <- mask$weightMatrix * popWeights
  popWeightsMasked <- sweep(popWeightsMasked, 2, colSums(popWeightsMasked), "/")
  return(popWeightsMasked)
}

aggregateYearlyAverageMonth <- function(filePath, aggregationMask, population, cellArea, funs) {
  
  data <- loadGlobalDailyNc(filePath)
  
  stopifnot(
    data$lon == aggregationMask$lon,
    data$lon == population$lon)
  stopifnot(
    data$lat == aggregationMask$lat,
    data$lat == population$lat)
  
  # get area weights
  areaWeightsMasked <- aggregationMask$weightMatrix * cellArea$weightMatrix
  areaWeightsMasked[areaWeightsMasked<0]  <- NA
  areaWeightsMasked <- sweep(areaWeightsMasked, 2, colSums(areaWeightsMasked), "/")
  
  pt <- proc.time()
  values <- sapply(
    seq_len(nrow(data$months)),
   function(i) {
      metaMonth <- data$months[i, ]
      
      cat(str_glue("{i}/{nrow(data$months)}: {metaMonth$month}/{metaMonth$year}"), "\n")
      
      popWeightsMasked <- getMaskedPopWeights(
        population, 
        year = metaMonth$year, 
        mask = aggregationMask)
      
      dataMonth <- applyPeriod(
        metaMonth$startIndex, 
        metaMonth$nDays, 
        data$nc, 
        funs, 
        list(
          grid = aggregationMask$weightMatrix, 
          area = areaWeightsMasked,
          pop = popWeightsMasked))
      
      return(dataMonth)
    },
    simplify = "array"
  )
  print(proc.time() - pt)
  
  RNetCDF::close.nc(data$nc)
  
  monthlyData <-
    as_tibble(values, rownames="iso") |> 
    pivot_longer(-iso) |> 
    separate(name, into = c("weight", "fun", "index"), convert=TRUE) |> 
    left_join(
      data$months |> rowid_to_column("index") |> select(year, month, index), 
      by = "index")
  
  yearlyData <-
    monthlyData |>
    summarise(value = mean(value), .by=c(iso, year, weight, fun)) |> 
    mutate(iso = substring(iso, 3))
  
  return(yearlyData)
}


aggregateYearly <- function(filePath, aggregationMask, population, cellArea, funs, isNaFun = \(x) !is.finite(x)) {
  
  cat("aggregateYearly()\n")
  
  cat("  load file\n")
  data <- loadGlobalDailyNc(filePath)
  
  stopifnot(
    data$lon == aggregationMask$lon,
    data$lon == population$lon,
    data$lon == cellArea$lon)
  stopifnot(
    data$lat == aggregationMask$lat,
    data$lat == population$lat,
    data$lat == cellArea$lat)
  
  # get area weights
  cat("  sweep\n")
  areaWeightsMasked <- sweep(aggregationMask$weightMatrix, 1, cellArea$weightMatrix, "*")
  areaWeightsMasked <- sweep(areaWeightsMasked, 2, colSums(areaWeightsMasked), "/")
  
  cat("  calculate values\n")
  pt <- proc.time()
  values <- sapply(
    seq_len(nrow(data$years)),
    function(i) {
      metaYear <- data$years[i, ]
      
      cat(str_glue("{i}/{nrow(data$years)}: {metaYear$year}"), "\n")
      
      popWeightsMasked <- getMaskedPopWeights(
        population, 
        year = metaYear$year, 
        mask = aggregationMask)
      
      dataYear <- applyPeriod(
        metaYear$startIndex, 
        metaYear$nDays, 
        data$nc, 
        funs, 
        list(
          grid = aggregationMask$weightMatrix, 
          area = areaWeightsMasked,
          pop = popWeightsMasked),
        isNaFun = isNaFun)
      
      return(dataYear)
    },
    simplify = "array"
  )
  print(proc.time() - pt)
  
  RNetCDF::close.nc(data$nc)
  
  cat("  rearrange values\n")
  yearlyData <-
    as_tibble(values, rownames="iso") |> 
    pivot_longer(-iso) |> 
    separate(name, into = c("weight", "fun", "index"), convert=TRUE) |> 
    left_join(
      data$years |> rowid_to_column("index") |> select(year, index), 
      by = "index") |> 
    mutate(iso = substring(iso, 3)) |> 
    select(-index) |> 
    relocate(iso, year, fun, weight, value) |> 
    arrange(iso, year, fun, weight)

  cat("  return\n")
  return(yearlyData)
}



