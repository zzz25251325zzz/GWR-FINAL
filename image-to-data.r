#install.packages("raster")
#install.packages("rgdal")
library(raster)
library(rgdal)

combineDataWithCoord <- function(dataFile, coordFile, gtColumn) {
  gt <- read.csv(coordFile)
  grx <- glob2rx("1/*")
  ind <- grep(grx, gt$time)
  jandata <- gt[ind,]
  data <- read.csv(dataFile)
  data[data == -9999] <- NA
  datenum <- sapply(strsplit(as.character(jandata$time), "/"), "[", 2)
  datenum <- as.double(datenum)
  jandata[, "datenum"] <- datenum
  subjan <- jandata[,c("time", "lat", "lon", "name", gtColumn, "datenum")]
  subjan[subjan=="-"] <- NA
  alldata <- merge(x = subjan, y = data, by.x = c("name", "datenum"), by.y = c("name", "day"), all = TRUE)
  return(alldata)
}

getDataFromImage <- function(folders, daynum) {
  for (f in folders) {
    im_name <- sprintf("%s/test_im_%03d.tif", f, daynum)
    im_image <- raster(im_name)
    im_data <- rasterToPoints(im_image, spatial = T)
    im_frame <- as.data.frame(im_data)
    colnames(im_frame)[1] = f
    assign(f, im_frame)
  }
  mergeall <- get(folders[1])
  for (f in folders[2:length(folders)]) {
    mergeall <- merge(mergeall, get(f), by = c("x", "y"), all = T)
  }
  return(mergeall)
}

makeNewImage <- function(formula, data, fitData, weight, adapt) {
    trainData.bw <- gwr.sel(formula, data, verbose = F, adapt = adapt, gweight = weight)
    if (adapt) {
        trainData.gwr <- gwr(formula, data, hatmatrix = T, gweight = weight, adapt = trainData.bw)
    } else {
        trainData.gwr <- gwr(formula, data, hatmatrix = T, gweight = weight, bandwidth = trainData.bw)        
    }
    fitData.gwr <- gwr(formula, data = data, bandwidth = trainData.gwr$bandwidth, predictions = T,
                       adapt = trainData.gwr$adapt, fit.points = fitData, se.fit = T,
                       fittedGWRobject = trainData.gwr)
    fitData$pred <- fitData.gwr$SDF$pred
    return(fitData)
}

convertToSpatial <- function(frame, lat, lon, projstr) {
  points <- cbind(frame[,lon], frame[,lat])
  colnames(points) <- c("LONG", "LAT")
  frame.data <- frame[setdiff(names(frame), c(lat, lon))]
  frame.sp <- SpatialPointsDataFrame(coords = points, 
                                     data = as.data.frame(frame.data),
                                     proj4string = projstr)
  return(frame.sp)
}

saveImage <- function(spdata, name) {
  gridded(spdata) <- T
  writeGDAL(spdata, name)
}
