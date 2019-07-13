library(spgwr)
library(Metrics)

setwd("~/Desktop/GWR")
source("image-to-data.r")
source("metrics.r")

tempdata <- combineDataWithCoord("Temp_data.csv", "obs_20140101_20141231.csv", "temp")
tempdata$average <- rowMeans(tempdata[, 7:11], na.rm = TRUE)
tempsub <- tempdata[, c("name", "datenum", "lat", "lon", "temp", "average")]
tempsub <- tempsub[complete.cases(tempsub), ]

alldays <- unique(tempsub$datenum)
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs")
tempdata.sp <- convertToSpatial(tempsub, "lat", "lon", proj)

tempFormula <- temp ~ average
numFolds <- 5
weightsFun <- c(gwr.Gauss, gwr.bisquare, gwr.tricube)
names(weightsFun) <- c("Gauss", "Bisquare", "Tricube")
adaptOpt <- c(T, F)
names(adaptOpt) <- c("Adaptive", "Fixed")

set.seed(1)
resultAll <- matrix(list(), 3, 2)
rownames(resultAll) <- names(weightsFun)
colnames(resultAll) <- names(adaptOpt)
for (i in rownames(resultAll)) {
  for (j in colnames(resultAll)) {
    resultAll[[i, j]] <- array(dim = c(3, length(alldays)))
  }
}

for (i in alldays) {
    partData <- tempdata.sp[tempdata.sp@data$datenum == i, ]
    print(i)
    for (w in names(weightsFun)) {
        for (a in names(adaptOpt)) {
            pref <- sprintf("%s & %s", a, w)
            err <- list(NA, NA, NA)
            while (any(is.na(err))) {
                err <- gwrCV(tempFormula, partData, weightsFun[[w]], adaptOpt[[a]], numFolds, "temp")
                if (any(is.na(err))) { print("NA") }
            }
            resultAll[[w, a]][,i] <- unlist(err, use.names = F)
            str <- sprintf("%s & %.3f & %.3f & %.3f \\", pref, err[1], err[2], err[3])
            print(str)
        }
    }
}

print("All days")
for (w in names(weightsFun)) {
  for (a in names(adaptOpt)) {
    pref <- sprintf("%s & %s", a, w)
    aver <- rowMeans(resultAll[[w, a]])
    #var <- rowVar(resultAll[[w, a]])
    str <- sprintf("%s & %.3f & %.3f & %.3f \\", pref, aver[1], aver[2], aver[3] * 100)
    #str <- sprintf("Var %s & %.3f & %.3f %.3f \\", pref, var[1], var[2], var[3])
    print(str)
    png(pref)
    plotKFoldError(resultAll[[w, a]], pref)
    dev.off()
  }
}

errorsAll <- array(dim = c(3, length(alldays)))
for (i in alldays) {
  partData <- tempdata.sp[tempdata.sp@data$datenum == i, ]
  gwrmodel.bw <- gwr.sel(tempFormula, partData, verbose = F, adapt = T)
  gwrmodel.gwr <- gwr(tempFormula, data = partData, adapt = gwrmodel.bw, hatmatrix = T)
  errF <- errorCal(partData@data[, "temp"], gwrmodel.gwr$SDF$pred)
  errorsAll[, i] <- unlist(errF, use.names = F)
}

setwd("~/Desktop/GWR/Temp/")
folderList <- c("MOD06", "MOD07", "MYD06", "MYD07", "VIIRS")
extractDate <- sample(huddata$datenum, 1)
extractedData <- getDataFromImage(folderList, extractDate)
extractedData$averagep <- rowMeans(extractedData[, 3:6], na.rm = TRUE)
extractedsub <- extractedData[, c("y", "x", "averagep")]
extractedsub <- extractedsub[complete.cases(extractedsub), ]
extractedsub.sp <- convertToSpatial(extractedsub, "y", "x", proj)

predictedData <- makeNewImage(hudForms[["Pressure"]], huddata.sp[huddata.sp@data$datenum == extractDate,], extractedsub.sp, gwr.Gauss, T)
predictedData$averagep <- NULL
saveImage(predictedData, "revised_hud_ex_28.tif")
