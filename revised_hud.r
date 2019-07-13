library(spgwr)
library(Metrics)

setwd("~/Desktop/GWR")
source("image-to-data.r")
source("metrics.r")

pressuredata <- combineDataWithCoord("Pressure_data.csv", "obs_20140101_20141231.csv", "hud")
vapordata <- combineDataWithCoord("WaterVapor_data.csv", "obs_20140101_20141231.csv", "hud")

pressuredata$averagep <- rowMeans(pressuredata[, 7:10], na.rm = TRUE)
pressuresub <- pressuredata[, c("name", "datenum", "lat", "lon", "hud", "averagep")]
vapordata$averagev <- rowMeans(vapordata[, 7:10], na.rm = TRUE)
vaporsub <- vapordata[, c("name", "datenum", "lat", "lon", "hud", "averagev")]

alldays <- unique(vaporsub$datenum)
huddata <- merge(pressuresub, vaporsub, by = c("name", "datenum", "lat", "lon", "hud"), all = T)
huddata <- huddata[complete.cases(huddata), ]

proj <- CRS("+proj=longlat +datum=WGS84 +no_defs")
huddata.sp <- convertToSpatial(huddata, "lat", "lon", proj)

hudFormula <- hud ~ averagev
numFolds <- 5
weightsFun <- c(gwr.Gauss, gwr.bisquare, gwr.tricube)
names(weightsFun) <- c("Gauss", "Bisquare", "Tricube")
adaptOpt <- c(T, F)
names(adaptOpt) <- c("Adaptive", "Fixed")

set.seed(1)
resultforms <- vector("list", 3)
hudForms <- c(hud ~ averagev, hud ~ averagep, hud ~ averagev + averagep)
names(hudForms) <- c("Vapor", "Pressure", "Both")
names(resultforms) <- names(hudForms)

for (f in names(hudForms)) {
  resultforms[[f]] <- matrix(list(), 3, 2)
  rownames(resultforms[[f]]) <- names(weightsFun)
  colnames(resultforms[[f]]) <- names(adaptOpt)
  for (i in rownames(resultAll)) {
    for (j in colnames(resultAll)) {
      resultforms[[f]][[i, j]] <- array(dim = c(3, length(alldays)))
    }
  }
  
  for (i in alldays) {
    partData <- huddata.sp[huddata.sp@data$datenum == i,]
    print(i)
    for (w in names(weightsFun)) {
      for (a in names(adaptOpt)) {
        pref <- sprintf("%s & %s & %s", f, a, w)
        err <- list(NA, NA, NA)
        while (any(is.na(err))) {
          err <- gwrCV(hudForms[[f]], partData, weightsFun[[w]], adaptOpt[[a]], numFolds, "hud")
          if (any(is.na(err))) { print("NA") }
        }
        str <- sprintf("%s & %.3f & %.3f & %.3f \\", pref, err[1], err[2], err[3])
        print(str)
        resultforms[[f]][[w, a]][,i] <- unlist(err, use.names = F)
      }
    }    
  }
}

print("All days")
for (f in names(hudForms)) {
  for (w in names(weightsFun)) {
    for (a in names(adaptOpt)) {
      pref <- sprintf("%s & %s & %s", f, a, w)
      aver <- rowMeans(resultforms[[f]][[w, a]])
      #var <- rowVar(resultforms[[f]][[w, a]])
      str <- sprintf("%s & %.3f & %.3f & %.3f \\", pref, aver[1], aver[2], aver[3] * 100)
      #str <- sprintf("Var %s & %.3f & %.3f %.3f \\", pref, var[1], var[2], var[3])
      print(str)
      png(pref)
      plotKFoldError(resultforms[[f]][[w, a]], pref)
      dev.off()
    }
  }
}

huderrorsAll <- array(dim = c(3, length(alldays)))
for (i in alldays) {
  partData <- huddata.sp[huddata.sp@data$datenum == i, ]
  gwrmodel.bw <- gwr.sel(hudForms[["Pressure"]], partData, verbose = F, adapt = T)
  gwrmodel.gwr <- gwr(hudForms[["Pressure"]], data = partData, adapt = gwrmodel.bw, hatmatrix = T)
  errF <- errorCal(partData@data[, "hud"], gwrmodel.gwr$SDF$pred)
  huderrorsAll[, i] <- unlist(errF, use.names = F)
}

setwd("~/Desktop/GWR/Pressure")
folderList <- c("MOD06", "MOD07", "MYD06", "MYD07")
extractDate <- sample(huddata$datenum, 1)
extractedData <- getDataFromImage(folderList, extractDate)
extractedData$averagep <- rowMeans(extractedData[, 3:6], na.rm = TRUE)
extractedsub <- extractedData[, c("y", "x", "averagep")]
extractedsub <- extractedsub[complete.cases(extractedsub), ]
extractedsub.sp <- convertToSpatial(extractedsub, "y", "x", proj)

predictedData <- makeNewImage(hudForms[["Pressure"]], huddata.sp[huddata.sp@data$datenum == extractDate,], extractedsub.sp, gwr.Gauss, T)
predictedData$averagep <- NULL
saveImage(predictedData, "revised_hud_ex_28.tif")
