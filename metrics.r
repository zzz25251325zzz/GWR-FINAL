library(Metrics)
library(Hmisc)
library(dplyr)
library(caret)

errorCal <- function(actual, predicted, prefix = NULL) {
    rmse <- rmse(actual, predicted)
    #for single variable regression model. Not adjusted
    err <- actual - predicted
    rss <- sum(err * err)
    meanActual <- mean(actual)
    totalVar <- sum((actual - meanActual) * (actual - meanActual))
    r2 <- 1 - rss / totalVar
    #mean relative error: precision or accuracy? 
    mre <- mean(ape(actual, predicted))
    if (!is.null(prefix)) {
        dat <- sprintf("%s & %.3f & %.3f & %.3f \\", prefix, rmse, mre, r2)
        print(dat)
    }
    return(list(rmse, r2, mre))
}

gwrLeaveOneOut <- function(formula, data, weight, adapt, prefix, col) {
    sum.rmse <- 0
    sum.mre <- 0
    for (i in 1:nrow(data)) {
        testData <- data[i,]
        trainData <- data[-i,]
        trainData.bw <- gwr.sel(formula, data = trainData, verbose = F, gweight = weight, adapt = adapt)
        if (adapt) {
            trainData.gwr <- gwr(formula, data = trainData, adapt = trainData.bw, hatmatrix = T, gweight = weight)
        } else {
            trainData.gwr <- gwr(formula, data = trainData, bandwidth = trainData.bw, hatmatrix = T, gweight = weight)
        }
        testData.gwr <- gwr(formula, data = trainData, predictions = T,
                            fit.points = testData, se.fit = T,
                            fittedGWRobject = trainData.gwr,
                            gweight = weight, adapt = trainData.gwr$adapt,
                            bandwidth = trainData.gwr$bandwidth)
        errors <- errorCal(testData@data[, col], testData.gwr$SDF$pred)
        sum.rmse <- sum.rmse + as.numeric(errors[1])
        sum.mre <- sum.mre + as.numeric(errors[3])
    }
    rmse <- sum.rmse / nrow(data)
    mre <- sum.mre / nrow(data)
    dat <- sprintf("%s & %.3f & %.3f \\", prefix, rmse, mre)
    print(dat)
}

gwrCV <- function(formula, data, weight, adapt, fold, col) {
    #simple random division of data into folds
    #data <- data[sample(nrow(data)),]
    #folds <- cut(seq(1, nrow(data)), breaks = fold, labels = F)

    all.rmse <- vector("numeric", fold)
    all.r2 <- vector("numeric", fold)
    all.mre <- vector("numeric", fold)
    #or fold with stratified latitude
    folds <- createFolds(data@coords[,2], k=fold)
    
    for (i in 1:fold) {
        #testIndices <- which(folds == i, arr.ind = T)
        testIndices <- folds[[i]]
        testData <- data[testIndices,]
        trainData <- data[-testIndices,]

        trainData.bw <- gwr.sel(formula, data = trainData, verbose = F, gweight = weight, adapt = adapt)
        if (adapt) {
            trainData.gwr <- gwr(formula, data = trainData, adapt = trainData.bw, hatmatrix = T, gweight = weight)
        } else {
            trainData.gwr <- gwr(formula, data = trainData, bandwidth = trainData.bw, hatmatrix = T, gweight = weight)
        }
        #errorCal(trainData@data[, col], trainData.gwr$SDF$pred, "train")
        testData.gwr <- gwr(formula, data = trainData, predictions = T,
                            fit.points = testData, se.fit = T,
                            fittedGWRobject = trainData.gwr,
                            gweight = weight, adapt = trainData.gwr$adapt,
                            bandwidth = trainData.gwr$bandwidth)
        errors <- errorCal(testData@data[, col], testData.gwr$SDF$pred)
        if (any(is.na(errors))) {
            bw <- ifelse(is.null(trainData.gwr$bandwidth), 0.0, trainData.bw)
            ad <- ifelse(is.null(trainData.gwr$adapt), 0.0, trainData.bw)            
            message <- sprintf("Error is NA in fold %d bandwidth %.3f adapt %.3f", i, bw, ad)
            print(message)
        }
        all.rmse[i] <- as.numeric(errors[1])
        all.r2[i] <- as.numeric(errors[2])
        all.mre[i] <- as.numeric(errors[3])
    }
    
    rmse <- mean(all.rmse)
    r2 <- median(all.r2)
    mre <- mean(all.mre)
    return(list(rmse, r2, mre))
}

plotKFoldError <- function(data, name) {
    par(mfrow = c(3, 1))
    plot(data[1,], xlab = "Day", ylab = "RMSE")
    mtext(name)
    plot(data[2,], xlab = "Day", ylab = "R2")
    plot(data[3,], xlab = "Day", ylab = "RE")    
}

plotFullError <- function(actual, predicted) {
    residual <- actual - predicted    
    par(mfrow=c(1,2))
    plot(predicted, residual)
    qqplot(predicted, residual)
}

plotOneDay <- function(data, day, form, col) {
  partData <- data[data@data$datenum == day, ]
  model.bw <- gwr.sel(form, partData, verbose = F, adapt = T)
  model.gwr <- gwr(form, data = partData, adapt = model.bw, hatmatrix = T)
  plotFullError(partData@data[,col], model.gwr$SDF$pred)
}

plotErrorOneDay <- function(formula, weight, adapt, data, days, col) {
    d <- sample(days, 1)
    print(d)
    sub <- data[data@data$datenum == d,]
    sub.bw <- gwr.sel(formula, sub, verbose = F, adapt = adapt)
    if (adapt) {
        sub.gwr <- gwr(formula, data = sub, adapt = sub.bw, gweight = weight)
    } else {
        sub.gwr <- gwr(formula, data = sub, bandwidth = sub.bw, gweight = weight)
    }
    predicted <- sub.gwr$SDF$pred
    actual <- sub@data[, col]
    residual <- actual - predicted
    sub$res <- residual
    spplot(sub, "res")
}

rowVar <- function(x) {
    rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

devErrorCal <- function(formula, train, bw, object, weight, adapt, data, days, col) {
    sum.rmse <- 0
    sum.mae <- 0
    sum.mape <- 0
    count <- 0
    for (d in days) {
        sub <- data[data@data$datenum == d,]
        sub.gwr <- gwr(formula, data = train, bandwidth = bw, predictions = T,
                       fit.points = sub, se.fit=T, fittedGWRobject = object,
                       gweight = weight, adapt = adapt)
        actual <- sub@data[, col]
        predicted <- sub.gwr$SDF$pred
        sub.rmse <- rmse(actual, predicted)
        sub.mae <- mean(ae(actual, predicted))
        sub.mape <- mean(ape(actual, predicted))
        count <- count + 1
        sum.rmse <- sum.rmse + sub.rmse
        sum.mae <- sum.mae + sub.mae
        sum.mape <- sum.mape + sub.mape
    }
    print(sum.rmse / count)
    print(sum.mae / count)
    print(sum.mape / count)
}

plotDev <- function(formula, train, bw, object, weight, adapt, data, days, col) {
    d <- sample(days, 1)
    sub <- data[data@data$datenum == d,]
    sub.gwr <- gwr(formula, data = train, bandwidth = bw, predictions = T,
                   fit.points = sub, se.fit=T, fittedGWRobject = object,
                   gweight = weight, adapt = adapt)
    actual <- sub@data[, col]
    predicted <- sub.gwr$SDF$pred
    residual <- actual - predicted
    par(mfrow=c(1,2))
    plot(predicted, residual)
    qqplot(predicted, residual)
}

