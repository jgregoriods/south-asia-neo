suppressPackageStartupMessages({
    library(gdistance)
    #library(rcarbon)
    library(rgdal)
    library(parallel)
    library(viridis)
})

ORIGIN <- c(42.45, 36.37)
START <- 11748
DATES <- read.csv("sites/dates100.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- CRS("+init=epsg:4326")

#coast <- ne_download(scale=10, type="coastline", category="physical")
coast <- readOGR("layers/ocean.shp")

#cal <- calibrate(DATES$C14, DATES$SD, verbose=FALSE)
#DATES$bp <- medCal(cal)

TEMP_RAW <- raster("layers/temp_cold.tif")
PREC_RAW <- raster("layers/prec.tif")


normRaster <- function(x) {
    return ((x - min(values(x), na.rm=TRUE)) / (max(values(x), na.rm=TRUE) - min(values(x), na.rm=TRUE)))
}


# Transform and scale
#ELE <- normRaster((ELE_RAW)^(1/3))
TEMP <- normRaster(TEMP_RAW)

# Remove outliers and 0 (for log transform)
PREC_RAW[values(PREC_RAW > quantile(PREC_RAW, .99))] <- quantile(PREC_RAW, .99)
#PREC <- normRaster(log(PREC_RAW + 1))
PREC <- normRaster(PREC_RAW)


simulateDispersal <- function(costRaster, origin, date) {
    tr <- transition(costRaster, function(x) 1 / mean(x), 16)
    tr <- geoCorrection(tr)
    ac <- accCost(tr, origin)
    ac <- ac / 1000
    ac[values(ac) == Inf] <- NA
    simDates <- date - ac
    return(simDates)
}


compareDates <- function(simRaster, dates) {
    dates$simbp <- extract(simRaster, dates)
    dates$dist <- spDistsN1(dates, ORIGIN, longlat=TRUE)

    model <- lm(simbp~poly(dist, 2), data=dates)
    x <- min(dates$dist, na.rm=T):max(dates$dist, na.rm=T)
    y <- predict(model, newdata=data.frame(dist=x))

    plot(dates$dist, dates$bp, xlab="Distance from origin (km)", ylab="cal BP", pch=20)
    points(dates$dist, dates$simbp, col=4, pch=20)
    lines(x, y, col=4)
}


testModel <- function(costRaster, sites=DATES, origin=ORIGIN, date=START) {
    simDates <- simulateDispersal(costRaster, origin, date)
    gc() 
    # Score
    sites$simbp <- extract(simDates, sites)
    sites <- sites[!is.na(sites$simbp),]
    rmse <- sqrt(sum((sites$simbp - sites$bp)^2) / nrow(sites))
    return(rmse)
}

# Genetic Algorithm

crossover <- function(x, y) {
    i <- sample(1:(length(x)-1), 1)
    return(c(x[1:i], y[(i+1):length(y)]))
}

mutate <- function(x) {
    i <- sample(1:length(x), 1)
    x[i] <- x[i] + rnorm(1)
    return(x)
}

numGenomes <- 100
numParents <- 50
numElite <- 10
mutationRate <- 0.1
numIter <- 20

# Initialize genomes
genomes <- as.data.frame(matrix(nrow=numGenomes, ncol=4))
for (i in 1:numGenomes) {
    genomes[i,] <- c(rnorm(3), Inf)
}

maxScores <- c()
avgScores <- c()

ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterEvalQ(cl, library("gdistance"))
clusterExport(cl, varlist=c("TEMP", "PREC", "ORIGIN", "START", "DATES",
                            "testModel", "simulateDispersal"),
                            envir=environment())

cat(paste("\nRunning genetic algorithm on", ncores,
          "parallel workers.\nThis may take a while...\n"))
pb <- txtProgressBar(min=0, max=numIter, style=3)
setTxtProgressBar(pb, 0)
for (iter in 1:numIter) {
    genomeList <- split(genomes, seq(nrow(genomes)))

    res <- parLapply(cl, genomeList, function(x) {
        if (x[1,4] != Inf) {
            return (x[1,4])
        } else {
            cost <- (TEMP * x[1,1]) + (PREC * x[1,2]) + x[1,3]
            if (min(values(cost), na.rm=T) <= 0) {
                return(Inf)
            } else {
                score <- testModel(cost)
                gc()
                return(score)
            }
        }
    })

    genomes[,4] <- unlist(res)

    avgScores[iter] <- mean(genomes[,4][!is.infinite(genomes[,4])])

    elite <- genomes[order(genomes[,4]),][1:numElite,]
    parents <- genomes[order(genomes[,4]),][1:numParents,]
    parents <- parents[order(as.numeric(rownames(parents))),]

    maxScores[iter] <- elite[1,4]

    children <- as.data.frame(matrix(nrow=numGenomes - numElite, ncol=4))
    j <- 1
    while (j <= numGenomes - numElite) {
        for (i in seq(1, numParents - 1, 2)) {
            if (j > numGenomes - numElite) {
                break
            }
            parent1 <- as.numeric(parents[i,])[1:3]
            parent2 <- as.numeric(parents[i+1,])[1:3]
            child <- crossover(parent1, parent2)
            if (runif(1) < mutationRate) {
                child <- mutate(child)
            }
            children[j,] <- c(child, Inf)
            j <- j+1
        }
    }
    genomes <- rbind(elite, children)
    rownames(genomes) <- sample(1:numGenomes)
    setTxtProgressBar(pb, iter)
}
close(pb)

stopCluster(cl)


plotMap <- function(r) {
    ext <- extent(r)
    plot(r, xlim=c(ext[1], ext[2]), ylim=c(ext[3], ext[4]), col=viridis_pal()(255),
         main="Simulated arrival times (cal BP)", xlab="Long", ylab="Lat")
    minDate <- floor(min(values(r), na.rm=T) / 1000) * 1000
    contour(r, add=T, levels=seq(minDate, max(values(r), na.rm=T), 1000))
    plot(coast, add=T, lwd=1.5)
}


plotSpeed <- function(r) {
    ext <- extent(r)
    plot(r, xlim=c(ext[1], ext[2]), ylim=c(ext[3], ext[4]), col=viridis_pal()(255),
         main="Simulated speed (km/yr)", xlab="Long", ylab="Lat")
    contour(r, add=T, levels=seq(0, max(values(r), na.rm=T), 0.1))
    plot(coast, add=T, lwd=1.5)
}


best <- as.numeric(genomes[1,])
costRaster <- ((TEMP * best[1]) + (PREC * best[2]) + best[3])
simDates <- simulateDispersal(costRaster, ORIGIN, START)

cat(paste("Best score:", best[4], "\n"))
write.csv(genomes, "results.csv")

par(mfrow=c(2,2))
plot(avgScores, col="blue", main="RMSE", ylim=c(min(maxScores), max(avgScores)),
     pch=0, type="b", ylab="RMSE", xlab="Generation")
lines(maxScores, col="red", pch=0, type="b")
legend("topright", legend = c("Avg score", "Best score"),
       col = c("blue", "red"), lty = 1:1)


compareDates(simDates, DATES)

plotSpeed(1/costRaster)
plotMap(simDates)
