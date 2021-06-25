suppressPackageStartupMessages({
    library(gdistance)
    library(rcarbon)
    library(parallel)
})

ORIGIN <- c(42.45, 36.37)
START <- 11748
DATES <- read.csv("sites/dates.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- CRS("+init=epsg:4326")

cal <- calibrate(DATES$C14, DATES$SD, verbose=FALSE)
DATES$bp <- medCal(cal)

ELE_RAW <- raster("layers/temp_cold.tif")
PREC_RAW <- raster("layers/prec.tif")


normRaster <- function(x) {
    return ((x - min(values(x), na.rm=TRUE)) / (max(values(x), na.rm=TRUE) - min(values(x), na.rm=TRUE)))
}


# Transform and scale
#ELE <- normRaster((ELE_RAW)^(1/3))
ELE <- normRaster(ELE_RAW)

# Remove outliers and 0 (for log transform)
maxVal <- quantile(PREC_RAW, .99)
PREC_RAW[values(PREC_RAW > maxVal)] <- maxVal
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
    plot(dates$dist, dates$bp)
    points(dates$dist, dates$simbp, col="red")
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
numElite <- 5
mutationRate <- 0.2
numIter <- 20

# Initialize genomes
genomes <- as.data.frame(matrix(nrow=numGenomes, ncol=5))
for (i in 1:numGenomes) {
    genomes[i,] <- c(rnorm(4, sd=3), Inf)
}

maxScores <- c()
avgScores <- c()

cat("Running genetic algorithm. This may take a while...\n")
pb <- txtProgressBar(min=0, max=numIter, style=3)
setTxtProgressBar(pb, 0)
for (iter in 1:numIter) {
    genomeList <- split(genomes, seq(nrow(genomes)))

    #ncores <- detectCores() - 1
    ncores <- 8
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library("gdistance"))
    clusterExport(cl, varlist=c("ELE", "PREC", "ORIGIN", "START", "DATES",
                                "testModel", "simulateDispersal"),
                                envir=environment())

    res <- parLapply(cl, genomeList, function(x) {
        if (x[1,5] != Inf) {
            return (x[1,5])
        } else {
            costEle <- (ELE * x[1,1]) + x[1,2]
            costPrc <- (PREC * x[1,3]) + x[1,4]
            if (min(values(costEle+costPrc), na.rm=T) <= 0) {
                return(Inf)
            } else {
                score <- testModel(costEle+costPrc)
                gc()
                return(score)
            }
        }
    })

    stopCluster(cl)
    genomes[,5] <- unlist(res)

    avgScores[iter] <- mean(genomes[,5][!is.infinite(genomes[,5])])

    elite <- genomes[order(genomes[,5]),][1:numElite,]
    parents <- genomes[order(genomes[,5]),][1:numParents,]
    parents <- parents[order(as.numeric(rownames(parents))),]

    maxScores[iter] <- elite[1,5]

    children <- as.data.frame(matrix(nrow=numGenomes - numElite, ncol=5))
    j <- 1
    while (j <= numGenomes - numElite) {
        for (i in seq(1, numParents - 1, 2)) {
            if (j > numGenomes - numElite) {
                break
            }
            parent1 <- as.numeric(parents[i,])[1:4]
            parent2 <- as.numeric(parents[i+1,])[1:4]
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

par(mfrow=c(2,2))
plot(avgScores, type="l", col="blue", main="RMSE", ylim=c(0, max(avgScores)))
lines(maxScores, col="red")

best <- as.numeric(genomes[1,])
costRaster <- (((ELE * best[1]) + best[2]) + ((PREC * best[3]) + best[4]))
simDates <- simulateDispersal(costRaster, ORIGIN, START)

compareDates(simDates, DATES)

plot(1/costRaster, main="Speed (km/yr)")
plot(simDates, main="Arrival (cal BP)")
