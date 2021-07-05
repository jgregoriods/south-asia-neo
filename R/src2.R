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

par(mfrow=c(2,2))
plot(avgScores, type="l", col="blue", main="RMSE", ylim=c(0, max(avgScores)))
lines(maxScores, col="red")

best <- as.numeric(genomes[1,])
costRaster <- ((TEMP * best[1]) + (PREC * best[2]) + best[3])
simDates <- simulateDispersal(costRaster, ORIGIN, START)

compareDates(simDates, DATES)

plot(1/costRaster, main="Speed (km/yr)")
plot(simDates, main="Arrival (cal BP)")
