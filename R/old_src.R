library(gdistance)
library(rcarbon)
library(parallel)

ORIGIN <- c(42.45, 36.37)
START <- 11748
DATES <- read.csv("sites/dates.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- CRS("+init=epsg:4326")

cal <- calibrate(DATES$C14, DATES$SD, verbose=FALSE)
DATES$bp <- medCal(cal)

ELE <- raster("layers/ele.tif")
PREC <- raster("layers/prec.tif")


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


qtl <- function(x, newVals=1:4) {
    vals <- quantile(x, na.rm=T)
    m <- c(vals[1], vals[2], newVals[1],
           vals[2], vals[3], newVals[2],
           vals[3], vals[4], newVals[3],
           vals[4], vals[5], newVals[4])
    rclmat <- matrix(m, ncol=3, byrow=T)
    rc <- reclassify(x, rclmat)
    rc[values(rc) == vals[1]] <- newVals[1]
    return(rc)
}


################################################################################

crossover <- function(x, y) {
    i <- sample(1:(length(x)-1), 1)
    return(c(x[1:i], y[(i+1):length(y)]))
}

mutate <- function(x) {
    i <- sample(1:length(x), 1)
    x[i] <- abs(x[i] + rnorm(1))
    return(x)
}

numGenomes <- 100
numParents <- 20
numElite <- 5
mutationRate <- 0.2
numIter <- 50

# Initialize genomes
genomes <- as.data.frame(matrix(nrow=numGenomes, ncol=9))
for (i in 1:numGenomes) {
    genomes[i,] <- c(abs(1 + rnorm(8)), Inf)
}

maxScores <- c()
avgScores <- c()

cat("Running genetic algorithm. This may take a while...\n")
#pb <- txtProgressBar(min=0, max=numIter, style=3)
#setTxtProgressBar(pb, 0)
for (iter in 1:numIter) {
    genomeList <- split(genomes, seq(nrow(genomes)))

    ncores <- detectCores() - 1
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library("gdistance"))
    clusterExport(cl, varlist=c("qtl", "ELE", "PREC", "ORIGIN", "START", "DATES",
                                "testModel", "simulateDispersal"),
                                envir=environment())

    #for (i in 1:numGenomes) {
    #    if (genomes[i,9] == Inf) {
    #        costEle <- qtl(ELE, as.numeric(genomes[i,][1:4]))
    #        costPrc <- qtl(PREC, as.numeric(genomes[i,][5:8]))
    #        genomes[i,9] <- testModel(costEle+costPrc)
    #        gc()
    #    }
    #}

    res <- parLapply(cl, genomeList, function(x) {
        if (x[1,9] != Inf) {
            return (x[1,9])
        } else {
            costEle <- qtl(ELE, as.numeric(x[1,1:4]))
            costPrc <- qtl(PREC, as.numeric(x[1,5:8]))
            score <- testModel(costEle+costPrc)
            gc()
            return(score)
        }
    })

    stopCluster(cl)
    genomes[,9] <- unlist(res)

    avgScores[iter] <- mean(genomes[,9])

    elite <- genomes[order(genomes[,9]),][1:numElite,]
    parents <- genomes[order(genomes[,9]),][1:numParents,]
    parents <- parents[order(as.numeric(rownames(parents))),]

    maxScores[iter] <- elite[1,9]

    children <- as.data.frame(matrix(nrow=numGenomes - numElite, ncol=9))
    j <- 1
    while (j <= numGenomes - numElite) {
        for (i in seq(1, numParents, 2)) {
            if (j > numGenomes - numElite) {
                break
            }
            parent1 <- as.numeric(parents[i,])[1:8]
            parent2 <- as.numeric(parents[i+1,])[1:8]
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
    #setTxtProgressBar(pb, iter)
}
#close(pb)

plot(avgScores, type="l", col="blue")
lines(maxScores, col="red")

best <- genomes[1,]
costEle <- qtl(ELE, as.numeric(best[1:4]))
costPrc <- qtl(PREC, as.numeric(best[5:8]))
simDates <- simulateDispersal(costEle+costPrc, ORIGIN, START)
plot(simDates)
