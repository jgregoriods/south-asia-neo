suppressPackageStartupMessages({
    library(gdistance)
    library(rgdal)
    library(parallel)
    library(viridis)
    library(RColorBrewer)
    library(pals)
    library(rasterVis)
    library(smatr)
})

set.seed(123)

ORIGIN <- c(42.45, 36.37)
START <- 11748
DATES <- read.csv("sites/dates100.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- CRS("+init=epsg:4326")

BIOMES <- raster("layers/biomes25.tif")

simulateDispersal <- function(costRaster, origin, date) {
    tr <- transition(costRaster, function(x) 1 / mean(x), 16)
    tr <- geoCorrection(tr)
    ac <- accCost(tr, origin)
    ac <- ac / 1000
    ac[values(ac) == Inf] <- NA
    simDates <- date - ac
    res <- list("dates"=simDates, "dists"=ac)
    return(res)
}

compareDates <- function(simRaster, dates) {
    dates$simbp <- extract(simRaster$dates, dates)
    dates$dist <- spDistsN1(dates, ORIGIN, longlat=TRUE)
    plot(dates$dist, dates$bp, xlab="Distance from origin (km)",
         ylab="Age (cal BP)", pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5)
}

testModel <- function(costRaster, sites=DATES, origin=ORIGIN, date=START) {
    simDates <- simulateDispersal(costRaster, origin, date)
    gc() 
    sites$dist <- extract(simDates$dists, sites)
    rma <- sma(bp~dist, data=sites, robust=T)
    return(rma$r2[[1]])
}

###########################################################
# Genetic Algorithm
###########################################################

crossover <- function(x, y) {
    i <- sample(1:(length(x)-1), 1)
    return(c(x[1:i], y[(i+1):length(y)]))
}

mutate <- function(x) {
    i <- sample(1:length(x), 1)
    x[i] <- x[i] + rnorm(1)
    return(x)
}

GA <- function(numGenes, numGenomes, numParents, numElite, mutationRate,
               numIter, cores=NULL) {
    # Initialize genomes
    genomes <- as.data.frame(matrix(nrow=numGenomes, ncol=numGenes+1))
    for (i in 1:numGenomes) {
        genomes[i,] <- c(rnorm(numGenes, mean=1), 0)
    }

    maxScores <- c()
    avgScores <- c()

    if (is.null(cores)) {
        ncores <- detectCores() - 1
    } else {
        ncores <- cores
    }

    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library("gdistance"))
    clusterEvalQ(cl, library("smatr"))
    clusterExport(cl, varlist=c("BIOMES", "ORIGIN", "START", "DATES",
                                "reclassRaster", "testModel",
                                "simulateDispersal"), envir=environment())

    cat(paste("\nRunning genetic algorithm on", ncores,
            "parallel workers.\nThis may take a while...\n"))

    for (iter in 1:numIter) {
        genomeList <- split(genomes, seq(nrow(genomes)))

        res <- parLapply(cl, genomeList, function(x) {
            if (x[1,numGenes+1] != 0) {
                return (x[1,numGenes+1])
            } else {
                if (min(x[1,1:numGenes]) <= 0) {
                    return(0)
                } else {
                    reclass_matrix <- cbind(1:6, as.numeric(x[1,1:numGenes]))
                    cost <- reclassify(BIOMES, reclass_matrix)
                    score <- testModel(cost)
                    gc()
                    return(score)
                }
            }
        })

        genomes[,numGenes+1] <- unlist(res)

        avgScores[iter] <- mean(genomes[,numGenes+1][genomes[,numGenes+1] > 0])

        elite <- genomes[order(genomes[,numGenes+1], decreasing=T),][1:numElite,]
        parents <- genomes[order(genomes[,numGenes+1], decreasing=T),][1:numParents,]
        parents <- parents[order(as.numeric(rownames(parents)), decreasing=T),]

        maxScores[iter] <- elite[1,numGenes+1]

        children <- as.data.frame(matrix(nrow=numGenomes-numElite,
                                         ncol=numGenes+1))
        j <- 1
        while (j <= numGenomes - numElite) {
            for (i in seq(1, numParents - 1, 2)) {
                if (j > numGenomes - numElite) {
                    break
                }
                parent1 <- as.numeric(parents[i,])[1:numGenes]
                parent2 <- as.numeric(parents[i+1,])[1:numGenes]
                child <- crossover(parent1, parent2)
                if (runif(1) < mutationRate) {
                    child <- mutate(child)
                }
                children[j,] <- c(child, 0)
                j <- j+1
            }
        }
        genomes <- rbind(elite, children)
        rownames(genomes) <- sample(1:numGenomes)
        cat("\r", floor(iter / numIter*100), "%\tAvg score:", avgScores[iter])
    }
    cat("\n")

    stopCluster(cl)

    return (list(genomes=genomes, maxScores=maxScores, avgScores=avgScores))
}

main <- function() {

    numGenes <- 6

    numGenomes <- 500
    numParents <- 200
    numElite <- 50
    mutationRate <- 0.2
    numIter <- 20

    start_time <- Sys.time()
    res <- GA(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, cores=10)
    total_time <- Sys.time() - start_time
    cat("Completed in",total_time[[1]],attributes(total_time)$units,"\n")

    if (F) {
    best <- as.numeric(res$genomes[1,])
    reclass_matrix <- cbind(1:6, best[1:numGenes])
    costRaster <- reclassify(BIOMES, reclass_matrix)

    simDates <- simulateDispersal(costRaster, ORIGIN, START)
    simDates.r <- crop(simDates$dates, extent(DATES))

    cat(paste("Best score:", best[numGenes+1], "\n"))
    write.csv(res$genomes, "results.csv")

    par(mfrow=c(2,2))
    plot(res$avgScores, col="blue", main="RMSE",
         ylim=c(min(res$maxScores), max(res$avgScores)),
         pch=0, type="b", ylab="RMSE", xlab="Generation")
    lines(res$maxScores, col="red", pch=0, type="b")
    legend("topright", legend = c("Avg score", "Best score"),
        col = c("blue", "red"), lty = 1:1)
    }

    save(res, file="ga0811.RData")
}
