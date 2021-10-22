suppressPackageStartupMessages({
    library(gdistance)
    #library(rcarbon)
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
#ORIGIN <- c(35.5511, 32.62)
#START <- 11537
DATES <- read.csv("sites/dates.csv")
#DATES <- read.csv("sites/merged.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- CRS("+init=epsg:4326")

#coast <- ne_download(scale=10, type="coastline", category="physical")
coast <- readOGR("layers/ocean.shp")

#cal <- calibrate(DATES$C14, DATES$SD, verbose=FALSE)
#DATES$bp <- medCal(cal)

BIOMES <- raster("layers/newBiomes_final.tif")


normRaster <- function(x) {
    return ((x - min(values(x), na.rm=TRUE)) / (max(values(x), na.rm=TRUE) - min(values(x), na.rm=TRUE)))
}


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
    #dates$dist <- spDistsN1(dates, ORIGIN, longlat=TRUE)
    dates$dist <- extract(simRaster$dists, dates)

    #model <- lm(simbp~poly(dist, 2), data=dates)
    #x <- min(dates$dist, na.rm=T):max(dates$dist, na.rm=T)
    #y <- predict(model, newdata=data.frame(dist=x))

    plot(dates$dist, dates$bp, xlab="Distance from origin (km)", ylab="Age (cal BP)", pch=20, cex=1.5,
         cex.axis=1.5, cex.lab=1.5)
    #points(dates$dist, dates$simbp, col=4, pch=3, cex=1.5)
    #lines(x, y, col=4)
    #legend("topright", legend=c("14C dates", "simulated"), col=c("black", 4),
    #                            pch=c(20, 3), cex=c(1.5,1.5), y.intersp=1.5, box.lty=0)
}

testModel <- function(costRaster, sites=DATES, origin=ORIGIN, date=START) {
    simDates <- simulateDispersal(costRaster, origin, date)
    gc() 
    # Score
    #sites$simbp <- extract(simDates, sites)
    sites$dist <- extract(simDates$dists, sites)
    #sites <- sites[!is.na(sites$simbp),]
    #rmse <- sqrt(sum((sites$simbp - sites$bp)^2) / nrow(sites))
    #return(rmse)
    rma <- sma(bp~dist, data=sites, robust=T)
    return(rma$r2[[1]])
}

reclassRaster <- function(r, vals) {
    r.new <- r
    codes <- 1:6
    #codes <- c(1,4,8,10,13,100)
    #codes <- c(1,4,8,10,13)
    #codes <- c(1,4,9,13)
    #codes <- c(1,2,4,5,6,9,11,12,13,18,19,21)
    for (i in 1:length(codes)) {
        r.new[values(r) == codes[i]] <- vals[i]
    }
    return (r.new)
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

GA <- function(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, cores=NULL) {
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
    clusterExport(cl, varlist=c("BIOMES", "ORIGIN", "START", "DATES", "reclassRaster",
                                "testModel", "simulateDispersal"),
                                envir=environment())

    cat(paste("\nRunning genetic algorithm on", ncores,
            "parallel workers.\nThis may take a while...\n"))
    #pb <- txtProgressBar(min=0, max=numIter, style=3)
    #setTxtProgressBar(pb, 0)
    for (iter in 1:numIter) {
        genomeList <- split(genomes, seq(nrow(genomes)))

        res <- parLapply(cl, genomeList, function(x) {
            if (x[1,numGenes+1] != 0) {
                return (x[1,numGenes+1])
            } else {
                if (min(x[1,1:numGenes]) <= 0) {
                    return(0)
                } else {
                    cost <- reclassRaster(BIOMES, as.numeric(x[1,1:numGenes]))
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

        children <- as.data.frame(matrix(nrow=numGenomes - numElite, ncol=numGenes+1))
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
        #setTxtProgressBar(pb, iter)
        cat("\r", floor(iter / numIter * 100), "%\tAvg score:", avgScores[iter])
    }
    cat("\n")
    #close(pb)

    stopCluster(cl)

    return (list(genomes=genomes, maxScores=maxScores, avgScores=avgScores))
}


plotMap <- function(r) {
    ext <- extent(r)
    plot(r, xlim=c(ext[1], ext[2]), ylim=c(ext[3], ext[4]), col=viridis_pal()(255),
         main="Simulated arrival times", xlab="Longitude", ylab="Latitude",
         legend.args=list(text="Age (yr BP)", side=2, font=1))
    minDate <- floor(min(values(r), na.rm=T) / 1000) * 1000
    contour(r, add=T, levels=seq(minDate, max(values(r), na.rm=T), 1000), col="white", lwd=1)
    plot(coast, add=T, lwd=1, col="lightgrey")
}


plotSpeed <- function(r) {
    ext <- extent(r)
    plot(r, xlim=c(ext[1], ext[2]), ylim=c(ext[3], ext[4]), col=viridis_pal()(255),
         main="Simulated speed (km/yr)", xlab="Long", ylab="Lat")
    #contour(r, add=T, levels=seq(0, max(values(r), na.rm=T), 0.1))
    plot(coast, add=T, lwd=1.5)
}

main <- function() {

    numGenes <- 6

    numGenomes <- 500
    numParents <- 200
    numElite <- 50
    mutationRate <- 0.1
    numIter <- 20

    res <- GA(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, cores=6)

    best <- as.numeric(res$genomes[1,])
    costRaster <- reclassRaster(BIOMES, best[1:numGenes])
    simDates <- simulateDispersal(costRaster, ORIGIN, START)
    simDates <- crop(simDates, extent(DATES))

    cat(paste("Best score:", best[numGenes+1], "\n"))
    write.csv(res$genomes, "results.csv")

    par(mfrow=c(2,2))
    plot(res$avgScores, col="blue", main="RMSE", ylim=c(min(res$maxScores), max(res$avgScores)),
        pch=0, type="b", ylab="RMSE", xlab="Generation")
    lines(res$maxScores, col="red", pch=0, type="b")
    legend("topright", legend = c("Avg score", "Best score"),
        col = c("blue", "red"), lty = 1:1)

    save(res, file="ga.RData")

}