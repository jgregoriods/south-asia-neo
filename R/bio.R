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

BIOMES <- raster("layers/biomesIndus.tif")


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


reclassRaster <- function(r, vals) {
    r.new <- r
    codes <- c(1,4,8,10,13,100)
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
        genomes[i,] <- c(rnorm(numGenes, mean=1), Inf)
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
            if (x[1,numGenes+1] != Inf) {
                return (x[1,numGenes+1])
            } else {
                if (min(x[1,1:numGenes]) <= 0) {
                    return(Inf)
                } else {
                    cost <- reclassRaster(BIOMES, as.numeric(x[1,1:numGenes]))
                    score <- testModel(cost)
                    gc()
                    return(score)
                }
            }
        })

        genomes[,numGenes+1] <- unlist(res)

        avgScores[iter] <- mean(genomes[,numGenes+1][!is.infinite(genomes[,numGenes+1])])

        elite <- genomes[order(genomes[,numGenes+1]),][1:numElite,]
        parents <- genomes[order(genomes[,numGenes+1]),][1:numParents,]
        parents <- parents[order(as.numeric(rownames(parents))),]

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
                children[j,] <- c(child, Inf)
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
         main="Simulated arrival times (cal BP)", xlab="Long", ylab="Lat")
    minDate <- floor(min(values(r), na.rm=T) / 1000) * 1000
    contour(r, add=T, levels=seq(minDate, max(values(r), na.rm=T), 1000))
    plot(coast, add=T, lwd=1.5)
}


plotSpeed <- function(r) {
    ext <- extent(r)
    plot(r, xlim=c(ext[1], ext[2]), ylim=c(ext[3], ext[4]), col=viridis_pal()(255),
         main="Simulated speed (km/yr)", xlab="Long", ylab="Lat")
    #contour(r, add=T, levels=seq(0, max(values(r), na.rm=T), 0.1))
    plot(coast, add=T, lwd=1.5)
}


numGenes <- 6

numGenomes <- 500
numParents <- 200
numElite <- 50
mutationRate <- 0.1
numIter <- 20

res <- GA(numGenes, numGenomes, numParents, numElite, mutationRate, numIter)

best <- as.numeric(res$genomes[1,])
costRaster <- reclassRaster(BIOMES, best[1:numGenes])
simDates <- simulateDispersal(costRaster, ORIGIN, START)

cat(paste("Best score:", best[numGenes+1], "\n"))
write.csv(res$genomes, "results1407.csv")

par(mfrow=c(2,2))
plot(res$avgScores, col="blue", main="RMSE", ylim=c(min(res$maxScores), max(res$avgScores)),
     pch=0, type="b", ylab="RMSE", xlab="Generation")
lines(res$maxScores, col="red", pch=0, type="b")
legend("topright", legend = c("Avg score", "Best score"),
       col = c("blue", "red"), lty = 1:1)


compareDates(simDates, DATES)

plotSpeed(1/costRaster)
plotMap(simDates)

save(res, file="ga1407.RData")

# scp jgregorio@marvin.s.upf.edu:/homes/users/jgregorio/south-asia-neo/Rplots.pdf Rplots.pdf
# scp jgregorio@marvin.s.upf.edu:/homes/users/jgregorio/SouthAsiaNeo/results.csv results1.csv