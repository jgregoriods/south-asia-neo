suppressPackageStartupMessages({
    library(gdistance)
    library(rgdal)
    library(parallel)
    library(rasterVis)
    library(rgrass7)
    library(raster)
    library(sp)
})

use_sp()

set.seed(123)

ORIGIN <- c(43.5, 36.33)    # M'lefaat
START <- 12855
DATES <- read.csv("sites/dates100.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- CRS("+init=epsg:4326")

BIOMES <- raster("layers/biomes25.tif")
BIOMES.sp <- as(BIOMES, "SpatialPixelsDataFrame")
names(BIOMES.sp@data) <- c("val")

setGRASS <- function(RASTER, RES) {
    writeRAST(RASTER, "GRID", zcol="val", flags=c("quiet", "overwrite"))
    execGRASS("g.region", flags=c("d"), res=as.character(RES), raster="GRID")
}

initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T)
execGRASS('g.proj', flags=c('c'), georef='layers/biomes25.tif')
execGRASS('r.import', input='layers/biomes25.tif', output='DEM')
execGRASS('g.region', raster='DEM')

#setGRASS(BIOMES.sp, 0.25)

costSurface <- function(coords) {
    execGRASS("r.cost", flags=c("k", "overwrite", "quiet"), input="GRID", start_coordinates=coords, output="COST", max_cost=0)
    cost <- raster(readRAST("COST"))
    return(cost)
}

#simulateDispersal <- function(costRaster, origin, start_date) {
simulateDispersal <- function(origin, start_date) {
    #costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
    #names(costSPDF@data) <- c("val")
#    tr <- transition(costRaster, function(x) 1 / mean(x), 16)
#    tr <- geoCorrection(tr)
#    ac <- accCost(tr, origin)
#    ac <- ac / 1000
#    ac[values(ac) == Inf] <- NA
    cs <- costSurface(origin)
    ac <- cs * 0.25 * 111
    simDates <- start_date - ac
    return(simDates)
}

compareDates <- function(simRaster, dates) {
    dates$simbp <- extract(simRaster, dates)
    dates$dist <- spDistsN1(dates, ORIGIN, longlat=TRUE)
    plot(dates$dist, dates$bp, xlab="Distance from origin (km)",
         ylab="Age (cal BP)", pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5)
    points(dates$dist, dates$simbp)
}

testModel <- function(costRaster, sites=DATES, origin=ORIGIN, start_date=START) {
    costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
    names(costSPDF@data) <- c("val")
    setGRASS(costSPDF, 0.25)
    #simDates <- simulateDispersal(costRaster, origin, start_date)
    simDates <- simulateDispersal(origin, start_date)
    gc() 
    sites$simbp <- extract(simDates, sites)
    sites <- sites[!is.na(sites$simbp),]
    rmse <- sqrt(sum((sites$simbp - sites$bp)^2) / nrow(sites))
    return(rmse)
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
    clusterEvalQ(cl, library("sp"))
    clusterEvalQ(cl, library("raster"))
    clusterEvalQ(cl, library("rgrass7"))
    clusterExport(cl, varlist=c("BIOMES", "BIOMES.sp", "ORIGIN", "START", "DATES",
                                "testModel", "simulateDispersal", "costSurface", "setGRASS"),
                                envir=environment())
    clusterEvalQ(cl, use_sp())
    clusterEvalQ(cl, initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T))
    clusterEvalQ(cl, execGRASS('g.proj', flags=c('c'), georef='layers/biomes25.tif'))
    clusterEvalQ(cl, execGRASS('r.import', input='layers/biomes25.tif', output='DEM'))
    clusterEvalQ(cl, execGRASS('g.region', raster='DEM'))

    #clusterEvalQ(cl, setGRASS(BIOMES.sp, 0.25))

    cat(paste("\nRunning genetic algorithm on", ncores,
            "parallel workers.\nThis may take a while...\n"))

    for (iter in 1:numIter) {
        genomeList <- split(genomes, seq(nrow(genomes)))

        res <- parLapply(cl, genomeList, function(x) {
            if (x[1,numGenes+1] != Inf) {
                return (x[1,numGenes+1])
            } else {
                if (min(x[1,1:numGenes]) <= 0) {
                    return(Inf)
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

        avgScores[iter] <- mean(genomes[,numGenes+1][genomes[,numGenes+1] < Inf])

        elite <- genomes[order(genomes[,numGenes+1]),][1:numElite,]
        parents <- genomes[order(genomes[,numGenes+1]),][1:numParents,]
        parents <- parents[order(as.numeric(rownames(parents))),]

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
                children[j,] <- c(child, Inf)
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
    costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
    names(costSPDF@data) <- c("val")

    setGRASS(costSPDF, 0.25)

    simDates <- simulateDispersal(ORIGIN, START)
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

    save(res, file="ga1111.RData")
}

#main()