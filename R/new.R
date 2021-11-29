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

scaleRaster <- function(x) {return ((x - min(values(x), na.rm=TRUE)) / (max(values(x), na.rm=TRUE) - min(values(x), na.rm=TRUE)))}

WGS <- CRS("+init=epsg:4326")
ALBERS <- CRS("+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#RASTER_FILE <- "layers/biomes__.tif"
#RASTER_FILE <- "layers/biomes_b_rec_m.tif"

ELEV_FILE <- "new_layers/ele_m.tif"
TEMP_FILE <- "new_layers/bio8_m.tif"
PREC_FILE <- "new_layers/bio19_m.tif"

ELEV <- scaleRaster(raster(ELEV_FILE))
TEMP <- scaleRaster(raster(TEMP_FILE))
PREC <- scaleRaster(raster(PREC_FILE))
RASTERS <- stack(ELEV, TEMP, PREC)

#ORIGIN <- c(43.5, 36.33)    # M'lefaat
#START <- 12855
#DATES <- read.csv("sites/dates100.csv")
DATES <- read.csv("sites/dates_.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- WGS
DATES.m <- spTransform(DATES, ALBERS)

ORIGIN_LL <- DATES[DATES$Site=="Dhra",]@coords
ORIGIN <- as.numeric(DATES.m[DATES.m$Site=="Dhra",]@coords)
START <- DATES.m[DATES.m$Site=="Dhra",]$bp

#ORIGIN <- as.numeric(DATES.m[DATES.m$Site=="M'lefaat",]@coords)
#START <- DATES.m[DATES.m$Site=="M'lefaat",]$bp

#BIOMES <- raster("layers/biomes25.tif")
#BIOMES <- raster("layers/newBiomes_final.tif")
#BIOMES.m <- projectRaster(BIOMES, res=25000, crs=ALBERS, method="ngb")
#writeRaster(BIOMES.m, "layers/biomes.tif", overwrite=T)
#BIOMES <- raster("layers/biomes.tif")
#BIOMES <- raster(RASTER_FILE)
#BIOMES.sp <- as(BIOMES, "SpatialPixelsDataFrame")
ELEV.sp <- as(ELEV, "SpatialPixelsDataFrame")
names(ELEV.sp@data) <- c("val")

setGRASS <- function(raster, res) {
    writeRAST(raster, "GRID", zcol="val", flags=c("quiet", "overwrite"))
    execGRASS("g.region", flags=c("d"), res=as.character(res), raster="GRID")
}

initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T)
execGRASS('g.proj', flags=c('c'), georef=ELEV_FILE)
execGRASS('r.import', input=ELEV_FILE, output='DEM')
execGRASS('g.region', raster='DEM')

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
    simDates <- start_date - (cs * 25)
    return(simDates)
}

compareDates <- function(simRaster, dates) {
    dates$simbp <- extract(simRaster, dates)
    dates$dist <- spDistsN1(dates, ORIGIN_LL, longlat=T)
    plot(dates$dist, dates$bp, xlab="Distance from origin (km)",
         ylab="Age (cal BP)", pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5)
    points(dates$dist, dates$simbp)
}
compareDates(simDates.ll, DATES)

testModel <- function(costRaster, sites=DATES, origin=ORIGIN, start_date=START) {
    costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
    names(costSPDF@data) <- c("val")
    setGRASS(costSPDF, 25000)
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
    clusterExport(cl, varlist=c("ELEV_FILE", "RASTERS", "ELEV.sp", "ORIGIN", "START", "DATES",
                                "testModel", "simulateDispersal", "costSurface", "setGRASS"),
                                envir=environment())
    clusterEvalQ(cl, use_sp())
    clusterEvalQ(cl, initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T))
    clusterEvalQ(cl, execGRASS('g.proj', flags=c('c'), georef=ELEV_FILE))
    clusterEvalQ(cl, execGRASS('r.import', input=ELEV_FILE, output='DEM'))
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
                    #reclass_matrix <- cbind(1:numGenes, as.numeric(x[1,1:numGenes]))
                    #cost <- reclassify(BIOMES, reclass_matrix)
                    cost <- sum(RASTERS * as.numeric(x[1,1:numGenes]))
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
    numGenes <- 3
    numGenomes <- 50
    numParents <- 25
    numElite <- 5
    mutationRate <- 0.2
    numIter <- 5

    start_time <- Sys.time()
    res <- GA(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, cores=10)
    total_time <- Sys.time() - start_time
    cat("Completed in",total_time[[1]],attributes(total_time)$units,"\n")

    best <- as.numeric(res$genomes[1,])
    #reclass_matrix <- cbind(1:numGenes, best[1:numGenes])
    #costRaster <- reclassify(BIOMES, reclass_matrix)
    cost <- sum(RASTERS * best[1:numGenes])
    costSPDF <- as(cost, "SpatialPixelsDataFrame")
    names(costSPDF@data) <- c("val")

    setGRASS(costSPDF, 25000)

    simDates <- simulateDispersal(ORIGIN, START)
    #simDates.r <- crop(simDates, extent(DATES.m))
    simDates.ll <- projectRaster(simDates, res=0.25, crs=WGS)
    simDates.r <- crop(simDates.ll, extent(DATES))

    cat(paste("Best score:", best[numGenes+1], "\n"))

    timeStamp <- round(as.numeric(Sys.time()))

    speeds <- data.frame("region"=1:numGenes, "speed"=1/best[1:numGenes])
    write.csv(speeds, file=paste("results/res", timeStamp, ".csv", sep=""))
    save(res, simDates.ll, file=paste("results/ga", timeStamp, ".RData", sep=""))
    writeRaster(simDates.r, paste("results/sim", timeStamp, ".tif", sep=""))

    plot(simDates.ll)
}

#main()