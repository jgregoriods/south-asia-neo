suppressPackageStartupMessages({
    library(gdistance)
    library(rgdal)
    library(parallel)
    library(rasterVis)
    library(rgrass7)
    library(raster)
    library(sp)
    library(ggplot2)
})

use_sp()

set.seed(100)

WGS <- CRS("+init=epsg:4326")
ALBERS <- CRS("+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#RASTER_FILE <- "layers/biomes__.tif"
RASTER_FILE <- "layers/biomes.tif"

#ORIGIN <- c(43.5, 36.33)    # M'lefaat
#START <- 12855
#DATES <- read.csv("sites/dates100.csv")
DATES <- read.csv("dates/dates.csv")
coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- WGS
DATES.m <- spTransform(DATES, ALBERS)

ORIGIN_SITE <- "Mureybet"
ORIGIN <- as.numeric(DATES.m[DATES.m$Site==ORIGIN_SITE,]@coords)
START <- DATES.m[DATES.m$Site==ORIGIN_SITE,]$bp

#ORIGIN <- as.numeric(DATES.m[DATES.m$Site=="M'lefaat",]@coords)
#START <- DATES.m[DATES.m$Site=="M'lefaat",]$bp

#BIOMES <- raster("layers/biomes25.tif")
#BIOMES <- raster("layers/newBiomes_final.tif")
#BIOMES.m <- projectRaster(BIOMES, res=25000, crs=ALBERS, method="ngb")
#writeRaster(BIOMES.m, "layers/biomes.tif", overwrite=T)
#BIOMES <- raster("layers/biomes.tif")
BIOMES <- raster(RASTER_FILE)
BIOMES.sp <- as(BIOMES, "SpatialPixelsDataFrame")
names(BIOMES.sp@data) <- c("val")

setGRASS <- function(RASTER, RES) {
    writeRAST(RASTER, "GRID", zcol="val", flags=c("quiet", "overwrite"))
    execGRASS("g.region", flags=c("d"), res=as.character(RES), raster="GRID")
}

initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T)
execGRASS('g.proj', flags=c('c'), georef=RASTER_FILE)
execGRASS('r.import', input=RASTER_FILE, output='DEM')
execGRASS('g.region', raster='DEM')

#setGRASS(BIOMES.sp, 25000)

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

compareDates <- function(simRaster, dates, origin, rmse) {
    #dates$simbp <- extract(simRaster, dates)
    dist <- spDistsN1(dates, origin) / 1000
    merged <- rbind(data.frame("dist"=dist, "age"=dates$bp, "type"="C14"),
                    data.frame("dist"=dist, "age"=extract(simRaster, dates), "type"="Simulated"))
    ggplot(merged) +
        geom_point(aes(x=dist, y=age, color=type, shape=type), size=2) +
        annotate(geom="text", x=max(dist, na.rm=T) / 3, y=(max(merged$age, na.rm=T) + min(merged$age, na.rm=T)) / 3, label=paste("RMSE=", round(rmse), sep=""), color=4) +
        #geom_point(aes(x=dist, y=simbp, fill="Simulated"), shape=3, size=2, col=4) +
        scale_color_manual(values=c("black", 4)) +
        scale_shape_manual(values=c(1, 3)) +
        labs(x="Distance from origin (km)", y="Age (yr BP)") +
        theme_classic() +
        theme(legend.position=c(0.85,0.95), legend.title=element_blank(), axis.text=element_text(color="black"))
    #plot(dist, dates$bp, xlab="Distance from origin (km)",
    #     ylab="Age (cal BP)", pch=20, cex=1.5, cex.axis=1.5, cex.lab=1.5)
    #points(dist, simbp)
}

#compareDates(simDates, DATES.m, ORIGIN, res$maxScores[20])

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
        #genomes[i,] <- c(rnorm(numGenes, mean=1), Inf)
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
    clusterExport(cl, varlist=c("RASTER_FILE", "BIOMES", "BIOMES.sp", "ORIGIN", "START", "DATES",
                                "testModel", "simulateDispersal", "costSurface", "setGRASS"),
                                envir=environment())
    clusterEvalQ(cl, use_sp())
    clusterEvalQ(cl, initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T))
    clusterEvalQ(cl, execGRASS('g.proj', flags=c('c'), georef=RASTER_FILE))
    clusterEvalQ(cl, execGRASS('r.import', input=RASTER_FILE, output='DEM'))
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
                    reclass_matrix <- cbind(1:numGenes, as.numeric(x[1,1:numGenes]))
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
    numGenes <- 12
    numGenomes <- 500
    numParents <- 250
    numElite <- 50
    mutationRate <- 0.2
    numIter <- 20

    start_time <- Sys.time()
    res <- GA(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, cores=10)
    total_time <- Sys.time() - start_time
    cat("Completed in",total_time[[1]],attributes(total_time)$units,"\n")

    best <- as.numeric(res$genomes[1,])
    reclass_matrix <- cbind(1:numGenes, best[1:numGenes])
    costRaster <- reclassify(BIOMES, reclass_matrix)
    costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
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
    save(res, simDates, file=paste("results/ga", timeStamp, ".RData", sep=""))
    writeRaster(simDates.r, paste("results/sim", timeStamp, ".tif", sep=""))

    plot(simDates.ll)
}

#main()