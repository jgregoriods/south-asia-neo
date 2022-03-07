source("src/src.R")


main <- function() {
    numGenes <- 12
    numGenomes <- 500
    numParents <- 250
    numElite <- 50
    mutationRate <- 0.2
    numIter <- 20

    startTime <- Sys.time()
    res <- GA(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, cores=10)
    totalTime <- Sys.time() - startTime
    cat("Completed in",totalTime[[1]],attributes(totalTime)$units,"\n")

    best <- as.numeric(res$genomes[1,])
    reclass_matrix <- cbind(1:numGenes, best[1:numGenes])
    costRaster <- reclassify(BIOMES, reclass_matrix)
    costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
    names(costSPDF@data) <- c("val")

    setGRASS(costSPDF, 25000)

    simDates <- simulateDispersal(ORIGIN, START)
    simDates.wgs <- projectRaster(simDates, res=0.25, crs=WGS)
    simDates.r <- crop(simDates.wgs, extent(DATES))

    cat(paste("Best score:", best[numGenes+1], "\n"))

    timeStamp <- round(as.numeric(Sys.time()))

    speeds <- data.frame("region"=1:numGenes, "speed"=1/best[1:numGenes])
    write.csv(speeds, file=paste("results/res", timeStamp, ".csv", sep=""))
    save(res, simDates, file=paste("results/ga", timeStamp, ".RData", sep=""))
    writeRaster(simDates.r, paste("results/sim", timeStamp, ".tif", sep=""))
}


main()