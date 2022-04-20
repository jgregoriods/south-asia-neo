source("src/src.R")

# Parameters for the genetic algorithm. We will be using 10 terrain classes,
# a population of 500, 250 parents to be transmitted to the next generation,
# 50 elite (preserved without mutation), 2% chance of mutation and 20 generations.
numGenes <- 10
numGenomes <- 500
numParents <- 250
numElite <- 50
mutationRate <- 0.2
numIter <- 20


# Testing the sites of Mureybet and Dhra as potential origins
for (originSite in c("Mureybet", "Dhra")) {
    # See src for documentation of genetic algorithm
    res <- GA(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, originSite, cores=10)

    # Get the values from the best model and generate a friction surface
    best <- as.numeric(res$genomes[1,])
    reclassMatrix <- cbind(1:numGenes, best[1:numGenes])
    costRaster <- reclassify(BIOMES, reclassMatrix)
    costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
    names(costSPDF@data) <- c("val")

    # Load cost surface in GRASS with a resolution of 25 km. See src for
    # documentation of the GRASS interface functions.
    setGRASS(costSPDF, 25000)

    # Calculate the speeds (inverse of the cost)
    speedRaster <- 1 / costRaster
    speedRaster.wgs <- projectRaster(speedRaster, res=0.25, crs=WGS, method="ngb")
    speedRaster.cropped <- crop(speedRaster.wgs, extent(DATES))

    # Coordinates and start yr BP of the chosen origin site
    ORIGIN <- as.numeric(DATES.m[DATES.m$Site==originSite,]@coords)
    START <- DATES.m[DATES.m$Site==originSite,]$bp

    # Simulate arrival times. See src for function documentation
    simDates <- simulateDispersal(ORIGIN, START)
    simDates.wgs <- projectRaster(simDates, res=0.25, crs=WGS)
    simDates.r <- crop(simDates.wgs, extent(DATES))

    cat(paste("Best score:", best[numGenes+1], "\n"))

    # Write files
    # Rasters used in QGIS to produce the maps in the paper
    speeds <- data.frame("region"=1:numGenes, "speed"=1/best[1:numGenes])
    writeRaster(speedRaster.cropped, paste("results/speed_", originSite, ".tif", sep=""), overwrite=T)
    write.csv(speeds, file=paste("results/res_", originSite, ".csv", sep=""))
    save(res, simDates, file=paste("results/ga_", originSite, ".RData", sep=""))
    writeRaster(simDates.r, paste("results/sim_", originSite, ".tif", sep=""), overwrite=T)

    # Fig 3 and S1 Fig
    print(plotDates(simDates, DATES.m, ORIGIN, best[11]))
    try(dev.print(tiff, paste("results/plot_", originSite, ".tiff", sep=""), res=300, width=1280))
    try(dev.off())

    # S3-S4 Fig
    par(bg="white")
    plotGA(res, paste("Origin:", originSite))
    try(dev.print(tiff, paste("results/scores_", originSite, ".tiff", sep=""), res=300, width=2000))
    try(dev.off())
}
