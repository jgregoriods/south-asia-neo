library(parallel)
library(ggplot2)
library(raster)
library(rgrass7)


# to use sp objects with rgrass7
use_sp()

# for reproducibility
set.seed(99)

# define projections
WGS <- CRS("+init=epsg:4326")
ALBERS <- CRS("+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# read files
RASTER_FILE <- "layers/biomes.tif"
DATES <- read.csv("dates/dates.csv")

coordinates(DATES) <- ~Longitude+Latitude
proj4string(DATES) <- WGS
DATES.m <- spTransform(DATES, ALBERS)

BIOMES <- raster(RASTER_FILE)
BIOMES.sp <- as(BIOMES, "SpatialPixelsDataFrame")
names(BIOMES.sp@data) <- c("val")

# Import the biome file as a raster map "DEM" and use it to set the
# projection and region in GRASS
initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T)
execGRASS('g.proj', flags=c('c'), georef=RASTER_FILE)
execGRASS('r.import', input=RASTER_FILE, output='DEM')
execGRASS('g.region', raster='DEM')


# -----------------------------------------------------------------------------
# Functions for calculating cost surfaces and simulated arrival times
# -----------------------------------------------------------------------------

#' Set the current geographic region in GRASS. A raster map "GRID" is created
#' in GRASS from the specified raster object.
#' @param RASTER A raster object to calculate the region from.
#' @param RES The raster cell size.
#' @return Void.
#'
setGRASS <- function(RASTER, RES) {
    writeRAST(RASTER, "GRID", zcol="val", flags=c("quiet", "overwrite"))
    execGRASS("g.region", flags=c("d"), res=as.character(RES), raster="GRID")
}

#' Calculate a cost surface using the specified origin coordinates.
#' @param coords A numeric vector with X,Y coordinates.
#' @return A raster object with the cost surface.
#'
costSurface <- function(coords) {
    # A raster map "GRID" with a friction surface must have been created
    # using the setGRASS function.
    execGRASS("r.cost", flags=c("k", "overwrite", "quiet"), input="GRID",
              start_coordinates=coords, output="COST", max_cost=0)
    cost <- raster(readRAST("COST"))
    return(cost)
}

#' Simulate arrival times based on a cost surface, point of origin and start date.
#' @param origin A numeric vector with X,Y coordinates.
#' @param startDate The start year BP.
#' @return A raster object with simulated arrival times in years BP.
#' 
simulateDispersal <- function(origin, startDate) {
    cost <- costSurface(origin)
    # Must be multiplied by 25, which is our raster resolution in km
    simDates <- startDate - (cost * 25)
    return(simDates)
}


# -----------------------------------------------------------------------------
# Functions for running the genetic algorithm
# -----------------------------------------------------------------------------

#' Perform crossover of two "genomes". The two are split at a random index
#' and their halves are combined.
#' @param x A numeric vector. The first parent.
#' @param y A numeric vector. The second parent.
#' @return A numeric vector encoding the new genome after crossover.
#'
crossover <- function(x, y) {
    i <- sample(1:(length(x)-1), 1)
    return(c(x[1:i], y[(i+1):length(y)]))
}

#' Mutate a "genome". A random "gene" is selected and its value is shifted
#' by a value taken from a random normal distribution.
#' @param x A numeric vector encoding the genome.
#' @return The mutated numeric vector.
#'
mutate <- function(x) {
    i <- sample(1:length(x), 1)
    x[i] <- x[i] + rnorm(1)
    return(x)
}

#' Evaluate how closely the simulated arrival times using a given
#' cost surface model approach the radiocarbon dates.
#' @param costRaster A raster object with the assigned costs.
#' @param sites A data frame with site coordinates and dates.
#' @param origin A vector of X,Y coordinates.
#' @param startDate Start of the dispersal in years BP.
#' @return The RMSE between simulated and real dates.
#'
testModel <- function(costRaster, sites, origin, startDate) {
    costSPDF <- as(costRaster, "SpatialPixelsDataFrame")
    names(costSPDF@data) <- c("val")
    # We are using a cell size of 25 km
    setGRASS(costSPDF, 25000)
    simDates <- simulateDispersal(origin, startDate)
    gc()
    sites$simbp <- extract(simDates, sites)
    sites <- sites[!is.na(sites$simbp),]
    rmse <- sqrt(sum((sites$simbp - sites$bp)^2) / nrow(sites))
    return(rmse)
}

#' Main function for running the genetic algorithm. A number of models are
#' randomly initialized, evaluated and the best are selected to be transmitted
#' to the next iteration. New models are bred with random crossovers and mutation.
#' @param numGenes Number of genes, i.e. terrain classes.
#' @param numGenomes Numer of genomes, i.e. models to be tested.
#' @param numParents Number of genomes to be preserved every generation.
#' @param numElite Number of genomes that will be preserved without mutation.
#' @param mutationRate Chance of mutations happening.
#' @param numIter Number of generations (iterations) of the algorithm.
#' @param cores Number of cores for parallel processing.
#' @param originSite A string. Name of a site in the DATES data frame to be
#' used as origin of the dispersal.
#' @return A list of vectors with the following attributes: genomes, containing
#' all the genome vectors; maxScores, containing the best score per generation;
#' and avgScores, with the average scores per generation.
#'
GA <- function(numGenes, numGenomes, numParents, numElite, mutationRate, numIter, originSite, cores=NULL) {
    startTime <- Sys.time()

    # Coordinates of the origin site and start yr BP
    ORIGIN <- as.numeric(DATES.m[DATES.m$Site==originSite,]@coords)
    START <- DATES.m[DATES.m$Site==originSite,]$bp

    # Initialize genomes with values from a random normal distribution
    # with mean = 1. The last value in the genome is its score
    genomes <- as.data.frame(matrix(nrow=numGenomes, ncol=numGenes+1))
    for (i in 1:numGenomes) {
        # Initialize scores to Inf. This is the worst possible score,
        # as the RMSE will be used for evaluation
        genomes[i,] <- c(rnorm(numGenes, mean=1), Inf)
    }

    # To store the max and average scores every generation
    maxScores <- c()
    avgScores <- c()

    # Use all cores but one if the parameter is not specified
    if (is.null(cores)) {
        ncores <- detectCores() - 1
    } else {
        ncores <- cores
    }

    # Export necessary packages and objects to the cluster
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library("sp"))
    clusterEvalQ(cl, library("raster"))
    clusterEvalQ(cl, library("rgrass7"))
    clusterExport(cl, varlist=c("RASTER_FILE", "BIOMES", "BIOMES.sp", "ORIGIN",
                                "START", "DATES", "testModel", "simulateDispersal",
                                "costSurface", "setGRASS"), envir=environment())

    # GRASS must be initialized in every core. The biome file is imported
    # as a raster map "DEM" and used to speficy the projection and region.
    clusterEvalQ(cl, use_sp())
    clusterEvalQ(cl, initGRASS("/usr/lib/grass78", home=tempdir(), mapset="PERMANENT", override=T))
    clusterEvalQ(cl, execGRASS('g.proj', flags=c('c'), georef=RASTER_FILE))
    clusterEvalQ(cl, execGRASS('r.import', input=RASTER_FILE, output='DEM'))
    clusterEvalQ(cl, execGRASS('g.region', raster='DEM'))

    cat(paste("\nRunning genetic algorithm on", ncores,
              "parallel workers.\nThis may take a while...\n"))

    for (iter in 1:numIter) {
        genomeList <- split(genomes, seq(nrow(genomes)))

        res <- parLapply(cl, genomeList, function(x) {
            # Preserve the scores of the models that were passed without mutation
            if (x[1,numGenes+1] != Inf) {
                return (x[1,numGenes+1])
            } else {
                # Assign worst score (Inf) to genomes containing negative values
                if (min(x[1,1:numGenes]) <= 0) {
                    return(Inf)
                } else {
                    # Otherwise, reclassify the biome raster using the genome
                    # values, simulate the arrival times and evaluate
                    reclassMatrix <- cbind(1:numGenes, as.numeric(x[1,1:numGenes]))
                    cost <- reclassify(BIOMES, reclassMatrix)
                    score <- testModel(cost, DATES, ORIGIN, START)
                    gc()
                    return(score)
                }
            }
        })

        genomes[,numGenes+1] <- unlist(res)

        avgScores[iter] <- mean(genomes[,numGenes+1][genomes[,numGenes+1] < Inf])

        # Order genomes from best to worst score and select the elite and parents
        # to be transmitted to the next generation
        elite <- genomes[order(genomes[,numGenes+1]),][1:numElite,]
        parents <- genomes[order(genomes[,numGenes+1]),][1:numParents,]
        parents <- parents[order(as.numeric(rownames(parents))),]

        # best score
        maxScores[iter] <- elite[1,numGenes+1]

        # crossover between parents to generate as many new genomes as
        # necessary to keep population size
        children <- as.data.frame(matrix(nrow=numGenomes-numElite, ncol=numGenes+1))
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

    totalTime <- Sys.time() - startTime
    cat("Completed in", totalTime[[1]], attributes(totalTime)$units, "\n")

    return (list(genomes=genomes, maxScores=maxScores, avgScores=avgScores))
}


# -----------------------------------------------------------------------------
# Functions for plotting the results
# -----------------------------------------------------------------------------

plotDates <- function(simRaster, dates, origin, rmse) {
    zagrosBegin <-  spDistsN1(matrix(c(5200000, 3700000), nrow=1), origin) / 1000
    zagrosEnd <-  spDistsN1(matrix(c(6500000, 2500000), nrow=1), origin) / 1000

    indusBegin <-  spDistsN1(matrix(c(8000000, 2100000), nrow=1), origin) / 1000
    indusEnd <-  spDistsN1(matrix(c(9100000, 1900000), nrow=1), origin) / 1000

    dist <- spDistsN1(dates, origin) / 1000

    merged <- rbind(data.frame("dist"=dist, "age"=dates$bp, "type"="C14"),
                    data.frame("dist"=dist, "age"=extract(simRaster, dates),
                    "type"="Simulated"))

    ggplot(merged) +
        geom_rect(aes(xmin=zagrosBegin, xmax=zagrosEnd, ymin=min(age, na.rm=T), ymax=max(age, na.rm=T)), fill="grey95") +
        geom_rect(aes(xmin=indusBegin, xmax=indusEnd, ymin=min(age, na.rm=T), ymax=max(age, na.rm=T)), fill="grey95") +

        geom_point(aes(x=dist, y=age, color=type, shape=type), size=2) +
        annotate(geom="text", x=mean(c(zagrosBegin, zagrosEnd)), y=5000, label="Zagros", color="grey") +
        annotate(geom="text", x=mean(c(indusBegin, indusEnd)), y=8000, label="Indus", color="grey") +
        scale_color_manual(values=c("black", 4)) +
        scale_shape_manual(values=c(1, 3)) +
        labs(x="Distance from origin (km)", y="Age (yr BP)") +
        ggtitle(paste("RMSE=", round(rmse), sep="")) +
        theme_classic() +
        theme(legend.position=c(0.88,0.95), legend.title=element_blank(),
              legend.background = element_blank(),
              axis.text=element_text(color="black"))
}

plotGA <- function(x, title) {
    plot(1:20, x$avgScores, type="b", col="blue", xlab="Generation", ylab="RMSE", main=title)
    lines(1:20, x$maxScores, type="b", col="red")
    legend("topright", legend=c("Average score", "Best score"), col=c("blue", "red"), lty=c(1, 1))
}
