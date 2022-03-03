library(raster)
library(rcarbon)
library(dplyr)

filterDates <- function(sites, radius) {
    x <- c(colnames(coordinates(sites)))[1]
    y <- c(colnames(coordinates(sites)))[2]

    clusters <- zerodist(sites, zero = radius, unique.ID = TRUE)

    sites$clusterID <- clusters
    sites.df <- as.data.frame(sites)
    sites.max <- as.data.frame(sites.df %>% group_by(clusterID) %>%
                               top_n(1, bp))

    xy <- cbind(sites.max[[x]], sites.max[[y]])

    sites.max.spdf <- SpatialPointsDataFrame(xy, sites.max)
    proj4string(sites.max.spdf) <- proj4string(sites)

    return(sites.max.spdf)
}

ffilterDates <- function(sites, radius) {
	return (filterDates(filterDates(sites, radius), radius))
}