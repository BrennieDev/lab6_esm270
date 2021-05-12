library(here)
library(tidyverse)
library(raster)


## Function to output a binary raster based on a user-given quantile (default is top 20%) ###
reclassify_topx <- function(rast,quant=0.8) {
  topx <- quantile(rast,quant) #find the 80% quantile of the raster values
  maxVal <- cellStats(rast,max) #find the maximum
  rcl <- c(-Inf,topx,0,
           topx,maxVal,1) # reclassify matrix (see help file for ?reclassify)
  out <- reclassify(rast,rcl=rcl)
  print(paste("max val:", maxVal))
  print(topx)
  return(out) # returns the new binary raster
}

## Function to output the intersected hotspots raster
calculate_hotspots <- function(richness_raster, impact_raster, save_with = NULL) {
  
  # Get top 25% of richness values
  richness_raster[richness_raster <= 0] <- NA
  richness_top_25 <- reclassify_topx(richness_raster, quant=0.75)
  
  # Get top 25% of impact values
  impact_raster[impact_raster <= 0] <- NA
  impact_top_25 <- reclassify_topx(impact_raster, quant=0.75)
  
  # Resample richness to higher resolution
  richness_top_25 <- resample(richness_top_25,impact_raster,method='ngb',progress='text')
  
  # Overlay rasters to get regions that are in both top 25% of impact and top 25% richness
  hotspots <- overlay(richness_top_25,impact_top_25,fun=function(x,y){x*y})
  
  # Set any values of 0 to NA
  hotspots <- reclassify(hotspots, c(-Inf, 0, NA), right=TRUE)
  
  if (!is.null(save_with)) {
    writeRaster(richness_top_25, paste(save_with,"_richness_25.tif"), overwrite = TRUE)
    writeRaster(impact_top_25, paste(save_with,"_impact_25.tif"), overwrite = TRUE)
  }
  return(hotspots)
}

### Sebastes Analysis ###

# Set our probability threshold for a species to be considered "present"
presence_threshold = 0.5

# Get list of all the sebastes files
files <- list.files(path=here("sebastes"), pattern="*.tif", full.names=TRUE, recursive=FALSE)

# Read in first sebastes raster
rockfish_raster <- raster(files[1]) %>% reclassify(rcl = c(-Inf, presence_threshold, 0, presence_threshold,1,1))

# Read in each of the other sebastes rasters
for (i in c(1:length(files))) {
  # Reclassify as present (1) or not present (0) based on 0.75 probability threshold
  temp <- raster(files[i]) %>% reclassify(rcl = c(-Inf, presence_threshold, 0, presence_threshold,1,1))
  # Also set any NAs to 0
  temp[is.na(temp[])] <- 0 
  
  # Add this raster to the cumulative sum of other sebastes presence rasters
  rockfish_raster <- sum(rockfish_raster, temp)
}

#writeRaster(rockfish_raster, "rockfish_richness.tif", overwrite = TRUE)

# Calculate rockfish+fishing hotspots
files <- list.files(path=here("lab6"), pattern="*.tif", full.names=TRUE, recursive=FALSE)

# Read in first fishing raster
fishing_raster <- raster(files[1])
fishing_raster[is.na(fishing_raster[])] <- 0
fishing_raster <- reclassify(fishing_raster, c(-Inf, 0, 0), right=TRUE)

# Read in each of the other commercial fishing rasters
for (i in c(1:length(files))) {
  temp <- raster(files[i])
  # Also set any NAs to 0
  temp[is.na(temp[])] <- 0 
  temp <- reclassify(temp, c(-Inf, 0, 0), right=TRUE)
  
  
  # Add this raster to the cumulative sum of other commercial fishing threat rasters
  fishing_raster <- sum(fishing_raster, temp)
}

writeRaster(fishing_raster, "all_fishing_threat.tif", overwrite = TRUE)

# Calculate rockfish+all fishing hotspots
#threats_rf <- raster(here("fishing_hotspots_all.tif"))
rockfish_fishing_hotspots <- calculate_hotspots(rockfish_raster, fishing_raster, save_with = "all_fishing")
writeRaster(rockfish_fishing_hotspots, "fishing_hotspots_all.tif", overwrite = TRUE)

# Read in first commercial fishing raster (now excluding destructive)
no_d_fishing_raster <- raster(files[2])
no_d_fishing_raster[is.na(no_d_fishing_raster[])] <- 0
no_d_fishing_raster <- reclassify(no_d_fishing_raster, c(-Inf, 0, 0), right=TRUE)

# Read in each of the other commercial fishing rasters
for (i in c(3:length(files))) {
  temp <- raster(files[i])
  # Also set any NAs to 0
  temp[is.na(temp[])] <- 0 
  temp <- reclassify(temp, c(-Inf, 0, 0), right=TRUE)
  
  
  # Add this raster to the cumulative sum of other commercial fishing threat rasters
  no_d_fishing_raster <- sum(no_d_fishing_raster, temp)
}

writeRaster(no_d_fishing_raster, "no_d_fishing_threat.tif", overwrite = TRUE)


# Calculate rockfish+no more destructive fishing hotspots
rockfish_no_d_hotspots <- calculate_hotspots(rockfish_raster, no_d_fishing_raster, save_with = "no_d")
writeRaster(rockfish_no_d_hotspots, "fishing_hotspots_no_d.tif", overwrite = TRUE)

hotspot_intersect <- sum(rockfish_no_d_hotspots, rockfish_fishing_hotspots)
writeRaster(hotspot_intersect, "fishing_hotspots_intersect.tif", overwrite = TRUE)


