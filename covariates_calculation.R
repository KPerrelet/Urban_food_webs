################################################################################
############################ covariates calculation ############################
################################################################################

library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(sf)
library(exactextractr)
library(landscapemetrics)
library(units)
library(lidaRtRee)
library(lidR)

#This function calculates the shortest distance to the forest
calculate_distance_forest <- function(Site){ 
  df <- data.frame() # Initialize empty dataframe to store distances
  for (i in c(1:dim(Site)[1])){ # Loop through each site
    distmin <- min(st_distance(Site[i, ], HabitatMap.sf[HabitatMap.sf$lrtyp1text=="Wälder", ])) # Calculate minimum distance to forest
    df<-rbind(df, distmin) # Append distance to dataframe
  }
  return(df)
}

calculate_distance_pond <- function(Site){ 
  df <- data.frame() # Initialize empty dataframe to store distances
  for (i in c(1:dim(Site)[1])){ # Loop through each site
    dist <- (st_distance(Site[i, ], HabitatMap.sf[HabitatMap.sf$layer == "gsz_stillgewaesserinventar copy", ])) # Calculate minimum distance to ponds
    NoZero <- data.frame() 
    for (i in c(1:length(dist))){
      if (dist[i]>units::set_units(0,m)) NoZero<-rbind(NoZero,dist[i]) # Remove zero distances
    }
    distmin <- min(NoZero) # Determine shortest nonzero distance
    df<-rbind(df, distmin) # Append to dataframe
  }
  return(df)
}

# Load data
sites.sf <- st_read("ponds.shp")
HabitatMap.sf <- st_read(paste0(path, "R/Data/HabitatMap/data/HabitatMapPonds.shp"))

# Rename second column for better readability (otherwise truncated)
colnames(sites.sf)[2] <- "natural_banks"
HabitatMap.sf[HabitatMap.sf$layer == "gsz_stillgewaesserinventar copy", 6] <- "Gewässer"

# Classify habitat types into "Blue", "Green", and "Grey"
Matrix <- c("Gebäude-, Verkehrs- und Industrieflächen", "Gewässer")
HabitatMap.sf$Habitat <- 
  ifelse(HabitatMap.sf$lrtyp1text %in% Matrix, ifelse(HabitatMap.sf$lrtyp1text=="Gewässer", "Blue", "Grey"), "Green")

# Rasterization of the habitat map
grd <- rast(ext(HabitatMap.sf), resolution = 1, crs = "EPSG:2056")
HabitatMatrix.rast <- terra::rasterize(vect(HabitatMap.sf), grd, field = "Habitat")
Blue.rast <- terra::rasterize(vect(HabitatMap.sf[HabitatMap.sf$layer == "gsz_stillgewaesserinventar copy", ]), grd)

# 1) LULC covariates
# Create 50m buffer around sites
local.sf <- st_buffer(sites.sf$geometry, dist = 50, unit = 'meter', ratio = 0.1)
# Retrieve the area and fraction covered by blue, green, and grey surfaces
local.lulc <- cbind(as.data.frame(t(exact_extract(HabitatMatrix.rast, local.sf,  function(value, cov_frac) {table(value)}))),
                    exact_extract(HabitatMatrix.rast, local.sf, fun = "frac"))

# Create 500m buffer around sites
landscape.sf <- st_buffer(sites.sf$geometry, dist = 500, unit = 'meter')
landscape.lulc <- cbind(as.data.frame(t(exact_extract(HabitatMatrix.rast, landscape.sf,  function(value, cov_frac) {table(value)}))),
                        exact_extract(HabitatMatrix.rast, landscape.sf, fun = "frac"))

# Assign column names
colnames(local.lulc) <- c("blue50", "green50", "grey50", "blue50_frac", "green50_frac", "grey50_frac")
colnames(landscape.lulc) <- c("blue500", "green500", "grey500", "blue500_frac", "green500_frac", "grey500_frac")

# Combine extracted LULC metrics with site data
lulc.df <- cbind(sites.sf$site_no,
                 local.lulc,
                 landscape.lulc
)
colnames(lulc.df)[1] <- "site_no"
sites.sf <- merge(sites.sf, lulc.df, by = "site_no")

# 2) Connectivity
# Calculate and add distance to the forest and to the closest pond
distance_forest <- calculate_distance_forest(as.data.frame(sites.sf$geometry))
distance_pond <- calculate_distance_pond(as.data.frame(sites.sf$geometry))
sites.sf$distance_forest <- distance_forest[, 1]
sites.sf$distance_pond <- distance_pond[, 1]

# Compute terrestrial connectivity using landscape cohesion index
terrestial_connectivity <- sample_lsm(HabitatMatrix.rast, 
                                      y = st_centroid(sites.sf$geometry), 
                                      plot_id = sites.sf$site_no,
                                      what = c("lsm_c_cohesion"), # Selecting cohesion metric
                                      shape = "circle", # Defining sampling shape
                                      size = 500) # Sampling radius in meters

# Filter data to include only green patches (class = 1)
terrestial_connectivity <- terrestial_connectivity[(terrestial_connectivity$class==1),]
# Retain only relevant columns (site ID and cohesion values)
terrestial_connectivity <- terrestial_connectivity[, c(5:7)]
# Transform data from long to wide format and merge it with site information
terrestial_connectivity <- terrestial_connectivity %>%
  pivot_wider(names_from = metric, values_from = value)
colnames(terrestial_connectivity) <- c("site_no", "terrestrial_connectivity")
sites.sf <- merge(sites.sf, terrestial_connectivity, by = "site_no")

# Compute aquatic connectivity (average distance between ponds in a 500 m buffer)
aquatic_connectivity <- sample_lsm(Blue.rast, 
                                  y = st_centroid(sites.sf$geometry), 
                                  plot_id = as.factor(sites.sf$site_no),
                                  what = c("lsm_c_enn_mn"), # Mean nearest-neighbor distance
                                  shape ="circle",
                                  size = 500)

aquatic_connectivity <- aquatic_connectivity[, c(6:7)]
colnames(aquatic_connectivity) <- c("aquatic_connectivity", "site_no")
aquatic_connectivity[is.na(aquatic_connectivity$aquatic_connectivity), "aquatic_connectivity"] <- 500

sites.sf <- merge(sites.sf, aquatic_connectivity, by = "site_no")

# 3) Habitat diversity 
# Based on joint entropy (join), which represents the uncertainty in determining the category of a focus cell
# high entropy = small diversity of values in co-occurence matrix
# low entropy = high diversity of values in co-occurence matrix

# Rasterize the habitat map based on a finer level of details (62 habitat types for green spaces only)
HabitatDiversity.rast <- terra::rasterize(vect(HabitatMap.sf[HabitatMap.sf$Habitat == "Green", ]), 
                                          grd, 
                                          field = "lrtyp2text")

# Loop through each site to compute habitat diversity metrics
sites.sf$complementarity <- NA
for (site_i in 1:nrow(sites.sf)) {
  landscape.rast <- terra::crop(HabitatDiversity.rast, landscape.sf, mask = T)
  sites.sf[site_i, "complementarity"] <- lsm_l_joinent(landscape.rast)$value
}

# 4) Climate 
# The climate and NO2 data is available at https://data.stadt-zuerich.ch/dataset/geo_klimadaten

overwarming.rast <- rast(paste0(path, "/R/data/Climate/data/22_Local Overwarming Zurich.tif"))
# overwarming.rast <- rast("overwarming.tif")
no2.rast <- rast(paste0(path, "/R/data/Climate/data/44_Jahresmittel der NO2-Konzentration Stadt.tif"))
# no2.rast <- rast("no2.tif")

# Extract climate variables for each site based on centroid location
sites.sf$overwarming <- terra::extract(overwarming.rast, st_centroid(sites.sf), method = "bilinear")[, 2]
sites.sf$no2 <- terra::extract(no2.rast, st_centroid(sites.sf), method = "bilinear")[, 2]

# 4) Population density
# The population density data is available at https://www.stadt-zuerich.ch/geodaten/download/63

pop.sf <- st_read(paste0(path, "/R/Data/Population/Raumliche_Bevolkerungsstatistik_-OGD/BEVOELKERUNG_HA_F.shp"))
# pop.sf <- st_read("population.shp")

# Transform the data (shp file) to a raster as it is grid-based
pop.rast <- terra::rasterize(vect(pop.sf), grd, field = "PERS_N", background = 0)
pop.rast[pop.rast < 0] <- 0

# Extract population density for each site
sites.sf$population <- terra::extract(pop.rast, st_centroid(sites.sf), method = "simple")[, 2]

# 5) LiDAR data
# The LiDAR data is available at https://www.swisstopo.admin.ch/en/lidar-data-swisstopo
# This code snippet is based on the lidRtree tutorial, available at https://lidar.pages.mia.inra.fr/lidaRtRee/articles/forest.structure.metrics.html

# Set up parallel processing
future::plan("multisession", workers = 4L)
options(future.rng.onMisuse = "ignore")

# Load LiDAR data and 
cata <- lidR::catalog("D:/LIDAR_ZH/")
# cata <- lidR::catalog("lidar_foler/")
# Set projection
lidR::projection(cata) <- 2056
# Filter for the city of Zurich
crop_extent <- HabitatMap.sf
crop_extend500 <- st_buffer(crop_extent, 500)
cata <- catalog_intersect(cata, crop_extend500) 

# Define LiDAR processing parameters
resolution <- 5 # Output resolution in meters
res_chm <- 0.2 # Canopy height model resolution
buffer_size <- 20 # Processing buffer for accurate edge calculations
points_max_h <- 60 # Maximum valid height for vegetation points
class_points <- c(2,3, 4, 5) # Vegetation and ground classification
class_ground <- 2 # Ground class

# Fonction to computed raster statistics from multi-scale smoothing
raster_stats_rectangle <- function(x) {
  data.frame(
    CHM_sd = sd(x$smoothed_image_0[x$smoothed_image_0 > 0]), #Standard deviation CHM
    CHM_mean = mean(x$smoothed_image_0[x$smoothed_image_0 > 0]), #Mean CHM
    CHM_mean_tree = mean(x$smoothed_image_0[x$smoothed_image_0 >= 2.5]), #Mean CHM > 2.5m 
    CHM_mean_grasses = mean(x$smoothed_image_0[x$smoothed_image_0 < 1 & x$smoothed_image_0 > 0]), #Mean CHM < 1m 
    CHM_mean_shrubs = mean(x$smoothed_image_0[x$smoothed_image_0 >= 1 & x$smoothed_image_0 < 2.5]),#Mean CHM between 1m and 2.5m
    CHM_Cover = sum(x$smoothed_image_0 > 0)*mf_rectangle, #Proportion of surface covered by vegetation
    CHM_Cover_tree = sum(x$smoothed_image_0 >= 2.5)*mf_rectangle, #Proportion of surface covered by trees
    CHM_Cover_grasses = sum(x$smoothed_image_0 < 1 & x$smoothed_image_0 > 0)*mf_rectangle, #Proportion of surface covered by grasses
    CHM_Cover_shrubs = sum(x$smoothed_image_0 >= 1 & x$smoothed_image_0 < 2.5)*mf_rectangle #Proportion of surface covered by shrubs
  )
}

# Define normalization factor
mf_rectangle <- 100 / (resolution / res_chm)^2

# Define buffer size
buffer_pond <- 50

# Initialize empty dataframe for results
lidar.df <- data.frame()

# Process each site in the study area
for (i in c(1:nrow(sites.sf))) {
  print(i) # Print progress indicator
  
  # Clip LiDAR data to a buffered site polygon
  Site_buffer <- try(lidR::clip_polygon(
    cata, 
    st_coordinates(st_buffer(sites.sf[i,], buffer_pond))[, 1],
    st_coordinates(st_buffer(sites.sf[i,], buffer_pond))[, 2]))

  # Exclude industrial, urban, and water areas
  Site_buffer <- merge_spatial(Site_buffer, 
                               crop_extent[!c(crop_extent$lrtyp1text == "Gebäude-, Verkehrs- und Industrieflächen" | crop_extent$lrtyp1text == "Gewässer"), ], attribute = "in_poly")
  
  # Filter only vegetation points within the valid height range
  Site_buffer <- lidR::filter_poi(Site_buffer,is.element(Classification, class_points) & Z <= points_max_h)
  
  # Classify low vegetation
  Site_buffer$Classification[Site_buffer$in_poly] <- LASLOWVEGETATION
  
  # Normalize heights; ensure no negative heights
  Site_buffer$Z[Site_buffer$Z < 0] <- 0
  
  # Generate a Canopy Height Model (CHM) raster
  chm <- lidR::rasterize_canopy(Site_buffer, res = res_chm, algorithm = lidR::p2r(),
                                pkg = "terra")
  chm[chm < 0 | chm > points_max_h] <- 0
  
  # Detect individual tree tops
  segms <- lidaRtRee::tree_segmentation(chm, hmin = 0.9)
  # Extract tree structure metrics
  segms.df <- lidaRtRee::tree_extraction(segms$filled_dem, segms$local_maxima, segms$segments_id)

  # Filter valid vegetation
  vegetation <- segms.df[segms.df$h > 1, ]
  
  if (nrow(vegetation) > 0) {
    # Compute vegetation metrics at the defined resolution
    metrics_all <- lidaRtRee::raster_metrics(
      vegetation[,-1],
      res = resolution,
      fun = function(x) {
        lidaRtRee::std_tree_metrics(x, resolution^2 / 10000)
      },
      output = "dataframe"
    )
    if (nrow(metrics_all>1))
      metrics_all <-rbind(metrics_all, colMeans(metrics_all, na.rm = TRUE)) # Aggregate data
  }
  colnames(metrics_all) <- c("X",
                             "Y",
                             "vegetation_meanH",
                             "vegetation_height_sd",
                             "vegetation_giniH",
                             "vegetation_density",
                             "vegetationInf10_density",
                             "vegetationSup10_density",
                             "vegetationSup20_density",
                             "vegetationSup30_density",
                             "vegetation_meanCrownSurface",
                             "vegetation_meanCrownVolume",
                             "vegetationCanopy_meanH")
  # Define result structure
  site.res <- data.frame(sites.sf[i,]$site_no)
  colnames(site.res) <- "site_no"
  if (nrow(vegetation)>0) site.res <- cbind (site.res,metrics_all[nrow(metrics_all),])
  
  # Remove unnecessary columns
  remove <- c("X",
              "Y",
              "vegetation_meanH",
              "vegetation_giniH",
              "vegetationInf10_density",
              "vegetationSup10_density",
              "vegetationSup20_density",
              "vegetationSup30_density",
              "vegetation_meanCrownSurface",
              "vegetation_meanCrownVolume",
              "vegetationCanopy_meanH")
  site.res <- site.res[, -which(names(site.res) %in% remove)]
  
  lidar.df <- rbind(lidar.df, site.res) # Append results
}

# Rename columns based on buffer size
colnames(lidar.df)[2:ncol(lidar.df)] <- paste(colnames(lidar.df)[2:ncol(lidar.df)], buffer_pond, sep = "_")

# Merge results with site data
sites.sf <- merge(sites.sf, lidar.df, by = "site_no")

# Remove columns
sites.sf$geometry <- NULL
remove_covariates <- c("green500", "green50_frac", "green500_frac", 
                       "blue50", "blue500", "blue50_frac", "blue500_frac", 
                       "grey50", "grey500", "grey50_frac", 
                       "geometry")

sites.sf <- sites.sf %>%
  select(-any_of(remove_covariates))

# Save results
write.csv(sites.sf, "covariates2.csv", row.names = F)

