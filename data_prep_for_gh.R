################################################################################
##### This scripts prepares and calculates a list of predictive variables ######
################################################################################
library(sf)
library(terra)
library(raster)

Ponds_CH2.shp <- st_read(paste0(path, "R/Data/Sites_CH2/Ponds_CH2_metadata.shp"))
Ponds_CH2.shp$sample <- paste("Water", Ponds_CH2.shp$Extrctn_n, sep = "")
colnames(Ponds_CH2.shp)[which(colnames(Ponds_CH2.shp) == "Sit_N_y" )] <- "Site_No"
water.sf <- Ponds_CH2.shp[Ponds_CH2.shp$Replct_N == 1, ]
water.sf <- water.sf %>% arrange(water.sf$Site_No)
water.sf <- water.sf[, c(1:9, 11, 17, 20:26, 28)]
colnames(water.sf) <- c("Uniq_cd", "Ownership", "Maintenance", "Species", "Quality", "Permanence", "FishGSZ", "State", "Address", "Smplng_d", "Site_No", "Smplng_t", 
                        "Banks_vege", "Banks_rocks", "Banks_concrete", "Banks_earth", "Marcophyte", "Walls", "Extraction_responsible","geometry")

water.sf <- water.sf %>%
  relocate(Site_No, Uniq_cd, Address, Smplng_d, Smplng_t, 
           Ownership, Maintenance, Species, Quality, Permanence, FishGSZ, State,  
           Banks_vege, Banks_rocks, Banks_concrete, Banks_earth, Marcophyte, Walls, Extraction_responsible,
           geometry)

water.sf$natural_banks <- 1 - water.sf$Banks_concrete
water.sf$marcophyte <- water.sf$Marcophyte

water.sf <- water.sf %>%
  select(c("Site_No", "natural_banks", "marcophyte", "geometry"))

write_sf(water.sf, "ponds.shp")

