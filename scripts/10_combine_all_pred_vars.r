# Script to create a raster stack for sampling associated with ensemble modeling
# Jenny Hansen
# 25 September 2025

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(terra)

# Import data -------------------------------------------------------------

# from script 03:
rast_stack <- rast("raster/predictor_stack.tif")

# from script 04:
ndvi_summer <- rast("raster/ndvi_summer.tif")
names(ndvi_summer) <- "ndvi_summer"

# from script 05:
forest_nb <- rast("raster/neighbor_prct_forest.tif")
open_nb <- rast("raster/neighbor_prct_open.tif")
imperv_nb <- rast("raster/neighbor_prct_imperv.tif")  

# from script 06:
coldest_temp <- rast("raster/coldest_winter_temp.tif")

# from script 07:
enebolig_density <- rast("raster/enebolig_density.tif")
names(enebolig_density) <- "enebolig_density"

# from script 08:
avfall <- rast("raster/distance_to_avfall.tif")

# from script 09b:
edge_density <- rast("raster/forest_edge_density.tif")


# from Andrew (can be used if desired):
ndvi_sd <- rast("raster/NDVIstdev.tif")
ndvi_max <- rast("raster/NDVImax.tif")
ndwi_median <- rast("raster/NDWImedian.tif")
prct_industrial <- rast("raster/percentage_industrial.tif")
prct_roads <- rast("raster/percentage_roads.tif")
prct_water <- rast("raster/percentage_water.tif")
prct_buildings <- rast("raster/percentage_buildings.tif")
prct_forest <- rast("raster/percentage_forest.tif")
prct_ag <- rast("raster/percentage_agriculture.tif")

# Prepare for stacking ----------------------------------------------------

# fix extent issue
ndvi_sd      <- resample(ndvi_sd, rast_stack)
ndvi_max     <- resample(ndvi_max, rast_stack)
ndwi_median  <- resample(ndwi_median, rast_stack)

# fix CRS issue
epsg25833 <- "EPSG:25833"

# apply to all rasters
rasters <- list(rast_stack, ndvi_summer, forest_nb, open_nb, imperv_nb,
                coldest_temp, enebolig_density, avfall, edge_density,
                ndvi_sd, ndvi_max, ndwi_median,
                prct_industrial, prct_roads, prct_water,
                prct_buildings, prct_forest, prct_ag)

rasters <- lapply(rasters, \(r) { crs(r) <- epsg25833; r })

# this straggler!
crs(prct_industrial) <- "EPSG:25833"

# double check!
for (nm in names(rasters)) {
  crs(rasters[[nm]]) <- "EPSG:25833"
}

all_rasters <- do.call(c, rasters)

# additional NA handling

template <- rast_stack[[1]]   # first layer, 250m resolution
landmask <- !is.na(template)

density_layers   <- c("edge_density", "enebolig_density")
distance_layers  <- c("distance_to_public_road", "distance_to_private_road",
                      "distance_to_railway", "distance_to_housing",
                      "distance_to_industrial_area", "distance_to_garden_center",
                      "distance_to_lumberyard", "distance_to_avfall")
continuous_layers <- setdiff(names(all_rasters), c(density_layers, distance_layers))

# clean/alignment function
clean_layer <- function(r, name, template, landmask) {
  r <- resample(r, template, method = "bilinear")
  crs(r) <- crs(template)
  r <- mask(r, landmask)
  
  if (name %in% density_layers) {
    r[is.na(r)] <- 0
  } else if (name %in% distance_layers) {
    r <- focal(r, w = 3, fun = mean, na.policy = "only", na.rm = TRUE)
  }
  
  names(r) <- name  
  return(r)
}


all_cleaned <- lapply(names(all_rasters), function(nm) {
  clean_layer(all_rasters[[nm]], nm, template, landmask)
})
all_rasters <- rast(all_cleaned)


# export stack
writeRaster(all_rasters, "raster/complete_prediction_stack.tif",
            overwrite = TRUE)
