# Script to import 250 m SSB rutenett and assign sampling grids to ruter
# Jenny Hansen
# 20 August 2024

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(mapview)
library(terra)
library(rgrass)

# connect to GRASS mapset
NinaR::grassConnect()

# Import data -------------------------------------------------------------

# created in 00_explore_data.r
samp_grid <- st_read("vector/cleaned_grid.geojson") %>% 
  group_by(name) %>% 
  slice(1) %>%  
  ungroup() %>% 
  st_make_valid()

samp_grid_cents <- samp_grid %>% 
  st_centroid()

# study are based on 250 m rutenett from SSB (provided by Ida)
rutenett_250 <- st_read("vector/ssb250_studyarea.shp") 

# Assign sampling grid to 250 m rutenett ----------------------------------

filtered_rn <- rutenett_250 %>% 
  st_filter(samp_grid, .predicate = st_contains) %>% 
  select(SSBID)

samp_grid_250 <- filtered_rn %>% 
  st_join(samp_grid, join = st_intersects)
mapview(samp_grid_250, col.regions = "orange") + 
  mapview(samp_grid, col.regions = "dodgerblue")


# Rasterize pred_nett -----------------------------------------------------

template <- rast(vect(rutenett_250), res = 250)
pred_nett_raster <- terra::rasterize(vect(rutenett_250), template,
                                     'SSBID')
plot(pred_nett_raster)


# Rasterize sampling grid -------------------------------------------------

samp_grid_raster <- terra::rasterize(vect(samp_grid_250), template,
                                     'SSBID')
plot(samp_grid_raster)

# Export files ------------------------------------------------------------

st_write(samp_grid_cents, "vector/sampling_grid_centroids.geojson")
writeRaster(pred_nett_raster, "raster/prediction_area_grid_raster.tif")
writeRaster(samp_grid_raster, "raster/sampling_grid_raster.tif")


# set correct crs for working with GRASS
gdal_utils(util = "translate",
           source = "raster/prediction_area_grid_raster.tif",
           destination = "raster/sa_grid_250m_fixed.tif",
           options = c("-a_srs", "EPSG:25833"))

sa_raster_fixed <- rast("raster/sa_grid_250m_fixed.tif")

# remove previous GRASS mask
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")

# write sampling grid to grass mapset
write_RAST(sa_raster_fixed, "tv_sa_250_grid")
