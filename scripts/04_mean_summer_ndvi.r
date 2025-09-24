# Script to calculate median NDVI over summer months in study area
# Jenny Hansen
# 17 September 20205

# working in TidligVarsling project


# Load required libraries -------------------------------------------------

library(terra)
library(rgrass)

# connect to GRASS mapset
NinaR::grassConnect()

# Setup GRASS environment -------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid", flags = "p")
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")


# Import data -------------------------------------------------------------

# blank raster over SA (raster created in 02_rasterize_ssb_grid_over_sa.r)
template_rast <- rast("raster/sa_grid_250m_fixed.tif")

# NDVI was created in the following GEE script:
# https://code.earthengine.google.com/093e3a97ed46dbbdc5f69f099321ebe0
ndvi_rast <- rast("raster/ndvi_summer.tif")
ndvi_rast <- project(ndvi_rast, template_rast)

# Align with GRASS rasters ------------------------------------------------

write_RAST(ndvi_rast, "mean_ndvi_tv_study")

ndvi_rast <- read_RAST("mean_ndvi_tv_study")

cols <- c("brown", "yellow", "darkgreen")
plot(ndvi_rast, col = colorRampPalette(cols)(100))

# overwrite original raster
writeRaster(ndvi_rast, "raster/ndvi_summer.tif", overwrite = TRUE)
