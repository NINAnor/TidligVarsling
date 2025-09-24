# script to calculate forest edge landscape metrics in grid cells
# Jenny Hansen
# 24 September 2025

# working in TidligVarsling projects

# Load required libraries -------------------------------------------------

library(terra)
library(rgrass)

# connect to GRASS mapset
NinaR::grassConnect()

# Setup GRASS environment -------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid", res = "3", flags = c("p"))
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")

# Import data -------------------------------------------------------------

# was created in HELPER_rasterize_ar5.r
ar5 <- rast("raster/ar5_raster_3m.tif")
ar5 <- project(ar5, "EPSG:25833")

# write to grass mapset
write_RAST(ar5, "ar5_raster_3m")

# Prep/configuration ------------------------------------------------------

# NB: an attempt was made to calculate edge density using the
# landscapemetrics package, but it repeatedly crashed. Going with
# the grass equivalent instead:
# https://grass.osgeo.org/grass-stable/manuals/r.li.edgedensity.html

# create forest mask
execGRASS("r.mapcalc",
          expression = "forest_mask = if(ar5_raster_3m == 30, 1, null())",
          flags = "overwrite")

# I manually created a config file using g.gui.lisetup with
# forest_mask as the map area, moving window, square, and 
# cell height and width = 83

# Get edge density --------------------------------------------------------

# get edge density per 250 m grid cell
execGRASS("r.li.edgedensity",
          parameters = list(
            input  = "forest_mask",
            config = "forest_edge_250m",
            output = "forest_edgedens_250m",
            patch_type = "1"),
          flags = "overwrite")

# Read into R, export -----------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid")

forest_edgedens <- read_RAST("forest_edgedens_250m")
plot(forest_edgedens)

writeRaster(forest_edgedens, "raster/forest_edge_density.tif")


