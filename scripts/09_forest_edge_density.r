# trying out another alternative to getting an edge density using
# lsmetrics and coarser resolution raster
# Jenny Hansen
# 25 September 2025

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(terra)
library(rgrass)

# connect to GRASS mapset
NinaR::grassConnect()

# Setup GRASS environment -------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid", res = "10", flags = c("p"))
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")

# Import data -------------------------------------------------------------

# was created in HELPER_rasterize_ar5.r
ar5 <- read_RAST("ar5_tv_study_rast_10m")
names(ar5) <- "ar5_arealtype"

# Create binary raster ----------------------------------------------------

# create forest mask
execGRASS("r.mapcalc",
          expression = "forest_mask_10m = if(ar5_tv_study_rast_10m == 30, 1, null())",
          flags = "overwrite")

# remove session mask (if active)
try(execGRASS("r.mask", flags = "r"), silent = TRUE)

# Alternative edge density calculation ------------------------------------

# convert NULL to 0
execGRASS("r.mapcalc",
          expression = "forest_mask_10m_0 = if(isnull(forest_mask_10m), 0, forest_mask_10m)",
          flags = "overwrite")

execGRASS("r.neighbors",
          input = "forest_mask_10m_0",
          output = "forest_mask_mean10m",
          method = "average",
          size = 3,
          flags = "overwrite")

# extract edges
execGRASS("r.mapcalc",
          expression = "forest_edge_pix10m = if(forest_mask_mean10m > 0 && forest_mask_mean10m < 1, 1, null())",
          flags = "overwrite")

execGRASS("r.mapcalc",
          expression = "forest_edge_pix10m_0 = if(isnull(forest_edge_pix10m), 0, forest_edge_pix10m)",
          flags = "overwrite")


# convert to length
execGRASS("r.mapcalc",
          expression = "forest_edge_len10m = forest_edge_pix10m_0 * 10",
          flags = "overwrite")

# aggregate
execGRASS("g.region", raster = "tv_sa_250_grid")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")
execGRASS("r.resamp.stats",
          input = "forest_edge_len10m",
          output = "forest_edge_len_250m",
          method = "sum",
          flags = "overwrite")
execGRASS("r.mapcalc",
          expression = "forest_edge_density_250m = forest_edge_len_250m / 6.25",
          flags = "overwrite")

# Read into R and inspect -------------------------------------------------

forest_mask   <- read_RAST("forest_mask_10m")

edge_density <- read_RAST("forest_edge_density_250m")
plot(edge_density)
names(edge_density) <- "forest_edge_density"

forest_edge <- read_RAST("forest_edge_pix10m")

plot(forest_mask, col = "darkgreen", legend = FALSE)
plot(forest_edge, col = "red", add = TRUE, legend = FALSE)
plot(forest_edge)

execGRASS("r.univar", map = "forest_edge_density_250m")

# these values look reasonable

# write to raster folder
writeRaster(edge_density, "raster/forest_edge_density.tif", overwrite = TRUE)
