# Moving window analysis to get percent cover of neighboring grids
# Jenny Hansen
# 17 September 2025

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(terra)
library(rgrass)

# connect to GRASS mapset
NinaR::grassConnect()

# Setup GRASS environment -------------------------------------------------

# set initial res to 10; change when aggregating
execGRASS("g.region", raster = "tv_sa_250_grid", flags = "p", res = "10")
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")


# Import data -------------------------------------------------------------

# impervious density downloaded from
# https://land.copernicus.eu/en/products/high-resolution-layer-imperviousness/imperviousness-density-2018
# clipped to study extent in QGIS
imp <- rast("raster/impervious_clipped.tif")

# convert from 0-255 to 0-100 (red to integer)
imp_vals <- as.integer(values(imp))

# replace nodata code (255) with NA
imp_vals[imp_vals == 255] <- NA

# update raster
values(imp) <- imp_vals

writeRaster(imp, "raster/impervious_gray.tif", overwrite=TRUE)
write_RAST(imp, "impervious_rast", flags="overwrite")

# Percent cover forest/open -----------------------------------------------

# create raster that contains percent cover inside each cell

# forest mask at 10 m
execGRASS("r.mapcalc",
          expression = "forest_mask = if(ar5_tv_study_rast_10m == 30, 1, 0)",
          flags = c("overwrite"))

# open landscapes mask @ 10 m
execGRASS("r.mapcalc",
          expression = "open_mask = if(ar5_tv_study_rast_10m == 22 || ar5_tv_study_rast_10m == 23 || ar5_tv_study_rast_10m == 50, 1, 0)",
          flags = c("overwrite"))

# change back to 250 m
execGRASS("g.region", raster = "tv_sa_250_grid")

# forest percent cover @ 250 m
execGRASS("r.resamp.stats",
          input = "forest_mask",
          output = "forest_cover_tv_study",
          method = "average",
          flags = c("overwrite"))

# open landscape percent cover @ 250 m 
execGRASS("r.resamp.stats",
          input = "open_mask",
          output = "open_cover_tv_study",
          method = "average",
          flags = c("overwrite"))


# impervious surface cover @ 250 m
execGRASS("r.resamp.stats",
          input  = "impervious_rast",
          output = "impervious_cover_tv_study",
          method = "average",
          flags  = "overwrite")


# Moving window analysis --------------------------------------------------

# forest cover neighborhood mean
execGRASS("r.neighbors",
          input = "forest_cover_tv_study",
          output = "forest_cover_3x3",
          method = "average",
          size = 3,
          flags = "overwrite")

# remove 'focal cell' from calculation
execGRASS("r.mapcalc",
          expression = "forest_cover_8nb = ((forest_cover_3x3 * 9) - forest_cover_tv_study) / 8.0",
          flags = "overwrite")

# convert to percentage
execGRASS("r.mapcalc",
          expression = "forest_cover_8nb_pct = forest_cover_8nb * 100",
          flags = "overwrite")

# open land cover neighborhood mean
execGRASS("r.neighbors",
          input = "open_cover_tv_study",
          output = "open_cover_3x3",
          method = "average",
          size = 3,
          flags = "overwrite")

execGRASS("r.mapcalc",
          expression = "open_cover_8nb = ((open_cover_3x3 * 9) - open_cover_tv_study) / 8.0",
          flags = "overwrite")

execGRASS("r.mapcalc",
          expression = "open_cover_8nb_pct = open_cover_8nb * 100",
          flags = "overwrite")

# impervious cover neighborhood mean
execGRASS("r.neighbors",
          input = "impervious_cover_tv_study",
          output = "impervious_cover_3x3",
          method = "average",
          size = 3,
          flags = "overwrite")

# remove focal cell (8-neighbors only)
execGRASS("r.mapcalc",
          expression = "impervious_cover_8nb = ((impervious_cover_3x3 * 9) - impervious_cover_tv_study) / 8.0",
          flags = "overwrite")

# rename to match other rasters
execGRASS("g.rename",
          raster = "impervious_cover_8nb,impervious_cover_8nb_pct",
          flags = "overwrite")

# read into session
forest_nb <- read_RAST("forest_cover_8nb_pct")
open_nb   <- read_RAST("open_cover_8nb_pct")
imperv_nb <- read_RAST("impervious_cover_8nb_pct")

# assign intuitive names
names(forest_nb) <- "neighbor_prct_forest"
names(open_nb)   <- "neighbor_prct_open"
names(imperv_nb) <- "neighbor_prct_impervious"

plot(forest_nb)
plot(open_nb)
plot(imperv_nb)

# write to raster folder
writeRaster(forest_nb, "raster/neighbor_prct_forest.tif", overwrite = TRUE)
writeRaster(open_nb, "raster/neighbor_prct_open.tif", overwrite = TRUE)
writeRaster(imperv_nb, "raster/neighbor_prct_imperv.tif", overwrite = TRUE)
