# Script to get coldest temperatures from the coldest months
# during the survey period
# Jenny Hansen
# 18 September 2025

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

# seNorge data
years <- c(2019:2024)
dir <- "/data/R/GeoSpatialData/Meteorology/Norway_SeNorge2018_v22.09/Original/rr_tm"
pattern <- "seNorge2018_"
file_ext <- ".nc"


# build paths & read in only minimum temp layers 
files <- file.path(dir, paste0(pattern, years, file_ext))
rasters_tn <- lapply(files, function(f) {
  vnames <- terra::sds(f)
  vnames[["tn"]]
})


# Get coldest temps in the winter months ----------------------------------

# function to extract coldest daily minimum 
get_winter_min <- function(r) {

  dts <- terra::time(r)
  months <- as.numeric(format(dts, "%m"))
  
  # select only Dec:Feb
  winter <- r[[which(months %in% c(12, 1, 2))]]
  
  # pixel-wise minimum (coldest winter night)
  winter_min <- terra::app(winter, min, na.rm = TRUE)
  return(winter_min)
}

# apply to all years
winter_min_list <- lapply(rasters_tn, get_winter_min)
plot(winter_min_list[[1]])

# Get average across years ------------------------------------------------

winter_min_stack <- terra::rast(winter_min_list)
winter_min_avg   <- terra::app(winter_min_stack, mean, na.rm = TRUE)
plot(winter_min_avg)

winter_min <- project(winter_min_avg, template_rast, method = "bilinear")

write_RAST(winter_min, "coldest_winter_temp_tv_study",
           flags = "overwrite")


# Read in with GRASS ------------------------------------------------------

# NB: original resolution of seNorge data is 1000 x 1000
# here it is 250 x 250

winter_min <- read_RAST("coldest_winter_temp_tv_study")
plot(winter_min)
names(winter_min) <- "coldest_winter_temperature"

# write to raster folder
writeRaster(winter_min, "raster/coldest_winter_temp.tif",
            overwrite = TRUE)

