# script to calculate forest edge landscape metrics in grid cells
# Jenny Hansen
# 23 September 2025

# working in TidligVarsling projects

# Load required libraries -------------------------------------------------

library(terra)
library(landscapemetrics)
library(rgrass)
library(sf)
library(future.apply)

# connect to GRASS mapset
NinaR::grassConnect()

# Setup GRASS environment -------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid", flags = "p")
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")

# Import data -------------------------------------------------------------

# was created in HELPER_rasterize_ar5.r
ar5 <- rast("raster/ar5_raster_3m.tif")

# vector of SSB grid shapefiles
grid_files <- list.files(
  "/data/R/GeoSpatialData/Population_demography/Norway_SSB/Original/ssb_250m",
  pattern = "\\.shp$",
  full.names = TRUE
)

# read and merge into a single sf object
grid_250 <- grid_files %>%
  lapply(st_read, quiet = TRUE) %>%
  do.call(rbind, .)

# assign CRS (it was missing)
st_crs(grid_250) <- crs(ar5)

# Prepare data ------------------------------------------------------------

# convert raster extent to sf polygon
ar5_extent <- st_as_sfc(st_bbox(ar5))

# clip
grid_250_clip <- grid_250 %>% 
  st_filter(ar5_extent, .predicate = st_intersects) %>% 
  st_intersection(ar5_extent)

# quick check
plot(grid_250_clip[,2])

# Calculate metrics -------------------------------------------------------

# NB: I tried running once, but quit after 19 hours of runtime; trying 
# batch/parallel processing instead

plan(multisession, workers = 10)   

batches <- split(grid_250_clip, (seq_len(nrow(grid_250_clip))-1) %/% 10000)

results <- future_lapply(
  batches,
  function(batch) {
    sample_lsm(
      landscape = ar5,
      y = batch,
      what = "lsm_c_te",
      classes = 30,
      unique_id = "SSBID"
    )
  }
)

forest_edge <- do.call(rbind, results)

# add results to grid_250

forest_edge <- forest_edge %>%
  dplyr::select(plot_id, value) %>%
  dplyr::rename(SSBID = plot_id, forest_edge_m = value)

grid_250_clip <- grid_250_clip %>%
  left_join(forest_edge, by = "SSBID")

# write to GRASS
write_VECT(grid_250_clip, "forest_edge_tv_study")


# Rasterize forest edge data ----------------------------------------------

execGRASS("v.to.rast",
          input = "forest_edge_tv_study",
          output = "forest_edge_250m",
          use = "attr",
          attribute_column = "forest_edge_m",
          type = "area",
          flags = "overwrite")



