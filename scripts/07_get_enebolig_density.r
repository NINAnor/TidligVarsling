# Script to extract 'enebolig' from FKB-bygning and calculate density of 
# enebloig in each 250 m grid cell in the study area
# Jenny Hansen
# 08 August 2025

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(mapview)
library(duckdb)
library(glue)

library(rgrass)

# connect to GRASS mapset
NinaR::grassConnect()

# Setup GRASS environment -------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid", flags = "p")
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")


# Import data -------------------------------------------------------------

# 250 m grid for projection
template_rast <- rast("raster/ndvi_summer.tif")

# bbox of the study area
study_bbox <- st_as_sfc(st_bbox(template_rast))

# Set up db ---------------------------------------------------------------

con <- dbConnect(duckdb())
dbExecute(con, "INSTALL spatial FROM core_nightly; LOAD spatial;")

gdb_path <- "vector/FKB_bygning.gdb"

# Register the FKB layer
dbExecute(con, glue("
  CREATE TABLE bygning AS
  SELECT * FROM ST_Read('{gdb_path}', layer='fkb_bygning_omrade')
"))

# Define bbox from your template raster
bbox <- st_bbox(template_rast)

# Query only enebolig types (111–113) inside bbox
query <- glue("
  SELECT
    *,
    ST_AsWKB(ST_Force2D(SHAPE)) AS geom_wkb
  FROM bygning
  WHERE bygningstype IN (111,112,113)
    AND ST_Intersects(
          SHAPE,
          ST_MakeEnvelope({bbox['xmin']},{bbox['ymin']},{bbox['xmax']},{bbox['ymax']})
        )
")

res <- dbSendQuery(con, query)
df <- dbFetch(res)
dbClearResult(res)

# Convert to sf
geom_sfc <- st_as_sfc(df$geom_wkb, EWKB = TRUE, crs = 25833)

enebolig <- df %>%
  select(bygningstype) %>%
  mutate(geometry = geom_sfc) %>%
  st_sf()


# Create enebolig density raster ------------------------------------------

# create centroids
enebolig_pts <- st_point_on_surface(enebolig) 

# make SpatVector for terra
enebolig_vect <- vect(enebolig_pts)

# Rasterize count per cell
house_density <- rasterize(
  enebolig_vect,
  template_rast,
  fun = "count",     # counts number of points per cell
  background = 0     # cells with no houses = 0
)

# make pretty
density_rast <- mask(house_density, template_rast)
plot(density_rast)

writeRaster(density_rast, "raster/enebolig_density.tif",
            overwrite = TRUE)


# Compare with matrikkelen ------------------------------------------------


gdb_path <- "vector/matrikkelen.gdb"

# Register the FKB layer
dbExecute(con, glue("
  CREATE TABLE bygning AS
  SELECT * FROM ST_Read('{gdb_path}', layer='bygning')
"))

# Define bbox from your template raster
bbox <- st_bbox(template_rast)

# Query only enebolig types (111–113) inside bbox
query <- glue("
  SELECT
    bygningstype,
    ST_AsWKB(ST_Force2D(SHAPE)) AS geom_wkb
  FROM bygning
  WHERE bygningstype IN (111,112,113)
    AND ST_Intersects(
          SHAPE,
          ST_MakeEnvelope({bbox['xmin']},{bbox['ymin']},{bbox['xmax']},{bbox['ymax']})
        )
")
res <- dbSendQuery(con, query)
df <- dbFetch(res)
dbClearResult(res)

# Convert to sf
geom_sfc <- st_as_sfc(df$geom_wkb, EWKB = TRUE, crs = 25833)

enebolig <- df %>%
  select(-geom_wkb) %>%
  mutate(geometry = geom_sfc) %>%
  st_sf()


# Create enebolig density raster ------------------------------------------

# create centroids
enebolig_pts <- st_point_on_surface(enebolig) 

# make SpatVector for terra
enebolig_vect <- vect(enebolig_pts)

# Rasterize count per cell
house_density <- rasterize(
  enebolig_vect,
  template_rast,
  fun = "count",     # counts number of points per cell
  background = 0     # cells with no houses = 0
)

# make pretty
density_rast_alt <- mask(house_density, template_rast)
plot(density_rast_alt)


# correlatoin between two datasets
all.equal(ext(density_rast), ext(density_rast_alt))
all.equal(res(density_rast), res(density_rast_alt))
all.equal(crs(density_rast), crs(density_rast_alt))


vals <- cbind(values(density_rast), values(density_rast_alt))
colnames(vals) <- c("rast1", "rast2")

# Drop NA rows
vals <- na.omit(vals)

# Pearson correlation
correlation <- cor(vals[, "rast1"], vals[, "rast2"], method = "pearson")
correlation # 0.999

# these are the same data


