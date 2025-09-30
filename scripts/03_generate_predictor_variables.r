# script to generate predictor variables for ensemble models
# Jenny Hansen
# 15 Sept. 2025

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(mapview)
library(terra)
library(rgrass)
library(osmdata)
library(googleway)
library(dbscan)

# connect to GRASS mapset
NinaR::grassConnect()

# Set API key (from .Renviron)
set_key(Sys.getenv("GOOGLE_API_KEY"))

# Setup GRASS environment -------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid", flags = "p")
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")


# Import data -------------------------------------------------------------

pot_nye <- st_read("vector/pot_nye_fremmede.geojson")

# blank raster over SA (raster created in 02_rasterize_ssb_grid_over_sa.r)
template_rast <- rast("raster/sa_grid_250m_fixed.tif")

# centroids of municipalities
ostlandet_centers <- read.csv("data/ostlandet_centers.csv")


# Import rasters from GRASS -----------------------------------------------

# read in rasters from SAM project folders

# distance to privat veier
map_names <- execGRASS("g.list", type = "raster", 
                       mapset = "p_sam_transport_urban,p_sam_industry,p_sam_tourism",
                       pattern = "*euclidean*",
                       flags = "m",
                       intern = TRUE)

# vector of relevant names
names_vec <- c("roads_public_nose_sum_zoi_nearest_euclidean@p_sam_transport_urban",
               "roads_private_nose_sum_zoi_nearest_euclidean@p_sam_transport_urban",
               "railways_nose_sum_zoi_nearest_euclidean@p_sam_transport_urban",
               "houses_nose_zoi_nearest_euclidean@p_sam_transport_urban" )

rast_stack <- read_RAST(names_vec)


# Calculate dist to industrial --------------------------------------------

# log into PostGIS database
NinaR::postgreSQLConnect(
  host = "gisdata-db.nina.no",
  dbname = "gisdata",
  username = "postgjest",
  password = "gjestpost"
)

#---
# DBI::dbSendQuery(con, "DROP TABLE sam_tmp.tmp_industrial_buildings_N50_no;")
# DBI::dbSendQuery(con, "DROP TABLE sam_env.industrial_buildings_no;")
# industrial buildings in Norway for energy supply in Norway, from N50
map <- "n50_2013_utm33n.n50_2013_byggoganlegg_points"
conditions <- c("(objtype) = 'Bygning'",
                "(byggtyp_nbr) IN ('211','212','214','216','219','221','223','229')")
# "(byggtyp_nbr) IN ('211', '212', '214', '216', '219', '221', '223', '229')") 
# 211: Factory building, building for industrial series production.
# 212: Workshop building, Building for special production or repair
# 214: Building for treatment plant, i.a. sewage pumping station
# 216: Building for water supply, i.a. pumping station.
# 219: Other industrial building, or building that is closely connected to / serves such and similar building (s)
# 221: Power station
# 223: Transformer drive
# 229: Other energy supply building, or building that is closely connected to / serves such and similar building (s)
query1 <- paste0(
  "SELECT * FROM ",
  map,
  " WHERE ", 
  paste(conditions, collapse = " AND "))
query1

indust <- st_read(con, query = query1)

# write vector to grass
write_VECT(vect(indust), "industrial_buildings", flags = c("overwrite"))
execGRASS("v.clean", input = "industrial_buildings", 
          output = "industrial_buildings_clean", 
          tool = "break", flags = c("overwrite"))

# rasterize
execGRASS("v.to.rast", input = "industrial_buildings_clean",
          output = "industrial_buildings_rast",
          use = "cat", type = "point", flags = c("overwrite"))

# create distance-to raster
execGRASS("r.grow.distance",
          input    = "industrial_buildings_rast",
          distance = "industrial_buildings_rast_zoi_nearest_euclidean",
          metric   = "euclidean",
          flags    = c("overwrite","quiet"))

indust_dist <- read_RAST("industrial_buildings_rast_zoi_nearest_euclidean")

# add to raster stack
rast_stack <- c(rast_stack, indust_dist)

# Create distance to garden centers ---------------------------------------

# import garden centers from OSM
norway_bbox <- c(4.93959, 58.07888, 31.29342, 71.18548) 

query <- opq(bbox = norway_bbox) %>%
  add_osm_feature(key = "shop", value = "garden_centre")

osm_data <- osmdata_sf(query)

# some garden centers are represented with both points
# and polygons; select only unique points
hs_pts <- osm_data$osm_points


# check with Google Places API
# function to loop over kommuner
get_places <- function(lat, lng, kommune, api_key) {
  res <- google_places(
    location = c(lat, lng),
    radius = 50000,  
    keyword = "building materials",
    key = api_key
  )
  
  if (is.null(res$results)) return(NULL)
  
  df <- data.frame(
    kommune = kommune,
    name = res$results$name,
    address = res$results$vicinity,
    lat = res$results$geometry$location$lat,
    lng = res$results$geometry$location$lng,
    place_id = res$results$place_id,
    stringsAsFactors = FALSE
  )
  
  df
}

search_terms <- c("Garden center", "Garden centre")
results_list <- list()
counter <- 1
api_key <- Sys.getenv("GOOGLE_API_KEY")

for (i in seq_len(nrow(ostlandet_centers))) {
  for (term in search_terms) {
    
    term_human <- term
    term_api <- URLencode(enc2utf8(term), reserved = TRUE)
    
    cat("Searching:", term_human, "in", ostlandet_centers$kommune[i], "\n")
    
    res <- google_places(
      location = c(ostlandet_centers$lat[i], ostlandet_centers$lng[i]),
      radius   = 50000,
      keyword  = term_api,
      key      = api_key
    )
    
    if (length(res$results) > 0) {
      df <- data.frame(
        kommune      = ostlandet_centers$kommune[i],
        search_term  = term_human,
        search_term_encoded = term_api,
        name         = res$results$name,
        address      = res$results$vicinity,
        lat          = res$results$geometry$location$lat,
        lng          = res$results$geometry$location$lng,
        place_id     = res$results$place_id,
        stringsAsFactors = FALSE
      )
      results_list[[counter]] <- df
      counter <- counter + 1
    }
    
    Sys.sleep(1)
  }
}

places_df_all <- do.call(rbind, results_list)

places_df_unique <- places_df_all[!duplicated(places_df_all$place_id), ] %>% 
  group_by(address) %>%
  slice_max(order_by = place_id, n = 1, with_ties = FALSE) %>%
  ungroup() # n = 249

places_sf <- st_as_sf(places_df_unique, coords = c("lng", "lat"), crs = 4326)
places_sf <- places_sf %>% mutate(source = "google") %>% select(source)

# add to OSM data
hs_pts <- hs_pts %>% mutate(source = "osm") %>% select(source)

final_garden_centers <- rbind(hs_pts, places_sf) %>% 
  st_transform(25833)

# remove duplicates
coords <- st_coordinates(final_garden_centers)

# cluster within 100 m
cl <- dbscan(coords, eps = 250, minPts = 1)

final_garden_centers$cluster <- cl$cluster

# keep first in each cluster
final_garden_centers_nodup <- final_garden_centers %>%
  group_by(cluster) %>%
  slice(1) %>%
  ungroup() %>%
  select(-cluster)

# make vect
garden_centers_vect <- final_garden_centers_nodup %>% 
  st_cast("POINT") %>% 
  vect()

write_VECT(garden_centers_vect, "garden_centers",
           flags = "overwrite")

# rasterize the vector
execGRASS("v.to.rast", input = "garden_centers",
          output = "garden_centers_rast",
          use = "cat", type = "point", flags = c("overwrite"))

# compute distance 
execGRASS("r.grow.distance",
          input    = "garden_centers_rast",
          distance = "garden_centers_rast_zoi_nearest_euclidean",
          metric   = "euclidean",
          flags    = c("overwrite","quiet"))
gc_dist <- read_RAST("garden_centers_rast_zoi_nearest_euclidean")

# add to raster stack
rast_stack <- c(rast_stack, gc_dist)


# Distance to lumberyards -------------------------------------------------

# use Google Places API
search_terms <- c("lumber", "trelast", "byggevarer", "byggvarehus", "sagbruk")
results_list <- list()
counter <- 1
api_key <- Sys.getenv("GOOGLE_API_KEY")

for (i in seq_len(nrow(ostlandet_centers))) {
  for (term in search_terms) {
    
    term_human <- term
    term_api <- URLencode(enc2utf8(term), reserved = TRUE)
    
    cat("Searching:", term_human, "in", ostlandet_centers$kommune[i], "\n")
    
    res <- google_places(
      location = c(ostlandet_centers$lat[i], ostlandet_centers$lng[i]),
      radius   = 50000,
      keyword  = term_api,
      key      = api_key
    )
    
    if (length(res$results) > 0) {
      df <- data.frame(
        kommune      = ostlandet_centers$kommune[i],
        search_term  = term_human,
        search_term_encoded = term_api,
        name         = res$results$name,
        address      = res$results$vicinity,
        lat          = res$results$geometry$location$lat,
        lng          = res$results$geometry$location$lng,
        place_id     = res$results$place_id,
        stringsAsFactors = FALSE
      )
      results_list[[counter]] <- df
      counter <- counter + 1
    }
    
    Sys.sleep(1)
  }
}

places_df_all <- do.call(rbind, results_list)

places_df_unique <- places_df_all[!duplicated(places_df_all$place_id), ] %>% 
  group_by(address) %>%
  slice_max(order_by = place_id, n = 1, with_ties = FALSE) %>%
  ungroup() # n =  551

places_sf <- st_as_sf(places_df_unique, coords = c("lng", "lat"), crs = 4326)

# remove locations after 'ground truthing' on street view & sat images
# remove Omsland sag
# remove Bark-Trefelling
# remove Helge Støvern Hiåskogen Sagbruk (couldn't verify)
# remove Sørum sag (couldn't verify)
# remove Krogsrud Sag AS 
# remove Hustveit sag (could not verify)
# remove Solerød sag Hans Christian Bjerkest (could not verify)
# remove Stokstad sag og høvel (could not verify)

remove <- c("Omslag sag", "Bark-Trefelling", "Helge Støvern Hiåskogen Sagbruk",
            "Sørum sag", "Krogsrud Sag AS", "Hustveit sag",
            "Solerød sag Hans Christian Bjerkest", "Stokstad sag og høvel")
clean_places <- places_sf %>% filter(!name %in% remove) %>% 
  st_transform(crs(template_rast))

write_VECT(vect(clean_places), "lumberyards", flags = "overwrite")

# rasterize the vector
execGRASS("v.to.rast", input = "lumberyards",
          output = "lumberyards_rast",
          use = "cat", type = "point", flags = c("overwrite"))

# compute distance 
execGRASS("r.grow.distance",
          input    = "lumberyards_rast",
          distance = "lumberyards_rast_zoi_nearest_euclidean",
          metric   = "euclidean",
          flags    = c("overwrite","quiet"))
lumber_dist <- read_RAST("lumberyards_rast_zoi_nearest_euclidean")

# add to raster stack
rast_stack <- c(rast_stack, lumber_dist)

# rename rasters
names(rast_stack) <- c("distance_to_public_road", "distance_to_private_road",
                       "distance_to_railway", "distance_to_housing",
                       "distance_to_industrial_area", 
                       "distance_to_garden_center",
                       "distance_to_lumberyard")

# write to raster folder
writeRaster(rast_stack, "raster/predictor_stack.tif", overwrite = TRUE)
