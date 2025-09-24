# Script to import avfallsstasjoner from Google Places API
# add to Ida's OSM vector
# Jenny Hansen
# 23 September 2025

# working in TidligVarsling project

Sys.setlocale("LC_CTYPE", "nb_NO.UTF-8")


# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(mapview)
library(googleway)
library(terra)
library(stringr)

# Set API key (from .Renviron)
set_key(Sys.getenv("GOOGLE_API_KEY"))

# Import data -------------------------------------------------------------

# raster for masking
template_rast <- rast("raster/ndvi_summer.tif")

# centroids of municipalities
ostlandet_centers <- read.csv("data/ostlandet_centers.csv")

# Ida's OSM data
osm_locs <- st_read("vector/recycling_pros6.geojson")


# Loop through kommuner ---------------------------------------------------

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


# Locate avfallsstasjoner -------------------------------------------------

# identified after poking around Google Maps
search_terms <- c(
  "avfallsstasjon",
  "gjenvinningsstasjon",
  "miljøstasjon",
  "returpunkt",
  "Recycling center",
  "Waste management service",
  "Garbage dump",
  "Garbage collection service"
)

results_list <- list()
api_key <- Sys.getenv("GOOGLE_API_KEY")
counter <- 1


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


# combine & remove duplicates

places_df_all <- do.call(rbind, results_list)

places_df_unique <- places_df_all %>%
  filter(!duplicated(place_id)) %>%
  group_by(address) %>%
  slice_max(order_by = place_id, n = 1, with_ties = FALSE) %>%
  ungroup()

places_sf <- st_as_sf(places_df_unique, coords = c("lng", "lat"), crs = 4326)

places_sf <- places_sf %>%
  mutate(
    across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8", sub = ""))
  )

mapview(places_sf)



# Remove bad matches ------------------------------------------------------

# based on the names in the places_sf df

keep_patterns <- c(
  "gjenvinn",       # gjenvinning, gjenvinningsstasjon
  "gjenbruks?",     # gjenbruk, gjenbruksstasjon
  "avfall",         # avfall, avfallsplass, avfallsdunk
  "s\u00F8ppel",    # søppel, søppelplass
  "milj\u00F8",     # miljø, miljøtorg, miljøstasjon
  "retur",          # returpunkt
  "omlast",         
  "renovasj",       
  "mottak",         
  "deponi",         
  "t\u00F8mmestasjon" # tømmestasjon
)


drop_patterns <- c(
  "Butikk","Comfort","Jernia","Elkjøp","Megaflis","Taxi",
  "Sushi","Massage","Frisør","Hotel","Extra","REMA","Kiwi",
  "Bygg","Service","Rør","AS","Shop","Mat","Restaurant"
)

places_filtered <- places_sf %>%
  mutate(name_lower = str_to_lower(name)) %>%
  filter(
    str_detect(name_lower, str_c(keep_patterns, collapse = "|")) &
      !str_detect(name, str_c(drop_patterns, collapse = "|"))
  ) %>%
  select(-name_lower) %>% 
  st_transform(25833)

mapview(places_filtered)


# Combine with Ida's data -------------------------------------------------

osm_pts <- st_centroid(osm_locs)

osm_pts <- osm_pts %>%
  mutate(source = "osm") %>%
  select(source, geometry)

places_pts <- places_filtered %>%
  mutate(source = "google_places") %>%
  select(source, geometry)

# combine
all_locs <- rbind(osm_pts, places_pts)


# remove duplicates
all_locs <- all_locs %>% 
  distinct(geometry, .keep_all = TRUE)


# Rasterize points to SSB grid --------------------------------------------

# use SSB template to map presence values (1) and generate a distance to
# nearest feature raster

presence_rast <- rasterize(vect(all_locs), template_rast, 
                           field = 1, background = NA)
plot(presence_rast, col = "darkgreen")

# distance-to-nearest-feature
dist_rast <- distance(presence_rast)
plot(dist_rast)

dist_masked <- mask(dist_rast, template_rast)
plot(dist_masked)

names(dist_masked) <- "distance_to_avfall"

# write to tif
writeRaster(dist_masked, "raster/distance_to_avfall.tif", overwrite = TRUE)
