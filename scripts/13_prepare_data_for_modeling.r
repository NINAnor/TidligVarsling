# Script to prepare data for modeling and check for collinearity
# Jenny Hansen
# 26 September 2025
# updated 08 Oct 2025 to include new survey data & get maxent bg pts
# updated 10 October 2025 to include national monitoring & artskart data
# update 04 Nov 2025 after receiving additional plant data

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(corrplot)
library(exactextractr)
library(spatstat.random)

# Import data -------------------------------------------------------------

# insect spp presence data, created in script 11
insect_data <- st_read("vector/insect_spp_presence_data.geojson")

# plant spp presence data, created in script 12
plant_data <- st_read("vector/plant_spp_presence_data.geojson") %>% 
  rename(name = lokname)

# predictor rasters, created in script 10
pred <- rast("raster/complete_prediction_stack.tif")

# bounding area for sampling (from Ida)
# study are based on 250 m rutenett from SSB (from Ida)
#aoi <- st_read("vector/ssb250_studyarea.shp") 

# ar5 for masking
#ar5 <- st_read("vector/ar5_pred_area.gpkg", layer = "ar5")


# Fix location mismatch ---------------------------------------------------

# after hours of troubleshooting, I realized that the points for insects
# and plants differ- assigning insect locations to plants here; they are
# more specific
plant_data <- plant_data %>%
  select(-geometry) %>%
  left_join(st_drop_geometry(insect_data[, c("name", "geometry")]), 
            by = "name") %>%
  st_as_sf()


# Mask fjord from aoi -----------------------------------------------------

# NB: commented out after initial run; very time consuming

# union aoi first
# aoi_union <- aoi %>%
#   st_union() %>%
#   st_as_sf() %>%
#   mutate(id = 1)
# 
# # use ar5 extent to mask non-land areas in aoi- I have masked them out
# # in the polygons
# ar5_union <- ar5 %>%
#   st_make_valid() %>%
#   st_union() %>%
#   st_as_sf() %>%
#   mutate(id = 1)
# 
# # remove area from aoi
# aoi_land <- st_intersection(aoi_union, ar5_union)
# 
# mapview::mapview(aoi_land)
# st_write(aoi_land, "vector/bg_sampling_area.geojson")

# Create exclusion zone ---------------------------------------------------

aoi_land <- st_read("vector/bg_sampling_area.geojson")

# this is to prevent randomly-sampled pseudoabsences from being
# drawn near presence locations
ist_buffer <- insect_data %>% st_buffer(500)
plt_buffer <- plant_data %>% st_buffer(500)

# union buffered areas
ist_exclusion_zone <- st_union(ist_buffer)
plt_exclusion_zone <- st_union(plt_buffer)

# remove exclusion area from aoi
ist_aoi_background <- st_difference(aoi_land, ist_exclusion_zone)
plt_aoi_background <- st_difference(aoi_land, plt_exclusion_zone)

# Generate background points ----------------------------------------------

set.seed(42)

# NB: these bg points are for GAM, BRT, and RF only
# using a 1:3 sample design (1:1 for plants)

# insect bg pts
insect_pa <- st_sample(ist_aoi_background, size = 531)
insect_pa <- st_as_sf(insect_pa)

# plants bg pts
plants_pa <- st_sample(plt_aoi_background, size = 3712)
plants_pa <- st_as_sf(plants_pa)


# Join presence & pseudoabsence data --------------------------------------

# create a source column in the original data
insect_obs <- insect_data %>%
  mutate(source = "presence")

# assign '0' as species richness in the PA data
insect_abs <- insect_pa %>%
  st_as_sf() %>%
  rename(geometry = x) %>% 
  mutate(species_richness = 0,
         source = "pseudoabsence",
         name = paste0("pseudo_", row_number()))

# combine into a single set
insect_full <- bind_rows(insect_obs, insect_abs)

# repeat for plants
plant_obs <- plant_data %>%
  mutate(source = "presence")

plant_abs <- plants_pa %>%
  st_as_sf() %>%
  rename(geometry = x) %>% 
  mutate(species_richness = 0,
         source = "pseudoabsence",
         name = paste0("pseudo_", row_number()))

plant_full <- bind_rows(plant_obs, plant_abs)


# Extract predictor values ------------------------------------------------

# NB: using exact_extract instead of terra because three of the
# sample grid points are right on the fjord edge and get NA
# values; this mitigates that problem

insect_full <- insect_full %>%
  st_cast("POINT")  

plant_full <- plant_full %>%
  st_cast("POINT")

# combine results back with attributes
insect_pred <- exact_extract(pred, st_buffer(insect_full, 250), 
                             fun = "mean")%>% 
  rename_with(~ sub("^mean\\.", "", .x)) %>% 
  bind_cols(st_drop_geometry(insect_full)) %>% 
  select(species_richness:source, 
         distance_to_public_road:percentage_agriculture)

plant_pred <- exact_extract(pred, st_buffer(plant_full, 250), 
                            fun = "mean") %>% 
  rename_with(~ sub("^mean\\.", "", .x)) %>% 
  bind_cols(st_drop_geometry(plant_full)) %>% 
  select(species_richness:source, 
         distance_to_public_road:percentage_agriculture)


# Check for completeness --------------------------------------------------

# isolate only predictor vars for insects
insect_vars <- insect_pred %>%
  select(where(is.numeric)) %>%    
  select(-species_richness)        
sum(is.na(insect_vars)) # 0

# same for plants
plant_vars <- plant_pred %>%
  select(where(is.numeric)) %>%
  select(-species_richness)
sum(is.na(plant_vars)) # 128

# NB: 128, too many to replace as before
# after visualizing, it is clear they are all on the mask edge

# identify complete rows
complete_idx <- complete.cases(plant_vars)

# keep only complete
plant_vars  <- plant_pred[complete_idx, ]
plant_full  <- plant_full[complete_idx, ]


# Create background data for MaxEnt (presence–background) -----------------

# using spatstat random to ensure pts are >= 1000 m apart, as in Ida's method

# after many NAs, I see the extent of aoi is slightly larger that the 
# prediction raster; fixing below

pred_complete <- app(pred, fun = function(x) all(!is.na(x)))

valid_area_vect <- as.polygons(pred_complete) %>%
  st_as_sf() %>%
  st_make_valid()

aoi_valid <- st_intersection(aoi_land, valid_area_vect)

aoi_fixed <- aoi_valid %>%
  st_collection_extract("POLYGON") %>%
  st_union() %>%
  st_cast("MULTIPOLYGON") %>%
  st_sf()

aoi_spat_valid <- as.owin(st_geometry(aoi_fixed))

# sample with min dist = 1000 m
bg_insects_ppp <- rSSI(r = 1000, n = 10000, win = aoi_spat_valid)
bg_plants_ppp  <- rSSI(r = 1000, n = 10000, win = aoi_spat_valid)

# convert  to sf
bg_insects_sf <- st_as_sf(as.data.frame(bg_insects_ppp), 
                          coords = c("x", "y"), 
                  crs = st_crs(aoi_land))
bg_plant_sf <- st_as_sf(as.data.frame(bg_plants_ppp), 
                        coords = c("x", "y"), 
                          crs = st_crs(aoi_land))


# add presence to existing insect and plant occurrence data
insect_pres <- insect_data %>%
  mutate(present = 1,
         source = "presence")

plant_pres <- plant_data %>%
  mutate(present = 1,
         source = "presence")

# add metadata to background points
bg_insects_sf <- bg_insects_sf %>%
  mutate(present = 0,
         source = "background",
         name = paste0("bg_insect_", row_number()))

bg_plants_sf <- bg_plant_sf %>%
  mutate(present = 0,
         source = "background",
         name = paste0("bg_plant_", row_number()))

# combine presences and background for each group
insect_pb <- bind_rows(insect_pres, bg_insects_sf)
plant_pb  <- bind_rows(plant_pres,  bg_plants_sf)

# cast to point
insect_pb <- insect_pb %>%
  st_cast("POINT") %>%
  st_sf()  

plant_pb <- plant_pb %>%
  st_cast("POINT") %>%
  st_sf()


# handle annoying edge issues
pred_masked <- mask(pred, pred_complete)

# extract
insect_vals <- extract(pred_masked, vect(insect_pb))
plant_vals <- extract(pred_masked, vect(plant_pb))

# join
insect_pb_pred <- insect_pb %>%
  dplyr::mutate(presence = ifelse(source == "presence", 1, 0)) %>%
  dplyr::select(presence) %>%                   
  dplyr::bind_cols(sf::st_drop_geometry(insect_vals[, -1])) 

plant_pb_pred <- plant_pb %>%
  dplyr::mutate(presence = ifelse(source == "presence", 1, 0)) %>%
  dplyr::select(presence) %>%
  dplyr::bind_cols(sf::st_drop_geometry(plant_vals[, -1]))


# check completeness
sum(is.na(insect_pb_pred)) # missing 39 values
sum(is.na(plant_pb_pred)) # missing 1726 values

# removing missing values; these are edge cases along fjord mask
complete_idx <- complete.cases(st_drop_geometry(plant_pb_pred))
plant_pb_pred <- plant_pb_pred[complete_idx, ]


# Write df to file --------------------------------------------------------

# join back response variables
insect_df <- cbind(insect_full %>% select(species_richness),
                   insect_vars) %>% st_drop_geometry()

plant_df <- cbind(plant_full %>% select(species_richness),
                  plant_vars) %>% st_drop_geometry()

write.csv(insect_df, "data/insect_model_df.csv", row.names = FALSE)
write.csv(plant_df, "data/plant_model_df.csv", row.names = FALSE)

write.csv(st_drop_geometry(insect_pb_pred), "data/insect_maxent_data.csv", 
          row.names = FALSE)
write.csv(st_drop_geometry(plant_pb_pred),  "data/plant_maxent_data.csv",  
          row.names = FALSE)


# Write sf to vector ------------------------------------------------------

insect_sf <- cbind(insect_full %>% select(species_richness),
                   insect_vars) 

plant_sf <- cbind(plant_full %>% select(species_richness),
                  plant_vars)

st_write(insect_sf, "vector/insect_model_sf.geojson")
st_write(plant_sf, "vector/plant_model_sf.geojson")

st_write(insect_pb_pred, "vector/insect_maxent_data.geojson")
st_write(plant_pb_pred,  "vector/plant_maxent_data.geojson")
