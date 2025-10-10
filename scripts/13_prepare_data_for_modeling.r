# Script to prepare data for modeling and check for collinearity
# Jenny Hansen
# 26 September 2025
# updated 08 Oct 2025 to include new survey data & get maxent bg pts

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
insect_data <- st_read("vector/insect_spp_presence_data_tv_only.geojson")

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
buffer <- insect_data %>% st_buffer(500)

# union buffered areas
exclusion_zone <- st_union(buffer)

# remove exclusion area from aoi
aoi_background <- st_difference(aoi_land, exclusion_zone)

# Generate background points ----------------------------------------------

set.seed(42)

# NB: these bg points are for GAM and BRT only
# I am going with an initial 1:3 sample design (higher for plants)

# insect bg pts
insect_pa <- st_sample(aoi_background, size = 225)
insect_pa <- st_as_sf(insect_pa)

# plants bg pts
plants_pa <- st_sample(aoi_background, size = 225)
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
  st_cast("POINT")  # keeps attributes by default

plant_full <- plant_full %>%
  st_cast("POINT")

# combine results back with attributes
insect_pred <- exact_extract(pred, st_buffer(insect_full, 250), fun = "mean")%>% 
  rename_with(~ sub("^mean\\.", "", .x)) %>% 
  bind_cols(st_drop_geometry(insect_full)) %>% 
  select(name:source, distance_to_public_road:percentage_agriculture)

plant_pred <- exact_extract(pred, st_buffer(plant_full, 250), fun = "mean") %>% 
  rename_with(~ sub("^mean\\.", "", .x)) %>% 
  bind_cols(st_drop_geometry(plant_full)) %>% 
  select(name:source, distance_to_public_road:percentage_agriculture)


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
sum(is.na(plant_vars)) # 1

# NB: there is 1 pseudoabsence NA- resmapling and replacing below

# Replace NA sample -------------------------------------------------------

# one plant pseudoabsence has an NA; need to resample/replace it

# identify na
bad_idx <- which(!complete.cases(plant_vars))

# remove it
bad_point <- plant_full[bad_idx, ]   # 198 pseudo_154 
plant_full <- plant_full[-bad_idx, ]
plant_vars <- plant_vars[-bad_idx, ]

# sample new pseudoabsence
new_point <- st_sample(aoi_background, size = 1) %>%
  st_as_sf() %>%
  st_set_geometry("x") %>%             
  rename(geometry = x) %>%  # rename to 'geometry'           
  mutate(source = "pseudoabsence",
         species_richness = 0,
         name = "pseudo_resample") %>%
  select(name, species_richness, source, geometry)

# extract pred vars for new point
new_pred <- exact_extract(pred, st_buffer(new_point, 250), 
                          fun = "mean") %>%
  rename_with(~ sub("^mean\\.", "", .x)) %>%
  bind_cols(st_drop_geometry(new_point)) %>% 
  select(-c(name:source))


# add back to complete dataset
plant_full <- rbind(plant_full, new_point)
plant_vars <- bind_rows(plant_vars, new_pred)

sum(is.na(plant_vars)) # 0    


# Create background data for MaxEnt (presence–background) -----------------

# using spatstat random to ensure pts are >= 1000 m apart, as in Ida's method

# after many NAs, I see the extent of aoi is slightly larger that the 
# prediction raster; clipping below
aoi_clip <- st_as_sfc(st_bbox(pred)) %>%
  st_sf(crs = st_crs(pred))

# intersect aoi with extent of pred rasters
aoi_clipped <- st_intersection(aoi_land, aoi_clip)

# fix geometry to work with spatstat
aoi_fixed <- aoi_clipped %>% 
  st_make_valid() %>% 
  st_collection_extract("POLYGON") %>% 
  st_union() %>%                 
  st_cast("MULTIPOLYGON") %>%    
  st_sf()

# convert to a spatstat window (necessary for enforcing dist requirement)
aoi_spat <- as.owin(st_geometry(aoi_fixed))

# sample with min dist = 1000 m
bg_insects_ppp <- rSSI(r = 1000, n = 10000, win = aoi_spat)
bg_plants_ppp  <- rSSI(r = 1000, n = 10000, win = aoi_spat)

# convert  to sf
bg_insects_sf <- st_as_sf(as.data.frame(bg_insects_ppp), coords = c("x", "y"), 
                  crs = st_crs(aoi_land))
bg_plant_sf <- st_as_sf(as.data.frame(bg_plants_ppp), coords = c("x", "y"), 
                          crs = st_crs(aoi_land))


# add presence to existing insect and plant occurrence data
insect_pres <- insect_data %>%
  mutate(present = 1,
         source = "presence")

plant_pres <- plant_data %>%
  mutate(present = 1,
         source = "presence")

# Add metadata to background points
bg_insects_sf <- bg_insects_sf %>%
  mutate(present = 0,
         source = "background",
         name = paste0("bg_insect_", row_number()))

bg_plants_sf <- bg_plant_sf %>%
  mutate(present = 0,
         source = "background",
         name = paste0("bg_plant_", row_number()))

# Combine presences and background for each group
insect_pb <- bind_rows(insect_pres, bg_insects_sf)
plant_pb  <- bind_rows(plant_pres,  bg_plants_sf)

# cast to point
insect_pb <- insect_pb %>%
  st_cast("POINT") %>%
  st_sf()  # ensure geometry column is pure POINT

plant_pb <- plant_pb %>%
  st_cast("POINT") %>%
  st_sf()



# handle annoying edge issues
pred_complete <- app(pred, fun = function(x) all(!is.na(x)))
pred_masked <- mask(pred, pred_complete)

# extract
insect_vals <- extract(pred_masked, vect(insect_pb))
plant_vals <- extract(pred_masked, vect(plant_pb))

# join
insect_pb_pred <- bind_cols(insect_pb, insect_vals[,-1])
insect_pb_pred <- insect_pb_pred %>%
  mutate(presence = ifelse(source == "presence", 1, 0))

plant_pb_pred <- bind_cols(plant_pb, plant_vals[,-1])
plant_pb_pred <- plant_pb_pred %>%
  mutate(presence = ifelse(source == "presence", 1, 0))


# check completeness
sum(is.na(insect_pb_pred)) # missing 21 values
sum(is.na(plant_pb_pred)) # missing 13 values


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
