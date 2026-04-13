# Rute utvalg for TidligVarsling project                    
# Ida M. Mienna and Jenny Hansen
# 24 March 2026

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(mapview)
library(kableExtra)

# Import data -------------------------------------------------------------

# from script 02
rutenett_250 <- st_read("vector/ssb250_studyarea.shp") 

# prediction area grid w/SSBID 250 m rutenett (script 02)
ssb_grid <- rast("raster/sa_grid_250m_fixed.tif")

# Boosted Regression Trees selected as highest performing for plants (script 15)
plant_pred <- rast("raster/plant_prob_brt.tif")

# RandomForest selected as highest performing for insects (script 15)
insect_pred <- rast("raster/insect_prob_rf_expanded.tif")

# sampling regions (drawn in QGIS)
samp_regions <- st_read("vector/sampling_regions.gpkg", 
                        layer = "sampling_regions")


# Create combined probability dataset -------------------------------------

combined_stack <- c(plant_pred, insect_pred, ssb_grid)
names(combined_stack) <- c("plant_pred", "insect_pred", "SSBID")

combined_df <- as.data.frame(combined_stack, na.rm = FALSE)

combined_sf <- rutenett_250 %>%
  select(SSBID, geometry) %>%
  left_join(combined_df, by = "SSBID") %>%
  st_as_sf(crs = 25833)


# Sample n = 50 each region -----------------------------------------------

# join top grid with sampling_regions, select 50 highest ranking
# ruter in each region
top_grid_regions <- st_join(combined_sf, samp_regions, left = FALSE)

# select top 50 per region based on insect probability
rute_utvalg_regions <- top_grid_regions %>%
  filter(!is.na(insect_pred)) %>%
  arrange(desc(insect_pred)) %>%
  group_by(regions) %>%
  slice_head(n = 50) %>%
  ungroup()
  
mapview(rute_utvalg_regions, zcol = "insect_pred") 

# export
st_write(rute_utvalg_regions, "vector/utvalgte_ruter_nytt.geojson")


# Selection of low-quality grids ------------------------------------------

low_prob_grids <- top_grid_regions %>% 
  filter(!is.na(insect_pred),
         insect_pred <= 0.1) %>% 
  arrange(desc(insect_pred)) %>%
  group_by(regions) %>%
  slice_head(n = 50) %>%
  ungroup()

mapview(low_prob_grids, zcol = "insect_pred")

# export
st_write(low_prob_grids, "vector/utvalgte_ruter_lav_sannsynlighet.geojson")

# Summary tables ----------------------------------------------------------

region_summary <- rute_utvalg_regions %>%
  st_drop_geometry() %>%
  group_by(regions) %>%
  summarise(
    n_ruter = n(),
    
    insect_mean = mean(insect_pred, na.rm = TRUE),
    insect_median = median(insect_pred, na.rm = TRUE),
    insect_min = min(insect_pred, na.rm = TRUE),
    insect_prop_high = mean(insect_pred > 0.7, na.rm = TRUE),
    
    plant_mean = mean(plant_pred, na.rm = TRUE),
    plant_median = median(plant_pred, na.rm = TRUE),
    plant_min = min(plant_pred, na.rm = TRUE),
    plant_prop_high = mean(plant_pred > 0.7, na.rm = TRUE)
  ) %>%
  arrange(desc(insect_mean))

region_summary %>%
  kable(digits = 3, caption = "Sampling summary per region") %>%
  kable_styling(full_width = FALSE)

low_plant_summary <- rute_utvalg_regions %>%
  st_drop_geometry() %>%
  group_by(regions) %>%
  summarise(
    n_low_plant = sum(plant_pred <= 0.7, na.rm = TRUE),
    prop_low_plant = mean(plant_pred <= 0.7, na.rm = TRUE)
  )

low_plant_summary %>% 
  kable(digits = 3, caption = "Low plant summary per region") %>%
  kable_styling(full_width = FALSE)
