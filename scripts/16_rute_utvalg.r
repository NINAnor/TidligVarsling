# Rute utvalg for TidligVarsling project                    
# Written by Ida M. Mienna
# Created 07.11.2025
# gently modified by Jenny Hansen 18 Feb 2026

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(mapview)
library(ggplot2)
library(purrr)

# Import data -------------------------------------------------------------

# from script 02
rutenett_250 <- st_read("vector/ssb250_studyarea.shp") 

# Boosted Regression Trees selected as highest performing for plants (script 15)
plant_pred <- rast("raster/plant_prob_brt.tif")

# RandomForest selected as highest performing for insects (script 15)
insect_pred <- rast("raster/insect_prob_rf_expanded.tif")

# prediction area grid w/SSBID 250 m rutenett (script 02)
ssb_grid <- rast("raster/sa_grid_250m_fixed.tif")

# combine to stack & assign names
plant_stack <- c(plant_pred, ssb_grid)
names(plant_stack) <- c("plant_pred", "ssbid")
insect_stack <- c(insect_pred, ssb_grid)
names(insect_stack) <- c("insect_pred", "ssbid")

# Rank by predicted values ------------------------------------------------

### Plants
plant_ranked <- plant_stack$plant_pred
values(plant_ranked) <- rank(values(plant_stack$plant_pred), na.last = "keep", 
                             ties.method = "average")

# get raster values
plant_vals <- values(plant_stack$plant_pred)

# identify top 500 cells (by predicted value)
plant_top_idx <- order(plant_vals, decreasing = TRUE)[1:500]

# get XY coordinates for those cells
plant_xy <- xyFromCell(plant_ranked, plant_top_idx)

# convert to SpatVector points
top_points_plants <- vect(plant_xy, crs = crs(plant_ranked))  

# add predicted values
top_points_plants$rank <- 1:500
top_points_plants$predicted_value <- plant_vals[plant_top_idx]
plot(top_points_plants)

# add SSBID
top_points_plants$SSBID <- extract(plant_stack$ssbid, top_points_plants)$ssbid
top_points_plants

# join with SSB rutenett
ssb_ru_plants <-rutenett_250 %>% 
  left_join(st_drop_geometry(st_as_sf(top_points_plants)),
    by = "SSBID")
mapview(ssb_ru_plants %>% filter(!is.na(rank)), zcol = "predicted_value")

# export
top_points_plants %>% 
  writeVector("vector/plant_ruteutvalg_points.geojson")
ssb_ru_plants %>% filter(!is.na(rank)) %>% 
  select(SSBID, rank, predicted_value) %>% 
  st_write("vector/plant_ruteutvalg_rutenett.geojson")


### Insects
insect_ranked <- insect_stack$insect_pred
values(insect_ranked) <- rank(values(insect_stack$insect_pred), na.last = "keep", 
                              ties.method = "average")

# get raster values
insect_vals <- values(insect_stack$insect_pred)

# identify top 500 cells (by predicted value)
insect_top_idx <- order(insect_vals, decreasing = TRUE)[1:500]

# get XY coordinates for those cells
insect_xy <- xyFromCell(insect_ranked, insect_top_idx)

# convert to SpatVector points
top_points_insects <- vect(insect_xy, crs = crs(insect_ranked))  

# add predicted values
top_points_insects$rank <- 1:500
top_points_insects$predicted_value <- insect_vals[insect_top_idx]

# add SSBID
top_points_insects$SSBID <- extract(insect_stack$ssbid, top_points_insects)$ssbid
top_points_insects

# join with SSB rutenett
ssb_ru_insects <-rutenett_250 %>% 
  left_join(st_drop_geometry(st_as_sf(top_points_insects)),
            by = "SSBID")
mapview(ssb_ru_insects %>% filter(!is.na(rank)), zcol = "predicted_value")

# export
top_points_insects %>% 
  writeVector("vector/insect_ruteutvalg_points.geojson")
ssb_ru_insects %>% filter(!is.na(rank)) %>% 
  select(SSBID, rank, predicted_value) %>% 
  st_write("vector/insect_ruteutvalg_rutenett.geojson")


# Combine ranks and find best combined ------------------------------------

# drop unnecessary columns & filter out NA vals
top_grid_plants <- ssb_ru_plants %>% 
  select(-c(RSIZE, ROW, COL, XCOOR, YCOOR)) %>% 
  filter(!is.na(rank))

top_grid_insects <- ssb_ru_insects %>% 
  select(-c(RSIZE, ROW, COL, XCOOR, YCOOR)) %>% 
  filter(!is.na(rank))

# weigh plants higher than insects
top_grid_plants <- top_grid_plants %>%
  mutate(rank_plants500 = 501 - rank) %>%
  rename(rank_plants1 = rank,
         predicted_value_plants = predicted_value) %>%
  as.data.frame()

top_grid_insects <- top_grid_insects %>%
  mutate(rank_insects500 = 501 - rank) %>%
  rename(rank_insects1 = rank,
         predicted_value_insects = predicted_value) %>%
  as.data.frame()

# combine
top_grid <- top_grid_plants %>%
  full_join(top_grid_insects, by ="geometry") %>%
  mutate(SSBID = coalesce(SSBID.x, SSBID.y)) %>% 
  select(-c(SSBID.x, SSBID.y)) %>% 
  # set NA to 0
  mutate(across(2:4, ~ tidyr::replace_na(.x, 0))) %>%
  # get combined rank
  mutate(combined_score = 0.75 * rank_plants500 + 0.25 * rank_insects500,
         combined_rank  = rank(-combined_score, ties.method = "average")) %>%
  st_as_sf(crs = 25833) %>% 
  select(SSBID, combined_rank, combined_score, rank_plants1:rank_insects500)

top_grid200 <- top_grid %>%
  filter(combined_rank < 201)

top_grid500 <- top_grid %>% 
  slice_min(combined_rank, n = 500)


mapview(top_grid, zcol="combined_rank")
mapview(top_grid200, zcol = "combined_rank")

top_grid500 %>% st_write("vector/combined_ruteutvalg_top_grid500.geojson")
