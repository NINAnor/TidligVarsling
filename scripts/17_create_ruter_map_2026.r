# script to collect and name all the ruter to be sampled for the 
# TidligVarsling project in 2026
# Jenny Hansen
# 30 April 2026

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(mapview)
library(terra)

# Import data -------------------------------------------------------------

# municipalities (for naming)
muni <- st_read("/data/R/GeoSpatialData/AdministrativeUnits/Norway_AdministrativeUnits/Original/Norway_Municipalities/Basisdata_0000_Norge_25833_Kommuner_FGDB/Basisdata_0000_Norge_25833_Kommuner_FGDB.gdb",
                layer = "kommune")

# ruter files were received from Rannveig on 29 April 2026
ruter_high <- st_read("vector/utvalgte_ruter_AE_HS.geojson")
ruter_low <- st_read("vector/utvalgte_ruter_lav_sannsynlighet_annotated.geojson")
ruter_2018 <- st_read("vector/Ruter_tidvars_2018.shp") |> 
  st_transform(25833) |> 
  st_zm()

# SSB rute to assign to 2018 data
ssb <- rast("raster/sa_grid_250m_fixed.tif")

# Filter and combine ------------------------------------------------------

select_high <- ruter_high |> 
  filter(Vurdert == "OK" | Vurdert == "Ok") |> 
  mutate(source = "hoy_sannsynlighet") |> 
  select(SSBID, source)

select_low <- ruter_low |> 
  filter(utvalg == "ok") |> 
  mutate(source = "lav_sannsynlighet") |> 
  select(SSBID, source)

select_2018 <- ruter_2018 |> 
  filter(name %in% c("1006", "1007", "189", "1001", "48")) 

# assign SSBID to 2018 ruter
cents <- st_centroid(select_2018)

# extract SSBID at centroid location
cents$SSBID <- extract(ssb, vect(cents))[,2]

# join SSB data to 2018 ruter, assign source
select_2018_ssb <- select_2018 |> 
  left_join(cents |> st_drop_geometry()) |> 
  mutate(source = "ruter_2018") |> 
  select(SSBID, source)

# combine datasets
ruter_combined <- rbind(select_high, select_low, select_2018_ssb)
mapview(ruter_combined, zcol = "source")

# drop høy sannsynlightet rute that overlaps 189 from 2018
ruter_combined <- ruter_combined |> 
  filter(!(SSBID == 22675006569000 & source == "hoy_sannsynlighet"))

# Assign unique ruter names -----------------------------------------------

# intersect ruter with municipalities and assign names
ruter_muni <- ruter_combined |> 
  st_join(muni |> select(kommunenavn)) |> 
  arrange(kommunenavn)

# it's clear one of these rute is overlapping two municipalities- need to fix
ruter_muni <- ruter_combined |> 
  st_intersection(muni |>  select(kommunenavn)) |> 
  mutate(area = st_area(geometry)) |> 
  group_by(SSBID) |> 
  slice_max(area, n = 1) |> # assign to municipality with greatest overlap
  ungroup()

# assign unique number for each ruter per municipality
ruter_numbered <- ruter_muni |> 
  group_by(kommunenavn) |> 
  mutate(id_num = row_number()) |> 
  ungroup()

# assign ruter ID
ruter_2026 <- ruter_numbered |> 
  tidyr::unite(ruter_id, c("kommunenavn", "id_num"), sep = "_") |> 
  select(-area)

# export 
st_write(ruter_2026, "vector/ruter_2026.geojson")
