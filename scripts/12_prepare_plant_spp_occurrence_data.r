# Script to prepare plant data for modeling
# Jenny Hansen
# 26 September 2025
# updated 03 November 2025 for additional plant data

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(readxl)
library(mapview)
library(tidyr)

# Import data -------------------------------------------------------------

# SSB grid over prediction area
ssb <- st_read("vector/ssb250_studyarea.shp")

# grid data for joining to plant data
samp_grid <- st_read("/data/R/Prosjekter/15821000_tidlig_oppdagelse_og_varsling_av_fremmede_arter/Ida/Data/GIS/cleaned_grid_18_25.geojson")

# plant data from Ida
plants <- read_excel("/data/R/Prosjekter/15821000_tidlig_oppdagelse_og_varsling_av_fremmede_arter/Ida/Data/Training_data/tidvars_18_25_plant_training.xlsx") %>%
  rename(lokname = "locality_v3")


# Make plants spatial -----------------------------------------------------

plants_sf_field <- plants %>%
  filter(coltype == "field") %>%
  left_join(samp_grid %>% st_drop_geometry(), by = c("lokname" = "name")) %>%
  left_join(samp_grid %>% select(name, geometry),  by = c("lokname" = "name")) %>%
  st_as_sf()

plants_sf_artskart <- plants %>%
  filter(coltype == "artskart") %>%
  st_as_sf(coords=c("decimalLongitude","decimalLatitude"), crs = 25833) 

plants_sf_artskart <- ssb %>%
  st_join(plants_sf_artskart, join = st_intersects) %>%
  filter(!is.na(species)) %>%
  select(species, year, lokname, coltype, SSBID)

plants_sf <- plants_sf_field %>%
  bind_rows(plants_sf_artskart) %>%
  # Fill lokname in Artskart too if geometry is the same
  group_by(geometry) %>%
  fill(lokname, .direction = "downup") %>%
  ungroup() %>%
  # Add SSB id to lokname for the ones that have NA
  mutate(lokname=ifelse(is.na(lokname),
                        SSBID,
                        lokname)) %>%
  select(-SSBID)


# Get richness per sampling location --------------------------------------

# get richness per site while accounting for variable sampling effort
plants_richness <- plants_sf %>% 
  group_by(lokname, year) %>%
  summarise(unique_species = n_distinct(species), .groups = "drop") %>%
  group_by(lokname) %>%
  summarise(
    species_richness = mean(unique_species, na.rm = TRUE),
    n_years = n_distinct(year),
    geometry = st_centroid(st_union(geometry)), 
    .groups = "drop"
  ) %>% 
  select(-n_years)

mapview(plants_richness, zcol = "species_richness")

# export richness data
st_write(plants_richness, "vector/plant_spp_presence_data.geojson")
#st_write(plants_richness, "/data/R/Prosjekter/15821000_tidlig_oppdagelse_og_varsling_av_fremmede_arter/Ida/Data/Training_data/plant_spp_presence_data.geojson")
