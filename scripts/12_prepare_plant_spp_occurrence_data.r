# Script to prepare plant data for modeling
# Jenny Hansen
# 26 September 2025

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(readxl)
library(mapview)

# Import data -------------------------------------------------------------

# SSB grid over prediction area
ssb <- st_read("vector/ssb250_studyarea.shp")

# grid data for joining to plant data
samp_grid <- st_read("vector/cleaned_grid.geojson") %>% 
  select(name) %>% 
  group_by(name) %>% 
  slice(1) %>% 
  ungroup()

# plant data from Ida
plants <- read_excel("/data/R/Prosjekter/15821000_tidlig_oppdagelse_og_varsling_av_fremmede_arter/Ida/Data/Fielddata/Processed/Tidvars_planter_2018-2023_traits.xlsx") %>%
  rename(lokname = "locality_v3") %>%
  filter(year != 2018,
         !Fremmedartsstatus %in% c("Etablert per år 1800", "Ikke fremmed"),
         Tvilsom != 1,
         Potensiell_ny_art == 1) 


# Make plants spatial -----------------------------------------------------

plants_sf <- plants %>%
  left_join(samp_grid %>% st_drop_geometry(), by = c("lokname" = "name")) %>%
  left_join(samp_grid %>% select(name, geometry),  by = c("lokname" = "name")) %>%
  st_as_sf()


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
