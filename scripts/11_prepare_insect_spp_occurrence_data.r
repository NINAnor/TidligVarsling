# Script to import insect data from Artskart and filter for potential new 
# invasive species to augment data collected by TidligVarsling surveys
# Jenny Hansen
# 25 September 2025
# updated 08 October 2025 after receiving new survey data

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(mapview)
library(rgbif)
library(readr)
library(lubridate)
library(terra)
library(stringr)


# Import data -------------------------------------------------------------

# SSB grid over prediction area
ssb <- st_read("vector/ssb250_studyarea.shp")

# tv insect species
insects <- st_read("vector/cleaned_insect_data.geojson") %>% 
  filter(final_sublist %in% c("potensiell_fremmede", 
                              "potensiell_fremmede_with_norsk_GBIF_occurrences", 
                              "not_in_GBIF",
                              "d\xf8rstokkart")) %>% 
  mutate(final_species_clean = str_replace_all(final_species, "_", " "))

# potential new sf
pnf <- st_read("vector/pot_nye_fremmede.geojson")

# National insect monitoring data downloaded from 
#https://www.gbif.org/occurrence/download?dataset_key=19fe96b0-0cf3-4a2e-90a5-7c1c19ac94ee
nim <- read_delim("vector/nim/occurrence.csv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)

# raster to calculate extent
r <- rast("raster/ndvi_summer.tif") %>%
  project("EPSG:4326")

# geometry to use as filter
wkt <- st_as_text(st_as_sfc(st_bbox(r)))
wkt

# download insect occurrences
insecta_key <- name_backbone(name = "Insecta")$usageKey
insecta_key

download_key <- occ_download(
  pred_and(
    pred("taxonKey", insecta_key),
    pred("geometry", wkt),
    pred("hasCoordinate", TRUE),
    pred("year", "1950,2025")  
  ),
  format = "SIMPLE_CSV"
)
download_key

occ_download_wait('0015076-250920141307145')

occ <- occ_download_get('0015076-250920141307145') %>%
  occ_download_import()


# Make spatial, filter ----------------------------------------------------

occ_select_cleaned <- occ %>%
  select(gbifID, occurrenceStatus, species, individualCount,
         x = decimalLongitude, y = decimalLatitude,
         coord_precision = coordinateUncertaintyInMeters, eventDate,
         basisOfRecord, institutionCode, catalogNumber, identifiedBy) %>%
  mutate(parsed_date = as.Date(eventDate, format = "%Y-%m-%d"),
         species_clean = word(species, 1, 2),  # keep only genus + species
         species_clean = str_trim(species_clean)) %>% 
  filter(!is.na(parsed_date),  
         occurrenceStatus != "ABSENT") %>%
  select(-parsed_date)


# filter for species occurring in TV study
occ_invasive <- occ_select_cleaned %>%
  filter(species_clean %in% insects$final_species_clean)

# make spatial 
occ_sf <- occ_invasive %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(25833)

mapview(occ_sf, col.regions = "dodgerblue") + 
  mapview(insects, col.regions = "orange")

# filter nim data
nim_invasive <- nim %>%
  select(gbifID, occurrenceStatus, species, individualCount,
         x = decimalLongitude, y = decimalLatitude,
         coord_precision = coordinateUncertaintyInMeters, eventDate,
         basisOfRecord, institutionCode, catalogNumber, identifiedBy) %>%
  filter(species %in% insects$final_species_clean,
         occurrenceStatus != "ABSENT")

# make spatial
nim_sf <- nim_invasive %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_intersection(st_as_sfc(st_bbox(r))) %>% 
  st_transform(25833) 

mapview(nim_sf)

mapview(occ_sf, col.regions = "dodgerblue") + 
  mapview(nim_sf, col.regions = "orange")

# combine occ_sf and nim_sf

# match columns before rbind
nim_sf2 <- nim_sf %>%
  mutate(
    institutionCode = as.character(institutionCode),
    catalogNumber   = as.character(catalogNumber),
    identifiedBy    = as.character(identifiedBy),
    species_clean   = species  # create to match occ_sf
  ) %>%
  select(names(occ_sf))  # reorder/align columns

all_sf <- rbind(occ_sf, nim_sf2)


# Remove duplicates -------------------------------------------------------

all_sf_unique <- all_sf %>%
  distinct(gbifID, .keep_all = TRUE)
mapview(all_sf_unique)

# Aggregate to 250 m grid -------------------------------------------------

# join with SSB grid
occ_with_grid <- st_join(all_sf_unique, ssb)
mapview(occ_with_grid)

# summarize spp. richness per overlapping grid
richness <- occ_with_grid %>%
  st_drop_geometry() %>%
  group_by(SSBID) %>%  # replace with your actual polygon ID column
  summarise(species_richness = n_distinct(species))

grid_richness <- ssb %>%
  left_join(richness, by = "SSBID") %>%
  mutate(species_richness = tidyr::replace_na(species_richness, 0))


# Combine TV data with other data -----------------------------------------

pnf_with_id <- st_join(pnf, ssb["SSBID"], join = st_within)

pnf_clean <- pnf_with_id %>%
  st_drop_geometry() %>%
  select(SSBID, mean_unique_species)

# combine and take highest value for the grid (if values differ)
spp_combined <- grid_richness %>%
  left_join(pnf_clean, by = "SSBID") %>%
  mutate(species_richness = pmax(species_richness, mean_unique_species, 
                                 na.rm = TRUE)) %>%
  select(-mean_unique_species)


# Centroid pts for sampling -----------------------------------------------

spp_pts <- spp_combined %>% st_centroid() %>% 
  filter(species_richness > 0)
mapview(spp_pts, zcol = "species_richness")

# export
st_write(spp_pts, "vector/insect_spp_presence_data.geojson")

# export only TV data
pnf %>% select(name, species_richness = mean_unique_species, geometry) %>% 
  st_write("vector/insect_spp_presence_data_tv_only.geojson")
