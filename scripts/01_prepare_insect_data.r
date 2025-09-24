# Script to get spatial data for insect groups
# Jenny Hansen
# 01 July 2024

# working in TidligVarsling project

Sys.setlocale("LC_CTYPE", "nb_NO.UTF-8")

# Load required libraries -------------------------------------------------

library(readr)
library(dplyr)
library(sf)
library(mapview)

# Import data -------------------------------------------------------------

# data came from Marie Davey (email from 08 May 2024)
insect_data <- 
  read.delim("data/tidlig_varsling_longform_combined_data_all_years_08052024.txt")

select_insect <- insect_data %>% 
  select(year, sample_id, sampling_name, sampling_number, UTM_east, UTM_north, 
         final_phylum, final_class,final_order, final_family, final_genus, 
         final_species, final_classification_level, final_sublist, ID_confidence, 
         Fremmedartsstatus, Etableringsklasse, extraction_type, 
         locality_selection_method) %>% 
  filter(extraction_type == "eDNA_lysed",
         ID_confidence %in% c("MODERATE", "HIGH"))


# Make spatial ------------------------------------------------------------

insect_sf <- select_insect %>% 
  st_as_sf(coords = c("UTM_east", "UTM_north"),
           crs = 25833)

# from 00_explore_grid_data.r
grid_locs <- st_read("vector/cleaned_grid.geojson") %>%
  select(name) %>% 
  group_by(name) %>% 
  slice(1) %>% 
  ungroup() %>% 
  st_buffer(10)

# add name column from grid_locs
insect_named <- st_join(insect_sf, grid_locs, join = st_intersects) %>% 
  filter(!is.na(name))

# write clean & filtered spatial data
st_write(insect_named, "vector/cleaned_insect_data.geojson")

# save filtered data for future work
write_csv(insect_named %>% st_drop_geometry(), "data/filtered_insect_data.csv")

# Aggregate by sampling location ------------------------------------------

# group & get centroids for trap locations
# these are averaged for each rute to account for
# uneven sampling effort

generate_spatial_data <- function(data, group_filter, group_name) {
  
  group_sf <- data %>%
    filter(final_sublist %in% group_filter) %>%
    group_by(name, year) %>%
    summarise(unique_species = n_distinct(final_species), .groups = 'drop') %>%
    group_by(name) %>%
    summarise(mean_unique_species = mean(unique_species), .groups = 'drop') %>% 
    group_by(geometry = st_sfc(lapply(geometry, function(geom) {
      if (st_geometry_type(geom) == "MULTIPOINT") {
        st_centroid(geom)
      } else {
        geom
      }
    }))) %>%
    ungroup()
  
  st_crs(group_sf) <- st_crs(data)
  
  # save the spatial data
  st_write(group_sf, paste0("vector/", group_name, ".geojson"))
  
  return(group_sf)
}

# apply function to each group
fremmede_sf <- generate_spatial_data(insect_named, "fremmede", 
                                     "fremmede")

mapview(fremmede_sf, cex = "mean_unique_species", zcol = "mean_unique_species")

potensielt_nye_sf <- generate_spatial_data(insect_named, c("potensiell_fremmede", 
                                                           "d\xf8rstokkart"), 
                                           "pot_nye_fremmede")
potensielt_nye_sf <- potensielt_nye_sf %>% 
  rename(mean_unique_species = total_abd)
mapview(potensielt_nye_sf, cex = "mean_unique_species", 
        zcol = "mean_unique_species")


overlooked_sf <- generate_spatial_data(insect_named, "overlooked", 
                                       "overlooked")
mapview(overlooked_sf, cex = "mean_unique_species", 
        zcol = "mean_unique_species")
