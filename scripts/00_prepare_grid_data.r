# Clean and export TidligVarsling grid data
# Jenny Hansen
# 10 January 2024

# working in TidligVarsling project

Sys.setlocale("LC_CTYPE", "nb_NO.UTF-8")

# Load required libraries -------------------------------------------------

library(dplyr)
library(readxl)
library(sf)
library(mapview)
library(stringr)

# Import data -------------------------------------------------------------

# data from Rannveig email 19 December 2023

ruter_2019 <- st_read("vector/Ruter_utplasserte_tidl_varsl_2019_250m.shp") %>% 
  st_cast("POLYGON") %>% 
  select(name, type = Utvalg) %>% 
  mutate(year = 2019)
ruter_2020 <- st_read("vector/Ruter_tidvars_2020.shp") %>% 
  st_cast("POLYGON") %>% 
  st_transform(st_crs(ruter_2019)) %>% 
  st_zm() %>% 
  select(name, type) %>% 
  mutate(year = 2020)
ruter_2021 <- st_read("vector/Ruter_Tidvars_2021a.shp") %>% 
  st_transform(st_crs(ruter_2019)) %>% 
  select(name, type = utvalgsmet) %>% 
  mutate(year = 2021)
ruter_2022 <- st_read("vector/Ruter 2022.shp") %>% 
  st_transform(st_crs(ruter_2019)) %>% 
  select(name, type = valgt) %>% 
  mutate(year = 2022,
         name = if_else(name == "Mjøndalen", "Mile", name)) # fix misnaming
ruter_2023 <- st_read("vector/Ruter_2023_Tidvars.shp") %>% 
  st_transform(st_crs(ruter_2019)) %>% 
  select(name = locality, type = habitat_ty) %>% 
  mutate(year = 2023)


# Clean and combine -------------------------------------------------------

ruter_2019 <- ruter_2019 %>% 
  mutate(name = str_replace(name, "^[^_]*_", ""))
ruter_2021 <- ruter_2021 %>% 
  mutate(type = if_else(type == "Subjektiv", "Manuelt", "Auto"))
ruter_2022 <- ruter_2022 %>% 
  mutate(type = if_else(type == "Tid_man", "Manuelt", "Auto"))
ruter_2023 <- ruter_2023 %>% 
  mutate(type = if_else(type == "Tid_man", "Manuelt", "Auto"))

ruter_list <- list(ruter_2019, ruter_2020, ruter_2021, ruter_2022, ruter_2023)

combined_sf <- do.call(rbind, ruter_list) %>% 
  mutate(year = as.factor(year)) %>% 
  arrange(name)

# naming follows notes in 'Ruter_alle_Tidlig varsling.xlsx' from Rannveig

# double check Konnerud & Konnerud Skole
combined_clean <- combined_sf %>% 
  mutate(name = case_when(name == "Borregaard" ~ "Borregård",
                          name == "Borregard" ~ "Borregård",
                          name == "Grytnes_ungdomsskole" ~ "Grytnes",
                          name == "Heroya_industri" ~ "Herøya_industri",
                          name == "Kjelsaas" ~ "Kjelsås",
                          name == "Lofteroed" ~ "Lofterød",
                          name == "Mjoendalen" ~ "Mjøndalen",
                          name == "Ora" ~ "Øra",
                          name == "Presteroedkilen" ~ "Presterødkilen",
                          name == "Raaholt" ~ "Råholt",
                          name == "Rolvsoey" ~ "Rolvsøy",
                          name == "Saetre" ~ "Sætre",
                          name == "Skjerkoya" ~ "Skjerkøya",
                          (name == "Slagentangen" & year == 2019) ~ "Slagentangen_I",
                          (name == "Slagentangen" & year == 2020) ~ "Slagentangen_II",
                          (name == "Slagentangen" & year == 2022) ~ "Slagentangen_III",
                          (name == "Slagentangen_I" & year == 2023) ~ "Slagentangen_II",
                          (name == "Solgaard" & year == 2019) ~ "Solgård_II",
                          (name == "Solgaard" & year == 2020) ~ "Solgård",
                          (name == "Solgård" & year == 2022) ~ "Solgård_II",
                          name == "Tofte_havn" ~ "Tofte",
                          name == "Tofte_I" ~ "Tofte",
                          name == "Tonsberg" ~ "Tønsberg",
                          TRUE ~ name)) 

unique(combined_clean$name)

# Map ---------------------------------------------------------------------

mapview(combined_clean, zcol = "year")

# Export ------------------------------------------------------------------

st_write(combined_clean, "vector/cleaned_grid.geojson")
