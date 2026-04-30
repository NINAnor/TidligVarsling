# script to get matrikkelen data for each of the 2026 ruter
# Jenny Hansen
# 30 April 2026

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(mapview)
library(stringr)
library(tidyr)
library(tibble)
library(writexl)

# Import data -------------------------------------------------------------

# ruter (from script 17)
ruter_2026 <- st_read("vector/ruter_2026.geojson")

# TEIG data 
# https://kartkatalog.geonorge.no/metadata/matrikkelen-eiendomskart-teig/74340c24-1c8a-4454-b813-bfe498e80f16
teig_viken <- st_read("/data/R/GeoSpatialData/CadastralParcels/Norway_Property/Processed/teig_Basisdata_30_Viken_25833_MatrikkelenEiendomskartTeig.shp")
teig_tv <- st_read("/data/R/GeoSpatialData/CadastralParcels/Norway_Property/Processed/teig_Basisdata_38_Telemark_og_Vestfold_25833_MatrikkelenEiendomskartTeig.shp")
teig <- rbind(teig_viken, teig_tv)


# Data munging teig -------------------------------------------------------

# create unique columns for gårsnummer (gnr) and bruksnummer (bnr)
teig_clean <- teig |>
  separate(matrikkeln, into = c("gnr", "rest"), # temp lump other numbers
            sep = "/", extra = "merge",  fill = "right")

# further split into bnr/fnr/snr
teig_clean <- teig_clean |>
  separate(rest,
           into = c("bnr_raw", "fnr", "snr"),
           sep = "/",
           fill = "right",   # handles missing fnr/snr
           extra = "merge")   # handles unexpected extra parts safely


# split cases where bnr has two numbers, e.g. 1-2
teig_expanded <- teig_clean |> 
  separate_rows(bnr_raw, sep = "[-,]") |> 
  rename(bnr = bnr_raw)

# create a property ID
teig_expanded <- teig_expanded |> 
  mutate(matrikkelen_id = paste(gnr, bnr, sep = "/"))

# Intersect with ruter ----------------------------------------------------

ruter_eiendom <- ruter_2026 |> 
  st_intersection(teig_expanded)

# remove excessive columns
ruter_eiendom <- ruter_eiendom |> 
  select(SSBID, source, ruter_id, matrikkelen_id,  gnr, bnr, fnr, snr, 
         kommunenavn = NAVN) |> 
  remove_rownames()

mapview(ruter_eiendom)

# Create a nice export table ----------------------------------------------

property_export <- ruter_eiendom |> 
  st_drop_geometry() |> 
  distinct(ruter_id, kommunenavn, gnr, bnr) |> 
  mutate(matrikkel = paste0(gnr, "/", bnr)) |> 
  arrange(ruter_id, gnr, bnr) |> 
  remove_rownames()

# export files
st_write(ruter_eiendom, "vector/ruter_2026_matrikkelen.geojson")
write_xlsx(property_export, "output/ruter_2026_matrikkelen.xlsx")
