# Script to rasterize AR5 data and extract values over prediction area raster
# Jenny Hansen
# 21 August 2024

# working in TidligVarsling project

# Load required libraries -------------------------------------------------

library(dplyr)
library(sf)
library(terra)
library(fasterize)
library(rgrass)

# connect to GRASS mapset
NinaR::grassConnect()

# Setup GRASS environment -------------------------------------------------

execGRASS("g.region", raster = "tv_sa_250_grid", flags = "p")
execGRASS("g.remove", type = "raster", name = "MASK", flags = "f")
execGRASS("r.mask", vector = "norway_limits_detailed@p_sam_tools")

# Import data -------------------------------------------------------------

# AR5 vector data was clipped to the prediction area polygon (with 5 km buffer)
# in QGIS and saved for work in R (faster)

ar5 <- st_read("vector/ar5_pred_area.gpkg",
               layer = "ar5")

# Add human-friendly codes ------------------------------------------------

ar5 <- ar5 %>% 
  mutate(landcover = case_when(arealtype == 21 ~ "Fulldyrka jord",
                               arealtype == 22 ~ "Overflatedyrka jord",
                               arealtype == 23 ~ "Innmarksbeite",
                               arealtype == 30 ~ "Skog",
                               arealtype == 60 ~ "Myr",
                               arealtype == 50 ~ "Aapen fastmark",
                               arealtype == 81 ~ "Ferskvann",
                               arealtype == 82 ~ "Hav",
                               arealtype == 11 ~ "Bebygd",
                               arealtype == 12 ~ "Samferdsel",
                               arealtype == 99 ~ "Unknown"))


# Map colors:
# "#FFD16E" # Fulldyrka jord
# "#FFFF4C" # Overflatedyrka jord
# "#FFFFAD" # Innmarksbeite
# "#9ECC73" # Skog
# "#D9D9D9" # Åpen fastmark
# "#D1D1FF" # Myr
# "#91E7FF" # Ferskvann
# "#CCFEFE" # Hav
# "#FCDBD6" # Bebygd
# "#B3784C" # Samferdsel
# "#FFFFFF" # Unknown 

pal <- c("#D9D9D9", "#FCDBD6", "#91E7FF",  "#FFD16E",  "#CCFEFE","#FFFFAD",
         "#D1D1FF", "#FFFF4C", "#B3784C", "#9ECC73")
#mapview(ar5, zcol = "landcover", col.regions = pal)

# this takes over an hour to write
write_VECT(vect(ar5), "ar5_tv_study")

# Rasterize AR5 -----------------------------------------------------------

template <- rast(vect(ar5), res = 3) 

# convert to a raster
template_raster <- raster::raster(template)

# rasterize (much faster than terra)
ar5_raster <- fasterize(ar5, template_raster, field = "arealtype")
plot(ar5_raster)
rm(ar5, template, template_raster)

#raster::writeRaster(ar5_raster, "raster/ar5_raster_3m.tif")

# convert back to terra

ar5_rast <- rast(ar5_raster)

# plot color dataframe
land_cover_colors <- data.frame(
  value = c(11, 12, 21, 22, 23, 30, 50, 60, 81, 82, 99),
  color = c("#FCDBD6",  # Bebygd
            "#B3784C",  # Samferdsel
            "#FFD16E",  # Fulldyrka jord
            "#FFFF4C",  # Overflatedyrka jord
            "#FFFFAD",  # Innmarksbeite
            "#9ECC73",  # Skog
            "#D9D9D9",  # Åpen fastmark
            "#D1D1FF",  # Myr
            "#91E7FF",  # Ferskvann
            "#CCFEFE",  # Hav
            "#FFFFFF"   # Unknown
  ))

plot(ar5_rast, col = land_cover_colors)

terra::writeRaster(ar5_rast, "raster/ar5_rast.tif")


# Rasterize in GRASS ------------------------------------------------------

# cleaning is necessary (for accurate topolgy) and also takes a long time :()
execGRASS("v.clean", input = "ar5_tv_study", 
          output = "ar5_tv_study_clean", 
          tool = "break", flags = c("overwrite"))

# rasterize at a 10 m resolution for better percent cover estimation
execGRASS("g.region", vector = "ar5_tv_study_clean", res = "10")

execGRASS("v.to.rast",
          input = "ar5_tv_study_clean",
          output = "ar5_tv_study_rast_10m",
          use = "attr",
          attribute_column = "arealtype",
          type = "area",
          flags = "overwrite")


ar5_rast <- read_RAST("ar5_tv_study_rast")
