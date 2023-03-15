# R packages used in 3135 module
# These can be manually installed from here in case of any issues on cloud

if (!require("pacman")) install.packages("pacman")
pkgs = 
  c("dplyr", "nycflights13", "haven", "lattice", "latticeExtra", # Introduction to R
    "raster", "sf", "readr", "rasterVis", "tmap", # Spatial data in R (dplyr also needed)
    "PrevMap", "ggmap", "Hmisc", "automap", "sp", "raster", "rgdal", "psych", "geoR", "rgeos", # Continuous risk (sf and tmap also needed)
    "raster", "sp", "mgcv", "nlme", # Ecological niche modelling
    "mapview", "leaflet", "ggplot2", "spdep", "spatialreg", "rmapshaper", # Discrete data in R (tmap, sf, dplyr also needed)
    "tidyverse", "readxl","fasterize","here", "janitor", # Extracting and managing spatial data (also, sf, leaflet, raster, tmap)
     #"IDSpatialStats", # archived from CRAN as depends on spatstat.core
    "truncnorm", "MASS", # tau extension (also needs dplyr)
    "spatstat", "readr", "sparr" # point process (also mapview and sparr)
  ) # End
pacman::p_load(pkgs, character.only = T)
