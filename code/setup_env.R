################################################################################
# Description: Load packages, set input/output directories, set plotting options
# 
# Author: Emily S Nightingale
# Date created: 25/02/2021
# 
################################################################################
################################################################################

# devtools::install_github('timcdlucas/INLAutils')
# install.packages("INLA",
#                  repos=c(getOption("repos"),
#                          INLA="https://inla.r-inla-download.org/R/stable"), 
#                  dep=TRUE)

# Packages required
packages <- c("tidyverse","lubridate","data.table","dtplyr","readxl","sf","INLA",
              "spdep","rgdal","ggspatial","patchwork","scales","here","rlist")

# Check and install if necessary
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Local data directory (for raw linelists)
# datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

# Set default theme for plotting
theme_set(theme_minimal())

# Load shapefiles
regions <- readRDS(here::here("data","LA_shp_wpops.rds")) %>%
  st_set_crs("+OSGB:1936 +units=m +no_defs") %>%
  filter(grepl("E", lad19cd))

border <- st_union(regions)
regions.df <- st_drop_geometry(regions)

# Specify measure to fit to (cases or deaths) and wave (1 or 2)
measure <- "deaths"
wave <- 1

# Source utility functions
list.files(here::here("code","utils"), full.names = TRUE) %>% purrr::walk(source)

