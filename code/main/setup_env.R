
library(tidyverse)
library(linelist)
library(INLA)
library(ggplot2)
library(rgdal)
library(spdep)
library(sf)
library(patchwork)
library(viridis)
library(here)
library(ggspatial)

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)
