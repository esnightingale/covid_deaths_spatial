################################################################################
# Description: Load packages, set input/output directories, set plotting options
# 
# Author: Emily S Nightingale
# Date created: 25/02/2021
# 
################################################################################
################################################################################

################################################################################
# SETUP
################################################################################

library(tidyverse)
library(lubridate)
library(linelist)
library(data.table)
library(dtplyr)
library(INLA)
library(inlabru)
library(ggplot2)
library(rgdal)
library(spdep)
library(sf)
library(ggspatial)
library(patchwork)
library(here)
library(scales)
# devtools::install_github('timcdlucas/INLAutils')

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

theme_set(theme_minimal())

