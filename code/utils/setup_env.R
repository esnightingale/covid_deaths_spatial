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
library(INLA)
library(inlabru)
library(ggplot2)
library(rgdal)
library(spdep)
library(sf)
library(patchwork)
library(here)
# devtools::install_github('timcdlucas/INLAutils')

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

theme_set(theme_minimal())

