################################################################################
# DATA SETUP - ALL
# 
# Take cleaned linelists of deaths and cases, calculate E over specified time
# period, aggregate by week and LTLA and merge all
# 
################################################################################

library(tidyverse)
library(lubridate)
library(sf)

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## LTLA covariates
covs <- readRDS(paste0(datadir,"covs.rds"))

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
regions.df <- sf::st_drop_geometry(regions)

list.files(here::here("code","utils"), full.names = TRUE) %>% purrr::walk(source)

# Define time period by week start
waves <- list(first = c(ymd("2020-01-01"), ymd("2020-06-30")),
              second = c(ymd("2020-07-01"), ymd("2020-12-02")))

## DEATHS ##
source(here::here("code","main","00_3_setup_deaths.R"))

deaths <- lapply(waves,
                 function(wave) setup_analysis_data(alldata, 
                                                    measure = "d",
                                                    start = wave[1], 
                                                    end = wave[2]))

deaths <- rlist::list.append(deaths, breaks = waves)

saveRDS(deaths, here::here("data","deaths.rds"))


## CASES ##
source(here::here("code","main","00_4_setup_cases.R"))

cases <- lapply(waves,
                function(wave) setup_analysis_data(alldata, 
                                                   measure = "c",
                                                   start = wave[1], 
                                                   end = wave[2]))

cases <- rlist::list.append(cases, breaks = waves)

saveRDS(cases, here::here("data","cases.rds"))

################################################################################
################################################################################
