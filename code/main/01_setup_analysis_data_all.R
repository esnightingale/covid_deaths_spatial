################################################################################
# DATA SETUP - ALL
# 
# Take cleaned linelists of deaths and cases, calculate E over specified time
# period, aggregate by week and LTLA and merge all
# 
################################################################################

library(tidyverse)
library(lubridate)

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## LTLA covariates
covs <- readRDS(paste0(datadir,"covs.rds"))

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))
regions.df <- st_drop_geometry(regions)

list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

# source(here::here("code","main","01_calculate_E.R"))
# source(here::here("code","main","01_setup_analysis_data.R"))


## DEATHS ##
source(here::here("code","main","00_setup_deaths.R"))

waves <- list(first = c(ymd("2020-01-01"), ymd("2020-06-30")),
              second = c(ymd("2020-07-01"), max(alldata$date, na.rm = T)))

deaths <- lapply(waves,
                 function(wave) setup_analysis_data(alldata, 
                                                    measure = "d",
                                                    start = wave[1], 
                                                    end = wave[2]))

deaths <- rlist::list.append(deaths, breaks = waves)

saveRDS(deaths, here::here("data","deaths.rds"))


## CASES ##
source(here::here("code","main","00_setup_cases.R"))

cases <- lapply(waves,
                function(wave) setup_analysis_data(alldata, 
                                                   measure = "c",
                                                   start = wave[1], 
                                                   end = wave[2]))

cases <- rlist::list.append(cases, breaks = waves)

saveRDS(cases, here::here("data","cases.rds"))


## Merge deaths and cases

first <- full_join(deaths[[1]], cases[[1]], 
                   by = c("w","week","la","lad19cd","lad19nm", "la_pop","geog","geography",
                          "area_km2","pop_dens","IMD","IMD_quint","prop_minority","prop_kw","w2","w3"),
                   suffix = c("_d","_c"))

second <- full_join(deaths[[2]], cases[[2]], 
                   by = c("w","week","la","lad19cd","lad19nm", "la_pop","geog","geography",
                          "area_km2","pop_dens","IMD","IMD_quint","prop_minority","prop_kw","w2","w3"),
                   suffix = c("_d","_c"))

merged <- list(first = first, second = second, breaks = waves)

saveRDS(merged, here::here("data","merged.rds"))

