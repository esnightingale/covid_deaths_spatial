################################################################################
# Description: Summarise and visualise death and case time series per LTLA. 
# Produce descriptive figures for paper.
# 
# Author: Emily S Nightingale
# Date created: 30/09/2020
# 
################################################################################
################################################################################

################################################################################
# SETUP
################################################################################

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

source(here::here("code","main","functions.R"))
source(here::here("code","main","plot_functions.R"))

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

# Figure output directory
# figdir <- "figures/first_wave"
# figdir <- "figures/second_wave"
figdir <- "figures/descriptive"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  st_set_crs("+OSGB:1936 +units=m +no_defs") %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
merged <- readRDS(here::here("data","merged.rds")) 

deaths <- readRDS(here::here("data","deaths.rds")) 
cases <- readRDS(here::here("data","cases.rds")) 

data.list <- list(deaths = deaths,cases = deaths)

theme_set(theme_minimal())

# ---------------------------------------------------------------------------- #

## COVARIATES ##

cov_names <- c("pop_dens", "IMD", "prop_minority","prop_kw")
deaths[[1]] %>%
  group_by(lad19cd) %>%
  summarise_at(vars(cov_names), .funs = base::mean) -> covs

regions %>%
  full_join(covs) %>%
  mutate(pop_dens = as.numeric(pop_dens)) -> regions_wcovs

map_dens <-
  basic_map(regions_wcovs, fill = "pop_dens") +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill = "", title = "Population per KM^2") +
  theme(plot.title = element_text(size=10))

map_imd <-
  basic_map(regions_wcovs, fill = "IMD") +
  labs(fill = "", title = "Index of Multiple Deprivation") +
  theme(plot.title = element_text(size=10))

map_mino <-
  basic_map(regions_wcovs, fill = "prop_minority") +
  labs(fill = "", title = "Proportion of black and \nminority ethnic population") +
  scale_fill_viridis_c(trans = "log10") +
  theme(plot.title = element_text(size=10))

map_kw <-
  basic_map(regions_wcovs, fill = "prop_kw") +
  labs(fill = "", title = "Proportion of population\n classified as key workers") +
  theme(plot.title = element_text(size=10))

png(here::here("figures","descriptive","map_covariates.png"), height = 1200, width = 1200, res = 150)
(map_dens + map_imd) /
  (map_mino + map_kw)
dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - by LTLA ##

titles <- matrix(c("COVID19-related deaths in England, by LTLA and week of death",
                   "COVID19-related deaths in England, by LTLA and week since first LTLA death",
                   "Confirmed COVID-19 cases in England, by LTLA and week of specimen",
                   "Confirmed COVID-19 cases in England, by LTLA and week since first LTLA case"),
                 nrow = 2, ncol = 2)

for (dat in seq_along(data.list)){
  for (wave in 1:2){
    
    print(sprintf("ts_%s_%s.png",names(data.list)[dat],wave))
    
    png(here::here(figdir,sprintf("ts_%s_%s.png",names(data.list)[dat],wave)), height = 800, width = 1200, res = 150)
    print(plot_la_ts(data.list[[dat]], wave, title1 = titles[1,dat], title2 = titles[dat,2]))
    dev.off()
  }
}

# ---------------------------------------------------------------------------- #

## MAP TOTALS ##

png(here::here(figdir,"map_totals.png"), height = 1500, width = 1500, res = 150)
  (plot_la_tot(cases,1, title = "Cases per 100,000") + plot_la_tot(cases,2, title = "Cases per 100,000")) / 
    (plot_la_tot(deaths,1, title = "Deaths per 100,000") + plot_la_tot(deaths,2, title = "Deaths per 100,000"))
dev.off()

# ---------------------------------------------------------------------------- #

## MAP EPIDEMIC TIMING ##

data.list <- list(cases = cases, deaths = deaths)

for (dat in seq_along(data.list)){
  for (wave in 1:2){
    
    print(sprintf("map_timing_%s_%s.png",names(data.list)[dat],wave))
    
    png(here::here(figdir,sprintf("map_timing_%s_%s.png",names(data.list)[dat],wave)), height = 800, width = 1200, res = 150)
    print(plot_epi_time(data.list[[dat]], wave, measure = names(data.list)[dat]))
    dev.off()
  }
}

# ---------------------------------------------------------------------------- #

## TIME SERIES - BY GEOGRAPHY ##

titles <- matrix(c("COVID19-related deaths in England, by geography and week of death",
                   "by week since first LTLA death",
                   "Confirmed COVID-19 cases in England, by geography and week of specimen",
                   "by week since first LTLA case"),
                 nrow = 2, ncol = 2)


for (dat in seq_along(data.list)){
  for (wave in 1:2){
    
    print(sprintf("ts_geog_%s_%s.png",names(data.list)[dat],wave))
    
    png(here::here(figdir,sprintf("ts_geog_%s_%s.png",names(data.list)[dat],wave)), height = 800, width = 1200, res = 150)
    print(plot_geog_ts(data.list[[dat]], wave, title1 = titles[1,dat], title2 = titles[2,dat]))
    dev.off()
  }
}


# png(here::here(figdir,"map_timing_ts_geog.png"), height = 1000, width = 1400, res = 150)
# (first_map + ts_geog_week ) / (peak_map + ts_geog_epiwk) + plot_annotation(tag_levels = 'A')
# dev.off()


# ---------------------------------------------------------------------------- #

## Compare case-death lag per LTLA ##

merged[[1]] %>%
  group_by(lad19cd, la_pop) %>%
  # mutate(n_c_roll = rollmean(n_c, 2),
  #        n_d_roll = rollmean(n_d, 2)) %>%
  summarise(peak_c = week[which.max(n_c/la_pop)],
            peak_d = week[which.max(n_d/la_pop)],
            lag_c_d = peak_d - peak_c) %>%
  right_join(regions) %>%
  basic_map(fill = "lag_c_d", rate1e5 = FALSE) +
  scale_fill_gradient2() +
  labs(title = "Lag between peak deaths and peak cases per 100,000",
       subtitle = paste(merged$breaks[[1]][1],"-",merged$breaks[[1]][2])) -> map_lag_rate_first


merged[[2]] %>%
  group_by(lad19cd, la_pop) %>%
  # mutate(n_c_roll = rollmean(n_c, 2),
  #        n_d_roll = rollmean(n_d, 2)) %>%
  summarise(peak_c = week[which.max(n_c/la_pop)],
            peak_d = week[which.max(n_d/la_pop)],
            lag_c_d = peak_d - peak_c) %>%
  right_join(regions) %>%
  basic_map(fill = "lag_c_d", rate1e5 = FALSE) +
  scale_fill_gradient2() +
  labs(title = "Lag between peak deaths and peak cases per 100,000",
       subtitle = paste(merged$breaks[[2]][1],"-",merged$breaks[[2]][2])) -> map_lag_rate_second

png(here::here(figdir,"map_case_death_lags.png"), height = 800, width = 1200, res = 150)
map_lag_rate_first + map_lag_rate_second
dev.off()
