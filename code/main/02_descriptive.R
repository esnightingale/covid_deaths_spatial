################################################################################
# Description: Summarise and visualise death time series per LTLA. Produce 
# descriptive figures for paper.
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

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

# Data version 
data <- "Deaths/dat_deaths_2020-01-01_2020-06-30.rds"
data2 <- "Deaths/dat_deaths_2020-07-01_2020-11-05.rds"
# data <- "Cases/dat_cases_2020-01-01_2020-05-30.rds" 

period <- unlist(strsplit(gsub(".rds","",data), split = "_"))[3:4]
period2 <- unlist(strsplit(gsub(".rds","",data2), split = "_"))[3:4]

# Figure output directory
# figdir <- "figures/first_wave"
# figdir <- "figures/second_wave"
figdir <- "figures"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  st_set_crs("+OSGB:1936 +units=m +no_defs") %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
dat <- readRDS(paste0(datadir, data)) 
dat2 <- readRDS(paste0(datadir, data2)) 

theme_set(theme_minimal())

# ---------------------------------------------------------------------------- #

## COVARIATES ##

# cov_names <- c("pop_dens", "IMD", "prop_minority","prop_kw")
# dat %>%
#   group_by(lad19cd) %>%
#   summarise_at(vars(cov_names), .funs = base::mean) -> covs
# 
# regions %>%
#   full_join(covs) -> regions_wcovsgit
# 
# map_dens <-
#   ggplot(regions_wcovs, aes(geometry = geometry, fill = pop_dens)) +
#   geom_sf() +
#   scale_fill_viridis_c() +
#   map_theme() +
#   scale_fill_viridis_c(trans = "log10") +
#   labs(fill = "", title = "Population per KM^2") +
#   theme(plot.title = element_text(size=10))
# 
# map_imd <-
#   ggplot(regions_wcovs, aes(geometry = geometry, fill = IMD)) +
#   geom_sf() +
#   scale_fill_viridis_c() +
#   map_theme() +
#   labs(fill = "", title = "Index of Multiple Deprivation") +
#   theme(plot.title = element_text(size=10))
# 
# map_mino <-
#   ggplot(regions_wcovs, aes(geometry = geometry, fill = prop_minority)) +
#   geom_sf() +
#   scale_fill_viridis_c() +
#   map_theme() +
#   labs(fill = "", title = "Proportion of black and \nminority ethnic population") +
#   scale_fill_viridis_c(trans = "log10") +
#   theme(plot.title = element_text(size=10))
# 
# map_kw <-
#   ggplot(regions_wcovs, aes(geometry = geometry, fill = prop_kw)) +
#   geom_sf() +
#   scale_fill_viridis_c() +
#   map_theme() +
#   labs(fill = "", title = "Proportion of population\n classified as key workers") +
#   theme(plot.title = element_text(size=10))
# 
# png(here::here("figures","map_covariates.png"), height = 1200, width = 1200, res = 150)
# (map_dens + map_imd) /
#   (map_mino + map_kw)
# dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - BY LTLA ##

dat %>%
  group_by(lad19cd) %>%
  mutate(peak = w[which.max(n)]) %>% 
  arrange(peak) %>%
  ggplot(aes(jitter(w), n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
  geom_line(alpha = 0.3) +
  labs(title = "COVID19-related deaths in England, by LTLA and week of death",
       subtitle = paste(period[1],"-",period[2]),
       x = "Calendar week",
       y = "Deaths per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  scale_colour_viridis_d(option = "cividis")  -> ts_by_wk

dat %>%
  group_by(lad19cd) %>%
  mutate(peak = wk_since_first[which.max(n)]) %>% 
  arrange(peak) %>%
  ggplot(aes(wk_since_first, n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
  geom_line(alpha = 0.3) +
  labs(title = "COVID19-related deaths in England, by LTLA and week since first LTLA death",
       x = "Weeks since first death",
       y = "Deaths per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  scale_colour_viridis_d(option = "cividis") -> ts_by_epiwk


png(here::here(figdir,"ts_ltla_peak.png"), height = 800, width = 1200, res = 150)
ts_by_wk / ts_by_epiwk
dev.off()

# ---------------------------------------------------------------------------- #

## MAP TOTALS ##


dat %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n)) %>% 
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE) +
  labs(title = "Total deaths per 100,000",
       subtitle = paste(period[1],"-",period[2])) + 
  map_theme() -> map_tot

dat2 %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n)) %>% 
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE) +
  labs(subtitle = paste(period2[1],"-",period2[2])) + 
  map_theme() -> map_tot2

dat %>%
  group_by(lad19cd) %>%
  summarise(n = unique(E)) %>% 
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE) + 
  map_theme() +
  labs(title = "Expected deaths per 100,000", subtitle = "according to age distribution of population") -> map_E

png(here::here(figdir,"map_totals.png"), height = 800, width = 1200, res = 150)
map_tot + map_E
dev.off()

# png(here::here(figdir,"map_totals_compwaves.png"), height = 800, width = 1200, res = 150)
pdf(here::here(figdir,"map_totals_compwaves.pdf"), height = 7, width = 12)
map_tot + map_tot2
dev.off()

# ---------------------------------------------------------------------------- #

## MAP EPIDEMIC TIMING ##

dat %>%
  group_by(lad19cd) %>%
  summarise(first = unique(first)) %>% 
  full_join(regions) %>%
  basic_map(fill = "first") +
  scale_fill_viridis(option = "plasma") +
  ggtitle("Week of first death")+ 
  map_theme()  -> first_map

dat %>%
  group_by(lad19cd) %>%
  mutate(peak = w[which.max(n)]) %>% 
  arrange(peak) %>%
  full_join(regions) %>%
  basic_map(fill = "peak") +
  scale_fill_viridis(option = "plasma") +
  ggtitle("Week of peak number of deaths") + 
  map_theme() -> peak_map

png(here::here(figdir,"map_timing.png"), height = 800, width = 1200, res = 150)
first_map + peak_map
dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - BY LTLA ##

dat %>%
  ggplot(aes(wod, n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  labs(title = "COVID19-related deaths in England, by week of death",
       x = "Calendar week",
       y = "Deaths per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) -> ts_by_wk2

dat %>%
  ggplot(aes(wk_since_first, n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  labs(title = "COVID19-related deaths in England, by week since first LA death",
       x = "Calendar week",
       y = "Deaths per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) -> ts_by_epiwk2

png(here::here(figdir,"ts_ltla.png"), height = 800, width = 1200, res = 150)
ts_by_wk2 / ts_by_epiwk2
dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - BY GEOGRAPHY ##

dat %>%
  group_by(w,wod,geography) %>%
  summarise(n = sum(n, na.rm= T),
            geog_pop = sum(la_pop)) %>%
  ggplot(aes(wod, n*1e5/geog_pop, group = geography, col = geography)) +
  geom_line() +
  geom_point() +
  labs(title = "COVID19-related deaths in England, by geography type and week of death",
       x = "Calendar week",
       y = "Deaths per 100,000",
       colour = "Geography type"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  theme(legend.position = c(0.9,0.7), legend.text=element_text(size=10),  legend.title=element_text(size=10)) -> ts_geog_week

dat %>%
  group_by(wk_since_first,geography) %>%
  summarise(n = sum(n, na.rm= T),
            geog_pop = sum(la_pop)) %>%
  ggplot(aes(wk_since_first, n*1e5/geog_pop, group = geography, col = geography)) +
  geom_line() +
  geom_point() +
  labs(title = "by week since first LA death",
       x = "Weeks since first observed death",
       y = "Deaths per 100,000",
       colour = "Geography type"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  theme(legend.position = "none") -> ts_geog_epiwk

png(here::here(figdir,"ts_geog.png"), height = 1000, width = 1200, res = 150)
ts_geog_week / ts_geog_epiwk
dev.off()


png(here::here(figdir,"map_timing_ts_geog.png"), height = 1000, width = 1400, res = 150)
(first_map + ts_geog_week ) / (peak_map + ts_geog_epiwk) + plot_annotation(tag_levels = 'A')
dev.off()




