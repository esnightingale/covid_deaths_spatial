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

source(here::here("code","functions.R"))

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

# Data version - spi-m linelists versus pipeline
# data <- "dat.rds" 
# data <- "dat_alt.rds" 
 
cases <- "Cases/dat_cases_2020-01-01_2020-06-30.rds" 
cases2 <- "Cases/dat_cases_2020-07-01_2020-11-05.rds" 
deaths <- "Deaths/dat_deaths_2020-01-01_2020-06-30.rds" 

period <- unlist(strsplit(gsub(".rds","",cases), split = "_"))[3:4]
period2 <- unlist(strsplit(gsub(".rds","",cases2), split = "_"))[3:4]

# period <- c(min(cases$wos), max(cases$wos))
# period2 <- c(min(cases2$wos), max(cases2$wos))

# Figure output directory
figdir <- "figures/Cases"
# figdir <- "figures/second_wave"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

# Cleaned linelist of cases
alldata <- readRDS(paste0(datadir, "Cases/alldata.rds"))

cases <- readRDS(paste0(datadir, cases)) 
cases2 <- readRDS(paste0(datadir, cases2)) 

# dat <- readRDS(paste0(datadir, paste0("Cases/", data))) %>%
#   filter(wor < ymd("2020-08-02")) %>% 
#   mutate(pop_dens = pop_dens_total) 

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
deaths <- readRDS(paste0(datadir, deaths)) 

# deaths <- readRDS(paste0(datadir, paste0("Deaths/", "dat_alt.rds"))) %>%
#   filter(wod < ymd("2020-08-02")) %>% 
#   mutate(pop_dens = pop_dens_total) 

theme_set(theme_minimal())

# ---------------------------------------------------------------------------- #
# Time series - total

png(here::here(figdir,"total_specdate_bypillar.png"), height = 600, width = 1000, res = 150)
alldata %>%
  mutate(pillar = factor(pillar, levels = c("pillar_1","pillar_2"), labels = c("Pillar 1","Pillar 2"))) %>%
  ggplot(aes(x = dos, fill = pillar)) +
  geom_bar() + 
  labs(x = "Date of specimen", y = "Count", fill = "") + 
  theme(legend.position = c(0.8, 0.9))
dev.off()

# ---------------------------------------------------------------------------- #


## TIME SERIES - BY LTLA ##

cases %>%
  group_by(lad19cd) %>%
  mutate(peak = w[which.max(n)]) %>% 
  arrange(peak) %>%
  ggplot(aes(jitter(w), n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
  geom_line(alpha = 0.3) +
  labs(title = "Cases of COVID-19 in England, by LTLA and week of reporting",
       x = "Calendar week",
       y = "Cases per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  scale_colour_viridis_d(option = "cividis")  -> ts_by_wk

cases %>%
  group_by(lad19cd) %>%
  mutate(peak = wk_since_first[which.max(n)]) %>% 
  arrange(peak) %>%
  ggplot(aes(wk_since_first, n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
  geom_line(alpha = 0.3) +
  labs(title = "Cases of COVID-19 in England, by LTLA and week since first LTLA case",
       x = "Calendar week",
       y = "Cases per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  scale_colour_viridis_d(option = "cividis") -> ts_by_epiwk


png(here::here(figdir,"ts_ltla_peak.png"), height = 800, width = 1200, res = 150)
ts_by_wk / ts_by_epiwk
dev.off()

# ---------------------------------------------------------------------------- #

## MAP TOTALS ##

cases %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n),
            la_pop = unique(la_pop),
            rate = sum(n)*1e5/unique(la_pop)) %>% 
  full_join(regions) %>%
  basic_map(fill = "rate") +
  labs(title = "Total cases per 100,000", subtitle = paste0(period[1]," - ",period[2])) -> map_tot

cases2 %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n),
            la_pop = unique(la_pop),
            rate = sum(n)*1e5/unique(la_pop)) %>% 
  full_join(regions) %>%
  basic_map(fill = "rate") +
  labs(subtitle = paste0(period2[1]," - ",period2[2])) -> map_tot2

cases %>%
  group_by(lad19cd) %>%
  summarise(rate = unique(E)*1e5/unique(la_pop)) %>% 
  full_join(regions) %>%
  basic_map(fill = "rate") +
  labs(title = "Expected cases per 100,000", subtitle = "according to age distribution of population") -> map_E

png(here::here(figdir,"map_totals.png"), height = 800, width = 1200, res = 150)
map_tot + map_E
dev.off()

png(here::here(figdir,"map_totals_compwaves.png"), height = 800, width = 1200, res = 150)
# pdf(here::here(figdir,"map_totals_compwaves.pdf"), height = 7, width = 12)
map_tot + map_tot2
dev.off()

## Pillar 1 vs 2

# alldata %>% 
#   group_by(lad19cd, pillar) %>% 
#   summarise(n = sum(n)) -> by_pillar
# 
# by_pillar %>%
#   ggplot(aes(n, fill = pillar)) +
#   geom_histogram(position = "dodge")

alldata %>%
  filter(pillar == "pillar_1",
         dos >= period[1] & dos <= period[2]) %>%
  group_by(lad19cd) %>%
  summarise(n = n()) %>%   #  sum(n, na.rm = T)) %>%
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE) +
  labs(title = "Cases per 100,000", subtitle = "Pillar 1 testing") -> map_pillar1


alldata %>%
  filter(pillar == "pillar_2",
         dos >= period[1] & dos <= period[2]) %>%
  group_by(lad19cd) %>%
  summarise(n = n()) %>%   #  sum(n, na.rm = T)) %>%
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE) +
  labs(subtitle = "Pillar 2 testing") -> map_pillar2

## Deaths

deaths %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n)) %>% 
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE) +
  ggtitle("Deaths per 100,000") -> map_deaths

png(here::here(figdir,"map_deaths_cases.png"), height = 600, width = 1800, res = 150)
map_deaths + map_pillar1 + map_pillar2 +
  plot_annotation(caption = paste0("Data between",period[1]," - ",period[2]))
dev.off()

# ---------------------------------------------------------------------------- #

## MAP EPIDEMIC TIMING ##

deaths %>%
  group_by(lad19cd) %>%
  summarise(first_death = unique(first)-3) -> first_death

dat %>%
  group_by(lad19cd) %>%
  summarise(first_case = unique(first)) %>% 
  full_join(first_death) %>%
  pivot_longer(c("first_case", "first_death")) %>%
  full_join(regions) %>%
  basic_map(fill = "value") +
  scale_fill_viridis(option = "plasma") +
  facet_wrap(~name) +
  ggtitle("Week of first case versus first death - 3wks") -> first_case_death_map

dat %>%
  group_by(lad19cd) %>%
  mutate(first = min(w)) %>% 
  arrange(first) %>%
  full_join(regions) %>%
  basic_map(fill = "first") +
  scale_fill_viridis(option = "plasma") +
  ggtitle("Week of first case") -> first_map

dat %>%
  group_by(lad19cd) %>%
  mutate(peak = w[which.max(n)]) %>% 
  arrange(peak) %>%
  full_join(regions) %>%
  basic_map(fill = "peak") +
  scale_fill_viridis(option = "plasma") +
  ggtitle("Week of peak number of cases") -> peak_map

# png(here::here(figdir,"map_timing.png"), height = 800, width = 1200, res = 150)
# first_map + peak_map
# dev.off()

png(here::here(figdir,"map_timing_cases_deaths.png"), height = 800, width = 1200, res = 150)
first_case_death_map
dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - BY LTLA ##

dat %>%
  ggplot(aes(wor, n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  labs(title = "Cases of COVID-19 in England, by week of reporting",
       x = "Calendar week",
       y = "Cases per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) -> ts_by_wk2

dat %>%
  ggplot(aes(wk_since_first, n*1e5/la_pop, group = lad19cd, col = lad19cd)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  labs(title = "Cases of COVID-19 in England, by week since first LTLA case",
       x = "Calendar week",
       y = "Cases per 100,000"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) -> ts_by_epiwk2

png(here::here(figdir,"ts_ltla.png"), height = 800, width = 1200, res = 150)
ts_by_wk2 / ts_by_epiwk2
dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - BY GEOGRAPHY ##

dat %>%
  group_by(w,wor,geography) %>%
  summarise(n = sum(n, na.rm= T),
            geog_pop = sum(la_pop)) %>%
  ggplot(aes(wor, n*1e5/geog_pop, group = geography, col = geography)) +
  geom_line() +
  geom_point() +
  labs(title = "Cases of COVID-19 in England, by geography type and week of reporting",
       x = "Calendar week",
       y = "Cases per 100,000",
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
  labs(title = "by week since first LTLA case",
       x = "Weeks since first observed case",
       y = "Cases per 100,000",
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




