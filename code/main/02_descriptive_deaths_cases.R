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
library(lubridate)
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
library(gganimate)
library(data.table)

# source(here::here("code","main","functions.R"))
# source(here::here("code","main","plot_functions.R"))

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

linelist_deaths <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","deaths")))%>%
  mutate(month = lubridate::month(date, label = TRUE)) 
linelist_cases <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","cases"))) %>%
  mutate(month = lubridate::month(date, label = TRUE)) 


theme_set(theme_minimal())


# ---------------------------------------------------------------------------- #

## GEOGRAPHY ##
png(here::here("figures","descriptive","map_geog.png"), height = 800, width = 900, res = 150)
regions %>% basic_map(fill = "geography") + scale_fill_viridis_d()
dev.off()

# ---------------------------------------------------------------------------- #

## COVARIATES ##

# Median age

png(here::here("figures","descriptive","map_age.png"), height = 800, width = 900, res = 150)
print(
regions %>%
  basic_map(fill = "med_age") +
  scale_fill_viridis_c() +
  labs(fill = "", title = "Median age")
)
dev.off()

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

## TIME SERIES - total ##

linelist_cases %>%
  group_by(lad19cd, week, age_group) %>%
  tally() %>% 
  ggplot(aes(week, n, fill = age_group)) +
  geom_col() +
  labs(title = "Total test-confirmed cases in England",
       x = "",
       y = "",
       fill = "Age group"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  scale_fill_viridis_d(option = "plasma") -> ts_cases

linelist_deaths %>%
  group_by(lad19cd, week, age_group) %>%
  tally() %>% 
  ggplot(aes(week, n, fill = age_group)) +
  geom_col() +
  labs(title = "Total COVID-19-related deaths in England",
       x = "",
       y = "",
       fill = "Age group"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  scale_fill_viridis_d(option = "plasma") -> ts_deaths

png(here::here(figdir,"total_ts.png"), height = 1000, width = 1200, res = 150)
ts_cases / ts_deaths
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

## MAP OVER TIME ##

png(here::here(figdir,"map_month_cases.png"), height = 3000, width = 3000, res = 150)
plot_la_mth(linelist_cases,title = "Cases per 100,000")
dev.off()

png(here::here(figdir,"map_month_deaths.png"), height = 3000, width = 3000, res = 150)
plot_la_mth(linelist_deaths,title = "Deaths per 100,000")
dev.off()

linelist_cases %>%
  filter(date < ymd("2020-07-01")) %>% 
  plot_la_mth_animate(title = "Cases per 100,000",
                      file = "map_month_cases_1.gif",
                      path = here::here(figdir))

linelist_deaths %>%
  filter(date < ymd("2020-07-01")) %>%
  plot_la_mth_animate(title = "Deaths per 100,000",
                      file = "map_month_deaths_1.gif",
                      path = here::here(figdir))


linelist_cases %>%
  filter(date >= ymd("2020-07-01")) %>%
  plot_la_mth_animate(title = "Cases per 100,000",
                      file = "map_month_cases_2.gif",
                      path = here::here(figdir))

linelist_deaths %>%
  filter(date >= ymd("2020-07-01")) %>%
  plot_la_mth_animate(title = "Deaths per 100,000",
                      file = "map_month_deaths_2.gif",
                      path = here::here(figdir))

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

## Overall swab-death lag ##

summary(linelist_deaths$swab_death)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -180.00    3.00    6.00    8.24   11.00  219.00   33430 

summary(linelist_deaths$swab_death[linelist_deaths$swab_death >= 0 & linelist_deaths$swab_death <= 100])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    3.00    6.00    8.22   11.00  100.00   33430 
   
png(here::here("figures","descriptive","swab_death_lag.png"), height = 800, width = 1000, res = 150)
linelist_deaths %>% 
  filter(swab_death >= 0) %>%
  ggplot(aes(x = swab_death)) + 
  geom_histogram(bins = 40, fill = "steelblue") +
  labs(x = "Days from swab to death", y = "Density") + 
  xlim(c(0,100))
dev.off()

# ---------------------------------------------------------------------------- #

## Compare case-death lag per LTLA ##

pal <- RColorBrewer::brewer.pal(3,"RdYlBu")
border <- st_union(regions)
# 
# merged[[1]] %>%
#   group_by(lad19cd, la_pop) %>%
#   # mutate(n_c_roll = rollmean(n_c, 2),
#   #        n_d_roll = rollmean(n_d, 2)) %>%
#   summarise(peak_c = week[which.max(n_c/la_pop)],
#             peak_d = week[which.max(n_d/la_pop)],
#             lag_c_d = peak_d - peak_c) %>%
#   right_join(regions) %>%
#   basic_map(fill = "lag_c_d", rate1e5 = FALSE) +
#   scale_fill_gradient2(low = pal[1], mid = pal[2], high = pal[3]) +
#   labs(title = "Lag between peak deaths and peak cases per 100,000",
#        subtitle = paste(merged$breaks[[1]][1],"-",merged$breaks[[1]][2])) -> map_lag_rate_first
# 
# 
# merged[[2]] %>%
#   group_by(lad19cd, la_pop) %>%
#   # mutate(n_c_roll = rollmean(n_c, 2),
#   #        n_d_roll = rollmean(n_d, 2)) %>%
#   summarise(peak_c = week[which.max(n_c/la_pop)],
#             peak_d = week[which.max(n_d/la_pop)],
#             lag_c_d = peak_d - peak_c) %>%
#   right_join(regions) %>%
#   basic_map(fill = "lag_c_d", rate1e5 = FALSE) +
#   scale_fill_gradient2(low = pal[1], mid = pal[2], high = pal[3]) +
#   labs(title = "Lag between peak deaths and peak cases per 100,000",
#        subtitle = paste(merged$breaks[[2]][1],"-",merged$breaks[[2]][2])) -> map_lag_rate_second
# 
# png(here::here(figdir,"map_case_death_lags.png"), height = 1000, width = 1500, res = 150)
# map_lag_rate_first + map_lag_rate_second
# dev.off()


## By rolling week avg from linelist

linelist_deaths %>%
  group_by(lad19cd, ltla_name, date) %>%
  tally() -> daily_deaths

date_seq <- seq(from = min(daily_deaths$date), to = max(daily_deaths$date), by = "days")

vars <- names(select(daily_deaths, lad19cd, ltla_name))
daily_deaths <- as.data.table(daily_deaths)
all_dates <- daily_deaths[,.(date=date_seq),by = vars]

# Merge and fill count with 0:
setkey(daily_deaths, lad19cd, ltla_name, date)
setkey(all_dates, lad19cd, ltla_name, date)
daily_deaths <- daily_deaths[all_dates,roll=TRUE]


daily_deaths %>%
  as.data.frame() %>%
  mutate(n = replace_na(n, 0)) %>%
  right_join(dplyr::select(regions.df, lad19cd, la_pop)) %>%
  group_by(lad19cd, ltla_name, la_pop) %>%
  mutate(roll_mean = zoo::rollmean(n*1e5/la_pop, k = 14, align = "right", fill = NA)) -> daily_deaths

# daily_deaths %>%
# ggplot(aes(date, roll_mean, group = lad19cd)) +
#   geom_line(alpha = 0.2)


# --- #
 
linelist_cases %>%
  group_by(lad19cd, lad19nm, date) %>%
  tally() -> daily_cases

date_seq <- seq(from = min(daily_cases$date), to = max(daily_cases$date), by = "days")


vars <- names(select(daily_cases, lad19cd, lad19nm))
daily_cases <- as.data.table(daily_cases)
all_dates <- daily_cases[,.(date=date_seq),by = vars]

# Merge and fill count with 0:
setkey(daily_cases, lad19cd, lad19nm, date)
setkey(all_dates, lad19cd, lad19nm, date)
daily_cases <- daily_cases[all_dates,roll=TRUE]


daily_cases %>%
  as.data.frame() %>%
  mutate(n = replace_na(n, 0)) %>%
  right_join(dplyr::select(regions.df, lad19cd, la_pop)) %>%
  group_by(lad19cd, lad19nm, la_pop) %>%
  mutate(roll_mean = zoo::rollmean(n*1e5/la_pop, k = 14, align = "right", fill = NA)) -> daily_cases

# daily_cases %>%
#   ggplot(aes(date, roll_mean, group = lad19cd)) +
#   geom_line(alpha = 0.2)

daily_all <- daily_cases %>%
  full_join(daily_deaths, by = c("lad19cd","la_pop","date"), suffix = c("_d","_c")) %>%
  mutate(wave = as.numeric(date < ymd("2020-07-01")) + 1)

daily_all %>%
  group_by(wave) %>%
  group_split() -> daily_all_split

# Check cross-correlation
ccf(daily_all_split[[1]]$n_c, daily_all_split[[1]]$n_c)
ccf(daily_all_split[[2]]$n_c, daily_all_split[[2]]$n_c)


# Join deaths and cases
daily_deaths %>%
  full_join(daily_cases, by = c("lad19cd","la_pop","date"), suffix = c("_d","_c")) %>%
  filter(date < ymd("2020-07-01")) %>%
  group_by(lad19cd) %>%
  summarise(peak_d = date[which.max(n_d)],
            peak_c = date[which.max(n_c)],
            lag_c_d = peak_d - peak_c) -> daily_all_first
summary(as.numeric(daily_all_first$lag_c_d))

daily_deaths %>%
  full_join(daily_cases, by = c("lad19cd","la_pop","date"), suffix = c("_d","_c")) %>%
  filter(date >= ymd("2020-07-01")) %>%
  group_by(lad19cd) %>%
  summarise(peak_d = date[which.max(n_d)],
            peak_c = date[which.max(n_c)],
            lag_c_d = peak_d - peak_c) -> daily_all_second
summary(as.numeric(daily_all_second$lag_c_d))

samp <- sample(daily_all$lad19cd, 1)
daily_all %>%
  filter(lad19cd == samp) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = roll_mean_d), col = "red") +
  geom_line(aes(y = roll_mean_c), col = "blue") +
  labs(title = daily_all$lad19nm[daily_all$lad19cd == samp], 
       subtitle = "Deaths red, cases blue")



border <- st_union(regions)

regions %>%
  full_join(daily_all_first) %>%
  basic_map(fill = "lag_c_d") +
  geom_sf(data = border, aes(geometry = geometry), fill = NA) +
  scale_fill_gradient2() +
  labs(title = "Lag between peak deaths and peak cases per 100,000",
       subtitle = paste(merged$breaks[[1]][1],"-",merged$breaks[[1]][2])) -> lag_first

regions %>%
  full_join(daily_all_second) %>%
  basic_map(fill = "lag_c_d") +
  geom_sf(data = border, aes(geometry = geometry), fill = NA) +
scale_fill_gradient2() +
  labs(title = "Lag between peak deaths and peak cases per 100,000",
       subtitle = paste(merged$breaks[[2]][1],"-",merged$breaks[[2]][2])) -> lag_second

png(here::here(figdir,"map_case_death_lags.png"), height = 1000, width = 1500, res = 150)
lag_first + lag_second
dev.off()

