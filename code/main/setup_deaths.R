################################################################################
# DATA SETUP - PIPELINE
# 
# Call deaths linelist from data pipeline and clean - defining 10-yr age groups
# and LTLA names to match with shapefiles and pops.
# 
################################################################################

library(reportfactory)
library(tidyverse)
library(lubridate)
library(linelist)
library(cyphr)
library(INLA)
library(rgdal)
library(spdep)
library(sf)
library(here)

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

source(here::here("code","functions.R"))

## Shapefiles

regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds"))  %>%
  st_transform(crs = "+init=epsg:27700 +units=m +no_defs") %>%
  filter(grepl("E",lad19cd)) 
# regions$area_m2 = st_area(regions)

regions.df <- st_drop_geometry(regions)
head(regions.df)

## Adjacency matrix 
# nb <- poly2nb(regions)
# translate to INLA format and save
# nb2INLA("regions_eng.adj", nb)

# read in the INLA graph we just created
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

## LTLA covariates
covs <- readRDS(paste0(datadir,"covs.rds"))

# --------------------------------------------------------------------------------------------#
## Load/tidy linelist

# READ ONS LINELIST
path_to_factory <- "~/COVID-19/covid19_automation"
file_path <- file.path(path_to_factory, "data", "rds", "deaths_eng_latest.rds")
key <- cyphr::data_key(file.path(path_to_factory, "data"))
covid_deaths_raw <- cyphr::decrypt(readRDS(file_path), key)

# Filter to deaths in England with non-missing LTLA
covid_deaths <- filter(covid_deaths_raw, grepl("e",ltla_code)) %>% # filter to England. death_type28 == 1 indicates official covid-related death type
  mutate(lad19cd = toupper(ltla_code),
         ltla_name = gsub(" And "," and ",
                        gsub(" Of "," of ",
                        stringr::str_to_title(gsub("_"," ",ltla_name)))),
         lad19nm = ltla_name)

# saveRDS(covid_deaths_raw,paste0(death_dir,"covid_deaths_raw.rds"))
# saveRDS(covid_deaths_E,paste0(death_dir,"covid_deaths_E.rds"))

# Aggregate Aylsbury Vale, Chiltern, South Bucks and Wycombe into Buckinghamshire unitary authority (created April 2020)
covid_deaths$lad19nm[covid_deaths$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
covid_deaths$lad19cd[covid_deaths$ltla_name %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "E06000060"

# Join city of London into westminster
covid_deaths$lad19cd[covid_deaths$lad19cd == "E09000001"] <- "E09000033"
covid_deaths$lad19nm[covid_deaths$lad19cd == "E09000033"] <- "Westminster"

# Tidy variables and define age group
covid_deaths %>%
  rename(dod = date_death,
         dor = ons_date_report) %>%
  mutate_at(vars(dod, dor), lubridate::ymd) %>%
  filter(!is.na(age)) %>% #n = 9
  mutate(report_delay = as.integer(difftime(dor, dod, units = "days")),
         wod = lubridate::floor_date(dod, unit = "week"),
         age_group = as.character(cut(age, breaks = c(0,seq(10,90,10),120), right = FALSE))) %>%
  group_by(lad19cd) %>%
  mutate(ltla_first = min(dod, na.rm = T)) %>%
  ungroup() %>%
  arrange(age) %>%
  dplyr::select(lad19cd, ltla_name, ltla_first, wod, dod, dor, report_delay, age_group) -> alldata

alldata$age_group[alldata$age_group == "[90,120)"] <- "[90,NA)"
alldata$age_group <- as.factor(alldata$age_group)
summary(alldata$age_group)

saveRDS(alldata,paste0(datadir,"Deaths/alldata_alt.rds"))

# --------------------------------------------------------------------------------------------#

## Calculate standardised mortality ratios

# How many days observed in total?
n_days <- as.integer(max(alldata$dod, na.rm = T)-min(alldata$dod, na.rm = T))
n_wks <- as.numeric(difftime(max(alldata$dod, na.rm = T),min(alldata$dod, na.rm = T), units = "week"))

# Aggregate by age band
alldata %>%
  group_by(age_group) %>%
  count() -> tot_byage

# LA populations in 5-year age bands, to match shapefile
pops.long <- readRDS(paste0(datadir,"maps/pops.long.rds")) %>%
  filter(grepl("E", lad19cd)) %>%
  group_by(lad19cd, lad19nm, geography, la_pop, mean_age, med_age, age_group10) %>%
  summarise(la_age_pop = sum(la_age_pop)) %>%
  rename(age_group = age_group10) %>%
  filter(lad19cd %in% unique(regions.df$lad19cd)) 

# Join age-specific totals with populations and calculate rates
pops.long %>%
  group_by(age_group) %>%
  summarise(tot_age_pop = sum(la_age_pop)) %>%
  full_join(tot_byage) %>%
  mutate(age_sp_rate_tot = n/tot_age_pop,
         age_sp_rate_day = age_sp_rate_tot/n_days, # approx. per day
         age_sp_rate_wk = age_sp_rate_tot/n_wks) %>% # approx. per week 
  dplyr::select(-n) %>%
  ungroup() -> tot_byage


# Apply these rates to age-specific LA populations
age_pops <- pops.long %>%
  full_join(tot_byage) %>%
  arrange(age_group) %>%
  mutate(la_age_tot_E = la_age_pop*age_sp_rate_tot,
         la_age_day_E = la_age_pop*age_sp_rate_day,
         la_age_wk_E = la_age_pop*age_sp_rate_wk)

la_pops <- age_pops %>%
  group_by(lad19cd, lad19nm) %>%
  summarise(geography = unique(geography),
            la_pop = unique(la_pop),
            mean_age = unique(mean_age),
            med_age = unique(med_age),
            la_tot_E = sum(la_age_tot_E, na.rm = T),
            la_wk_E = sum(la_age_wk_E, na.rm = T),
            la_day_E = sum(la_age_day_E, na.rm = T))

# --------------------------------------------------------------------------------------------#
## Analysis data

source(here::here("code","setup_analysis_data.R"))

# First wave
setup_analysis_data(alldata, start = ymd("2020-01-01"), end = ymd("2020-06-30"))
# Second wave
setup_analysis_data(alldata, start = ymd("2020-08-01"), end = ymd("2020-10-31"))

