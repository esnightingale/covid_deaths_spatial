################################################################################
# DATA SETUP - PIPELINE
# 
# Call deaths linelist from data pipeline, clean and aggregate to LTLA/week.
# Calculate age-adjusted, expected deaths per LTLA and add LTLA-level covariates
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

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

source(here::here("code","functions.R"))

## Shapefiles

regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds"))  %>%
  st_transform(crs = "+init=epsg:27700 +units=m +no_defs") %>%
  filter(grepl("E",lad19cd)) 
regions$area_m2 = st_area(regions)

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
  arrange(age) %>%
  dplyr::select(lad19cd, ltla_name, wod, dod, dor, report_delay, age_group) -> alldata

alldata$age_group[alldata$age_group == "[90,120)"] <- "[90,NA)"
alldata$age_group <- as.factor(alldata$age_group)
summary(alldata$age_group)


# Include deaths between 2020-03-01 and 2020-08-30
# alldata <- filter(alldata, dod %within% interval(ymd("2020-03-01"),ymd("2020-08-30")))

saveRDS(alldata,paste0(datadir,"Deaths/alldata_alt.rds"))

# --------------------------------------------------------------------------------------------#

## Standardised mortality ratios

# How many days observed in total?
n_days <- as.integer(max(alldata$dod)-min(alldata$dod)) #200
n_wks <- as.numeric(difftime(max(alldata$dod),min(alldata$dod), units = "week")) # 28.57

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

## Aggregate linelist

# Aggregate by age group, week and LA
alldata %>%
  group_by(lad19cd, wod, age_group) %>%
  count() %>%
  left_join(age_pops) %>% 
  ungroup() %>%
  mutate(SMR_wk = n/la_age_wk_E,
         SMR_wk_grp = cut(SMR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_age

# Aggregate by day and LA
alldata %>%
  group_by(lad19cd, dod) %>%
  count() %>%
  left_join(la_pops) %>%
  ungroup() %>%
  mutate(day = as.numeric(dod)- min(as.numeric(dod)),
         SMR_day = n/la_day_E) -> d_agg_day

# Aggregate by week and LA
alldata %>%
  group_by(lad19cd, wod) %>%
  count() %>%
  left_join(la_pops) %>% View()
  ungroup() %>%
  mutate(w = as.integer(lubridate::week(wod)),
         SMR_wk = n/la_wk_E,
         SMR_wk_grp = cut(SMR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_wk


# Aggregate by LA overall
alldata %>%
  group_by(lad19cd) %>%
  count(name = "n_total") %>%
  left_join(la_pops) %>%
  mutate(SMR = n_total/la_tot_E) %>%
  ungroup() %>%
  mutate(SMR_grp = cut(SMR, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_tot


# --------------------------------------------------------------------------------------------#


# Analysis data - total

dat <- d_agg_wk %>%
  full_join(d_agg_tot) %>% 
  left_join(dplyr::select(regions, lad19cd, lad19nm, area_m2)) %>% 
  left_join(covs) %>% 
  rename(E = la_tot_E,
         E_wk = la_wk_E) %>% 
  group_by(lad19cd) %>%
  mutate(first = min(w, na.rm = T),
         date_first = min(wod, na.rm = T)) %>%
  ungroup() %>%
  mutate(first_overall = min(w, na.rm = T),
         date_first_overall = min(wod, na.rm = T)) %>%
  mutate(n = replace_na(n, 0),
         geog = as.numeric(as.factor(geography)),
         wk_since_first = w - first,
         wk_first_la_overall = first - first_overall,
         area_km2 = as.numeric(1e-6*area_m2),
         pop_dens_total = la_pop/area_km2,
         pop_dens_dev = la_pop/(total_area_km2*prop_dev_land)) %>% # calculate population density per hectare of developed land
  dplyr::select(n, E, E_wk, w, wod, lad19cd, lad19nm, la_pop, geog, geography, total_area_km2, prop_dev_land, pop_dens_total, pop_dens_dev,
                first, wk_since_first, first_overall, wk_first_la_overall, IMD, prop_minority, prop_kw) %>%
  mutate(w2 = w, w3 = w, SMR = n/E) 

dat$la <- dat %>%
  group_by(lad19cd) %>%
  group_indices()

saveRDS(dat, paste0(datadir,"Deaths/dat_alt.rds"))
