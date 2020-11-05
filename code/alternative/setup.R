################################################################################
# DATA SETUP - ONS raw
# 
# Read deaths linelist from ONS weekly records, clean and aggregate to LTLA/week.
# Calculate age-adjusted, expected deaths per LTLA and add LTLA-level covariates
# 
################################################################################

library(tidyverse)
library(linelist)
library(INLA)
library(spdep)
library(sf)

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

source(here("code","functions.R"))

## Shapefiles

regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  st_transform(crs = "+init=epsg:27700 +units=m +no_defs") %>%
  filter(grepl("E",lad19cd)) 
regions$area_m2 = st_area(regions)

regions.df <- st_drop_geometry(regions)
head(regions.df)

## Adjacency matrix 
nb <- poly2nb(regions)
# translate to INLA format and save
nb2INLA("regions_eng.adj", nb)

# read in the INLA graph we just created
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

## LTLA covariates
covs <- readRDS(paste0(datadir,"covs.rds"))

## Death data
files <- list("deaths_wk11_19.xlsx", "SPI_M_Week_20_Final.xlsx", "SPI_M_Week_21_Final.xlsx","ONS deaths Week 22 Final.xlsx")

read_deaths <- function(file){
  
  path <- paste0(datadir,"Deaths/",file)
  
  ment <- grep("mention",tolower(readxl::excel_sheets(path)))
  und <- grep("underlying",tolower(readxl::excel_sheets(path)))
  
  deaths <- bind_rows(mutate(readxl::read_xlsx(path, sheet = ment), type = "mentioned"),
                      mutate(readxl::read_xlsx(path, sheet = und), type = "underlying")) 
  names(deaths) <- tolower(names(deaths))
  names(deaths)[grep("age", tolower(names(deaths)))] <- "agegroup"
  
  return(deaths)
}

covid_deaths_raw <- bind_rows(lapply(files, read_deaths)) 

# Add earlier deaths
early <- read_deaths("ons_deaths.xlsx") %>%
  filter(lubridate::week(lubridate::ymd(dod)) < 11) 
covid_deaths_raw <- bind_rows(early, covid_deaths_raw)  

# Filter to deaths in England
covid_deaths_E <- filter(covid_deaths_raw, grepl("E",laua)) 

saveRDS(covid_deaths_raw,paste0(datadir,"Deaths/covid_deaths_raw.rds"))
saveRDS(covid_deaths_E,paste0(datadir,"Deaths/covid_deaths_E.rds"))

# Filter to those with non-missing LTLA #& type == "underlying"
covid_deaths <- filter(covid_deaths_E,laua != "NULL")

# Aggregate Aylsbury Vale, Chiltern, South Bucks and Wycombe into Buckinghamshire unitary authority (created April 2020)
covid_deaths$lad19nm[covid_deaths$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
covid_deaths$laua[covid_deaths$lad19nm == "Buckinghamshire"] <- "E06000060"

# Join city of London into westminster
covid_deaths$laua[covid_deaths$laua %in% c("E09000001")] <- "E09000033"
covid_deaths$lad19nm[covid_deaths$laua == "E09000033"] <- "Westminster"

# covid_deaths$Age[covid_deaths$Age == "<1"] <- "0 to 1"
covid_deaths$agegroup[covid_deaths$agegroup == "<1"] <- "0 to 1"

covid_deaths %>%
  rename(lad19cd = laua) %>%
  mutate_at(vars(dod, dor), lubridate::ymd) %>%
  mutate(report_delay = as.integer(difftime(dor, dod, units = "days")),
         wod = lubridate::floor_date(dod, unit = "week")) %>%
  separate(agegroup, into = c("agemin","agemax"), sep = " to ", remove = F, convert = T) %>%
  mutate(agemin = as.numeric(substr(agemin,1,2)),
         age_group = paste0("[",agemin,",",agemax,"]")) %>%
  arrange(agemin) -> alldata

alldata$age_group[alldata$agemin >= 90] <- "[90,NA)"
alldata$age_group[alldata$agemin %in% c(0,1)] <- "[0,5)"

# Add decade age groups
alldata <- alldata %>%
  mutate(age_group10 = case_when(agemin %in% c(90,95) ~ "[90,NA)",
                                 agemin %in% c(80,85) ~ "[80,90)",
                                 agemin %in% c(70,75) ~ "[70,80)",
                                 agemin %in% c(60,65) ~ "[60,70)",
                                 agemin %in% c(50,55) ~ "[50,60)",
                                 agemin %in% c(40,45) ~ "[40,50)",
                                 agemin %in% c(30,35) ~ "[30,40)",
                                 agemin %in% c(20,25) ~ "[20,30)",
                                 agemin %in% c(10,15) ~ "[10,20)",
                                 agemin %in% c(0,1,5) ~ "[0,10)")) %>%
  mutate_at(vars("age_group","age_group10"), as.factor)

# Redefine levels for place of death
summary(as.factor(alldata$pod))
#                    Care home                    Elsewhere                         Home                      Hospice                     Hospital Other communal establishment Other Communal Establishment 
#                         5885                          109                         1307                          300                        19590                           85                           20 

alldata <- 
  alldata %>%
  mutate(pod2 = case_when(pod %in% c("Hospital", "Care home","Home") ~ pod,
                          !pod %in% c("Hospital", "Care home","Home") ~ "Elsewhere"))

alldata$pod2 <- factor(alldata$pod2, levels = c("Hospital", "Care home","Home","Elsewhere"))
summary(alldata$pod2)
 # Hospital Care home      Home Elsewhere 
 #    19590      5885      1307       514 
 

## Define which age group to use

alldata <- alldata %>%
  dplyr::select(-age_group) %>%
  rename(age_group = age_group10) %>%
  ungroup()

saveRDS(alldata,paste0(datadir,"Deaths/alldata.rds"))


## Standardised mortality ratios

# How many days observed in total?
n_days <- as.integer(max(alldata$dod)-min(alldata$dod)) #53
n_wks <- as.numeric(difftime(max(alldata$dod),min(alldata$dod), units = "week")) # 7.57

# Aggregate by age band
alldata %>%
  group_by(age_group) %>%
  summarise(n = sum(deaths, na.rm = T)) -> tot_byage


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


# Aggregate by age group, week and LA
alldata %>%
  group_by(lad19cd, lad19nm, wod, age_group) %>%
  summarise(n = sum(deaths, na.rm = T)) %>%
  left_join(age_pops) %>% 
  ungroup() %>%
  mutate(SMR_wk = n/la_age_wk_E,
         SMR_wk_grp = cut(SMR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_age

# Aggregate by day and LA
alldata %>%
  group_by(lad19cd, lad19nm, dod) %>%
  summarise(n = sum(deaths, na.rm = T)) %>%
  left_join(la_pops) %>%
  ungroup() %>%
  mutate(day = as.numeric(dod)- min(as.numeric(dod)),
         SMR_day = n/la_day_E) -> d_agg_day

# Aggregate by week and LA
alldata %>%
  group_by(lad19cd, lad19nm, wod) %>%
  summarise(n = sum(deaths, na.rm = T)) %>%
  left_join(la_pops) %>%
  ungroup() %>%
  mutate(w = as.integer(lubridate::week(wod)),
         SMR_wk = n/la_wk_E,
         SMR_wk_grp = cut(SMR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_wk


# Aggregate by LA overall
alldata %>%
  group_by(lad19cd, lad19nm) %>%
  summarise(n_total = sum(deaths, na.rm = T)) %>%
  left_join(la_pops) %>%
  mutate(SMR = n_total/la_tot_E) %>%
  ungroup() %>%
  mutate(SMR_grp = cut(SMR, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_tot


# Analysis data - total

dat <- d_agg_wk %>%
  full_join(d_agg_tot) %>% 
  left_join(covs) %>%
  # left_join(imd) %>% 
  # left_join(ethn) %>%
  # left_join(kw) %>%
  # left_join(landuse) %>%
  left_join(dplyr::select(regions, lad19cd, area_m2)) %>% 
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

saveRDS(dat, paste0(datadir,"Deaths/dat.rds"))

