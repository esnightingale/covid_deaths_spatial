################################################################################
# DATA SETUP - CASES
# 
# Call cases linelist from data pipeline, clean and aggregate to LTLA/week.
# 
################################################################################

library(tidyverse)

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## LTLA covariates
covs <- readRDS(paste0(datadir,"covs.rds"))

# --------------------------------------------------------------------------------------------#
## Load/tidy linelist

path_to_factory <- "~/COVID-19/covid19_automation"
file_path <- file.path(path_to_factory, "data", "rds", "linelist_latest.rds")
key <- cyphr::data_key(file.path(path_to_factory, "data"))
linelist_raw <- cyphr::decrypt(readRDS(file_path), key)

# Filter to deaths in England with non-missing LTLA
linelist <- filter(linelist_raw, grepl("e",ltla_code)) %>% # filter to England. death_type28 == 1 indicates official covid-related death type
  mutate(lad19cd = toupper(ltla_code),
         ltla_name = gsub(" And "," and ",
                          gsub(" Of "," of ",
                               stringr::str_to_title(gsub("_"," ",ltla_name)))),
         lad19nm = ltla_name)


# Aggregate Aylsbury Vale, Chiltern, South Bucks and Wycombe into Buckinghamshire unitary authority (created April 2020)
linelist$lad19nm[linelist$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
linelist$lad19cd[linelist$ltla_name %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "E06000060"

# Join city of London into westminster
linelist$lad19cd[linelist$lad19cd == "E09000001"] <- "E09000033"
linelist$lad19nm[linelist$lad19nm == "City of London"] <- "Westminster"

# Exclude Isles of Scilly
linelist <- filter(linelist, lad19nm != "Isles of Scilly")

# Tidy variables and define age group
linelist %>%
  rename(date_report = date_lab_report,
         date = date_specimen) %>%
  mutate_at(vars(date_report, date), lubridate::ymd) %>%
  filter(!is.na(age), !is.na(date)) %>% 
  mutate(result_delay = as.integer(difftime(date_report, date, units = "days")),
         week = lubridate::floor_date(date, unit = "week", week_start = 3),
         age_group = as.character(cut(age, breaks = c(0,seq(10,90,10),120), right = FALSE))) %>%
  group_by(lad19cd) %>%
  mutate(ltla_first = min(date, na.rm = T)) %>%
  ungroup() %>%
  arrange(age) %>%
  dplyr::select(lad19cd, lad19nm, week, date, date_report, ltla_first, result_delay, age_group, pillar) -> alldata

alldata$age_group[alldata$age_group == "[90,120)"] <- "[90,NA)"
alldata$age_group <- as.factor(alldata$age_group)
summary(alldata$age_group)

saveRDS(alldata,paste0(datadir,"linelist_cases.rds"))

# --------------------------------------------------------------------------------------------#

# n_days <- as.integer(max(alldata$dor)-min(alldata$dor)) #200
# n_wks <- as.numeric(difftime(max(alldata$dor),min(alldata$dor), units = "week")) # 28.57
# 
# # Aggregate by age band
# alldata %>%
#   group_by(age_group) %>%
#   count() -> tot_byage
# 
# 
# # LA populations in 5-year age bands, to match shapefile
# pops.long <- readRDS(paste0(datadir,"populations/pops.long.rds")) %>%
#   filter(grepl("E", lad19cd)) %>%
#   group_by(lad19cd, lad19nm, geography, la_pop, mean_age, med_age, age_group10) %>%
#   summarise(la_age_pop = sum(la_age_pop)) %>%
#   rename(age_group = age_group10) %>%
#   filter(lad19cd %in% unique(regions.df$lad19cd)) 
# 
# # Join age-specific totals with populations and calculate rates
# pops.long %>%
#   group_by(age_group) %>%
#   summarise(tot_age_pop = sum(la_age_pop)) %>%
#   full_join(tot_byage) %>%
#   mutate(age_sp_rate_tot = n/tot_age_pop,
#          age_sp_rate_day = age_sp_rate_tot/n_days, # approx. per day
#          age_sp_rate_wk = age_sp_rate_tot/n_wks) %>% # approx. per week 
#   dplyr::select(-n) %>%
#   ungroup() -> tot_byage
# 
# 
# # Apply these rates to age-specific LA populations
# age_pops <- pops.long %>%
#   full_join(tot_byage) %>%
#   arrange(age_group) %>%
#   mutate(la_age_tot_E = la_age_pop*age_sp_rate_tot,
#          la_age_day_E = la_age_pop*age_sp_rate_day,
#          la_age_wk_E = la_age_pop*age_sp_rate_wk)
# 
# la_pops <- age_pops %>%
#   group_by(lad19cd, lad19nm) %>%
#   summarise(geography = unique(geography),
#             la_pop = unique(la_pop),
#             mean_age = unique(mean_age),
#             med_age = unique(med_age),
#             la_tot_E = sum(la_age_tot_E, na.rm = T),
#             la_wk_E = sum(la_age_wk_E, na.rm = T),
#             la_day_E = sum(la_age_day_E, na.rm = T))


# --------------------------------------------------------------------------------------------#
## Analysis data
# 
# source(here::here("code","main","setup_analysis_data_cases.R"))
# 
# # First wave
# setup_analysis_data_cases(alldata, start = ymd("2020-01-01"), end = ymd("2020-06-30"))
# # Second wave
# setup_analysis_data_cases(alldata, start = ymd("2020-07-01"), end = max(alldata$dos, na.rm = T))


# # --------------------------------------------------------------------------------------------#
# 
# ## Aggregate linelist
# 
# # Aggregate by age group, week and LA
# alldata %>%
#   group_by(lad19cd, wor, age_group) %>%
#   count() %>%
#   left_join(age_pops) %>% 
#   ungroup() %>%
#   mutate(SIR_wk = n/la_age_wk_E,
#          SIR_wk_grp = cut(SIR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_age
# 
# # Aggregate by day and LA
# # alldata %>%
# #   group_by(lad19cd, dor) %>%
# #   count() %>%
# #   left_join(la_pops) %>%
# #   ungroup() %>%
# #   mutate(day = as.numeric(dor)- min(as.numeric(dor)),
# #          SIR_day = n/la_day_E) -> d_agg_day
# 
# # Aggregate by week and LA
# alldata %>%
#   group_by(lad19cd, wor) %>%
#   count() %>%
#   left_join(la_pops) %>%
#   ungroup() %>%
#   mutate(w = as.integer(lubridate::week(wor)),
#          SIR_wk = n/la_wk_E,
#          SIR_wk_grp = cut(SIR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_wk
# 
# # Aggregate by LA overall
# alldata %>%
#   group_by(lad19cd) %>%
#   count(name = "n_total") %>%
#   left_join(la_pops) %>%
#   mutate(SIR = n_total/la_tot_E) %>%
#   ungroup() %>%
#   mutate(SIR_grp = cut(SIR, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_tot
# 
# 
# # --------------------------------------------------------------------------------------------#
# 
# ## Analysis data
# 
# dat <- d_agg_wk %>%
#   full_join(d_agg_tot) %>% 
#   left_join(dplyr::select(regions, lad19cd, lad19nm, area_m2)) %>% 
#   left_join(covs) %>% 
#   rename(E = la_tot_E,
#          E_wk = la_wk_E) %>% 
#   group_by(lad19cd) %>%
#   mutate(first = min(w, na.rm = T),
#          date_first = min(wor, na.rm = T)) %>%
#   ungroup() %>%
#   mutate(first_overall = min(w, na.rm = T),
#          date_first_overall = min(wor, na.rm = T)) %>%
#   mutate(n = replace_na(n, 0),
#          geog = as.numeric(as.factor(geography)),
#          wk_since_first = w - first,
#          wk_first_la_overall = first - first_overall,
#          area_km2 = as.numeric(1e-6*area_m2),
#          pop_dens_total = la_pop/area_km2,
#          pop_dens_dev = la_pop/(total_area_km2*prop_dev_land)) %>% # calculate population density per hectare of developed land
#   dplyr::select(n, E, E_wk, w, wor, lad19cd, lad19nm, la_pop, geog, geography, total_area_km2, prop_dev_land, pop_dens_total, pop_dens_dev,
#                 first, wk_since_first, first_overall, wk_first_la_overall, IMD, prop_minority, prop_kw) %>%
#   mutate(w2 = w, w3 = w, SIR = n/E) 
# 
# dat$la <- dat %>%
#   group_by(lad19cd) %>%
#   group_indices()
# 
# saveRDS(dat, paste0(datadir,"Cases/dat.rds"))

