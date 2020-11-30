################################################################################
# ANALYSIS DATA SETUP
# 
# Calculate age-adjusted, expected deaths per LTLA and add LTLA-level covariates
# 
################################################################################

# Function which takes alldata object from setup_deaths.R and a start/end date
# for analysis

setup_analysis_data <- function(alldata, 
                                start = ymd("2020-03-01"), 
                                end = ymd("2020-08-01")){

  
alldata <- alldata %>%
    filter(dod >= start & dod <= end)

# --------------------------------------------------------------------------------------------#

## Aggregate linelist

# Aggregate by age group, week and LA
# alldata %>%
#   group_by(lad19cd, wod, age_group) %>%
#   count() %>%
#   left_join(age_pops) %>% 
#   ungroup() %>%
#   mutate(SMR_wk = n/la_age_wk_E,
#          SMR_wk_grp = cut(SMR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_age

# Aggregate by day and LA
# alldata %>%
#   group_by(lad19cd, dod) %>%
#   count() %>%
#   ungroup() %>%
#   mutate(day = as.numeric(dod)- min(as.numeric(dod))) -> d_agg_day

# Aggregate by week and LA
alldata %>%
  group_by(lad19cd, wod, ltla_first) %>%
  count() %>%
  ungroup() %>%
  full_join(la_pops) %>%
  mutate(n = replace_na(n, 0),
         w = as.integer(lubridate::week(wod)),
         SIR_wk = n/la_wk_E,
         SIR_wk_grp = cut(SIR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_wk

# Aggregate by LA overall
alldata %>%
  group_by(lad19cd, ltla_first) %>%
  count(name = "n_total") %>%
  full_join(la_pops) %>%
  mutate(n_total = replace_na(n_total, 0),
         SIR = n_total/la_tot_E) %>%
  ungroup() %>%
  mutate(SIR_grp = cut(SIR, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_tot


# --------------------------------------------------------------------------------------------#

# Final analysis dataset - merged with total incidence, covariates, populations and shapefile areas

dat <- d_agg_wk %>%
  full_join(d_agg_tot) %>% 
  full_join(dplyr::select(regions, lad19cd, lad19nm, area_km2)) %>% 
  left_join(covs) %>% 
  rename(E = la_tot_E,
         E_wk = la_wk_E) %>%
  mutate(first = lubridate::week(lubridate::floor_date(ltla_first, unit = "week")),
         first_overall = min(first),
         date_first_overall = min(ltla_first, na.rm = T)) %>%
  mutate(n = replace_na(n, 0),
         geog = as.numeric(as.factor(geography)),
         wk_since_first = w - first,
         wk_first_la_overall = first - first_overall,
         pop_dens = la_pop/area_km2) %>% 
  dplyr::select(n, E, E_wk, w, wod, lad19cd, lad19nm, la_pop, geog, geography, area_km2, pop_dens,
                first, wk_since_first, first_overall, wk_first_la_overall, IMD, prop_minority, prop_kw) %>%
  mutate(w2 = w, w3 = w, SMR = n/E) 

# add numeric indices for each LTLA
dat$la <- dat %>%
  group_by(lad19cd) %>%
  group_indices()

saveRDS(dat, paste0(datadir,sprintf("Deaths/dat_deaths_%s_%s.rds",start,end)))

}
