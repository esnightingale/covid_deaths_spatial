################################################################################
# ANALYSIS DATA SETUP
# 
# Calculate age-adjusted, expected deaths per LTLA and add LTLA-level covariates
# 
################################################################################

# Function which takes alldata object from setup_deaths.R and a start/end date
# for analysis

setup_analysis_data <- function(alldata, 
                                measure,
                                start = ymd("2020-03-01"), 
                                end = ymd("2020-08-01")){

  
alldata_sub <- alldata %>%
    filter(date >= start & date <= end)

# --------------------------------------------------------------------------------------------#

## Calculate age-stratified expected counts

la_pops <- calc_E(alldata_sub)

# --------------------------------------------------------------------------------------------#

## Aggregate linelist by week and LA

alldata_sub %>%
  group_by(lad19cd, week, ltla_first) %>%
  count() %>%
  ungroup() %>%
  full_join(la_pops) %>% 
  mutate(n = replace_na(n, 0),
         w = as.integer(lubridate::week(week)),
         SIR_wk = n/la_wk_E,
         SIR_wk_grp = cut(SIR_wk, breaks = 5, include.lowest = T, ordered_result = T)) -> d_agg_wk

# Aggregate by LA overall
alldata_sub %>%
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
  dplyr::select(n, E, E_wk, w, week, lad19cd, lad19nm, la_pop, geog, geography, area_km2, pop_dens,
                first, wk_since_first, first_overall, wk_first_la_overall, IMD, prop_minority, prop_kw) %>%
  mutate(w2 = w, w3 = w, SIR = n/E) 

# add numeric indices for each LTLA
dat$la <- dat %>%
  group_by(lad19cd) %>%
  group_indices()

saveRDS(dat, paste0(datadir,sprintf("/by_week_la/dat_%s_%s_%s.rds",measure,start,end)))

return(dat)

}
