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
                                start = lubridate::ymd("2020-01-01"), 
                                end = lubridate::ymd("2020-06-30")){

# start = ymd("2020-01-01")
# end = ymd("2020-06-30")

week_seq <- lubridate::ymd(seq(start, end, by = "week"))
week_seq
  
alldata_sub <- alldata %>%
  dplyr::filter(between(date, start, end))
  
# --------------------------------------------------------------------------------------------#

## Calculate age-stratified expected counts, according to rates over specified 
## time period

la_pops <- calc_E(alldata_sub)

# --------------------------------------------------------------------------------------------#

# Set up data with LTLA chars for each week in period
ltla_first <- alldata_sub %>%
  dplyr::select(lad19cd, ltla_first) %>%
  unique()

expand <- data.frame(week = rep(week_seq, dplyr::n_distinct(alldata_sub$lad19cd)),
                     lad19cd = rep(unique(alldata_sub$lad19cd), each = length(week_seq))) %>%
  dplyr::full_join(ltla_first) %>%
  dplyr::mutate(w = as.integer(lubridate::week(week))) %>%
  dplyr::left_join(la_pops)

## Aggregate linelist by week and LA
alldata_sub %>%
  dplyr::group_by(lad19cd, week) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::full_join(expand) %>%
  dplyr::full_join(ltla_first) %>%
  dplyr::left_join(la_pops) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0),
         w = as.integer(lubridate::week(week))) %>%
  dplyr::arrange(lad19cd, week) -> d_agg_wk

# Aggregate by LA overall
alldata_sub %>%
  dplyr::group_by(lad19cd, ltla_first) %>%
  dplyr::count(name = "n_total") %>%
  dplyr::full_join(la_pops) %>%
  dplyr::mutate(n_total = tidyr::replace_na(n_total, 0)) %>%
  dplyr::ungroup() -> d_agg_tot


# --------------------------------------------------------------------------------------------#

# Final analysis dataset - merged with total incidence, covariates, populations and shapefile areas

dat <- d_agg_wk %>%
  dplyr::full_join(d_agg_tot) %>% 
  dplyr::full_join(dplyr::select(regions, lad19cd, lad19nm, area_km2)) %>% 
  dplyr::left_join(covs) %>% 
  dplyr::rename(E = la_tot_E,
         E_wk = la_wk_E,
         E_wk_unstrat = la_unstrat_wk_E) %>%
  dplyr::mutate(first = lubridate::floor_date(ltla_first, unit = "week", week_start = 3),
         first_overall = min(first),
         date_first_overall = min(ltla_first, na.rm = T)) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0),
         geog = as.numeric(as.factor(geography)),
         wk_since_first = as.numeric(lubridate::interval(first,week), "weeks"),
         wk_first_la_overall = as.numeric(lubridate::interval(first_overall,first), "weeks"),
         pop_dens = la_pop/area_km2,
         IMD_quint = cut(IMD, 5)) %>%  
  dplyr::select(n, E, E_wk, E_wk_unstrat, w, week, lad19cd, lad19nm, la_pop, geog, geography, area_km2, pop_dens,
                first, wk_since_first, first_overall, wk_first_la_overall, med_age, IMD, IMD_quint, prop_minority) %>%
  dplyr::mutate(w2 = w, w3 = w, SIR = n/E) 

# add numeric indices for each LTLA
dat$la <- dat %>%
  dplyr::group_by(lad19cd) %>%
  dplyr::group_indices()

# saveRDS(dat, here::here("data","aggregated",sprintf("/%s_%s_%s.rds",measure,start,end)))

return(dat)

}
