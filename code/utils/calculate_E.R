# --------------------------------------------------------------------------------------------#

## Calculate expected counts per LTLA given population size and total age-specific rates

calc_E <- function(alldata){

# How many days observed in total?
n_days <- as.integer(max(alldata$date, na.rm = T)-min(alldata$date, na.rm = T))
n_wks <- as.numeric(difftime(max(alldata$date, na.rm = T),min(alldata$date, na.rm = T), units = "week"))

# Aggregate by age band
alldata %>%
  group_by(age_group) %>%
  count() -> tot_byage

# LA populations in 10-year age bands, to match shapefile
pops.long <- readRDS(here::here("data","pops.long.rds")) %>%
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
  ungroup() %>%
  mutate(tot_rate = sum(n)/sum(tot_age_pop),
         tot_rate_wk = tot_rate/n_wks)  %>%
  dplyr::select(-n) -> tot_byage


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
            la_unstrat_wk_E = la_pop*unique(tot_rate_wk)) %>%
            # "la_tot_E_{measure}" := sum(la_age_tot_E, na.rm = T),
            # la_day_E = sum(la_age_day_E, na.rm = T),
            # "la_wk_E_{measure}" := sum(la_age_wk_E, na.rm = T))
   ungroup()

return(la_pops)

}
