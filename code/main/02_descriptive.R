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

figdir <- "figures/descriptive"

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
deaths <- readRDS(here::here("data","aggregated","deaths.rds")) 
cases <- readRDS(here::here("data","aggregated","cases.rds")) 

data.list <- list(deaths = deaths,cases = deaths)

regions <- readRDS(here::here("data","LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
regions.df <- sf::st_drop_geometry(regions)

################################################################################
# DESCRIPTIVE SUMMARIES/PLOTS
################################################################################

## MAP TOTALS ##

period <- cases[[3]][[1]]
cases[[1]] %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n, na.rm = TRUE)) %>% 
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE, scale = F) -> case_map
# labs(title = "Confirmed cases per 100,000",
#      subtitle = paste(period[1],"-",period[2])) + 

deaths[[1]] %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n, na.rm = TRUE)) %>% 
  full_join(regions) %>%
  basic_map(fill = "n", rate1e5 = TRUE, scale = F) -> death_map
# labs(title = "Deaths per 100,000") + 

png(here::here(figdir,"map_totals.png"), height = 1000, width = 1500, res = 150)
case_map + death_map
dev.off()

## TIME SERIES - BY GEOGRAPHY ##

period <- deaths$breaks[[1]]

deaths[[1]] %>%
  group_by(w,week,geography) %>%
  summarise(n = sum(n, na.rm= T),
            geog_pop = sum(la_pop)) %>%
  ggplot(aes(week, n*1e5/geog_pop, group = geography, col = geography)) +
  geom_line() +
  labs(subtitle = "COVID19-related deaths in England, by geography and week of death",
       x = "",
       y = "Rate per 100,000", 
       colour = "Geography type") +
  geom_vline(xintercept = ymd("2020-03-23"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-03-22"), y = 15, label = "National lockdown enforced", cex = 2, hjust = "right") +
  scale_x_date(date_minor_breaks = "1 week", 
               date_breaks = "1 month",
               date_labels = "%b",
               limits = period) +
  theme(legend.position = c(0.16,0.60), 
        legend.text=element_text(size=8),  
        legend.title=element_text(size=8)) -> ts_geog_deaths

cases[[1]] %>%
  group_by(w,week,geography) %>%
  summarise(n = sum(n, na.rm= T),
            geog_pop = sum(la_pop)) %>%
  ggplot() +
  geom_line(aes(week, n*1e5/geog_pop, group = geography, col = geography)) +
  geom_vline(xintercept = ymd("2020-03-12"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-03-13"), y = 40, label = "Community testing halted", cex = 2, hjust = "left",) +
  geom_vline(xintercept = ymd("2020-03-23"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-03-24"), y = 45, label = "National lockdown enforced", cex = 2, hjust = "left") +
  geom_vline(xintercept = ymd("2020-04-15"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-04-16"), y = 55, label = "P2 available to care home residents and staff", cex = 2, hjust = "left") +
  geom_vline(xintercept = ymd("2020-05-18"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-05-19"), y = 60, label = "P2 available to all symptomatic cases", cex = 2, hjust = "left") +
  labs(subtitle = "Confirmed COVID-19 cases in England, by geography and week of specimen",
       x = "Calendar week",
       y = "Rate per 100,000") +
  guides(col = "none") +
  scale_x_date(date_minor_breaks = "1 week", 
               date_breaks = "1 month",
               date_labels = "%b",
               limits = period) -> ts_geog_cases


png(here::here(figdir,"fig1.png"), height = 1200, width = 2000, res = 150)
(ts_geog_deaths | death_map) / (ts_geog_cases | case_map ) + plot_layout(widths = c(2,1))
dev.off()

# ---------------------------------------------------------------------------- #

## GEOGRAPHY ##

png(here::here(figdir,"map_geog.png"), height = 800, width = 900, res = 150)
regions %>% basic_map(fill = "geography") + scale_fill_viridis_d()
dev.off()


# ---------------------------------------------------------------------------- #

## COVARIATES ##

cov_names <- c("med_age","pop_dens", "IMD", "prop_minority")
deaths[[1]] %>%
  group_by(geography, lad19cd) %>%
  summarise_at(all_of(cov_names), .funs = base::mean) %>%
  full_join(dplyr::select(regions.df, lad19cd, med_age, prop_male_all)) %>%
  ungroup() -> covs

# Summarise covariates
get_quants <- function(var){ paste(round(quantile(var, p = c(0.25,0.5,0.75)),2), collapse = ", ")}

# By geography
covs %>% 
  group_by(geography) %>%
  summarise(across(c(med_age,IMD,prop_minority, prop_male_all), get_quants))
#   geography                 med_age         IMD                 prop_minority    prop_male_all   
# 1 London Borough            33, 34.5, 36    13.89, 20.4, 26.46  0.31, 0.39, 0.47 0.49, 0.5, 0.5  
# 2 Metropolitan District     35, 39, 41      21.4, 27.2, 30.99   0.04, 0.11, 0.19 0.49, 0.49, 0.5 
# 3 Non-metropolitan District 40, 43, 46      10.78, 13.77, 18.38 0.02, 0.04, 0.07 0.49, 0.49, 0.49
# 4 Unitary Authority         35.75, 39.5, 43 12.95, 19.14, 23.87 0.03, 0.06, 0.14 0.49, 0.5, 0.5 

# Overall
covs %>% 
  summarise(across(c(med_age,IMD,prop_minority, prop_male_all), get_quants))
#   med_age    IMD                 prop_minority    prop_male_all  
# 1 37, 41, 45 11.43, 16.11, 22.44 0.03, 0.05, 0.13 0.49, 0.49, 0.5

regions %>%
  full_join(covs) %>%
  mutate(pop_dens = as.numeric(pop_dens)) -> regions_wcovs

map_dens <-
  basic_map(regions_wcovs, fill = "pop_dens", scale = F) +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill = "", title = "Population density \n(per KM-squared)") +
  theme(plot.title = element_text(size=10))

map_pop <-
  basic_map(regions_wcovs, fill = "la_pop", scale = F) +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill = "", title = "Population size") +
  theme(plot.title = element_text(size=10))

map_imd <-
  basic_map(regions_wcovs, fill = "IMD", scale = F) +
  labs(fill = "", title = "Index of Multiple Deprivation \n(median score)") +
  theme(plot.title = element_text(size=10))

map_mino <-
  basic_map(regions_wcovs, fill = "prop_minority", scale = F) +
  labs(fill = "", title = "Proportion of minority \nethnicities in population") +
  scale_fill_viridis_c(trans = "log10") +
  theme(plot.title = element_text(size=10))

map_age <-
  basic_map(regions_wcovs, fill = "med_age", scale = F) +
  labs(fill = "", title = "Median age") +
  theme(plot.title = element_text(size=10))

map_sex <-
  basic_map(regions_wcovs, fill = "prop_male_all", scale = F) +
  labs(fill = "", title = "Proportion male") +
  theme(plot.title = element_text(size=10))

png(here::here(figdir,"map_covariates.png"), height = 2000, width = 2000, res = 300)
(map_age + map_pop) /
  (map_mino + map_imd)
dev.off()

png(here::here(figdir,"map_sex.png"), height = 1000, width = 1000, res = 300)
map_sex
dev.off()

png(here::here(figdir,"map_covariates3.png"), height = 600, width = 1800, res = 150)
(map_age + map_mino + map_imd)
dev.off()

# ---------------------------------------------------------------------------- #

deaths[[1]] %>% 
  group_by(lad19nm) %>%
  summarise(N = sum(n, na.rm = T),
            pop = mean(la_pop),
            rate = N*1e5/pop) -> death_rates
summary(death_rates$rate)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.34   71.42   88.84   90.56  112.06  196.34 

hist(death_rates$rate, breaks = 40)

cases[[1]] %>% 
  group_by(lad19nm) %>%
  summarise(N = sum(n, na.rm = T),
            pop = mean(la_pop),
            rate = N*1e5/pop) -> case_rates
summary(case_rates$rate)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 71.78  298.29  379.17  403.77  491.51 1039.74 

hist(case_rates$rate, breaks = 40)

################################################################################
################################################################################
