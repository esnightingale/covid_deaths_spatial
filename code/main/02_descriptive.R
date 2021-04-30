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
merged <- readRDS(here::here("data","merged.rds")) 

deaths <- readRDS(here::here("data","deaths.rds")) 
cases <- readRDS(here::here("data","cases.rds")) 

data.list <- list(deaths = deaths,cases = deaths)

linelist_deaths <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","deaths")))%>%
  mutate(month = lubridate::month(date, label = TRUE)) 
linelist_cases <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","cases"))) %>%
  mutate(month = lubridate::month(date, label = TRUE)) 

################################################################################
# DESCRIPTIVE SUMMARIES/PLOTS
################################################################################

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
  guides(col = FALSE) +
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

cov_names <- c("med_age","pop_dens", "IMD", "prop_minority", "prop_kw")
deaths[[1]] %>%
  group_by(geography, lad19cd) %>%
  summarise_at(all_of(cov_names), .funs = base::mean) %>%
  full_join(dplyr::select(regions.df, lad19cd, med_age)) %>%
  ungroup() -> covs

# Summarise covariates
get_quants <- function(var){ paste(round(quantile(var, p = c(0.25,0.5,0.75)),2), collapse = ", ")}

# By geography
covs %>% 
  group_by(geography) %>%
  summarise(across(c(med_age,IMD,prop_minority), get_quants))

# Overall
covs %>% 
  summarise(across(c(med_age,IMD,prop_minority), get_quants))

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
  labs(fill = "", title = "Proportion of black and \nminority ethnic population") +
  scale_fill_viridis_c(trans = "log10") +
  theme(plot.title = element_text(size=10))

map_age <-
  basic_map(regions_wcovs, fill = "med_age", scale = F) +
  labs(fill = "", title = "Median age") +
  theme(plot.title = element_text(size=10))

png(here::here(figdir,"map_covariates.png"), height = 1200, width = 1200, res = 150)
(map_age + map_pop) /
  (map_mino + map_imd)
dev.off()

png(here::here(figdir,"map_covariates3.png"), height = 600, width = 1800, res = 150)
(map_age + map_mino + map_imd)
dev.off()


# ---------------------------------------------------------------------------- #

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


# ---------------------------------------------------------------------------- #
## Swab-death lag ##

linelist_deaths %>% 
  mutate(postP2 = (date_swab >= ymd("2020-05-18"))) %>%
  group_by(postP2) %>%
  summarise(mean = mean(swab_death, na.rm = T),
            sd = sqrt(var(swab_death, na.rm = T)),
            quants = paste0(quantile(swab_death, c(0.05,0.25,0.5,0.75,0.95), na.rm = T),collapse = "-"))

summary(linelist_deaths$swab_death)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -180.00    3.00    6.00    8.24   11.00  219.00   33430 

summary(linelist_deaths$swab_death[linelist_deaths$swab_death >= 0 & linelist_deaths$swab_death <= 100])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    3.00    6.00    8.22   11.00  100.00   33430 

## By geography
linelist_deaths %>% 
  mutate(postP2 = (date_swab >= ymd("2020-05-18"))) %>%
  filter(!is.na(postP2)) %>%
  left_join(dplyr::select(regions.df, lad19cd, geography)) %>%
  group_by(geography, postP2) %>%
  summarise(mean = mean(swab_death, na.rm = T),
            sd = sqrt(var(swab_death, na.rm = T)),
            quants = paste0(quantile(swab_death, c(0.05,0.25,0.5,0.75,0.95), na.rm = T),collapse = "-"))

png(here::here("figures","descriptive","swab_death_lag.png"), height = 800, width = 1000, res = 150)
linelist_deaths %>% 
  filter(swab_death >= 0 & date_swab >= ymd("2020-05-18")) %>%
  ggplot(aes(x = swab_death)) + 
  geom_histogram(bins = 40, fill = "steelblue") +
  labs(x = "Days from swab to death", y = "Density") + 
  xlim(c(0,100))
dev.off()

png(here::here("figures","descriptive","swab_death_lag_bypillar.png"), height = 800, width = 1000, res = 150)
linelist_deaths %>% 
  filter(swab_death >= 0 & date_swab >= ymd("2020-05-18") & pillars %in% c("pillar1", "pillar2")) %>%
  ggplot(aes(x = swab_death)) + 
  geom_histogram(bins = 40, fill = "steelblue") +
  facet_wrap(~pillars) +
  labs(x = "Days from swab to death", y = "Density") + 
  xlim(c(0,100))
dev.off()

linelist_deaths %>%
  group_by(pillars) %>%
  summarise(mean = mean(swab_death, na.rm = T),
            sd = sqrt(var(swab_death, na.rm = T)),
            quants = paste0(quantile(swab_death, c(0.05,0.25,0.5,0.75,0.95), na.rm = T),collapse = "-"))

# ---------------------------------------------------------------------------- #
