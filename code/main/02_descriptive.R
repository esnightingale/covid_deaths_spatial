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

list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

# Figure output directory
figdir <- "figures/descriptive"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  st_set_crs("+OSGB:1936 +units=m +no_defs") %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
merged <- readRDS(here::here("data","merged.rds")) 

deaths <- readRDS(here::here("data","deaths.rds")) 
cases <- readRDS(here::here("data","cases.rds")) 

data.list <- list(deaths = deaths,cases = deaths)

linelist_deaths <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","deaths")))%>%
  mutate(month = lubridate::month(date, label = TRUE)) 
linelist_cases <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","cases"))) %>%
  mutate(month = lubridate::month(date, label = TRUE)) 

# ---------------------------------------------------------------------------- #

## GEOGRAPHY ##
png(here::here("figures","descriptive","map_geog.png"), height = 800, width = 900, res = 150)
regions %>% basic_map(fill = "geography") + scale_fill_viridis_d()
dev.off()

# ---------------------------------------------------------------------------- #

## COVARIATES ##

# Median age

png(here::here("figures","descriptive","map_age.png"), height = 800, width = 900, res = 150)
print(
regions %>%
  basic_map(fill = "med_age") +
  scale_fill_viridis_c() +
  labs(fill = "", title = "Median age")
)
dev.off()

cov_names <- c("pop_dens", "IMD", "prop_minority","prop_kw")
deaths[[1]] %>%
  group_by(lad19cd) %>%
  summarise_at(vars(cov_names), .funs = base::mean) -> covs

regions %>%
  full_join(covs) %>%
  mutate(pop_dens = as.numeric(pop_dens)) -> regions_wcovs

map_dens <-
  basic_map(regions_wcovs, fill = "pop_dens") +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill = "", title = "Population per KM^2") +
  theme(plot.title = element_text(size=10))

map_imd <-
  basic_map(regions_wcovs, fill = "IMD") +
  labs(fill = "", title = "Index of Multiple Deprivation") +
  theme(plot.title = element_text(size=10))

map_mino <-
  basic_map(regions_wcovs, fill = "prop_minority") +
  labs(fill = "", title = "Proportion of black and \nminority ethnic population") +
  scale_fill_viridis_c(trans = "log10") +
  theme(plot.title = element_text(size=10))

map_kw <-
  basic_map(regions_wcovs, fill = "prop_kw") +
  labs(fill = "", title = "Proportion of population\n classified as key workers") +
  theme(plot.title = element_text(size=10))

png(here::here("figures","descriptive","map_covariates.png"), height = 1200, width = 1200, res = 150)
(map_dens + map_imd) /
  (map_mino + map_kw)
dev.off()

# ---------------------------------------------------------------------------- #

## RATE BY COVARIATES ##

plot_covs <- function(cov, y_trans = "log10", x_trans = "identity", method = "lm"){
  
    deaths[[1]] %>%
      ggplot(aes(as.numeric(!!sym(cov)), n*1e5/la_pop)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = method) +
      scale_x_continuous(trans = x_trans) +
      scale_y_continuous(trans = y_trans) +
      labs(x = "", y = "Weekly rate per 100,000", subtitle = "Deaths") -> d_rate_imd
    
    cases[[1]] %>%
      ggplot(aes(as.numeric(!!sym(cov)), n*1e5/la_pop)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = method) +
      scale_x_continuous(trans = x_trans) +
      scale_y_continuous(trans = y_trans) +
      labs(x = cov, y = "Weekly rate per 100,000" , subtitle = "Cases") -> c_rate_imd
    
    deaths[[1]] %>%
      ggplot(aes(as.numeric(!!sym(cov)), n/E_wk)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = method) +
      scale_x_continuous(trans = x_trans) +
      scale_y_continuous(trans = y_trans) +
      labs(x = "", y = "Count / expected (age-adjusted)") -> d_sir_imd
    
    cases[[1]] %>%
      ggplot(aes(as.numeric(!!sym(cov)), n/E_wk)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = method) +
      scale_x_continuous(trans = x_trans) +
      scale_y_continuous(trans = y_trans) +
      labs(x = cov, y = "Count / expected (age-adjusted)") -> c_sir_imd

    p <- (d_rate_imd + d_sir_imd) / (c_rate_imd + c_sir_imd)

  return(p)

}

pdf(here::here(figdir,"case_death_bycovs.pdf"), height = 8, width = 10)
plot_covs("IMD")
plot_covs("prop_minority", x_trans = "log10")
plot_covs("prop_kw")
plot_covs("pop_dens", x_trans = "log10")
dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - total ##

linelist_cases %>%
  group_by(lad19cd, week, age_group) %>%
  tally() %>% 
  ggplot(aes(week, n, fill = age_group)) +
  geom_col() +
  labs(title = "Total test-confirmed cases in England",
       x = "",
       y = "",
       fill = "Age group"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  scale_fill_viridis_d(option = "plasma") +
  xlim(c(ymd("2020-01-01"),ymd("2020-06-30"))) -> ts_cases

linelist_deaths %>%
  group_by(lad19cd, week, age_group) %>%
  tally() %>% 
  ggplot(aes(week, n, fill = age_group)) +
  geom_col() +
  labs(title = "Total COVID-19-related deaths in England",
       x = "",
       y = "",
       fill = "Age group"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  scale_fill_viridis_d(option = "plasma") +
  xlim(c(ymd("2020-01-01"),ymd("2020-06-30"))) -> ts_deaths

png(here::here(figdir,"total_ts_firstwave.png"), height = 1000, width = 1200, res = 150)
ts_cases / ts_deaths
dev.off()


# ---------------------------------------------------------------------------- #

## MAP TOTALS ##

png(here::here(figdir,"map_totals.png"), height = 1500, width = 1500, res = 150)
  (plot_la_tot(cases,1, title = "Cases per 100,000") + plot_la_tot(cases,2, title = "Cases per 100,000")) / 
    (plot_la_tot(deaths,1, title = "Deaths per 100,000") + plot_la_tot(deaths,2, title = "Deaths per 100,000"))
dev.off()

# ---------------------------------------------------------------------------- #

## TIME SERIES - BY GEOGRAPHY ##

period <- deaths$breaks[[1]]

deaths[[1]] %>%
  group_by(w,week,geography) %>%
  summarise(n = sum(n, na.rm= T),
            geog_pop = sum(la_pop)) %>%
  ggplot(aes(week, n*1e5/geog_pop, group = geography, col = geography)) +
  geom_line() +
  labs(title = "COVID19-related deaths in England, by geography and week of death",
       subtitle = paste(period[1],"-",period[2]),
       x = "",
       y = "Rate per 100,000", 
       colour = "Geography type") +
  geom_vline(xintercept = ymd("2020-03-23"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-03-22"), y = 15, label = "National lockdown enforced", cex = 2, hjust = "right") +
  theme(legend.position = c(0.16,0.60), legend.text=element_text(size=8),  legend.title=element_text(size=8)) +
  xlim(period) -> ts_geog_deaths

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
  # geom_vline(xintercept = ymd("2020-03-23"), lty = "dashed") +
  # annotate("text", xmin = ymd("2020-03-15"), y = 52, label = "National lockdown enforced", cex = 2) +
  # geom_vline(xintercept = ymd("2020-04-01"), lty = "dashed") +
  # annotate("text", x = ymd("2020-04-02"), y = 50, label = "Pillar 2 testing piloted", cex = 2, hjust = "left") +
  geom_vline(xintercept = ymd("2020-04-15"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-04-16"), y = 55, label = "P2 available to care home residents and staff", cex = 2, hjust = "left") +
  geom_vline(xintercept = ymd("2020-05-18"), lty = "dashed", lwd = 0.2) +
  annotate("text", x = ymd("2020-05-19"), y = 60, label = "P2 available to all symptomatic cases", cex = 2, hjust = "left") +
  labs(title = "Confirmed COVID-19 cases in England, by geography and week of specimen",
       x = "Calendar week",
       y = "Rate per 100,000") +
  guides(col = FALSE) +
  xlim(period) -> ts_geog_cases

png(here::here(figdir,"ts_geog.png"), height = 1000, width = 1500, res = 150)
ts_geog_deaths / ts_geog_cases
dev.off()

# ---------------------------------------------------------------------------- #

## Death : case ratio
## Calculate only with cases confirmed after pillar 2 testing introduced

calc_ratio <- function(lag){
deaths[[1]] %>%
  mutate(week = week-lag) %>%
  inner_join(cases[[1]], by = c("lad19cd","week"), suffix = c("_d","_c")) %>%
  mutate(CFR_obs = n_d/n_c,
         period = factor(case_when(week < ymd("2020-05-18") ~ 1,
                                      week >= ymd("2020-05-18") ~ 2),
                            labels = c("2020-01-04 - 2020-05-17","2020-05-18 - 2020-06-28"))) %>%
  group_by(geography_d) %>%
  mutate(scale = median(CFR_obs)) %>%
  ungroup() -> ratio

print(summary(ratio$CFR_obs))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004008 0.090909 0.153846 0.232584 0.272727 2.000000 

## Total ratio, with lag
ratio %>%
  filter(period == "2020-05-18 - 2020-06-28") %>%
  summarise(n_d = sum(n_d, na.rm = T),
            n_c = sum(n_c, na.rm = T),
            CFR_obs = n_d/n_c) -> ratio_tot
print(ratio_tot)
#       n_d   n_c CFR_obs
#   1  3254 26756   0.122

## Ratio per LTLA
ratio %>%
  group_by(lad19cd, period) %>%
  summarise(n_d = sum(n_d, na.rm = T),
            n_c = sum(n_c, na.rm = T),
            CFR_obs = n_d/n_c) -> ratio_la

print(summary(ratio_la$CFR_obs))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02935 0.11765 0.16667 0.22368 0.26087 2.00000 


png(here::here("figures","compare",paste0("map_obs_cfr_",lag,".png")), height = 800, width = 1000, res = 150)
print(
regions %>%
  full_join(ratio_la) %>%
  basic_map(fill = "CFR_obs") +
  facet_wrap(~period) +
  labs(title = "Observed ratio of deaths to cases before/after expansion of Pillar 2 testing",
       subtitle = paste(lag,"day lag"),
       fill = "Ratio")
)
dev.off()

return(ratio)

}

ratio_7 <- calc_ratio(7)
ratio_14 <- calc_ratio(14)

ratio_la %>%
  group_by(geography, period) %>%
  summarise(scale = unique(scale)) -> geog_scale

geog_scale
# geography_d               scale
# <chr>                     <dbl>
#   1 London Borough            0.2  
# 2 Metropolitan District     0.167
# 3 Non-metropolitan District 0.231
# 4 Unitary Authority         0.189

# ---------------------------------------------------------------------------- #

## Overall swab-death lag ##

summary(linelist_deaths$swab_death)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -180.00    3.00    6.00    8.24   11.00  219.00   33430 

summary(linelist_deaths$swab_death[linelist_deaths$swab_death >= 0 & linelist_deaths$swab_death <= 100])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    3.00    6.00    8.22   11.00  100.00   33430 
   
png(here::here("figures","descriptive","swab_death_lag.png"), height = 800, width = 1000, res = 150)
linelist_deaths %>% 
  filter(swab_death >= 0) %>%
  ggplot(aes(x = swab_death)) + 
  geom_histogram(bins = 40, fill = "steelblue") +
  labs(x = "Days from swab to death", y = "Density") + 
  xlim(c(0,100))
dev.off()

# ---------------------------------------------------------------------------- #
