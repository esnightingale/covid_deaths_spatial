################################################################################
# Description: Summarise and visualise death and case linelists.
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

linelist_deaths <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","deaths")))%>%
  mutate(month = lubridate::month(date, label = TRUE)) 
linelist_cases <- readRDS(paste0(datadir, sprintf("linelist_%s.rds","cases"))) %>%
  mutate(month = lubridate::month(date, label = TRUE)) 

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


################################################################################
################################################################################
