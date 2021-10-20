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

linelist_deaths <- readRDS(here::here("data","raw", sprintf("linelist_%s.rds","deaths")))%>%
  mutate(month = lubridate::month(date, label = TRUE)) 
linelist_cases <- readRDS(here::here("data","raw", sprintf("linelist_%s.rds","cases"))) %>%
  mutate(month = lubridate::month(date, label = TRUE)) 

# ---------------------------------------------------------------------------- #
## Swab-death lag ##

## By time period
linelist_deaths %>% 
  mutate(postP2 = (date_swab >= ymd("2020-05-18"))) %>%
  group_by(postP2) %>%
  summarise(mean = mean(swab_death, na.rm = T),
            sd = sqrt(var(swab_death, na.rm = T)),
            quants = paste0(quantile(swab_death, c(0.05,0.25,0.5,0.75,0.95), na.rm = T),collapse = "-"))
# postP2   mean    sd quants        
# 1 FALSE    8.22  9.02 1-3-6-11-23   
# 2 TRUE     8.28  7.10 1-3-7-11-21 

summary(linelist_deaths$swab_death)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -180.00    3.00    6.00    8.24   11.00  219.00   33430 

summary(linelist_deaths$swab_death[linelist_deaths$swab_death >= 0 & linelist_deaths$swab_death <= 100])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    3.00    6.00    8.22   11.00  100.00   33430 

## By pillar
linelist_deaths %>% 
  group_by(pillars) %>%
  summarise(n = n(),
            mean = mean(swab_death, na.rm = T),
            sd = sqrt(var(swab_death, na.rm = T)),
            quants = paste0(quantile(swab_death, c(0.05,0.25,0.5,0.75,0.95), na.rm = T),collapse = "-"))
#   pillars                             n   mean    sd quants         
# 1 notifications_only_hpt            777 NaN    NA    NA-NA-NA-NA-NA 
# 2 notifications_only_hpt_nhse         6   2.2   1.30 0.4-2-3-3-3    
# 3 notifications_only_hpt_ons       1447 NaN    NA    NA-NA-NA-NA-NA 
# 4 notifications_only_hpt_ons_nhse    41   5.08  5.89 0-1-3-6-15.3   
# 5 notifications_only_nhse             9   3.22  3.15 -0.2-1-3-4-8   
# 6 notifications_only_ons          10722 NaN    NA    NA-NA-NA-NA-NA 
# 7 notifications_only_ons_nhse      2611   6.58  6.87 0-2-5-9-19     
# 8 pillar1                         52953   8.36  8.58 1-3-6-11-22    
# 9 pillar1_and_pillar2               132   8.91  7.40 1.8-3-6-12-23.4
# 10 pillar2                          5679   8.19  6.67 1-3-6-11.5-21.3

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
