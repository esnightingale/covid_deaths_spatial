################################################################################
# Description: Compare model output to cases
# 
# Author: Emily S Nightingale
# Date created: 30/09/2020
# 
################################################################################
################################################################################

################################################################################
# SETUP
################################################################################

library(tidyverse)
library(INLA)
library(spdep)
library(sf)
library(patchwork)
library(viridis)
library(INLAutils)
library(ggspatial)
# devtools::install_github('timcdlucas/INLAutils')

theme_set(theme_minimal())


# measure <- "deaths"
wave <- 1

source(here::here("code","main","functions.R"))

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
merged <- readRDS(here::here("data","merged.rds")) 

deaths <- readRDS(here::here("data","deaths.rds")) 
cases <- readRDS(here::here("data","cases.rds")) 

data.list <- list(deaths = deaths,cases = deaths)

period <- merged$breaks[[wave]]


# Fitted models and posterior samples
fits_d <- readRDS(file = here::here("output",
                                  sprintf("fits_deaths_%s.rds", wave)))
fit_final_d <- fits_d[["BYM_geog"]]

fits_c <- readRDS(file = here::here("output",
                                    sprintf("fits_cases_%s.rds", wave)))
fit_final_c <- fits_d[["BYM_geog"]]

samples_d <- readRDS(file = here::here("output",
                                     sprintf("samples_deaths_%s.rds", wave)))
samples_final_d <- samples_d[["BYM_geog"]]
samples_c <- readRDS(file = here::here("output",
                                       sprintf("samples_cases_%s.rds", wave)))
samples_final_c <- samples_c[["iid_geog"]]

# outputdir <- "deaths"
dat_d <- deaths[[wave]]
dat_c <- cases[[wave]]

################################################################################
# Compare fitted random effects to cases
################################################################################

## Temporal ##

get_fit_peak <- function(dat, samples, offset){
  
fitted <- bind_cols(lapply(samples, function(s) exp(s$latent[1:nrow(dat)]))) 

fitted %>%
  mutate(across(everything(), function(x) x*offset)) %>%
  rowwise() %>%
  mutate(med = median(c_across(everything())),
         lo = quantile(c_across(everything()), p = 0.05),
         hi = quantile(c_across(everything()), p = 0.95)) %>%
  ungroup() -> fitted

dat <- cbind(dat, fitted[,c("med","lo","hi")])

dat %>%
  group_by(lad19cd) %>%
  summarise(fit_peak = week[which.max(med)],
            fit_peak_n = max(med),
            obs_peak = week[which.max(n)],
            obs_peak_n = max(n, na.rm = T),
            peak_n_err = (obs_peak_n - fit_peak_n)/obs_peak_n,
            peak_err = obs_peak - fit_peak) -> dat_peak

return(list(dat_peak,dat))
}

peak_d <- get_fit_peak(dat_d, samples_final_d, offset = dat$E_wk)
peak_c <- get_fit_peak(dat_c, samples_final_c, offset = log(dat$la_pop))

peak <- peak_c[[1]] %>%
  full_join(peak_d[[1]], by = "lad19cd", suffix = c("_c","_d")) %>%
  mutate(lag = as.factor(difftime(fit_peak_c,fit_peak_d, units = "week")))

png(here::here("figures","compare","map_lag_case_deaths.png"), height = 800, width = 1200, res = 150)
regions %>%
  full_join(peak) %>%
  basic_map(fill = "lag") +
  labs(title = "Difference in weeks between fitted peak cases and peak deaths",
       subtitle = "2020-01-01 - 2020-06-30") +
  scale_fill_viridis_d()
dev.off()

png(here::here("figures","compare","map_err_peak.png"), height = 800, width = 1600, res = 150)
regions %>%
  full_join(peak) %>%
  basic_map(fill = "peak_n_err_c") +
  labs(title = "Percentage error between fit and observed peak count - cases",
       subtitle = "2020-01-01 - 2020-06-30") +
  scale_fill_gradient2() -> err_c

regions %>%
  full_join(peak) %>%
  basic_map(fill = "peak_n_err_d") +
  labs(title = "Deaths") +
  scale_fill_gradient2() -> err_d

err_c + err_d

dev.off()


dat <- dat_c %>%
  full_join(dat_d) %>% 

peak_d[[2]] %>%
  filter(la == 16) %>% #Leicester
  ggplot(aes(week, med, ymin = LQ, ymax = UQ)) +
  geom_ribbon(alpha = 0.1) +
  geom_line()
  
## Spatial ##

dat_d <- full_join(dat_d, bind_cols(la = 1:312, fit_spatial = fit_final_d$marginals.random$la$mean[1:312]))


regions %>%
  full_join(dat_d) %>%
  basic_map(fill = "fit_spatial")


View(sample_final_c[[1]]$latent)


################################################################################
# Plot sampled trajectories for LTLAs with different case-death lags
################################################################################


# la_samp <- sample(dat$la,9)
# la_list <- c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")
london <- unique(dat$lad19nm[dat$geography == "London Borough"])
la_list <- c("Leeds", "Bradford",sample(london,2))
la_samp <- unique(dat$la[dat$lad19nm %in% la_list])

get_preds <- function(sample){
  pred <- sample$latent[1:nrow(dat)]
}

get_rshp_preds <- function(dat, samples, la_samp){
preds <- bind_cols(lapply(samples, get_preds))

dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), preds) %>%
  pivot_longer(cols = -1:-8) %>%
  mutate(pred_n = exp(value)*E_wk)

return(dat_pred)

}

dat_pred_c <-get_rshp_preds(dat_c, samples_final_c, la_samp)
dat_pred_d <-get_rshp_preds(dat_d, samples_final_c, la_samp)

dat_pred %>%
  filter(la %in% la_samp) %>%
  ggplot() + 
  geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
  geom_point(aes(week, n)) + 
  facet_wrap(~lad19nm) +
  theme_minimal() -> p




