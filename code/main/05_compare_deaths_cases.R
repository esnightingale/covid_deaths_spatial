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

# source(here::here("code","main","functions.R"))
list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

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

# fits_c <- readRDS(file = here::here("output",
#                                     sprintf("fits_cases_%s.rds", wave)))
# fit_final_c <- fits_d[["BYM_geog"]]
# fit_final_c <- readRDS(file = here::here("output",
#                                              sprintf("fit_final_cases_%s.rds", wave)))


samples_d <- readRDS(file = here::here("output",
                                     sprintf("samples_deaths_%s.rds", wave)))
samples_final_d <- samples_d[["BYM_geog"]]
# samples_c <- readRDS(file = here::here("output",
#                                        sprintf("samples_cases_%s.rds", wave)))
# samples_final_c <- samples_c[["iid_geog"]]

# samples_final_c <- readRDS(file = here::here("output",
#                                              sprintf("samples_final_cases_%s.rds", wave)))

# outputdir <- "deaths"
dat_d <- deaths[[wave]]
dat_c <- cases[[wave]] %>%
  mutate(E_wk = log(la_pop))

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

peak_c <- get_fit_peak(dat_c, samples_final_c, offset = 1)
peak_d <- get_fit_peak(dat_d, samples_final_d, offset = dat_d$E_wk)

peak <- peak_c[[1]] %>%
  full_join(peak_d[[1]], by = "lad19cd", suffix = c("_c","_d")) %>%
  mutate(lag = as.factor(difftime(fit_peak_d,fit_peak_c, units = "week")))

png(here::here("figures","compare","map_lag_case_deaths.png"), height = 800, width = 1200, res = 150)
regions %>%
  full_join(peak) %>%
  basic_map(fill = "lag") +
  labs(title = "Difference in weeks between fitted peak cases and peak deaths",
       subtitle = "2020-01-01 - 2020-06-30") +
  scale_fill_viridis_d()
dev.off()

png(here::here("figures","compare","map_err_peak.png"), height = 800, width = 20000, res = 150)
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


# dat <- dat_c %>%
#   full_join(dat_d) %>% 
# 
# peak_d[[2]] %>%
#   filter(la == 16) %>% #Leicester
#   ggplot(aes(week, med, ymin = LQ, ymax = UQ)) +
#   geom_ribbon(alpha = 0.1) +
#   geom_line()
  
## Spatial ##

# dat_d <- full_join(dat_d, bind_cols(la = 1:312, fit_spatial = fit_final_d$marginals.random$la$mean[1:312]))


# regions %>%
#   full_join(dat_d) %>%
#   basic_map(fill = "fit_spatial")


# View(sample_final_c[[1]]$latent)


################################################################################
# Plot sampled trajectories for LTLAs with different case-death lags
################################################################################


# la_samp <- sample(dat$la,9)
# la_list <- c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")

# london <- unique(dat$lad19nm[dat$geography == "London Borough"])
# la_list <- c("Leeds", "Bradford",sample(london,2))
# la_samp <- unique(dat$la[dat$lad19nm %in% la_list])

get_preds <- function(sample,dat){
  pred <- sample$latent[1:nrow(dat)]
}

get_rshp_preds <- function(dat, samples, scale){
  
preds <- bind_cols(lapply(samples, get_preds, dat))


dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), preds) %>%
  pivot_longer(cols = -1:-8) 

# since log(pop) used as an offset rather than E, don't need to multiply
if (scale == "E"){
  dat_pred <- dat_pred %>%
    mutate(pred_n = exp(value)*E_wk)
}else{
  dat_pred <- dat_pred %>%
    mutate(pred_n = exp(value))
}

return(dat_pred)

}

dat_pred_c <-get_rshp_preds(dat_c, samples_final_c, scale = "offset") 
dat_pred_d <-get_rshp_preds(dat_d, samples_final_d, scale = "E") %>%
  left_join(peak) %>%
  mutate(week_lag = week - as.numeric(lag)*7)


# dat_pred <- dat_pred_c %>%
  # full_join(dat_pred_d, by = c("lad19cd","la_pop","week","name"), suffix = c("_c","_d"))
# 
# dat_pred %>%
#   filter(la %in% la_samp) %>%
#   ggplot() + 
#   geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
#   geom_point(aes(week, n)) + 
#   facet_wrap(~lad19nm) +
#   theme_minimal() -> p


################################################################################
# Plot back-shifted deaths versus observed cases
################################################################################

dat_pred_d %>%
  mutate(lag = round(rweibull(nrow(dat_pred_d),shape = 1.5,scale = 18)),
         week_onset = week - lag,
         week = lubridate::floor_date(week_onset, unit = "week")) -> dat_lag

dat_c %>%
  dplyr::select(lad19cd, week, n) %>%
  right_join(dat_lag, by = c("lad19cd","week"), suffix = c("_c","_d")) %>%
  mutate(CFR = n_d/n_c) -> dat_lag_cfr

# dat_lag %>%
#   full_join(dplyr::select(dat_c,lad19cd, week, n), by = c("lad19cd","week"), suffix = c("_d","_c")) %>%
#   mutate(CFR = n_d/n_c) -> dat_lag

summary(dat_lag$CFR)

dat_lag %>%
  group_by(week, name) %>%
  summarise(n_c = sum(n_c, na.rm = TRUE),
            n_d = sum(n_d, na.rm = TRUE),
            pred_n = sum(pred_n, na.rm = TRUE),
            pop = sum(la_pop), 
            CFR = mean(CFR, na.rm = TRUE)) %>%
  ggplot(aes(week, CFR, group = name)) +
  geom_line(col= "grey",alpha = 0.1) +
  ylim(c(0,50))


  # ggplot(aes(x = week)) +
  # geom_line(aes(y = pred_n, group = name), alpha = 0.1, col = "grey") +
  # geom_point(aes(y = n_c)) +
  # geom_point(aes(y = n_d), col = "steelblue4")


set.seed(101)
la_samp <- sample(dat_d$la,9)
# la_list <- c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")
# london <- unique(dat$lad19nm[dat$geography == "London Borough"])
# la_list <- c("Leeds", "Bradford",sample(london,2))
# la_samp <- unique(dat$la[dat$lad19nm %in% la_list])

png(here::here("figures","compare","death_vs_cases_sample_lag2.png"), height = 1000, width = 1500, res = 150)
print(
  dat_pred_c %>%
    filter(la %in% la_samp) %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    geom_line(data = filter(dat_pred_d,la %in% la_samp), aes(week-2, pred_n, group = name), alpha = 0.1, col = "steelblue3") +
    # geom_smooth(data = filter(dat_pred_d,la %in% la_samp), aes(week, n), se = TRUE, col = "steelblue") + 
    geom_point(data = filter(dat_pred_d,la %in% la_samp), aes(week-2, n), col = "steelblue4") + 
    facet_wrap(~lad19nm, scales = "free") +
    theme_minimal()
)
dev.off()

png(here::here("figures","compare","death_vs_cases_sample_imputedonset.png"), height = 1000, width = 1500, res = 150)
print(
  dat_c %>%
    filter(la %in% la_samp) %>%
    ggplot() + 
    # geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    geom_line(data = filter(dat_lag,la %in% la_samp), aes(week_onset, pred_n, group = name), alpha = 0.1, col = "steelblue3") +
    # geom_smooth(data = filter(dat_pred_d,la %in% la_samp), aes(week, n), se = TRUE, col = "steelblue") + 
    geom_point(data = filter(dat_lag,la %in% la_samp), aes(week_onset, n), col = "steelblue4") + 
    facet_wrap(~lad19nm, scales = "free") +
    theme_minimal()
)
dev.off()


## By geography
dat_pred_c_geog <- dat_pred_c %>%
  group_by(week, name, geography) %>%
  summarise(n = sum(n),
            pred_n = sum(pred_n))

dat_pred_d_geog <- dat_pred_d %>%
  group_by(week, geography, name) %>%
  summarise(n = sum(n),
            pred_n = sum(pred_n))

png(here::here("figures","compare","death_vs_cases_geog_nolag.png"), height = 1000, width = 1500, res = 150)

print(
  dat_pred_c_geog %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    geom_line(data = dat_pred_d_geog, aes(week, pred_n, group = name), alpha = 0.1, col = "steelblue") +
    # geom_smooth(data = dat_pred_c_geog, aes(week, n), se = TRUE, col = "steelblue") + 
    geom_point(data = dat_pred_d_geog, aes(week, n), col = "darkblue") +
    facet_wrap(~geography, scales = "free") +
    theme_minimal()
)
dev.off()


png(here::here("figures","compare","death_vs_cases_geog_lag2.png"), height = 1000, width = 1500, res = 150)

print(
  dat_pred_c_geog %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    geom_line(data = dat_pred_d_geog, aes(week-2, pred_n, group = name), alpha = 0.1, col = "steelblue") +
    # geom_smooth(data = dat_pred_c_geog, aes(week, n), se = TRUE, col = "steelblue") + 
    geom_point(data = dat_pred_d_geog, aes(week-2, n), col = "darkblue") +
    facet_wrap(~geography, scales = "free") +
    theme_minimal()
)
dev.off()



dat_pred_d_geog <- dat_pred_d %>%
  group_by(week_lag, geography, name) %>%
  summarise(n = sum(n),
            pred_n = sum(pred_n))

png(here::here("figures","compare","death_vs_cases_geog_wlag.png"), height = 1000, width = 1500, res = 150)

print(
  dat_pred_c_geog %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    geom_line(data = dat_pred_d_geog, aes(week_lag, pred_n, group = name), alpha = 0.1, col = "steelblue") +
    # geom_smooth(data = dat_pred_c_geog, aes(week, n), se = TRUE, col = "steelblue") + 
    geom_point(data = dat_pred_d_geog, aes(week_lag, n), col = "darkblue") +
    facet_wrap(~geography, scales = "free") +
    theme_minimal()
)
dev.off()

