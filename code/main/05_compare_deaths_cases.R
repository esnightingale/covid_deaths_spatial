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

# measure <- "deaths"
wave <- 1

list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))
border <- st_union(regions)

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
fit_final_c <- readRDS(file = here::here("output",
                                             sprintf("fit_final_cases_%s.rds", wave)))


samples_d <- readRDS(file = here::here("output",
                                     sprintf("samples_deaths_%s.rds", wave)))
samples_final_d <- samples_d[["BYM_geog"]]
# samples_c <- readRDS(file = here::here("output",
#                                        sprintf("samples_cases_%s.rds", wave)))
# samples_final_c <- samples_c[["iid_geog"]]

samples_final_c <- readRDS(file = here::here("output",
                                             sprintf("samples_final_cases_%s.rds", wave)))

# outputdir <- "deaths"
dat_d <- deaths[[wave]] %>%
  mutate(obs_peak = week[which.max(n)],
         obs_peak_n = max(n, na.rm = T))
dat_c <- cases[[wave]] %>%
  mutate(E_wk = E_wk_unstrat,
         obs_peak = week[which.max(n)],
         obs_peak_n = max(n, na.rm = T))

################################################################################
# Compare fitted random effects
################################################################################

## Total (structured + unstructured) spatial random effects (BYM):

bym_c <- data.frame(lad19cd = unique(dat_c$lad19cd), value = fit_final_c$summary.random$la[1:nrow(regions),"mean"], effect = "Total")
bym_d <- data.frame(lad19cd = unique(dat_d$lad19cd), value = fit_final_d$summary.random$la[1:nrow(regions),"mean"], effect = "Total")

regions %>%
  full_join(bym_d) %>%
  basic_map(fill = "value") +
  labs(subtitle = "Deaths") -> map_bym_d

regions %>%
  full_join(bym_c) %>%
  basic_map(fill = "value") +
  labs(subtitle = "Cases", title = "Total  random effect (BYM)")  -> map_bym_c


## Structured-only random effects (Besag):

besag_c <- data.frame(lad19cd = unique(dat_c$lad19cd), value = fit_final_c$summary.random$la[(nrow(regions)+1):(2*nrow(regions)),"mean"], effect = "Structured")
besag_d <- data.frame(lad19cd = unique(dat_d$lad19cd), value = fit_final_d$summary.random$la[(nrow(regions)+1):(2*nrow(regions)),"mean"], effect = "Structured")

regions %>%
  full_join(besag_d) %>%
  basic_map(fill = "value") +
  labs(subtitle = "Deaths") -> map_besag_d

regions %>%
  full_join(besag_c) %>%
  basic_map(fill = "value") +
  labs(subtitle = "Cases", subtitle = "Structured random effect") -> map_besag_c

png(here::here("figures","compare","case_death_spatial_re.png"), height = 800, width = 1000, res = 150)
(map_bym_c + map_bym_d ) / (map_besag_c + map_besag_d)
dev.off()

################################################################################
# Compare fitted peaks
################################################################################


get_fit_peak <- function(dat, samples, offset){
  
fitted <- bind_cols(lapply(samples, function(s) exp(s$latent[1:nrow(dat)]))) %>%
  mutate(across(everything(), function(x) x*offset)) %>%
  bind_cols(dat) %>%
  pivot_longer(1:1000,names_to = "sim") 
  
fitted %>%
  group_by(sim, lad19cd) %>%
  summarise(fit_peak = week[which.max(value)],
            fit_peak_n = max(value),
            obs_peak = unique(obs_peak),
            obs_peak_n= unique(obs_peak_n)
            ) %>%
  ungroup() -> dat_sim_peak

dat_sim_peak %>%
  group_by(lad19cd) %>%
  summarise(fit_peak = median(fit_peak),
            fit_peak_med = median(fit_peak_n),
            fit_peak_lo = quantile(fit_peak_n, 0.25),
            fit_peak_hi = quantile(fit_peak_n, 0.75),
            obs_peak = unique(obs_peak),
            obs_peak_n= unique(obs_peak_n),
            peak_n_err = (obs_peak_n - fit_peak_med)/obs_peak_n,
            peak_err = obs_peak - fit_peak
            ) %>%
  ungroup() -> dat_peak

return(list(dat_peak, dat_sim_peak, fitted))
}

peak_c <- get_fit_peak(dat_c, samples_final_c, offset = dat_c$E_wk)
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

png(here::here("figures","compare","map_err_peak.png"), height = 800, width = 2000, res = 150)
regions %>%
  full_join(peak) %>%
  basic_map(fill = "peak_n_err_c") +
  labs(title = "Percentage error between fit and observed peak count (2020-01-01 - 2020-06-30)",
       subtitle = "Cases") +
  scale_fill_gradient2() -> err_c

regions %>%
  full_join(peak) %>%
  basic_map(fill = "peak_n_err_d") +
  labs(subtitle = "Deaths") +
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


View(samples_final_c[[1]]$latent)

dat_d %>%
  filter(w %in% 22:25) %>% 
  group_by(lad19cd) %>%
  summarise(n = sum(n, na.rm = TRUE),
            la_pop = mean(la_pop),
            rate = n/la_pop) -> june_total_d
  
regions %>%
  full_join(june_total_d) %>%
  basic_map(fill = "n", rate1e5 = TRUE) 
  # labs(title = "Difference in weeks between fitted peak cases and peak deaths",
  #      subtitle = "2020-01-01 - 2020-06-30") 

################################################################################
# Plot sampled trajectories for LTLAs with different case-death lags
################################################################################

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

dat_pred_c <-get_rshp_preds(dat_c, samples_final_c, scale = "E") 
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

## By LTLA
set.seed(101)
la_samp <- sample(dat_d$la,9)

# london <- unique(dat$lad19nm[dat$geography == "London Borough"])
# la_list <- c("Leeds", "Bradford",sample(london,2))
# la_samp <- unique(dat$la[dat$lad19nm %in% la_list])

plot_one_la(124) # Ashford
plot_one_la() # Fylde

la_list <- c("Ashford", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")

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
    # geom_smooth(aes(week, n), se = TRUE, col = "grey", alpha = 0.5) +
    geom_line(data = dat_pred_d_geog, aes(week, pred_n, group = name), alpha = 0.1, col = "steelblue") +
    geom_point(data = dat_pred_d_geog, aes(week, n), col = "darkblue") +
    facet_wrap(~geography, scales = "free") +
    theme_minimal()
)
dev.off()


png(here::here("figures","compare","death_vs_cases_geog_lag7.png"), height = 1000, width = 1500, res = 150)

print(
  dat_pred_c_geog %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    # geom_smooth(aes(week, n), se = TRUE, col = "grey", alpha = 0.5) +
    geom_line(data = dat_pred_d_geog, aes(week-7, pred_n, group = name), alpha = 0.1, col = "steelblue") +
    geom_point(data = dat_pred_d_geog, aes(week-7, n), col = "darkblue") +
    facet_wrap(~geography, scales = "free") +
    theme_minimal()
)
dev.off()

dat_pred_c_geog %>%
  ggplot() + 
  geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
  geom_point(aes(week, n), col = "black") + 
  geom_line(data = dat_pred_d_geog, aes(week-2, pred_n, group = name), alpha = 0.1, col = "steelblue") +
  # geom_smooth(data = dat_pred_c_geog, aes(week, n), se = TRUE, col = "steelblue") + 
  geom_point(data = dat_pred_d_geog, aes(week-2, n), col = "darkblue") +
  facet_wrap(~geography, scales = "free") +
  theme_minimal()

dat_pred_d_geog2 <- dat_pred_d %>%
  group_by(week_lag, geography, name) %>%
  summarise(n = sum(n),
            pred_n = sum(pred_n))

png(here::here("figures","compare","death_vs_cases_geog_wlag.png"), height = 1000, width = 1500, res = 150)

print(
  dat_pred_c_geog %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    geom_line(data = dat_pred_d_geog2, aes(week_lag, pred_n, group = name), alpha = 0.1, col = "steelblue") +
    # geom_smooth(data = dat_pred_c_geog, aes(week, n), se = TRUE, col = "steelblue") + 
    geom_point(data = dat_pred_d_geog2, aes(week_lag, n), col = "darkblue") +
    facet_wrap(~geography, scales = "free") +
    theme_minimal()
)
dev.off()


################################################################################
# Calculate CFR from back-shifted deaths and plot over time
################################################################################

deaths[[1]] %>%
  mutate(week = week-14) %>%
  inner_join(filter(cases[[1]], week >= ymd("2020-05-18")), by = c("lad19cd","week"), suffix = c("_d","_c")) %>%
  mutate(CFR_obs = n_d/n_c) %>%
  group_by(geography_d) %>%
  mutate(scale_geog = median(CFR_obs)) %>%
  ungroup() -> ratio

scale = median(ratio$CFR_obs)

dat_pred_d %>%
  mutate(lag = 14,#round(rweibull(nrow(dat_pred_d),shape = 1.5,scale = 18)),
         week_orig = week,
         week = week - lag) %>% 
  group_by(lad19cd, week) %>%
  mutate(sim = row_number()) %>%
  ungroup() -> dat_pred_lag


dat_pred_c %>%
  group_by(lad19cd, week) %>%
  mutate(sim = row_number()) %>%
  ungroup() %>%
  dplyr::select(lad19cd, week, n, pred_n, sim) %>%
  inner_join(dat_pred_lag, by = c("lad19cd","week","sim"), suffix = c("_c","_d")) %>%
  mutate(CFR_obs = n_d/n_c,
         CFR_fit_d = pred_n_d/n_c,
         CFR_fit_d_c = pred_n_d/pred_n_c,
         pred_n_scale = pred_n_d/scale,
         n_scale = n_d/scale) -> dat_lag_cfr

ggplot(dat_lag_cfr, aes(CFR_fit_d)) +
  geom_histogram(fill = "white", col = "steelblue") +
  scale_x_continuous(trans = "log")

summary(dat_lag_cfr$CFR_fit_d)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    0.13    0.23    0.54    0.51   24.40  176000 


dat_lag_cfr %>%
  group_by(sim) %>%
  summarise(n_c = sum(n_c, na.rm = TRUE),
            n_d = sum(n_d, na.rm = TRUE),
            pred_n_c = sum(pred_n_c, na.rm = TRUE),
            pred_n_d = sum(pred_n_d, na.rm = TRUE),
            CFR_obs = n_d/n_c,
            CFR_fit_d = pred_n_d/n_c,
            CFR_fit_d2 = mean(CFR_fit_d),
            CFR_fit_d_c = pred_n_d/pred_n_c,
            CFR_fit_d_c2 = mean(CFR_fit_d_c)) %>% 
  ungroup() %>%
  ggplot() + 
  geom_histogram(aes(x = CFR_fit_d), alpha = 0.5, fill = "indianred") +
  geom_histogram(aes(x = CFR_fit_d_c), alpha = 0.5, fill = "steelblue") +
  geom_vline(aes(xintercept = CFR_obs), col = "forestgreen", lty = "dashed") +
  geom_vline(aes(xintercept = median(CFR_fit_d)), col = "red", lty = "dashed") +
  geom_vline(aes(xintercept = median(CFR_fit_d_c)), col = "blue", lty = "dashed") 

dat_lag_cfr %>%
  group_by(week, sim) %>%
  summarise(n_c = sum(n_c, na.rm = TRUE),
            n_d = sum(n_d, na.rm = TRUE),
            pred_n_c = sum(pred_n_c, na.rm = TRUE),
            pred_n_d = sum(pred_n_d, na.rm = TRUE),
            pop = sum(la_pop), 
            CFR_obs = n_d/n_c,
            CFR = pred_n_d/n_c) -> cfr_week_sim


## Overall
png(here::here("figures","compare","median_CFR_all_lag2.png"), height = 1000, width = 1500, res = 150)
dat_lag_cfr %>%
  group_by(week, sim) %>%
  summarise(CFR = median(CFR_fit_d, na.rm = TRUE),
            CFR_obs = median(CFR_obs, na.rm = TRUE)) %>% 
  ggplot(aes(week)) +
  geom_line(aes(y = CFR, group = sim), col= "forestgreen",alpha = 0.01) + 
  geom_line(aes(y = CFR_obs), col = "black") +
  xlim(c(ymd("2020-02-02", ymd("2020-06-30")))) +
  ylim(c(0,1)) +
  labs(title = "Observed deaths:cases (black) and simulations from predicted deaths:cases (green)",
       subtitle = "Median across all LTLAs")
dev.off()

## Overall
png(here::here("figures","compare","total_CFR_all_lag2.png"), height = 1000, width = 1500, res = 150)
dat_lag_cfr %>%
  group_by(week, sim) %>%
  summarise(CFR = sum(pred_n_d, na.rm = TRUE)/sum(n_c, na.rm = TRUE),
            CFR_obs = sum(n_d, na.rm = TRUE)/sum(n_c, na.rm = TRUE)) %>% 
  ggplot(aes(week)) +
  geom_line(aes(y = CFR, group = sim), col= "forestgreen",alpha = 0.01) + 
  geom_line(aes(y = CFR_obs), col = "black") +
  xlim(c(ymd("2020-02-02", ymd("2020-06-30")))) +
  ylim(c(0,1)) +
  labs(title = "Observed deaths:cases (black) and simulations from predicted deaths:cases (green)",
       subtitle = "Total ratio across all LTLAs")
dev.off()

## By geog
png(here::here("figures","compare","total_CFR_geog_lag2.png"), height = 1000, width = 1500, res = 150)
dat_lag_cfr %>%
  group_by(geography, week, sim) %>%
  summarise(n_c = sum(n_c, na.rm = TRUE),
            n_d = sum(n_d, na.rm = TRUE),
            pred_n_c = sum(pred_n_c, na.rm = TRUE),
            pred_n_d = sum(pred_n_d, na.rm = TRUE),
            pop = sum(la_pop), 
            CFR = pred_n_d/n_c,
            CFR_obs = n_d/n_c) %>%
  ggplot(aes(week)) +
  geom_line(aes(y = CFR, group = sim, col = geography),alpha = 0.01) + 
  geom_line(aes(y = CFR_obs, col = geography)) +
  facet_wrap(~geography) +
  ylim(c(0,1)) +
  xlim(c(ymd("2020-02-01"), ymd("2020-06-30"))) +
  guides(col = FALSE, fill = FALSE) +
  labs(title = "Observed deaths:cases (black) and simulations from predicted deaths:cases (green)",
       subtitle = "Total ratio across LTLAs in each geography")
dev.off()

# dat_lag_cfr %>%
#   group_by(week, sim) %>%
#   summarise(CFR = median(CFR_fit_d, na.rm = TRUE),
#             CFR_obs = mean(CFR_obs, na.rm = TRUE)) %>% 
#   group_by(week) %>%
#   summarise(CFR_lo = quantile(CFR, prob = 0.05, na.rm = TRUE),
#             CFR_hi = quantile(CFR, prob = 0.95, na.rm = TRUE),
#             CFR_obs = mean(CFR_obs, na.rm = TRUE)) %>% 
#   ggplot(aes(week)) +
#   geom_ribbon(aes(ymin = CFR_lo, ymax = CFR_hi), alpha = 0.3) +
#   geom_line(aes(y = CFR_obs)) +
#   ylim(c(0,1)) +
#   xlim(c(ymd("2020-02-01"), ymd("2020-06-30"))) +
#   labs(title = "Observed deaths:cases and simulated 90% quantile intervals from predicted deaths:cases")


png(here::here("figures","compare","median_CFR_geog_lag2.png"), height = 1000, width = 1500, res = 150)
dat_lag_cfr %>%
  group_by(geography, week, sim) %>%
  summarise(CFR = median(CFR_fit_d, na.rm = TRUE),
            CFR_obs = median(CFR_obs, na.rm = TRUE)) %>% 
  group_by(geography, week) %>%
  summarise(CFR_lo = quantile(CFR, prob = 0.05, na.rm = TRUE),
            CFR_hi = quantile(CFR, prob = 0.95, na.rm = TRUE),
            CFR_obs = median(CFR_obs, na.rm = TRUE)) %>% 
  ggplot(aes(week)) +
  geom_ribbon(aes(ymin = CFR_lo, ymax = CFR_hi, fill = geography), alpha = 0.3) +
  geom_line(aes(y = CFR_obs, col= geography)) +
  ylim(c(0,1)) +
  xlim(c(ymd("2020-02-01"), ymd("2020-06-30"))) +
  guides(fill = FALSE, col = FALSE) +
  facet_wrap(~geography) +
  labs(title = "Observed deaths:cases and simulated 90% quantile intervals from predicted deaths:cases",
       subtitle = "Median across LTLAs in each geography")
dev.off()

# dat_lag_cfr %>%
#   group_by(lad19cd, sim) %>%
#   summarise(pred_n_d = sum(pred_n_d, na.rm = TRUE)) %>%
#   group_by(lad19cd) %>%
#   summarise(pred_n_d = median(pred_n_d)/scale) -> med_pred_d_la
# 
# dat_c %>%
#   group_by(lad19cd) %>%
#   summarise(n_c = sum(n, na.rm = TRUE),
#             la_pop = unique(la_pop)) %>%
#   full_join(med_pred_d_la) -> med_pred_c_la
# 
# regions %>%
#   full_join(med_pred_c_la) %>%
#   basic_map(fill = "n_c", rate1e5 = TRUE) -> case_burden
# 
# regions %>%
#   full_join(med_pred_d_la) %>%
#   basic_map(fill = "pred_n_d", rate1e5 = TRUE) -> lag_pred_death_burden
# 
# case_burden + lag_pred_death_burden


################################################################################
# Inflate back-shifted deaths by baseline CFR and compare to cases
################################################################################


# dat_lag_cfr %>%
#   group_by(lad19cd, geography, sim) %>%
#   summarise(pred_n_d = sum(pred_n_d, na.rm = T),
#             n_d = sum(n_d, na.rm = T),
#             n_c = sum(n_c, na.rm = T),
#             CFR = n_d/n_c,
#             pred_n_scale = pred_n_d/CFR,
#             n_scale = n_d/CFR) -> dat_pred_la


dat_pred_d %>%
  group_by(geography, week, name) %>%
  summarise(pred_n_scale = sum(pred_n, na.rm = TRUE)/scale,
            n_scale = sum(n, na.rm = TRUE)/scale) %>%
  mutate(week = week - 14) %>% 
  inner_join(dat_pred_c_geog) -> dat_pred_geog_scale

png(here::here("figures","compare","death_inflate_geog_lag2wk.png"), height = 1000, width = 1500, res = 150)
print(
  dat_pred_geog_scale %>%
    ggplot(aes(week)) + 
    # geom_line(aes(week, pred_n, group = name), alpha = 0.01, col = "grey") +
    # geom_smooth(aes(week, n), se = TRUE, col = "grey", alpha = 0.5) +
    geom_line(aes(y = pred_n_scale, group = name), alpha = 0.01, col = "steelblue") +
    geom_point(aes(y = n_scale), col = "navy") +
    scale_x_date(limits = range(dat_pred_geog_scale$week)) +
    # geom_line(data = dat_pred_geog_scale, aes(week, n_scale), col = "navy", lty = "dashed") +
    geom_line(aes(y = n), col = "grey") +
    geom_point(aes(y = n), col = "black") + 
    facet_wrap(~geography, scales = "free") +
    theme_minimal()
)
dev.off()


dat_pred_d %>%
  mutate(pred_n_scale = pred_n/scale,
         n_scale = n/scale,
         week = week - 14) %>%
  select(name, lad19nm, week, pred_n_scale, n_scale) %>%
  inner_join(dat_pred_c) -> dat_pred_scale

png(here::here("figures","compare","death_inflate_lasamp_lag2wk.png"), height = 1000, width = 1500, res = 150)
print(
  dat_pred_scale %>%
    filter(la %in% la_samp) %>%
    ggplot(aes(week)) + 
    # geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_line(aes(y = pred_n_scale, group = name), alpha = 0.1, col = "steelblue") +
    geom_point(aes(y = n_scale), col = "navy") +
    # geom_line(aes(week-lag, n/scale), col = "steelblue", lty = "dashed") +
    geom_point(aes(y = n)) + 
    geom_line(aes(y = n), col = "grey") +
    facet_wrap(~lad19nm, scales = "free") +
    # scale_x_date(limits = range(dat_pred_geog_scale$week)) +
    theme_minimal()
)
dev.off()

## For each LTLA and posterior sample, calculate the total cases inferred from 
## back-shifted and upscaled deaths, and compare to total confirmed cases during 
## the same period. This percentage error suggests extent of underascertainent 
## across the country.
## 
## Based on: Predicted deaths per week, per LTLA, shifted back two weeks and scaled up by 0.2.
dat_lag_cfr %>%
  group_by(lad19nm, la, sim) %>%
  summarise(n_c = sum(n_c, na.rm = TRUE),
            n_d = sum(n_d, na.rm = TRUE),
            pred_n_d_scale = sum(pred_n_d/scale),
            pred_n_c = sum(pred_n_c),
            pred_n_d = sum(pred_n_d),
            E_tot = sum(E_wk),
            perc_diff = (pred_n_d_scale - n_c)*100/n_c,
            ratio = pred_n_d_scale/n_c) -> scale_deaths_sim


# Look into the variation of % difference by LTLA:
scale_deaths_sim %>%
  filter(la %in% la_samp) %>%
  ggplot(aes(ratio)) +
  geom_histogram(bins = 40) +
  geom_vline(aes(xintercept = 1)) +
  facet_wrap(~lad19nm)

summary(scale_deaths_sim$ratio)

# Map out the median differences
scale_deaths_sim %>%
  group_by(lad19nm) %>%
  summarise(n_c = unique(n_c),
            n_d = unique(n_d),
            pred_n_c = median(pred_n_c),
            pred_n_d = median(pred_n_d),
            pred_n_d_scale = median(pred_n_d_scale),
            perc_diff = median(perc_diff),
            ratio = median(ratio)) -> scale_deaths_tot # View()
  
summary(scale_deaths_tot$ratio)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5582  1.2991  1.5802  1.6106  1.8812  3.3708

png(here::here("figures","compare","underascertainment_map.png"), height = 1000, width = 1500, res = 150)
regions %>%
  full_join(scale_deaths_tot) %>%
  basic_map(fill = "ratio") + 
  scale_fill_gradient2(midpoint = 1) +
  labs(fill = "Ratio", title = "Relative difference between inferred total cases and observed test positives",
       subtitle = paste0("Predicted deaths back-dated by two weeks and scaled by a factor of ",round(1/scale,2)))
dev.off()

# Differences over time
scale_deaths_sim %>%
  group_by(week, sim) %>%
  summarise(n_c = unique(n_c),
            n_d = unique(n_d),
            pred_n_c = median(pred_n_c),
            pred_n_d = median(pred_n_d),
            pred_n_d_scale = median(pred_n_d_scale),
            perc_diff = median(perc_diff),
            ratio = median(ratio)) -> scale_deaths_time # View()

dat_lag_cfr %>%
  group_by(week, sim) %>%
  summarise(n_c = sum(n_c, na.rm = TRUE),
            pred_n_d_scale = sum(pred_n_d/scale),
            ratio = pred_n_d_scale/n_c) -> scale_deaths_time

png(here::here("figures","compare","underascertainment_timeline.png"), height = 1000, width = 1500, res = 150)
scale_deaths_time %>%
  ggplot(aes(week, ratio, group = sim)) +
  geom_line(alpha = 0.1, col = "forestgreen")
dev.off()

## Which LTLAs are most extreme?
# subset <- slice_max(scale_deaths_tot, abs(ratio - 1), n = 9)
subset <- slice_min(scale_deaths_tot, ratio, n = 6)
greatestdiff <- unique(dat_pred_c$la[dat_pred_c$lad19nm %in% subset$lad19nm])

png(here::here("figures","compare","death_inflate_under.png"), height = 1000, width = 1500, res = 150)
print(
  dat_pred_scale %>%
    filter(la %in% greatestdiff) %>%
    ggplot(aes(week)) + 
    geom_line(aes(y = pred_n_scale, group = name), alpha = 0.1, col = "steelblue") +
    geom_point(aes(y = n_scale), col = "navy") +
    geom_point(aes(y = n)) + 
    geom_line(aes(y = n), col = "grey") +
    facet_wrap(~lad19nm, scales = "free") +
    # scale_x_date(limits = range(dat_pred_geog_scale$week)) +
    theme_minimal()
)
dev.off()


dat_lag_cfr %>%
  group_by(lad19nm, la, sim) %>%
  summarise(n_c = sum(n_c, na.rm = TRUE),
            pred_n_d_scale = sum(pred_n_d/scale)) -> scale_deaths_sim


# Total number of cases missed
scale_deaths_sim %>%
  group_by(sim) %>%
  summarise(missed_cases = sum(pred_n_d_scale) - sum(n_c)) -> tot_missed_sim

summary(tot_missed_sim$missed_cases)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 103369  111154  112851  112829  114522  120251 

ggplot(tot_missed_sim, aes(missed_cases)) +
  geom_histogram(bins = 50)
  

