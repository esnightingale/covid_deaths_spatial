################################################################################
# Description: Summarise and visualise model output
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
list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

measure <- "deaths"
wave <- 1

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))
# border <- st_union(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data","expanded",paste0(measure,".rds")))

dat <- dat_all[[wave]] 
dat$n[dat$wk_since_first < 0] <- NA

n.la <- n_distinct(dat$la)
period <- dat_all$breaks[[wave]]

# Fitted models and posterior samples
fits <- readRDS(file = here::here("output","expanded_data",
                                  sprintf("fits_%s_%s.rds",measure, wave)))
fit_final <- fits[[6]]
samples <- readRDS(file = here::here("output","expanded_data",
                                     sprintf("samples_%s_%s.rds",measure, wave)))
samples_final <- samples[[6]]

cases <- readRDS(here::here("data","expanded","cases.rds"))

# ---------------------------------------------------------------------------- #
######################
#       PRIORS       #
######################

# Specify priors as before, based on real data

# SD over time
dat %>%
  group_by(w) %>%
  summarise(n = sum(n), E = sum(E_wk)) %>%
  mutate(rate = n/E) %>%
  pull(rate) %>%
  var() %>%
  sqrt() -> sd_time

# SD over space
dat %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n), E = unique(E_wk)) %>%
  mutate(rate = n/E) %>%
  pull(rate) %>%
  var() %>%
  sqrt() -> sd_space

# => Less variation over space than over time

# Undajusted for other covariates and structure so take these as upper bounds:
prior.prec.tp <- list(prec = list(prior = "pc.prec",
                                  param = c(sd_time/0.31, 0.01)))
prior.prec.sp <- list(prec = list(prior = "pc.prec",
                                  param = c(sd_space/0.31, 0.01)))

prior.prec.sp
prior.prec.tp

# ---------------------------------------------------------------------------- #

######################################
#   PREDICTION AT ALL TIME POINTS    #
######################################

# week_seq <- ymd(seq(min(dat$week), max(dat$week), by = "week"))
# 
# vars <- names(dplyr::select(dat, E:E_wk_unstrat, la, lad19cd:first, first_overall:prop_kw))
# dat.dt <- as.data.table(dat)
# 
# # Replicate per ltla (by vars are all values I want to copy down per date):
# all_dates <- dat.dt[,.(week=week_seq),by = vars]
# 
# # Merge and fill count with 0:
# setkey(dat.dt, week, E, E_wk, E_wk_unstrat, la, lad19cd, lad19nm, la_pop, geog, geography, area_km2, 
#        pop_dens, IMD, IMD_quint, prop_minority, prop_kw, first_overall, 
#        wk_first_la_overall, first)
# setkey(all_dates, week, E, E_wk, E_wk_unstrat, la, lad19cd, lad19nm, la_pop, geog, geography, area_km2, 
#        pop_dens, IMD, IMD_quint, prop_minority, prop_kw, first_overall, 
#        wk_first_la_overall, first)
# 
# dat.dt <- dat.dt[all_dates,roll=TRUE]
# 
# dat.dt %>%
#   as.data.frame() %>% 
#   mutate(w = as.integer(lubridate::week(week)),
#          w2 = w,
#          w3 = w,
#          wk_since_first = w - first) -> dat_expand
# 
# dat_expand$n[is.na(dat_expand$n)] <- 0
# dat_expand$wk_since_first[dat_expand$wk_since_first < 0] <- NA
# 
# f <- n ~ 
#   IMD_quint + prop_minority + 
#   f(w, model = "rw1",
#     hyper = list(prec = list("pc.prec", c(4.1, 0.01))),
#     # values = 1:26,
#     scale.model = T) +
#   f(wk_since_first, model = "rw2",
#     # values = -13:25,
#     hyper = list(prec = list("pc.prec", c(4.1, 0.01))),
#     replicate = geog,
#     scale.model = T) +
#   f(la, model = "bym2", graph = g,
#     scale.model = T,
#     constr = T,
#     hyper=list(
#       phi =list(param = c(0.5, 2/3)),
#       prec = list("pc.prec", c(35.1, 0.01))) 
#   )
# 
# fit_all <- inla(f,
#                  "zeroinflatedpoisson0",
#                  data = dat_expand,
#                  E = E_wk,
#                  control.compute=list(dic=TRUE, 
#                                       waic=TRUE, 
#                                       cpo = TRUE,
#                                       config = TRUE),
#                  control.predictor = list(compute = TRUE, link = 1),
#                  control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
#                  verbose = T)
# summary(fit_all)
# 
# # Samples are of *marginal* densities
# samples_all <- inla.posterior.sample(n = 1000, fit_all)
# 
# 
# saveRDS(list(fit = fit_all, samples = samples_all), file = here::here("output","expanded_data","fit_samples_alltps.rds"))


sims <- bind_cols(lapply(samples_final, get_preds, dat))
dat_sims <- dat %>%
  dplyr::select(geography, lad19cd, lad19nm, la, la_pop, week, E_wk) %>%
  bind_cols(sims) 

dat_sims <-  dat %>%
  dplyr::select(lad19cd, week, n) %>% 
  full_join(dat_sims) %>% 
  arrange(geography,lad19cd, week) %>%
  pivot_longer(cols = -1:-8) %>%
  mutate(pred = exp(value),
         pred_n = exp(value)*E_wk,
         n = replace_na(n, 0))


dat_sims %>%
  group_by(week, name) %>%
  summarise(pred = sum(pred_n),
            n = sum(n),
            pop = sum(la_pop)) %>%
  ggplot() + 
  geom_line(aes(week, pred*1e5/pop, group = name), alpha = 0.1, col = "grey") +
  geom_point(aes(week, n*1e5/pop)) + 
  labs(x = "Calendar week", y = "Rate per 100,000", title = "Total fit over time, by calendar week", subtitle = "Observed rates shown in black, with 1000 posterior samples in grey") -> plot_fit_time

dat_sims %>%
  group_by(week, name, geography) %>%
  summarise(pred = sum(pred_n),
            n = sum(n),
            pop = sum(la_pop)) %>%
  mutate(group = paste0(name, geography)) -> by_samp_geog
by_samp_geog %>%
  group_by(week, geography) %>%
  summarise(pred = mean(pred),
            n = mean(n),
            pop = mean(pop)) -> by_geog


ggplot(by_samp_geog) + 
  geom_line(aes(week, pred*1e5/pop, group = group, col = geography), alpha = 0.1) +
  geom_point(dat = by_geog, aes(week, n*1e5/pop, col = geography), pch = 21, fill = "white") + 
  labs(x = "Calendar week", y = "Rate per 100,000", title = "Total fit over time, by calendar week and geography", subtitle = "Observed rates shown in white, with 1000 posterior samples", col = "Geography") +
  theme(legend.position = c(0.2,0.7))  -> plot_fit_time_geog

png(here::here("figures",measure,"temp_fit_alltps.png"), height = 1000, width = 1500, res = 150)
plot_fit_time / plot_fit_time_geog
dev.off()


# ---------------------------------------------------------------------------- #

######################################
#  PREDICTION AT AVERAGE COVARIATES  #
######################################

dat %>%
  mutate(prop_minority = median(prop_minority),
         IMD_quint = "(21.3,29.2]",
         n = NA)  -> dat_avgcov

# Make pred data with average covariate values
dat_pred <- bind_rows(dat, dat_avgcov) %>%
  mutate(IMD_quint = factor(IMD_quint, levels = levels(dat$IMD_quint)))

f <- n ~ 
  IMD_quint + prop_minority + 
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = prior.prec.sp) 
  )

fit_pred <- INLA::inla(f,
                 "nbinomial",
                 data = dat_pred,
                 E = E_wk,
                 control.compute=list(dic=TRUE, 
                                      waic=TRUE, 
                                      cpo = TRUE,
                                      config = TRUE),
                 control.predictor = list(compute = TRUE, link = 1),
                 control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
                 verbose = T)
summary(fit_pred)

# Samples are of *marginal* densities
samples_pred <- inla.posterior.sample(n = 1000, fit_pred)

saveRDS(list(fit = fit_pred, samples = samples_pred), file = here::here("output","fit_samples_expanded_avgcov.rds"))

pred_avgcov <- readRDS(here::here("output","fit_samples_expanded_avgcov.rds"))
fit_pred <- pred_avgcov$fit
samples_pred <- pred_avgcov$samples

# ---------------------------------------------------------------------------- #
## Summarise samples by geography ##

sims <- bind_cols(lapply(samples_pred, get_preds, dat_pred))
idx <- -1:-nrow(dat) 

dat_sims <- dat_pred %>%
  dplyr::select(geography, lad19cd, lad19nm, la, la_pop, week, E_wk) %>%
  bind_cols(sims) %>%
  slice(idx)

dat_sims <-  dat %>%
  dplyr::select(lad19cd, week, n) %>% 
  full_join(dat_sims) %>% 
  arrange(geography,lad19cd, week) %>% 
  pivot_longer(cols = -1:-8) %>%
  mutate(pred = exp(value),
         pred_n = exp(value)*E_wk)

png(here::here("figures",measure,"expanded_data","pred_avgcov_geog.png"), height = 1000, width = 1200, res = 150)
print(
  dat_sims %>%
    group_by(week, name, geography) %>%
    summarise(pred_n = sum(pred_n),
              obs = sum(n, na.rm = TRUE),
              pop = sum(la_pop)) %>% 
    ggplot() + 
    geom_line(aes(week, pred_n*1e5/pop, group = name, col = geography), alpha = 0.05) +
    geom_point(aes(week, obs*1e5/pop)) +
    facet_wrap(~geography) + #
    labs(y = "Rate per 100,000", x = "Week", title = "Predicted COVID-19-related deaths at median covariate values") +
    guides(col = FALSE) +
    theme_minimal()
)
dev.off()

dat_sims %>%
  group_by(la, lad19cd, lad19nm, la_pop, geography, week) %>%
  summarise(q01 = quantile(pred_n, 0.01),
            q10 = quantile(pred_n, 0.1),
            q50 = quantile(pred_n, 0.5),
            q90 = quantile(pred_n, 0.9),
            q99 = quantile(pred_n, 0.99)) %>%
  ungroup() -> agg_sims

agg_sims %>%
  filter(la %in% 1:4) %>%
ggplot(aes(week, q50)) +
  geom_ribbon(aes(ymin = q01, ymax = q99), alpha = 0.5) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.5) +
  geom_line() +
  facet_wrap(~lad19nm)

write.csv(agg_sims, file = here::here("output","pred_quants_avgcov.csv"), row.names = F)


# ---------------------------------------------------------------------------- #
## Scale predictions by death:case ratio ##

lag <- 7
# Function: for each simulated trajectory, scale by a sampled CFR from the observed empirical distribution post-P2
# lag_rescale <- function(lag, ratio, sims, nsamp = 1){
#   
#   sims_scale <- sims[1,]
#   for (i in 1:nsamp){
#   sims %>%
#     group_by(name) %>%
#     mutate(scale = EnvStats::remp(1, ratio),
#            week = week - lag,
#            pred_n_scale = pred_n*scale) -> temp
#   sims_scale <- bind_rows(sims_scale, temp) 
#   }
#   return(sims_scale[-1,])
# }

lag_rescale <- function(scale, lag, sims){
  sims %>%
   mutate(week = week - lag,
          pred_n_scale = pred_n*scale) -> sims_scale
  return(sims_scale)
}

set.seed(101)
reconstruct <- function(sims, data, lag, scale_quants, plot_quants, denom_cutoff = 1, suffix = ""){

data %>%
  mutate(n = replace_na(n, 0),
         week = week-lag) %>%
  full_join(dplyr::select(cases[["first"]],lad19cd,geography,week,n), by = c("lad19cd","week","geography"), suffix = c("_d","_c")) %>%
  mutate(CFR_obs = n_c/n_d,
         period = case_when(week < ymd("2020-05-18") ~ 1,
                                   week >= ymd("2020-05-18") ~ 2),
         ID = row_number()) -> ratio

# Drop unstable ratios where denominator less than cutoff
ratio$CFR_obs[ratio$n_d < denom_cutoff | ratio$n_c < denom_cutoff] <- NA

ratio %>% 
  summarise(n_d = sum(n_d, na.rm = TRUE),
            n_c = sum(n_c, na.rm = TRUE),
            CFR_obs = median(CFR_obs, na.rm = TRUE)) %>%
  print()

ratio %>% 
  group_by(period) %>%
  summarise(n_d = sum(n_d, na.rm = TRUE),
            n_c = sum(n_c, na.rm = TRUE),
            CFR_obs = median(CFR_obs, na.rm = TRUE)) %>%
  print()

CFR_obs <- ratio$CFR_obs[!is.na(ratio$CFR_obs) & ratio$period == 2]
x <- 0:max(CFR_obs)
hist(CFR_obs, breaks = 100, prob = T)
lines(x, EnvStats::demp(x, CFR_obs), col = "red")

quants <- quantile(ratio$CFR_obs[ratio$period == 2], probs = scale_quants, na.rm = TRUE)
quants

scaled_sims <- bind_rows(lapply(quants, lag_rescale, lag = lag, sims = sims)) %>%
  group_by(lad19nm, week) %>%
  mutate(sim = row_number()) %>%
  ungroup() %>%
  rename(n_d = n) %>%
  left_join(dplyr::select(cases[["first"]], lad19nm, geography, week, n))

scaled_quants <- scaled_sims %>%
  group_by(lad19nm,la, week) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            obs = unique(n))

la_samp <- sample(scaled_quants$lad19nm, size =  4)
png(here::here("figures","compare","expanded",paste0("reconstr_lasamp_lag_",lag,"_",suffix,".png")), height = 1000, width = 1500, res = 150)
print(
  scaled_quants %>%
    filter(lad19nm %in% la_samp) %>% #View()
    ggplot(aes(week)) + 
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = med), col = "steelblue") +
    scale_x_date(limits = range(data$week)) +
    geom_point(aes(y = obs)) +
    facet_wrap(~lad19nm, scales = "free") +
    labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths - four sampled LTLAs",
         caption = paste0("Median and ",plot_quants[1]*100,"-",plot_quants[2]*100,
                           "% quantiles over 1000 posterior simulations, scaled by ",
                          # samples_cfr, 
                          # "samples from",
                          min(scale_quants)*100,"% to ", max(scale_quants)*100,"% quantiles of",
                          " the observed CFR distribution post-P2 expansion."
                          )
         ) +
    theme_minimal()
)
dev.off()

scaled_quants_geog <- scaled_sims %>%
  group_by(geography, week, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(geography, week) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1]),
            med = quantile(pred_n_scale, 0.5),
            high = quantile(pred_n_scale, plot_quants[2]),
            obs = unique(cases, na.rm = TRUE)) 

png(here::here("figures","compare","expanded",paste0("reconstr_geog_lag_",lag,"_",suffix,".png")), height = 1000, width = 1500, res = 150)
print(
  scaled_quants_geog %>%
    ggplot(aes(week)) + 
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = med), col = "steelblue") +
    scale_x_date(limits = range(data$week)) +
    geom_point(aes(y = obs)) +
    facet_wrap(~geography, scales = "free_y") +
    labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths - by geography type",
         caption = paste0("Median and ",plot_quants[1]*100,"-",plot_quants[2]*100,
                           "% quantiles over 1000 posterior simulations, scaled by ",
                          # samples_cfr, 
                          # "samples from",
                          min(scale_quants)*100,"% to ", max(scale_quants)*100, "% quantiles of",
                          " the observed CFR distribution post-P2 expansion.")) +
    theme_minimal()
)
dev.off()

return(scaled_sims)

}

scale_quants <- seq(0.1,0.9,0.1)
plot_quants <- c(0.1,0.9)

scaled_sims7 <- reconstruct(sims = dat_sims, data = dat, lag = 7, scale_quants = scale_quants, plot_quants = plot_quants, suffix = "avgcov")
scaled_sims14 <- reconstruct(sims = dat_sims, data = dat, lag = 14, scale_quants = scale_quants, plot_quants = plot_quants, suffix = "avgcov")


## Total cases per week - two lags
scaled_sims7 %>%
  group_by(week, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(week) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            cases = unique(cases, na.rm = TRUE)) %>%
  ggplot(aes(week)) + 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  scale_x_date() +
  geom_point(aes(y = cases)) +
  labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths across England",
       subtitle = "1 week lag") +
  theme_minimal() -> reconstr_total7

scaled_sims14 %>%
  group_by(week, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(week) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            cases = unique(cases, na.rm = TRUE)) %>%
  ggplot(aes(week)) + 
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  scale_x_date() +
  geom_point(aes(y = cases)) +
  labs(x = "",y = "", title = "",
       subtitle = "2 week lag",
       caption = paste0("Median and ",plot_quants[1]*100,"-",plot_quants[2]*100,
                        "% quantiles over 1000 posterior simulations, scaled by ",
                        # samples_cfr, 
                        "samples",
                        # min(scale_quants)*100,"% to ", max(scale_quants)*100, "% quantiles",
                        " from the observed CFR distribution post-P2 expansion.")) +
  theme_minimal() -> reconstr_total14

png(here::here("figures","compare","expanded","reconstr_totals.png"), height = 800, width = 2000, res = 150)
reconstr_total7 + reconstr_total14
dev.off()

## Total cases overall
scaled_sims14 %>%
  group_by(sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            observed = unique(cases, na.rm = TRUE)) -> totals
totals
#     low     med    high observed
# 107783. 334155. 884277.   229758

dplyr::select(totals, -observed) - totals$observed
#       low      med     high
# -121974.7 104396.6 654519.5

totals %>%
  mutate(across(-observed, function(x) (x - totals$observed)*100/totals$observed)) %>%
  as.data.frame()

## By geography
scaled_sims14 %>%
  group_by(geography, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(geography) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            observed = unique(cases, na.rm = TRUE)) %>%
  ungroup() -> geog_totals
geog_totals

geog_totals %>%
  mutate(across(c(-geography,-observed), function(x) (x - geog_totals$observed)*100/geog_totals$observed)) %>%
  as.data.frame()

# ---------------------------------------------------------------------------- #

# Map out the median differences
scaled_sims14 %>%
  group_by(lad19cd, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(lad19cd) %>%
  summarise(med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            cases = unique(cases, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(perc_diff = (med - cases)*100/cases) -> scale_by_la

summary(scale_by_la$perc_diff)

scale_by_la %>%
  summarise(low = quantile(perc_diff, plot_quants[1], na.rm = TRUE),
            med = quantile(perc_diff, 0.5, na.rm = TRUE),
            high = quantile(perc_diff, plot_quants[2], na.rm = TRUE)) %>%
  as.data.frame()
  



png(here::here("figures","compare","expanded","underascertainment_lag_14.png"), height = 1000, width = 1500, res = 150)
regions %>%
  full_join(scale_by_la) %>%
  basic_map(fill = "ratio") + 
  scale_fill_gradient2(midpoint = 0, trans = "log2") +
  labs(fill = "Ratio", title = "Relative difference between inferred total cases and observed test positives",
       subtitle = paste0("Predicted deaths back-dated by two weeks"))
dev.off()


# 
# 
# CFR <- ratio$CFR_obs[!is.na(ratio$CFR_obs)]
# hist(CFR, breaks = 40, prob = T)
# lines(dgamma(seq(0,10,0.1), shape = 4, scale = 1/8))
# lines(dexp(seq(0,10,0.1), rate = 7))
# 
# fitdistrplus::descdist(CFR)
# fitgam <- fitdistrplus::fitdist(log(CFR), distr = "gamma", lower = c(0, 0), start = list(shape = 3, rate = 4))
# 
# fitweib <- fitdistrplus::fitdist(ratio$CFR_obs, distr = "weibull", lower = c(0, 0), start = list(shape = 3, scale = 1/10))
# fitexp <- fitdistrplus::fitdist(ratio$CFR_obs, distr = "exp")
# fitbeta <- fitdistrplus::fitdist(CFR, distr = "beta", lower = c(0,0), start = list(shape1 = 0.5, shape2 = 0.5)) 
# 


# hist(ratio$CFR_obs, breaks = 100, prob = TRUE)
# lines(dgamma(seq(0,10,0.1), shape = fitgam$estimate[1], rate = fitgam$estimate[2]), col = "orange")
# lines(dgamma(seq(0,10,0.1), shape = 1, scale = 0.1), col = "red")
# lines(dweibull(seq(0,10,0.1), shape = fitweib$estimate[1], scale = fitweib$estimate[2]), col = "green")
# lines(dexp(seq(0,10,0.1), rate = fitexp$estimate[1]), col = "blue")


# ---------------------------------------------------------------------------- #
## Mean fitted value + HPD interval ##
# 
# get_hpd <- function(marg){
#   hpd <- as.data.frame(t(
#     # extract hpd for each value of the variable
#     sapply(marg, 
#            FUN = function(m) inla.hpdmarginal(0.98, m))))
#   hpd <- setNames(hpd, c("low","high"))
#   return(hpd)
# }
# 
# hpd.fit <- get_hpd(fit_pred$marginals.fitted.values)
# dat_pred <- bind_cols(dat_pred, hpd.fit) %>%
#   mutate(mean.fit = fit_pred$summary.fitted.values$mean)
# 
# 
# dat_pred %>%
#   slice(-1:-4499) %>%
#   filter(la %in% 1:4) %>%
#   ggplot(aes(week, mean.fit)) +
#   geom_ribbon(aes(ymin = low, ymax = high), fill = "red", alpha = 0.2) +
#   geom_line(col = "red") +
#   facet_wrap(~lad19nm)