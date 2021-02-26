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

list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

measure <- "deaths"
wave <- 1

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))
border <- st_union(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data",paste0(measure,".rds")))

dat <- dat_all[[wave]] 
n.la <- n_distinct(dat$la)
period <- dat_all$breaks[[wave]]

# Fitted models and posterior samples
fits <- readRDS(file = here::here("output",
                                  sprintf("fits_%s_%s.rds",measure, wave)))
fit_final <- fits[[6]]
samples <- readRDS(file = here::here("output",
                                     sprintf("samples_%s_%s.rds",measure, wave)))
samples_final <- samples[[6]]

cases <- readRDS(here::here("data","cases.rds"))

# ---------------------------------------------------------------------------- #
######################
#  PRIOR PREDICTION  #
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

f <- n ~ 
  IMD_quint + prop_minority + 
  f(w, model = "rw1",
    hyper = list(prec = list("pc.prec", c(4.1, 0.01))),
    values = 1:26,
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    values = 0:24,
    hyper = list(prec = list("pc.prec", c(4.1, 0.01))),
    replicate = geog,
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = list("pc.prec", c(35.1, 0.01))) 
  )

dat_na <- dat %>%
  mutate(n = NA)

fit_prior <- inla(f,
                  "nbinomial",
                  data = dat_na,
                  E = E_wk,
                  control.compute = list(cpo = TRUE,
                                         config = TRUE),
                  control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
                  verbose = T)
summary(fit_prior)

dat_na <- dat_na %>%
  mutate(fitted.value = fit_prior$summary.fitted.values$mean)

fit_prior.w <- data.frame(w = 1:26, fit.w = fit_prior$summary.random$w$mean) 
fit_prior.wk_since_first <- data.frame(geography = rep(unique(dat$geography), each = 25),
                                       wk_since_first = fit_prior$summary.random$wk_since_first$ID,
                                       fit.wk_since_first = fit_prior$summary.random$wk_since_first$mean) 
fit_prior.la <- data.frame(la = unique(dat$la),
                           comp = c(rep("v+u",n.la),
                                    rep("u",n.la)),
                           fit.la = fit_prior$summary.random$la$mean)

fit_prior.w <- data.frame(w = 1:26, fit.w = fit_prior$summary.random$w$mean) 

inla.hpdmarginal(0.98,fit_prior$marginals.random$w$index.1)

get_hpd <- function(marg){
  hpd <- as.data.frame(t(
    # extract hpd for each value of the variable
    sapply(marg, 
           FUN = function(m) inla.hpdmarginal(0.98, m))))
  hpd <- setNames(hpd, c("low","high"))
  return(hpd)
}

hpd.w <- get_hpd(fit_prior$marginals.random$w)
hpd.wk_since_first <- get_hpd(fit_prior$marginals.random$wk_since_first)

fit_prior.w <- bind_cols(data.frame(w = 1:26, fit.w = fit_prior$summary.random$w$mean),hpd.w) 
fit_prior.wk_since_first <- bind_cols(data.frame(geography = rep(unique(dat$geography), each = 25),
                                                 wk_since_first = fit_prior$summary.random$wk_since_first$ID,
                                                 fit.wk_since_first = fit_prior$summary.random$wk_since_first$mean),
                                      hpd.wk_since_first)

dat_na_pred <- dat_na %>%
  right_join(fit_prior.w) %>%
  left_join(fit_prior.wk_since_first) %>%
  left_join(fit_prior.la[1:n.la,]) 

## Calendar week

dat_na %>% 
  group_by(lad19cd) %>%
  summarise(E_wk = unique(E_wk)) %>%
  ungroup() %>% View
summarise(E_wk = sum(E_wk, na.rm = T))
dat_na %>%
  group_by(w) %>%
  summarise(n = sum(n, na.rm = T)) -> dat_w

rescale_w <- function(m){exp(m)*dat_w$E_wk}

sapply(fit_prior$marginals.random$w, 
       FUN = function(m) inla.emarginal(fun = rescale, m))

dat_na %>%
  group_by(w) %>%
  summarise(fit.w = sum(exp(fit.w)*E_wk)) %>%
  ggplot(aes(w, fit.w)) +
  geom_line(col = "grey", lwd = 1.2) +
  geom_line(data = dat_w, aes(y = n), col = "steelblue", lwd = 1.2) +
  labs(y = "Count per week", x = "Week", title = "Prior predicted trend by calendar week")

## Epidemic week
dat %>%
  group_by(wk_since_first, geography) %>%
  summarise(n = sum(n, na.rm = T)) -> dat_wk_geog
dat_na %>%
  group_by(wk_since_first, geography) %>%
  summarise(fit.wk_since_first = sum(exp(fit.wk_since_first)*E_wk)) %>%
  ggplot(aes(wk_since_first, fit.wk_since_first, col = geography)) +
  geom_line(lwd = 1.2, lty = "dashed") +
  geom_line(data = dat_wk_geog, aes(y = n), lwd = 1.2) +
  labs(y = "Count per week", x = "Weeks since first death", title = "Prior predicted trends by epidemic week") +
  theme(legend.position = c(0.8,0.8))

## LTLA
dat %>%
  group_by(la, lad19cd) %>%
  summarise(n = sum(n, na.rm = T)) -> dat_la
dat_na %>%
  group_by(la) %>%
  summarise(fit.la = sum(exp(fit.la)*E_wk)) %>%
  full_join(dat_la) %>%
  pivot_longer(c("fit.la","n"))-> dat_la

regions %>%
  full_join(dat_la) %>%
  basic_map(fill = "value", rate1e5 = TRUE) +
  facet_wrap(~name) + 
  labs(title = "Prior predicted count by LTLA", subtitle = "Adjusted for age-stratified expected count and prior predicted temporal trends") 


samples_prior <- inla.posterior.sample(fit_prior, n = 100)
sims <- bind_cols(lapply(samples_prior, get_preds))

pops <- dat %>%
  group_by(lad19cd) %>%
  summarise(la_pop = unique(la_pop))
tot_pop <- sum(pops$la_pop)

dat_sims <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), sims) %>%
  pivot_longer(cols = -1:-8) %>%
  mutate(pred_n = exp(value)*E_wk)


dat_sims %>%
  group_by(week, name) %>%
  summarise(pred = sum(pred_n),
            n = sum(n)) %>%
  group_by(week) %>%
  summarise(upper = quantile(pred, 0.75),
            lower = quantile(pred,0.25)) %>%
  ggplot() + 
  geom_ribbon(aes(week, ymin = lower, ymax = upper),  fill = "grey") 


# ---------------------------------------------------------------------------- #
######################################
#  PREDICTION AT AVERAGE COVARIATES  #
######################################

# Make pred data with average covariate values
dat_pred <- dat %>%
  mutate(prop_minority = median(prop_minority),
         IMD_quint = "(21.3,29.2]",
         E_wk = E_wk_unstrat,
         n = NA)
dat_pred <- bind_rows(dat, dat_pred)

f <- n ~ 
  IMD_quint + prop_minority + 
  f(w, model = "rw1",
    hyper = list(prec = list("pc.prec", c(4.1, 0.01))),
    values = 1:26,
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    values = 0:24,
    hyper = list(prec = list("pc.prec", c(4.1, 0.01))),
    replicate = geog,
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = list("pc.prec", c(35.1, 0.01))) 
  )

fit_pred <- inla(f,
                 "nbinomial",
                 data = dat_pred,
                 E = E_wk,
                 control.compute=list(dic=TRUE, 
                                      waic=TRUE, 
                                      cpo = TRUE,
                                      config = TRUE),
                 control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
                 verbose = T)
summary(fit_pred)

# Samples are of *marginal* densities
samples_pred <- inla.posterior.sample(n = 1000, fit_pred)


saveRDS(list(fit = fit_pred, samples = samples_pred), file = here::here("output","fit_samples_pred_avgcov.rds"))

# ---------------------------------------------------------------------------- #
## Mean fitted value + HPD interval ##

get_hpd <- function(marg){
  hpd <- as.data.frame(t(
    # extract hpd for each value of the variable
    sapply(marg, 
           FUN = function(m) inla.hpdmarginal(0.98, m))))
  hpd <- setNames(hpd, c("low","high"))
  return(hpd)
}

hpd.fit <- get_hpd(fit_pred$marginals.fitted.values)
dat_pred <- bind_cols(dat_pred, hpd.fit) %>%
  mutate(mean.fit = fit_pred$summary.fitted.values$mean)


dat_pred %>%
  slice(4500:8998) %>%
  filter(la %in% 1:4) %>%
  ggplot(aes(week, mean.fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "red", alpha = 0.2) +
  geom_line(col = "red") +
  facet_wrap(~lad19nm)

# ---------------------------------------------------------------------------- #
## Quantile interval across samples ##

sims <- bind_cols(lapply(samples_pred, get_preds, dat_pred))
idx <- (nrow(dat)+1):(2*nrow(dat))

dat_sims <- bind_cols(dplyr::select(dat_pred, geography, lad19cd, lad19nm, la, la_pop, week, E_wk), sims) %>%
  # Select only new data predictions
  slice(idx) %>%
  pivot_longer(cols = -1:-7) %>%
  mutate(pred = exp(value),
         pred_n = exp(value)*E_wk)

# png(here::here("figures",measure,"preds_geog.png"), height = 1000, width = 1200, res = 150)
print(
  dat_sims %>%
    group_by(week, name, geography) %>%
    summarise(pred = mean(pred),
              obs = sum(n)/sum(E_wk)) %>%
    ggplot() + 
    geom_line(aes(week, pred, group = name, col = geography), alpha = 0.1, col = "grey") +
    geom_point(aes(week, obs)) + 
    facet_wrap(~geography, ) +
    labs(y = "Predicted versus expected deaths", x = "Week") +
    theme_minimal()
)
# dev.off()


# ---------------------------------------------------------------------------- #
## Scale predictions by death:case ratio ##


lag_rescale <- function(lag,scale, dat){
  dat %>%
    mutate(week = week - lag,
           pred_n_scale = pred_n/scale) -> dat_pred_scale
  return(dat_pred_scale)
}


reconstruct <- function(dat_sims, lag, scale_quants, plot_quants){

dat %>%
  mutate(week = week-lag) %>%
  inner_join(dplyr::select(cases[[1]],lad19cd,geography,week,n), by = c("lad19cd","week","geography"), suffix = c("_d","_c")) %>%
  mutate(CFR_obs = n_d/n_c,
         period = case_when(week < ymd("2020-05-18") ~ 1,
                                   week >= ymd("2020-05-18") ~ 2)) %>%
  group_by(geography) %>%
  mutate(scale = median(CFR_obs)) %>%
  ungroup() -> ratio

# Drop unstable ratios where denominator less than 5 
# ratio$CFR_obs[ratio$n_c < 5] <- NA

ratio %>% 
  summarise(n_d = sum(n_d, na.rm = TRUE),
            n_c = sum(n_c, na.rm = TRUE)) %>%
  print()

hist(ratio$CFR_obs, breaks = 100)
descdist(ratio$CFR_obs[!is.na(ratio$CFR_obs)])

scale <- quantile(ratio$CFR_obs[ratio$period == 2], probs = scale_quants, na.rm = TRUE)

scaled_sims <- bind_rows(lapply(scale, lag_rescale, lag = lag, dat = dat_sims)) %>%
  group_by(lad19cd, week) %>%
  mutate(sim = row_number()) %>%
  ungroup() %>%
  left_join(dplyr::select(cases[[1]], lad19cd, geography, week, n))

scaled_quants <- scaled_sims %>%
  group_by(lad19nm,la, week) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            cases = unique(n, na.rm = TRUE))

png(here::here("figures","compare",paste0("reconstr_la1-4_lag",lag,".png")), height = 1000, width = 1500, res = 150)
print(
  scaled_quants %>%
    filter(la %in% 1:4) %>%
    ggplot(aes(week)) + 
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = med), col = "steelblue") +
    scale_x_date(limits = range(dat_pred_scale$week)) +
    geom_point(aes(y = cases)) +
    facet_wrap(~lad19nm, scales = "free") +
    labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths - four samples LTLAs",
         caption = paste0("Median and ",plot_quants[1]*100,"-",plot_quants[2]*100,
                           "% quantiles over 1000 posterior simulations, scaled by ",
                           min(scale_quants)*100,"% to ", max(scale_quants)*100,
                           "% quantiles of observed CFR.")) +
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
            cases = unique(cases, na.rm = TRUE)) 

png(here::here("figures","compare",paste0("reconstr_geog_lag",lag,".png")), height = 1000, width = 1500, res = 150)
print(
  scaled_quants_geog %>%
    ggplot(aes(week)) + 
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = med), col = "steelblue") +
    scale_x_date(limits = range(dat_sims$week)) +
    geom_point(aes(y = cases)) +
    facet_wrap(~geography, scales = "free_y") +
    labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths - by geography type",
         caption = paste0("Median and ",plot_quants[1]*100,"-",plot_quants[2]*100,
                           "% quantiles over 1000 posterior simulations, scaled by ",
                           min(scale_quants)*100,"% to ", max(scale_quants)*100,
                           "% quantiles of observed CFR.")) +
    theme_minimal()
)
dev.off()

return(scaled_sims)

}

scale_quants <- seq(0.1,0.9,0.1)
plot_quants <- c(0.1,0.9)

scaled_sims7 <- reconstruct(dat_sims, lag = 7, scale_quants, plot_quants)
scaled_sims14 <- reconstruct(dat_sims, lag = 14, scale_quants, plot_quants)


## Total cases
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
  scale_y_continuous(limits = c(0,160000)) +
  theme_minimal() -> reconstr_total7

## Total cases
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
                        min(scale_quants)*100,"% to ", max(scale_quants)*100,
                        "% quantiles of observed CFR.")) +
  scale_y_continuous(limits = c(0,160000)) +
  theme_minimal() -> reconstr_total14

png(here::here("figures","compare","reconstr_totals.png"), height = 800, width = 2000, res = 150)
reconstr_total7 + reconstr_total14
dev.off()


# ---------------------------------------------------------------------------- #
