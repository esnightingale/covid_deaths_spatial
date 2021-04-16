################################################################################
# Description: Fit and draw samples from prediction model, for fitted values
# and predicted with average covariates.
# 
# 
# Author: Emily S Nightingale
# Date created: 18/03/2021
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

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data",paste0(measure,".rds")))

dat <- dat_all[[wave]] 
dat$n[dat$wk_since_first < 0] <- NA

n.la <- n_distinct(dat$la)
period <- dat_all$breaks[[wave]]
nsims <- 1000

# Fitted models and posterior samples
fits <- readRDS(file = here::here("output",
                                  sprintf("fits_%s_%s.rds",measure, wave)))
fit_final <- fits[[6]]
samples <- readRDS(file = here::here("output",
                                     sprintf("samples_%s_%s.rds",measure, wave)))
samples_final <- samples[[6]]

# ---------------------------------------------------------------------------- #
######################
#       PRIORS       #
######################

# Specify priors as before, based on real data

# SD over time
dat %>%
  group_by(w) %>%
  summarise(n = sum(n, na.rm = TRUE), E = sum(E_wk)) %>%
  mutate(rate = n/E) %>%
  pull(rate) %>%
  var() %>%
  sqrt() -> sd_time

# SD over space
dat %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n, na.rm = TRUE), E = unique(E_wk)) %>%
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

# samples_final <- inla.posterior.sample(fit_final, n = nsims)

sims <- as.data.table(bind_cols(lapply(samples_final, get_preds, dat)))

setDT(dat)
dat_sims <-  dplyr::bind_cols(dat[,.(geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n)], 
                              sims)

dat_sims_long <- reshape2::melt(
  dat_sims, 
  id.vars = 1:9)

dat_sims_long <- dat_sims_long[, pred := exp(value)]
dat_sims_long <- dat_sims_long[, pred_n := exp(value)*E_wk]

dat_plot_tot <- dat_sims_long[,.(pred_n = sum(pred_n), 
                                  obs = sum(n, na.rm = T), 
                                  pop = sum(la_pop)), by = .(week, variable)] 

dat_plot_geog <- dat_sims_long[,.(pred_n = sum(pred_n), 
                                  obs = sum(n, na.rm = T), 
                                  pop = sum(la_pop)), by = .(week, variable, geography)] 
dat_plot_geog <- dat_plot_geog[, group := paste(variable, geography)]

ggplot(setDT(dat_plot_tot)) + 
  geom_line(aes(week, pred_n*1e5/pop, group = variable), alpha = 0.1, col = "grey") +
  geom_point(aes(week, obs*1e5/pop)) + 
  labs(x = "Calendar week", y = "Rate per 100,000", title = "Total fit over time, by calendar week", subtitle = paste0("Observed rates shown in black, with ", nsims, " posterior samples in grey")) -> plot_fit_time

ggplot(setDT(dat_plot_geog)) + 
  geom_line(aes(week, pred_n*1e5/pop, group = group, col = geography), alpha = 0.1) +
  geom_point(aes(week, obs*1e5/pop, col = geography), pch = 21, fill = "white") + 
  labs(x = "Calendar week", y = "Rate per 100,000", title = "Total fit over time, by calendar week and geography", subtitle = paste0("Observed rates shown in white, with ", nsims, " posterior samples"), col = "Geography") +
  theme(legend.position = c(0.2,0.7))  -> plot_fit_time_geog

png(here::here("figures",measure,"fit_total_geog.png"), height = 1000, width = 1500, res = 200)
plot_fit_time / plot_fit_time_geog
dev.off()

agg_sims <- dat_sims_long[,.(q05 = quantile(pred_n, 0.05),
                             q25 = quantile(pred_n, 0.25),
                             q50 = quantile(pred_n, 0.5),
                             q75 = quantile(pred_n, 0.75),
                             q95 = quantile(pred_n, 0.95),
                             obs = mean(n)),
                          by = .(la, lad19cd, lad19nm, la_pop, geography, week)]

calc_inc <- function(x) x*1e5/agg_sims$la_pop
agg_sims_plot <- mutate(agg_sims, across(q05:obs, calc_inc))

la_samp <- sample(dat$lad19cd, 6)
ggplot(as.data.frame(agg_sims), aes(x = week)) + #[lad19cd %in% la_samp]
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = q50), col = "steelblue") +
  geom_point(aes(y = obs), cex = 0.5) + 
  facet_wrap(~lad19nm, scales = "free") +
  labs(x = "Calendar week", 
       y = "Rate per 100,000", 
       title = "Total fit over time, by calendar week", 
       subtitle = paste0("Observed rates shown in black, with 50-90% quantiles over ", nsims, " posterior samples")) -> plot_fit_la

png(here::here("figures",measure,"fit_lasamp.png"), height = 1000, width = 1500, res = 200)
plot_fit_la
dev.off()

pdf(here::here("figures",measure,"fit_all_ltlas.pdf"), height = 25, width = 25)
plot_fit_la
dev.off()

# ---------------------------------------------------------------------------- #

######################################
#  PREDICTION AT AVERAGE COVARIATES  #
######################################

median(dat$IMD)
levels(dat$IMD_quint)

setDF(dat) %>%
  mutate(prop_minority = median(prop_minority),
         IMD_quint = "(12.3,20]",
         E_wk = E_wk_unstrat,
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
samples_pred <- inla.posterior.sample(n = nsims, fit_pred)

saveRDS(list(fit = fit_pred, samples = samples_pred, dat = dat_pred), file = here::here("output","fit_samples_avgcov.rds"))
