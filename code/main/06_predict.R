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

# Neighbourhood graph
g <- inla.read.graph(filename = here::here("data","regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data","aggregated",paste0(measure,".rds")))

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

priors <- set_priors()

# ---------------------------------------------------------------------------- #

######################################
#  PREDICTION AT AVERAGE COVARIATES  #
######################################

median(dat$IMD)
# 16.11325
levels(dat$IMD_quint)
# "(4.51,12.3]" "(12.3,20]"   "(20,27.7]"   "(27.7,35.4]" "(35.4,43.1]"

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
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]]),
    replicate = geog,
    scale.model = T) +
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]]),
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = priors[["space"]]) 
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

saveRDS(list(fit = fit_pred, samples = samples_pred, dat = dat_pred), file = here::here("output","predict","fit_samples_avgcov.rds"))

sims <- as.data.table(bind_cols(lapply(samples_pred, function(s) s$latent[(nrow(dat)+1):(2*nrow(dat))])))

# ----- #

setDT(dat)
dat_sims <-  dplyr::bind_cols(dat[,.(geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n)], 
                              sims)
dat_sims_long <- reshape2::melt(
  dat_sims, 
  id.vars = 1:9)
setDT(dat_sims_long)

dat_sims_long <- dat_sims_long[, pred_n := exp(value)*E_wk]

dat_sims_la <- dat_sims_long[,.(pred_n = sum(pred_n), 
                                 obs = sum(n, na.rm = T)), by = .(lad19nm, variable)] 
dat_plot_la <- dat_sims_la[,.(`Predicted (averaged LTLA covariates)` = median(pred_n), 
                                Observed = median(obs, na.rm = T)), by = .(lad19nm)] 
dat_plot_la <- reshape2::melt(
  dat_plot_la, 
  id.vars = 1)

regions %>%
  full_join(dat_plot_la) %>%
  basic_map(fill = "value", rate1e5 = TRUE, scale = FALSE, plot.border = FALSE) + 
  facet_wrap(~variable) +
  labs(title = "Predicted deaths per 100,000 with averaged LTLA covariate values", subtitle = "Median of 1000 posterior samples.", fill = "Incidence\nper 100,000") +
  theme(legend.position = c(0.05,0.5)) -> pred_map

png(here::here("figures","predict","pred_avgcov_la.png"), height = 1800, width = 2600, res = 300)
pred_map
dev.off()

######################################
#   PREDICTION WITHOUT COVARIATES    #
######################################

f <- n ~  
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]]),
    replicate = geog,
    scale.model = T) +
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]]),
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = priors[["space"]]) 
  )

fit_nocov <- INLA::inla(f,
                       "nbinomial",
                       data = dat,
                       E = E_wk,
                       control.compute=list(dic=TRUE, 
                                            waic=TRUE, 
                                            cpo = TRUE,
                                            config = TRUE),
                       control.predictor = list(compute = TRUE, link = 1),
                       control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
                       verbose = T)
summary(fit_nocov)

# Samples are of *marginal* densities
samples_nocov <- inla.posterior.sample(n = nsims, fit_nocov)

saveRDS(list(fit = fit_nocov, samples = samples_nocov, dat = dat), file = here::here("output","predict","fit_samples_nocov.rds"))

