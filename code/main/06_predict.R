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

