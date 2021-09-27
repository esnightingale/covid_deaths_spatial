################################################################################
# Description: Setup and run models, draw posterior samples
# 
# Author: Emily S Nightingale
# Date created: 30/09/2020
# 
################################################################################
################################################################################

################################################################################
# SETUP
################################################################################

measure <- "deaths"
wave <- 1
expected <- "E"
nsims <- 1000

# Neighbourhood graph
g <- inla.read.graph(filename = here::here("data","regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data","aggregated",paste0(measure,".rds")))

dat <- dat_all[[wave]]
period <- dat_all$breaks[[wave]]

dat$n[dat$wk_since_first < 0] <- NA

################################################################################
# PRIOR SPECIFICATION
################################################################################

# Dispersion in RR (n/E) - variance ~2x mean
var(dat$n/dat$E_wk, na.rm = TRUE)/mean(dat$n/dat$E_wk, na.rm = TRUE)

priors <- set_priors()

################################################################################
# FITTING
################################################################################

## Base model: No spatial effects, two temporal RWs, independent of geography
f_base <- n ~ 
  IMD_quint + prop_minority +
  # RW2 on epidemic week 
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]]),
    scale.model = T) +
  # RW1 on week overall
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]])) 

## No spatial effect, geography-dependent RW
f_base_geog <- n ~ 
  IMD_quint + prop_minority +
  # RW2 on epidemic week, replicated by geography
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]]),
    replicate = geog) +
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]])) 

## IID spatial
f_iid <- n ~ 
  IMD_quint + prop_minority + 
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]])) +
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]])) +
  f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = priors[["space"]]
    ))

## IID spatial with geography-dependent RW
f_iid_geog <- n ~ 
  IMD_quint + prop_minority + 
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]]),
    replicate = geog) +
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]])) +
  f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = priors[["space"]]
    ))

## BYM spatial
f_bym <- n ~ 
  IMD_quint + prop_minority + 
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]])) +
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]])) +
  f(la, model = "bym2", graph = g,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = priors[["space"]]) 
  )

## BYM spatial effect, geography-dependent RW
f_bym_geog <- n ~ 
  IMD_quint + prop_minority + 
  f(wk_since_first, model = "rw2",
    hyper = list(prec = priors[["time1"]]),
    replicate = geog) +
  f(w, model = "rw1",
    hyper = list(prec = priors[["time2"]])) +
  f(la, model = "bym2", graph = g,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = priors[["space"]]) 
  )

## Fit all models ##
formulae <- list(base = f_base, base_geog = f_base_geog, 
                 iid = f_iid, iid_geog = f_iid_geog, 
                 BYM = f_bym, BYM_geog = f_bym_geog) #, BYM_geog_nocovs = f_bym_geog_nocovs
fits <- lapply(formulae, fit_mod, dat, expected = expected)

# Check failed cpo values
lapply(fits, function(f) summary(f$cpo$failure))

saveRDS(fits, file = here::here("output",sprintf("fits_%s_%s.rds",measure, wave)))

# Just final model
# fit <- fit_mod(f_bym_geog, dat, expected = expected)
# saveRDS(fit, file = here::here("output",sprintf("fit_final_%s_%s.rds",measure, wave)))


# ## Draw posterior samples ##
samples <- lapply(fits, inla.posterior.sample,n = nsims)
saveRDS(samples, file =  here::here("output",sprintf("samples_%s_%s.rds",measure, wave)))

# Just final model
# samples <- inla.posterior.sample(fit, n = 1000)
# saveRDS(samples, file = here::here("output",sprintf("samples_final_%s_%s.rds",measure, wave)))
