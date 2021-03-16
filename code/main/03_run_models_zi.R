################################################################################
# Description: Run models
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

measure <- "deaths"
expected <- "E"
wave <- 1

list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data","expanded",paste0(measure,".rds")))

dat <- dat_all[[wave]]
period <- dat_all$breaks[[wave]]

# dat$n[dat$n == 0] <- NA

# dat$n[dat$wk_since_first < 0] <- NA

################################################################################
# PRIOR SPECIFICATION
################################################################################

# Dispersion in RR (n/E) - variance >2x mean
var(dat$n/dat$E_wk, na.rm = TRUE)/mean(dat$n/dat$E_wk, na.rm = TRUE)

# SD of RR over time
  dat %>%
    group_by(w) %>%
    summarise(n = sum(n, na.rm = TRUE), E = sum(E_wk)) %>%
    mutate(rate = n/E) %>%
    pull(rate) %>%
    var(na.rm = TRUE) %>%
    sqrt() -> sd_time
  
# SD of RR over space
  dat %>%
    group_by(lad19cd) %>%
    summarise(n = sum(n, na.rm = TRUE), E = unique(E_wk)) %>%
    mutate(rate = n/E) %>%
    pull(rate) %>%
    var(na.rm = TRUE) %>%
    sqrt() -> sd_space

# => Less variation over space than over time

# Undajusted for other covariates and structure so take these as upper bounds:
prior.prec.tp <- list(prec = list(prior = "pc.prec",
                               param = c(sd_time/0.31, 0.01)))
prior.prec.sp <- list(prec = list(prior = "pc.prec",
                               param = c(sd_space/0.31, 0.01),
                               initial = log(0.000001), fixed = TRUE))

# Sample and plot priors on SD scale
prior.samp.tp <- inla.pc.rprec(1000,u = sd_time/0.31, alpha = 0.01)
hist(1/sqrt(prior.samp.tp), breaks = 100)

prior.samp.sp <- inla.pc.rprec(1000,u = sd_space/0.31, alpha = 0.01)
hist(1/sqrt(prior.samp.sp), breaks = 100)

################################################################################
# FITTING
################################################################################


## Base model: No spatial effects, two temporal RWs, independent of geography
f_base <- n ~ 
  IMD_quint + prop_minority +
  # RW2 on epidemic week 
  f(wk_since_first, model = "rw2",
    # hyper = list(prec = prior.prec.tp),
    scale.model = T) +
  # RW1 on week overall
  f(w, model = "rw1",
    # hyper = list(prec = prior.prec.tp), 
    scale.model = T) 

## No spatial effect, geography-dependent RW
f_base_geog <- n ~ 
  IMD_quint + prop_minority +
  # RW2 on epidemic week, replicated by geography
  f(wk_since_first, model = "rw2",
    # hyper = list(prec = prior.prec.tp),
    replicate = geog,
    scale.model = T) +
  f(w, model = "rw1",
    # hyper = list(prec = prior.prec.tp), 
    scale.model = T) 

## IID spatial
f_iid <- n ~ 
  IMD_quint + prop_minority + 
  f(w, model = "rw1",
    # hyper = list(prec = prior.prec.tp),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # hyper = list(prec = prior.prec.tp),
    scale.model = T) +
 f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = prior.prec.sp
    ))

## IID spatial with geography-dependent RW
f_iid_geog <- n ~ 
  IMD_quint + prop_minority + 
  f(w, model = "rw1",
    # hyper = list(prec = prior.prec.tp),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # hyper = list(prec = prior.prec.tp),
    replicate = geog,
    scale.model = T) +
 f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = prior.prec.sp
    ))

## BYM spatial
f_bym <- n ~ 
  IMD_quint + prop_minority + 
  f(w, model = "rw1",
    # hyper = list(prec = prior.prec.tp),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # hyper = list(prec = prior.prec.tp),
    scale.model = T) +
 f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = prior.prec.sp) 
    )

## BYM spatial effect, geography-dependent RW
f_bym_geog <- n ~ 
  IMD_quint + prop_minority + 
  f(w, model = "rw1",
    # hyper = list(prec = prior.prec.tp),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # hyper = list(prec = prior.prec.tp),
    replicate = geog,
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = prior.prec.sp) 
  )

# ## Compare likelihoods on base model
likelihoods <- list(zip0 = "zeroinflatedpoisson0", zinb0 = "zeroinflatednbinomial0", zinb1 = "zeroinflatednbinomial1") #zip0 = "zeroinflatedpoisson0",  zinb0 = "zeroinflatednbinomial0",
fits_zi <- lapply(likelihoods, function(fam) fit_mod(f_bym_geog, dat0, expected, family = fam))
 

dat0 <- dat
dat0$n[is.na(dat0$n)] <- 0
fit_zip0 <- fit_mod(f_iid, dat0, expected, "zeroinflatedpoisson0")
fit_nb <- fit_mod(f_iid, dat0, expected, "nbinomial")
fit_zinb0 <- fit_mod(f_iid, dat0, expected, "zeroinflatednbinomial0")
 
# fits <- list(nb = fit_nb, nb2 = fit_nb2)
# 
# formulae <- list(rw1 = f_base_rw, rw2 = f_base_rw2, ar = f_base_ar) 
# fits <- lapply(formulae, fit_mod, dat2, expected = expected, family = "nbinomial")

# saveRDS(fits, file = here::here("output/expanded_data","fits_base_comp_llh.rds"))


# Check zero-inflation in final model
fit_zinb <- fit_mod(f_bym_geog, dat, expected = expected, family = "zeroinflatednbinomial0")
fit_zip <- fit_mod(f_bym_geog, dat, expected = expected, family = "zeroinflatedpoisson0")

model_comp(list(nb = fit_final, zip = fit_zip, zip1 = fit_zip1, zinb = fit_zinb))
# Model      DIC     WAIC logs diff_WAIC diff_DIC diff_logs
# 1    nb 24930.55 24946.20   NA      0.00     0.00        NA
# 2  zinb 25985.02 26001.61   NA   1055.41  1054.47        NA
# 3   zip 26691.32 26999.58   NA   2053.38  1760.77        NA

fits_zi <- list(nb = fit_final, zip = fit_zip, zip1 = fit_zip1, zinb = fit_zinb)
saveRDS(fits, file = here::here("output",sprintf("fits_final_zi_%s_%s.rds",measure, wave)))

# ## Draw posterior samples ##
samples <- lapply(fits, inla.posterior.sample,n = 1000)
saveRDS(samples, file =  here::here("output/expanded_data",sprintf("samples_%s_%s.rds",measure, wave)))

# Just final model
# samples <- inla.posterior.sample(fit, n = 1000)
# saveRDS(samples, file = here::here("output",sprintf("samples_final_%s_%s.rds",measure, wave)))
