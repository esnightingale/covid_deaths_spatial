################################################################################
# Description: Summarise and visualise death time series per LTLA. Produce 
# descriptive figures for paper.
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
library(lubridate)
library(INLA)
library(spdep)
library(sf)
library(here)

source(here::here("code","functions.R"))

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
dat <- readRDS(paste0(datadir, "Deaths/dat_deaths_2020-01-01_2020-06-30.rds"))

################################################################################
# PRIOR SPECIFICATION
################################################################################

# Dispersion in n
sqrt(var(dat$n))
mean(dat$n)

# Overall SD in SMR
sqrt(var(dat$SMR))

# SD over time
dat %>%
  group_by(w) %>%
  summarise(n = sum(n), E = sum(E_wk)) %>%
  mutate(SMR = n/E) %>%
  pull(SMR) %>%
  var() %>%
  sqrt() -> sd_time
sd_time

# SD over space
dat %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n), E = unique(E)) %>%
  mutate(SMR = n/E) %>%
  pull(SMR) %>%
  var() %>%
  sqrt() -> sd_space
sd_space

# => Less variation over space than over time

# Undajusted for other covariates and structure so take these as upper bounds:
prior.prec.tp <- list(prec = list(prior = "pc.prec",
                               param = c(sd_time/0.31, 0.01)))
prior.prec.sp <- list(prec = list(prior = "pc.prec",
                               param = c(sd_space/0.31, 0.01)))

# Plot priors
plot(inla.pc.dprec(seq(0,100,0.01),u = sd_time/0.31, alpha = 0.01))
plot(inla.pc.dprec(seq(0,100,0.01),u = sd_space/0.31, alpha = 0.01))



################################################################################
# FITTING
################################################################################


## Base model: No spatial effects, two temporal RWs, independent of geography
f_base <- n ~ 
   IMD + prop_minority + log(pop_dens) + 
   #tb_inc + #cv_mort + can_mort +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    hyper = list(prec = prior.prec.tp),
    scale.model = T)

## No spatial effect, geography-dependent RW
f_base_geog <- n ~ 
   IMD + prop_minority + log(pop_dens) + 
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)), 
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    scale.model = T) 


## IID spatial
f_iid <- n ~ 
   IMD + prop_minority + log(pop_dens) +
   #prop_kw +
   #tb_inc + #cv_mort + can_mort +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    hyper = list(prec = prior.prec.tp),
    scale.model = T) +
 f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = prior.prec.sp
    ))

## IID spatial with geography-dependent RW
f_iid_geog <- n ~ 
   IMD + prop_minority + log(pop_dens) + 
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)), 
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    scale.model = T) +
 f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = prior.prec.sp
    ))

## BYM spatial
f_bym <- n ~ 
   IMD + prop_minority + log(pop_dens) +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    hyper = list(prec = prior.prec.tp),
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
   IMD + prop_minority + log(pop_dens) + 
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),   
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    scale.model = T) +
 f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param =c(0.5, 2/3)),
      prec = prior.prec.sp) 
    )

## BYM spatial effect, geography-dependent RW, no covariates
f_bym_geog_nocovs <- n ~
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),   
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    scale.model = T) +
  f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param =c(0.5, 2/3)),
      prec = prior.prec.sp) 
  )


## Fit all models ##
formulae <- list(base = f_base, base_geog = f_base_geog, iid = f_iid, iid_geog = f_iid_geog, BYM = f_bym, BYM_geog = f_bym_geog, BYM_geog_nocovs = f_bym_geog_nocovs)
fits <- lapply(formulae, fit_mod, dat)

saveRDS(fits, file = "output/fits_final_1stwave.rds")

## Draw posterior samples ##
samples <- lapply(fits, inla.posterior.sample,n = 1000)
saveRDS(samples, file = "output/samples_final_1stwave.rds")

