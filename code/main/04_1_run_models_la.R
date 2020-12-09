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
library(lubridate)
library(INLA)
library(spdep)
library(sf)
library(here)

measure <- "deaths"
wave <- 2

source(here::here("code","main","functions.R"))

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(datadir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data",paste0(measure,".rds")))

dat <- dat_all[[wave]]
period <- dat_all$breaks[[wave]]

# dat <- readRDS(paste0(datadir, "Deaths/dat_deaths_2020-01-01_2020-06-30.rds"))


################################################################################
# PRIOR SPECIFICATION
################################################################################

# Dispersion in n
sqrt(var(dat$n))
mean(dat$n)

# Overall SD in SIR
sqrt(var(dat$SIR))

# SD over time
dat %>%
  group_by(w) %>%
  summarise(n = sum(n), E = sum(E_wk)) %>%
  mutate(SIR = n/E) %>%
  pull(SIR) %>%
  var() %>%
  sqrt() -> sd_time
sd_time

# SD over space
dat %>%
  group_by(lad19cd) %>%
  summarise(n = sum(n), E = unique(E)) %>%
  mutate(SIR = n/E) %>%
  pull(SIR) %>%
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


## Base model: No spatial effects, temporal RW, independent of geography
f_base <- n ~ 
   IMD + prop_minority + #log(pop_dens) + 
   #tb_inc + #cv_mort + can_mort +
  f(w, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),
    scale.model = T) 
  # f(wk_since_first, model = "rw2",
  #   values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
  #   hyper = list(prec = prior.prec.tp),
  #   scale.model = T)

## No spatial effect, geography-dependent RW
f_base_la <- n ~ 
   IMD + prop_minority + #log(pop_dens) + 
  # f(w, model = "rw1",
  #   hyper = list(prec = prior.prec.tp),
  #   values = seq(min(dat$w),max(dat$w)), 
  #   scale.model = T) +
  f(w, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = la,
    values = seq(min(dat$w),max(dat$w)),
    scale.model = T) 


## IID spatial
f_iid <- n ~ 
   IMD + prop_minority + # log(pop_dens) +
   #prop_kw +
   #tb_inc + #cv_mort + can_mort +
  f(w, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  # f(wk_since_first, model = "rw2",
  #   values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
  #   hyper = list(prec = prior.prec.tp),
  #   scale.model = T) +
 f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = prior.prec.sp
    ))

## IID spatial with sp-tp interaction
f_iid_int <- n ~ 
   IMD + prop_minority + #log(pop_dens) + 
 f(la, model = "iid", 
    constr = T,
   group = w,
   control.group = list(model = "rw2", values = seq(min(dat$w),max(dat$w))),
    hyper=list(
      prec = prior.prec.sp
    ))

## BYM spatial
f_bym <- n ~ 
   IMD + prop_minority + #log(pop_dens) +
  f(w, model = "rw2",
    values = seq(min(dat$w),max(dat$w)),
    hyper = list(prec = prior.prec.tp),
    scale.model = T) +
 f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    hyper=list(
      phi =list(param = c(0.5, 2/3)),
      prec = prior.prec.sp) 
    )

## BYM spatial effect, sptp interaction
f_bym_int <- n ~ 
   IMD + prop_minority + #log(pop_dens) + 
 f(la, model = "bym2", graph = g,
    scale.model = T,
    constr = T,
    group = w,
    control.group = list(model = "rw2", values = seq(min(dat$w),max(dat$w))),
    hyper=list(
      phi =list(param =c(0.5, 2/3)),
      prec = prior.prec.sp) 
    )

## BYM spatial effect, geography-dependent RW, no covariates
# f_bym_geog_nocovs <- n ~
#   f(w, model = "rw1",
#     hyper = list(prec = prior.prec.tp),
#     values = seq(min(dat$w),max(dat$w)),   
#     scale.model = T) +
#   f(wk_since_first, model = "rw2",
#     hyper = list(prec = prior.prec.tp),
#     replicate = geog,
#     values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
#     scale.model = T) +
#   f(la, model = "bym2", graph = g,
#     scale.model = T,
#     constr = T,
#     hyper=list(
#       phi =list(param =c(0.5, 2/3)),
#       prec = prior.prec.sp) 
#   )

## Fit all models ##
formulae <- list(base = f_base, base_la = f_base_la, iid = f_iid, iid_int = f_iid_int, BYM = f_bym, BYM_int = f_bym_int)
fits <- lapply(formulae, fit_mod, dat)

saveRDS(fits, file = here::here("output",sprintf("fits_%s_%s_altmods.rds",measure, wave)))

# ## Draw posterior samples ##
samples <- lapply(fits, inla.posterior.sample,n = 1000)
saveRDS(samples, file =  here::here("output",sprintf("samples_%s_%s_altmods.rds",measure, wave)))



# Fit interaction model


f_temp_by_la <- n ~ 
  # IMD + prop_minority + #log(pop_dens) + 
  f(w, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    values = seq(min(dat$w),max(dat$w)),
    scale.model = T,
    group = la, 
    control.group = list(model = "besag", 
                         graph = g, 
                         hyper=list(prec = prior.prec.sp),
                         scale.model = T)) 

fit_temp_by_la <- fit_mod(f_temp_by_la, dat_c)
saveRDS(fit_temp_by_la, here::here("output","fit_cases_1_temp_by_la.rds"))
saveRDS(samples_temp_by_la, here::here("output","samples_cases_1_temp_by_la.rds"))