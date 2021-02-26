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
dat_all <- readRDS(here::here("data",paste0(measure,".rds")))

dat <- dat_all[[wave]]
period <- dat_all$breaks[[wave]]

################################################################################
# PRIOR SPECIFICATION
################################################################################

# Dispersion in n
sqrt(var(dat$n))
mean(dat$n)

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
   IMD_quint + prop_minority + #log(pop_dens) +
   #tb_inc + #cv_mort + can_mort +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    # values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    hyper = list(prec = prior.prec.tp),
    scale.model = T)

## No spatial effect, geography-dependent RW
f_base_geog <- n ~ 
  IMD_quint + prop_minority + #log(pop_dens) +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    # values = seq(min(dat$w),max(dat$w)), 
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    # values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    scale.model = T) 


## IID spatial
f_iid <- n ~ 
  IMD_quint + prop_minority + #log(pop_dens) +
   #prop_kw +
   #tb_inc + #cv_mort + can_mort +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    # values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    hyper = list(prec = prior.prec.tp),
    scale.model = T) +
 f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = prior.prec.sp
    ))

## IID spatial with geography-dependent RW
f_iid_geog <- n ~ 
  IMD_quint + prop_minority + #log(pop_dens) +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    # values = seq(min(dat$w),max(dat$w)), 
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    hyper = list(prec = prior.prec.tp),
    replicate = geog,
    # values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
    scale.model = T) +
 f(la, model = "iid", 
    constr = T,
    hyper=list(
      prec = prior.prec.sp
    ))

## BYM spatial
f_bym <- n ~ 
  IMD_quint + prop_minority + #log(pop_dens) +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    # values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
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
  IMD_quint + prop_minority + #log(pop_dens) +
  f(w, model = "rw1",
    hyper = list(prec = prior.prec.tp),
    # values = seq(min(dat$w),max(dat$w)),
    scale.model = T) +
  f(wk_since_first, model = "rw2",
    # values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
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

# f_besag_la <- n ~ 
#   # IMD + prop_minority + log(pop_dens) +
#   # f(w, model = "rw1",
#   #   hyper = list(prec = prior.prec.tp),
#   #   values = seq(min(dat$w),max(dat$w)),
#   #   scale.model = T) +
#   f(wk_since_first, model = "rw2",
#     values = seq(min(dat$wk_since_first),max(dat$wk_since_first)),
#     replicate = la,
#     hyper = list(prec = prior.prec.tp),
#     scale.model = T) +
#  f(la, model = "besag", graph = g,
#     scale.model = T,
#     constr = T,
#     hyper=list(
#       # phi =list(param =c(0.5, 2/3)),
#       prec = prior.prec.sp) 
#     )

# f_iid <- n ~ 
#   f(la, model = "iid", 
#     constr = T,
#     group = w,
#     control.group = list(
#       model = "rw1",
#       hyper = list(prec = prior.prec.tp),
#       scale.model = T),
#     hyper=list(
#       prec = prior.prec.sp
#     ))

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
formulae <- list(base = f_base, base_geog = f_base_geog, iid = f_iid, iid_geog = f_iid_geog, BYM = f_bym, BYM_geog = f_bym_geog) #, BYM_geog_nocovs = f_bym_geog_nocovs
fits <- lapply(formulae, fit_mod, dat, expected = expected)

saveRDS(fits, file = here::here("output/expanded_data",sprintf("fits_%s_%s.rds",measure, wave)))

# Just final model
# fit <- fit_mod(f_bym_geog, dat, expected = expected)
# saveRDS(fit, file = here::here("output",sprintf("fit_final_%s_%s.rds",measure, wave)))


# ## Draw posterior samples ##
samples <- lapply(fits, inla.posterior.sample,n = 1000)
saveRDS(samples, file =  here::here("output/expanded_data",sprintf("samples_%s_%s.rds",measure, wave)))

# Just final model
# samples <- inla.posterior.sample(fit, n = 1000)
# saveRDS(samples, file = here::here("output",sprintf("samples_final_%s_%s.rds",measure, wave)))
