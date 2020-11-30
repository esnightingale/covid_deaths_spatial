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
library(INLA)
library(spdep)
library(sf)

source(here("code","functions.R"))

dir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## Shapefiles
regions <- readRDS(paste0(dir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(dir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
dat <- readRDS(paste0(dir, "Deaths/dat_alt.rds")) %>%
  # filter(wod > ymd("2020-04-20")) # four weeks post lockdown
  filter(wod > ymd("2020-08-02")) %>% # August onwards
  mutate(pop_dens = pop_dens_total) 

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


f_bym_geog <- n ~ 
  IMD + prop_minority + log(pop_dens) + prop_kw + 
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

fit <- fit_mod(f_bym_geog, dat)
saveRDS(fit, file = here("output","fit_secondwave.rds"))

## Draw posterior samples ##
sample <- inla.posterior.sample(fit, n = 1000)
saveRDS(sample, file = here("output","sample_secondwave.rds"))


################################################################################
# OUTPUT
################################################################################

summary(fit)

get_resid <- function(fit){
  dat %>% 
    mutate(mu = fit$summary.fitted.values$mean,
           sigma2 = mu*(1 + mu/fit$summary.hyperpar[1,"mean"]),
           resid = (SMR-mu)/sqrt(sigma2),
           SE = mu^2) %>%
    pull(resid)
}

resids <- get_resid(fit) 
dat_resid <- mutate(dat, resids = resids) 

png(here("figures","secondwave","map_resids.png"), height = 800, width = 1200, res = 150)
# tiff(filename = "./figures/final_altdata/map_resids.tif", height = 600, width = 1000)
dat_resid %>%
  group_by(lad19cd) %>%
  summarise(value = mean(resids)) %>%
  left_join(regions) %>%
  basic_map(fill = "value") +
  # scale_fill_viridis_c(trans = "log10") +
  labs(title = "Mean Squared Error per local authority")
dev.off()

