################################################################################
# Rerun final model using pipeline data from the beginning of the second wave
# and assess evidence for the same spatial structure.
################################################################################

# ---------------------------------------------------------------------------- #
# SETUP
# ---------------------------------------------------------------------------- #

library(reportfactory)
library(tidyverse)
library(linelist)
library(cyphr)
library(INLA)
library(ggplot2)
library(rgdal)
library(spdep)
library(sf)
library(ggpubr)
library(patchwork)
library(viridis)
library(ggspatial)
library(spatial)
library(splines)
library(MASS)
library(gstat)
library(mgcv)

source("./code/functions.R")

dir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## Shapefiles
regions <- readRDS(paste0(dir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

regions.df <- st_drop_geometry(regions)

# Neighbourhood graph
g <- inla.read.graph(filename = paste0(dir,"maps/regions_eng.adj"))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
dat <- readRDS(paste0(dir, "Deaths/dat_alt.rds")) %>%
  rename(pop_dens = pop_dens_total)

# ---------------------------------------------------------------------------- #
# DESCRIPTIVE
# ---------------------------------------------------------------------------- #

pdf("./figures/descriptive_pipeline_total.pdf", height = 7, width = 9)
dat %>%
  group_by(lad19cd) %>%
  mutate(peak = w[which.max(n)]) %>% 
  arrange(peak) %>%
  ggplot(aes(w, n, group = lad19cd, col = lad19cd)) +
  geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
  geom_line(alpha = 0.3) +
  labs(title = "Total COVID19-related deaths in England, by week of death"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  scale_colour_viridis_d(option = "cividis") +
  theme_minimal() -> tot_by_wk

dat %>%
  group_by(lad19cd) %>%
  mutate(peak = wk_since_first[which.max(n)]) %>% 
  arrange(peak) %>%
  ggplot(aes(wk_since_first, n, group = lad19cd, col = lad19cd)) +
  geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
  geom_line(alpha = 0.3) +
  labs(title = "Total COVID19-related deaths in England, by week since first LA death"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  scale_colour_viridis_d(option = "cividis") +
  theme_minimal() -> tot_by_epiwk

tot_by_wk 
tot_by_epiwk

dat %>%
  group_by(lad19cd) %>%
  summarise(first = unique(first)) %>% 
  full_join(regions) %>%
  basic_map(fill = "first") +
  scale_fill_viridis(option = "plasma") +
  ggtitle("Week of first death") -> first_map

dat %>%
  group_by(lad19cd) %>%
  mutate(peak = w[which.max(n)]) %>% 
  arrange(peak) %>%
  full_join(regions) %>%
  basic_map(fill = "peak") +
  scale_fill_viridis(option = "plasma") +
  ggtitle("Week of peak number of deaths") -> peak_map

first_map + peak_map

dat %>%
  ggplot(aes(w, n, group = lad19cd, col = lad19cd)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  labs(title = "Total COVID19-related deaths in England, by week of death"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  theme_minimal() -> E_by_wk

dat %>%
  ggplot(aes(wk_since_first, n, group = lad19cd, col = lad19cd)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  labs(title = "Total COVID19-related deaths in England, by week since first LA death"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  guides(col = F) +
  theme_minimal() -> E_by_epiwk

E_by_wk 
E_by_epiwk


dat %>%
  group_by(w,geography) %>%
  summarise(n = sum(n, na.rm= T)) %>%
  ggplot(aes(w, n, group = geography, col = geography)) +
  geom_line() +
  geom_point() +
  labs(title = "Total COVID19-related deaths in England, by geography type and week of death",
       x = "Calendar week",
       colour = "Geography type"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  theme_minimal() -> E_by_wk

dat %>%
  group_by(wk_since_first,geography) %>%
  summarise(n = sum(n, na.rm= T)) %>%
  ggplot(aes(wk_since_first, n, group = geography, col = geography)) +
  geom_line() +
  geom_point() +
  labs(title = "by week since first LA death",
       x = "Weeks since first observed death",
       colour = "Geography type"
       # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
  ) +
  theme_minimal() -> E_by_epiwk

E_by_wk / E_by_epiwk

dev.off()

# ---------------------------------------------------------------------------- #
# Restrict analysis period
# ---------------------------------------------------------------------------- #

# Filter data to after first wave



# ---------------------------------------------------------------------------- #
# PRIOR SETUP
# ---------------------------------------------------------------------------- #

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


# ---------------------------------------------------------------------------- #
# FIT MODEL
# ---------------------------------------------------------------------------- #

## BYM spatial effect, geography-dependent RW
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

fit <- inla(f_bym_geog,
              "nbinomial",
              data= dat,
              E = E_wk,
              # offset = log(la_age_pop),
              control.compute=list(dic=TRUE, 
                                   waic=TRUE, 
                                   cpo = TRUE,
                                   config = TRUE),
              control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
              verbose = T)

saveRDS(fit, file = "./fits/fits_final_pipelinetotal.rds")

## Draw posterior samples ##

samples <- inla.posterior.sample(fit,n = 1000)
saveRDS(samples, file = "./fits/samples_final_pipelinetotal.rds")


# ---------------------------------------------------------------------------- #
# Summarise error
# ---------------------------------------------------------------------------- #

pdf(paste0("./figures/final_altdata/fit_summ_",names(fits)[[fit]],".pdf"))
# pdf(paste0("./figures/0408/fit_summ_","repgeog",".pdf"))

summary(fits[[fit]]) %>% print()

pit_hist(fits[[fit]]) %>% print()


dat$fitted <- fits[[fit]]$summary.fitted.values[,"mean"]
dat$fitted_sd <- fits[[fit]]$summary.fitted.values[,"sd"]

get_error <- function(x){
  error <- (dat$n - (x*dat$E_wk))/sqrt(x*dat$E_wk)
}

# error <- inla.tmarginal(get_error, fit_main$marginals.fitted.values)

dat %>%
  mutate(mu = fitted*E_wk,
         P.resid = (n - mu)/sqrt(mu*(1 + (mu/fits[[fit]]$summary.hyperpar[1,"mean"]))),
         raw.err = (n - mu)) -> dat

print(
  dat%>%
    ggplot(aes(n, fitted*E_wk)) +
    geom_point(alpha = 0.1) + 
    geom_smooth() + 
    geom_hline(yintercept = 0)
) #%>% print()

print(
  dat%>%
    ggplot(aes(n, P.resid)) +
    geom_point(alpha = 0.1) + 
    geom_smooth() + 
    geom_hline(yintercept = 0, col = "red", lty = "dashed") +
    scale_x_continuous(trans = "log") +
    geom_hline(yintercept = 0)
) #%>% print()

print(
  dat %>%
    ggplot(aes(n, raw.err)) +
    geom_point(alpha = 0.1) + 
    geom_smooth() + 
    geom_hline(yintercept = 0, col = "red", lty = "dashed") +
    scale_x_continuous(trans = "log") +
    geom_hline(yintercept = 0)
) #%>% print()

print(
  dat%>%
    ggplot(aes(w, raw.err)) +
    geom_point(alpha = 0.1) + 
    # geom_jitter() + 
    geom_smooth() + 
    geom_hline(yintercept = 0)
) #%>% print()

print(
  dat%>%
    ggplot(aes(wk_since_first, raw.err)) +
    geom_point(alpha = 0.1) + 
    # geom_jitter() + 
    geom_smooth() + 
    geom_hline(yintercept = 0)
) #%>% print()

print(
  dat %>%
    group_by(lad19cd) %>%
    summarise(P.resid = mean(P.resid)) %>%
    left_join(regions) %>%
    basic_map(fill = "P.resid") +
    labs(title = "Mean standardised error (Pearson residual) per local authority")
) #%>% print()

-mean(log(fits[[fit]]$cpo$cpo))
