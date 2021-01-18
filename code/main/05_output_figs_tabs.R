################################################################################
# Description: Summarise and visualise model output
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
library(patchwork)
library(viridis)
library(INLAutils)
library(cowplot)
library(ggspatial)
library(ggforce)
# devtools::install_github('timcdlucas/INLAutils')

theme_set(theme_minimal())

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

measure <- "deaths"
wave <- 1

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data",paste0(measure,".rds")))

dat <- dat_all[[wave]] 
period <- dat_all$breaks[[wave]]

linelist <- readRDS(paste0(datadir, sprintf("linelist_%s.rds",measure))) 

# %>%
#   filter(wod < ymd("2020-08-02")) # filter to last week in July
# covid_deaths_raw <- readRDS(paste0(datadir, "Deaths/covid_deaths_raw.rds"))
# covid_deaths_E <- readRDS(paste0(datadir, "Deaths/covid_deaths_E.rds"))

weekrange <- seq(min(dat$w), max(dat$w))

# Fitted models and posterior samples
fits <- readRDS(file = here::here("output",
                                  sprintf("fits_%s_%s.rds",measure, wave)))
samples <- readRDS(file = here::here("output",
                                     sprintf("samples_%s_%s.rds",measure, wave)))

# ---------------------------------------------------------------------------- #
# MODEL COMPARISON TABLE

model_comp <- data.frame(WAIC = sapply(fits, function(m){return(m$waic$waic)}),
                         logs = sapply(fits, function(m){return(-mean(log(m$cpo$cpo)))})) %>% 
  rownames_to_column(var = "Model")

model_comp %>%
  mutate(diff_WAIC = WAIC - min(WAIC),
         diff_logs = logs - min(logs)) %>%
  arrange(diff_WAIC) 

# ---------------------------------------------------------------------------- #
# MODEL SUMMARY PLOTS

pdf(here::here("figures",measure,"fits.pdf"), height = 8, width = 12)
for (m in seq_along(fits)){
  plot.model(fits[[m]])
}
dev.off()

# ---------------------------------------------------------------------------- #
# MAP MODEL MSE

get_resid <- function(fit){
  dat %>% 
    mutate(mu = fit$summary.fitted.values$mean*E_wk,
           sigma2 = mu*(1 + mu/fit$summary.hyperpar[1,"mean"]),
           resid = (n-mu)/sqrt(sigma2),
           perc_err = (n-mu)*100/n,
           MSqE = (n-mu)^2,
           logs = log(fit$cpo$cpo)) 
}

resids <- lapply(fits, function(fit) get_resid(fit)$resid) %>%
  bind_cols() %>%
  setNames(unlist(names(fits)))

dat_resid <- bind_cols(dat, resids) 

png(here::here("figures",measure,"map_resids.png"), height = 800, width = 1200, res = 150)
# tiff(filename = "./figures/final_altdata/map_resids.tif", height = 600, width = 1000)
dat_resid %>%
  # dplyr::select(-BYM_geog_nocovs) %>%
  # pivot_longer(cols = base:BYM_geog) %>% 
  # mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
  pivot_longer(cols = base:BYM_geog_nocovs) %>% 
  mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog","BYM_geog_nocovs"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial","B + BYM spatial - covs"))) %>%
  group_by(lad19cd, name) %>%
  summarise(value = mean(value)) %>%
  left_join(regions) %>%
  basic_map(fill = "value") +
  scale_fill_gradient2(midpoint = 0)+
  # scale_fill_viridis_c(trans = "log2") +
  facet_wrap(~name) +
  labs(title = "Average pearson residual per local authority")
dev.off()


# Map log score per LTLA

logs <- lapply(fits, function(fit) get_resid(fit)$logs) %>%
  bind_cols() %>%
  setNames(unlist(names(fits)))

dat_logs <- bind_cols(dat, logs) 

png(here::here("figures",measure,"map_logs.png"), height = 1000, width = 1600, res = 150)
dat_logs %>%
  # dplyr::select(-BYM_geog_nocovs) %>%
  # pivot_longer(cols = base:BYM_geog) %>%
  pivot_longer(cols = base:BYM_geog_nocovs) %>% 
  group_by(lad19cd, name) %>%
  summarise(value = -mean(value)) %>%
  # mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>% 
  mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog","BYM_geog_nocovs"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial","B + BYM spatial - covs"))) %>%
  group_by(name) %>%
  mutate(name2 = as.factor(paste0(name," (",round(mean(value),2),")"))) %>% #View()
  left_join(regions) %>%
  basic_map(fill = "value") +
  # scale_fill_viridis_c(trans = "log2") +
  facet_wrap(~name2) +
  labs(title = "Average log score of prediction per local authority",
       subtitle = "Mean score over all LTLAs given in brackets")
dev.off()

# ---------------------------------------------------------------------------- #
# SUMMARISE POSTERIOR SAMPLES

nval <- nrow(dat)

for (s in seq_along(samples)){
  
  preds <- bind_cols(lapply(samples[[s]], get_preds))
  
  pdf(here::here("figures",measure,paste0("summ_post_", names(fits)[s],".pdf")), height = 8, width = 10)
  
  # par(mfrow = c(2,1))
  # hist(exp(as.matrix(preds))*dat$E_wk, breaks = 30, xlim = c(0,200), prob = T)
  # hist(dat$n, breaks = 30, xlim = c(0,200), prob = T)
  
  dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), preds) %>%
    pivot_longer(cols = -1:-8) %>%
    mutate(pred_n = exp(value)*E_wk)
  
  print(
  dat_pred %>%
    group_by(lad19cd, name) %>%
    # sum over weeks
    summarise(value = sum(pred_n)) %>%
    # average over samples
    group_by(lad19cd) %>%
    summarise(q50 = median(value),
              q25 = quantile(value, 0.25),
              q75 = quantile(value, 0.75)) %>% 
    pivot_longer(cols = contains("q")) %>% 
    left_join(regions) %>%
    basic_map(fill = "value", rate1e5 = TRUE) +
    facet_wrap(~name) +
    labs(title = "Predicted deaths per 100,000: Quantiles of 1000 posterior samples per local authority")
  )
  
  print(
    dat_pred %>%
      group_by(week, name) %>%
      summarise(pred_n = sum(pred_n, na.rm = T),
                n = sum(n, na.rm = T)) %>%
      ggplot() + 
      geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
      geom_point(aes(week, n)) + 
      theme_minimal()
  )
  
  print(
    dat_pred %>%
      group_by(week, name, geography) %>%
      summarise(pred_n = sum(pred_n),
                n = sum(n)) %>%
      ggplot() + 
      geom_line(aes(week, pred_n, group = name, col = geography), alpha = 0.1, col = "grey") +
      geom_point(aes(week, n)) + 
      facet_wrap(~geography) +
      theme_minimal()
  )
  
  # la_samp <- sample(dat$la,9)
  # la_list <- c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")
  la_list <- c("Leeds", "Bradford")
  la_samp <- unique(dat$la[dat$lad19nm %in% la_list | dat$geography == "London Borough"])
  print(
    dat_pred %>%
      filter(la %in% la_samp) %>%
      ggplot() + 
      geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
      geom_point(aes(week, n)) + 
      facet_wrap(~lad19nm) +
      theme_minimal()
  )
  
  
  plot_parm <- function(parm, opt = 1){
    
    if (opt == 1){ d <- fits[[s]]$marginals.fixed[[parm]]
    }else{ d <- fits[[s]]$marginals.hyperpar[[parm]] }
    
    print(
      ggplot(data.frame(inla.smarginal(d)), aes(x, y)) +
        geom_line() +
        geom_vline(xintercept = 0, col = "red", lty = "dashed") +
        labs(title = parm) +
        theme_bw()
    )
  }
  
  lapply(names(fits[[s]]$marginals.fixed), plot_parm, opt = 1)
  
  lapply(names(fits[[s]]$marginals.hyperpar), plot_parm, opt = 2)
  
  dev.off()
  
  print(s)
}

# ---------------------------------------------------------------------------- #
# Print sampled trajectories for each LTLA

preds <- bind_cols(lapply(samples[["BYM_geog"]], get_preds))

pdf(here::here("figures",measure,paste0("summ_post_final_allLTLA.pdf")), height = 8, width = 10)

dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), preds) %>%
  pivot_longer(cols = -1:-8) %>%
  mutate(pred_n = exp(value)*E_wk)

print(
dat_pred %>%
  group_by(lad19cd, name) %>%
  # sum over weeks
  summarise(value = sum(pred_n)) %>%
  # average over samples
  group_by(lad19cd) %>%
  summarise(q50 = median(value),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75)) %>% 
  pivot_longer(cols = contains("q")) %>% 
  left_join(regions) %>%
  basic_map(fill = "value", rate1e5 = TRUE) +
  facet_wrap(~name) +
  labs(title = "Predicted deaths per 100,000: Quantiles of 1000 posterior samples per local authority")
)


# la_samp <- sample(dat$la,9)
# la_list <- c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")
london <- unique(dat$lad19nm[dat$geography == "London Borough"])
la_list <- c("Leeds", "Bradford",sample(london,2))
la_samp <- unique(dat$la[dat$lad19nm %in% la_list])
print(
  dat_pred %>%
    filter(la %in% la_samp) %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n)) + 
    facet_wrap(~lad19nm) +
    theme_minimal()
)

dev.off()

# ---------------------------------------------------------------------------- #
# SUMMARISE COVARIATE EFFECTS

png(here::here("figures",measure,"covariates_final.png"), height = 800, width = 1000, res = 150)
# tiff(filename = "./figures/final_altdata/map_resids.tif", height = 600, width = 1000)

fits[["BYM_geog"]]$summary.fixed[-1,] %>% #rep_BYM
  rownames_to_column(var = "Covariate") %>%
  ggplot(aes(x = Covariate, y = mean,  ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_pointrange() +
  geom_hline(aes(yintercept = 0), lty = "dashed",col = "red") +
  labs(y = "Estimate") 
dev.off()


# ---------------------------------------------------------------------------- #
# Map fitted values over time

