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

figdir <- "figures/model comparison"
measure <- "deaths"
wave <- 1

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data","aggregated","deaths.rds"))

dat <- dat_all[[wave]] 
dat$n[dat$wk_since_first < 0] <- NA
period <- dat_all$breaks[[wave]]

dat_tot <- dat %>%
  group_by(la, lad19cd, lad19nm, geography) %>%
  summarise(n = sum(n, na.rm = T),
            E = unique(E))

weekrange <- seq(min(dat$w), max(dat$w))
weekseq <- seq(min(dat$week), max(dat$week), by = "week")

# Fitted models and posterior samples
fits <- readRDS(file = here::here("output",sprintf("fits_%s_%s.rds","deaths", wave)))
samples <- readRDS(file = here::here("output",sprintf("samples_%s_%s.rds","deaths", wave)))
nsims <- length(samples[[1]])

# ---------------------------------------------------------------------------- #
# CHECK CPO

pits <- lapply(fits, pit_hist, bins= 50)

png(here::here(figdir,"pit.png"), height = 1000, width = 1600, res = 150)
(pits[[1]] + pits[[2]] + pits[[3]])/(pits[[4]] + pits[[5]] + pits[[6]])
dev.off()

# ---------------------------------------------------------------------------- #
# MODEL COMPARISON TABLE

table <- model_comp(fits)
table

write.csv(table, here::here("output",paste0("model_comp.csv")), row.names = F)

# ---------------------------------------------------------------------------- #
# FIXED EFFECTS

coeffs <- lapply(fits, get_coeffs)
coeffs

saveRDS(coeffs, here::here("output",paste0("model_coeffs.rds")))
        
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

dat_resid %>%  
  pivot_longer(cols = base:BYM_geog) %>% 
  group_by(name) %>% 
  summarise(median = median(value, na.rm = T),
            mean = mean(value, na.rm = T)) -> avg_resid

png(here::here(figdir,"map_resids_6mods.png"), height = 1000, width = 1600, res = 150)
dat_resid %>%
  pivot_longer(cols = base:BYM_geog) %>% 
  group_by(lad19cd, name) %>%
  summarise(value = mean(value, na.rm = T)) %>%
  left_join(avg_resid) %>%
  mutate(name = factor(name, 
                       levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), 
                       labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
  group_by(name) %>%
  mutate(name2 = as.factor(paste0(name," (",round(mean,4),")"))) %>% # View()
  left_join(regions) %>%
  basic_map(fill = "value", scale = F, plot.border = T) +
  scale_fill_gradient2(midpoint = 0)+
  # scale_fill_viridis_c(trans = "log2") +
  facet_wrap(~name2) +
  theme(legend.position = c(0,0.5)) +
  labs(title = "Mean pearson residual per local authority", 
       subtitle = "Mean over all LTLAs given in brackets")
dev.off()


# Map log score per LTLA

logs <- lapply(fits, function(fit) get_resid(fit)$logs) %>%
  bind_cols() %>%
  setNames(unlist(names(fits)))

dat_logs <- bind_cols(dat, logs) 

png(here::here(figdir,"map_logs_6mods.png"), height = 1000, width = 1600, res = 150)
dat_logs %>%
  pivot_longer(cols = base:BYM_geog) %>% 
  group_by(lad19cd, name) %>%
  summarise(value = -mean(value, na.rm = T)) %>%
  mutate(name = factor(name, 
                       levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), 
                       labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
  group_by(name) %>%
  mutate(name2 = as.factor(paste0(name," (",round(mean(value),2),")"))) %>% #View()
  left_join(regions) %>%
  basic_map(fill = "value", scale = F) +
  facet_wrap(~name2) +
  theme(legend.position = c(0,0.5)) +
  labs(title = "Log score of predictions per local authority",
       subtitle = "Mean over all LTLAs given in brackets")
dev.off()

# ---------------------------------------------------------------------------- #
# SUMMARISE POSTERIOR SAMPLES - ALL MODELS

nval <- nrow(dat)
la_samp <- sample(dat_tot$la,9)

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

for (s in seq_along(samples)){
  
  preds <- bind_cols(lapply(samples[[s]], get_preds, dat))
  
  pdf(here::here("figures","posterior summaries",paste0("summ_post_", names(fits)[s],".pdf")), height = 8, width = 10)
  
  dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), 
                        preds) %>%
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
                q01 = quantile(value, 0.01),
                q99 = quantile(value, 0.99)) %>% 
      pivot_longer(cols = contains("q")) %>% 
      left_join(regions) %>%
      basic_map(fill = "value", rate1e5 = TRUE) +
      facet_wrap(~name) +
      labs(title = "Predicted deaths per 100,000: Quantiles of 1,000 posterior samples per local authority")
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

  plot_la_samp(dat_pred, la_samp)
  
  lapply(names(fits[[s]]$marginals.fixed), plot_parm, opt = 1)
  
  lapply(names(fits[[s]]$marginals.hyperpar), plot_parm, opt = 2)
  
  dev.off()
  
  print(s)
  
}


################################################################################
################################################################################

