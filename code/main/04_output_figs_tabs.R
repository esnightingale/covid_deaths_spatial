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
list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

measure <- "deaths"
wave <- 1

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))
border <- st_union(regions)

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data","expanded",paste0(measure,".rds")))

dat <- dat_all[[wave]] 
dat$n[dat$wk_since_first < 0] <- NA
period <- dat_all$breaks[[wave]]

dat_tot <- dat %>%
  group_by(la, lad19cd, lad19nm, geography) %>%
  summarise(n = sum(n, na.rm = T),
            E = unique(E))

linelist <- readRDS(paste0(datadir, sprintf("linelist_%s.rds",measure))) 

# covid_deaths_raw <- readRDS(paste0(datadir, "Deaths/covid_deaths_raw.rds"))
# covid_deaths_E <- readRDS(paste0(datadir, "Deaths/covid_deaths_E.rds"))

weekrange <- seq(min(dat$w), max(dat$w))
weekseq <- seq(min(dat$week), max(dat$week), by = "week")

# Fitted models and posterior samples
fits <- readRDS(file = here::here("output/expanded_data",
                                  sprintf("fits_%s_%s.rds",measure, wave)))
fit_final <- fits[[6]]
samples <- readRDS(file = here::here("output/expanded_data",
                                     sprintf("samples_%s_%s.rds",measure, wave)))
samples_final <- samples[[6]]


# ---------------------------------------------------------------------------- #
# CHECK CPO

pits <- lapply(fits, pit_hist, bins= 50)

png(here::here("figures",measure,"expanded_data","pit.png"), height = 1000, width = 1600, res = 150)
(pits[[1]] + pits[[2]] + pits[[3]])/(pits[[4]] + pits[[5]] + pits[[6]])
dev.off()

# ---------------------------------------------------------------------------- #
# MODEL COMPARISON TABLE

table <- model_comp(fits)
table

saveRDS(table, here::here("output/expanded_data",paste0("model_comp_",measure,".rds")))

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

dat_resid %>%  pivot_longer(cols = base:BYM_geog) %>% group_by(name) %>% summarise(median = median(value, na.rm = T), mean = mean(value, na.rm = T)) -> avg_resid
png(here::here("figures",measure,"expanded_data","map_resids_6mods_mean.png"), height = 1000, width = 1600, res = 150)
# tiff(filename = "./figures/final_altdata/map_resids.tif", height = 600, width = 1000)
dat_resid %>%
  pivot_longer(cols = base:BYM_geog) %>% 
  group_by(lad19cd, name) %>%
  summarise(value = mean(value, na.rm = T)) %>%
  left_join(avg_resid) %>%
  mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
  group_by(name) %>%
  mutate(name2 = as.factor(paste0(name," (",round(mean,4),")"))) %>% # View()
  left_join(regions) %>%
  basic_map(fill = "value", scale = F, plot.border = T) +
  scale_fill_gradient2(midpoint = 0)+
  # scale_fill_viridis_c(trans = "log2") +
  facet_wrap(~name2) +
  theme(legend.position = c(0,0.5)) +
  labs(title = "Mean pearson residual per local authority", subtitle = "Mean over all LTLAs given in brackets")
dev.off()


# Map log score per LTLA

logs <- lapply(fits, function(fit) get_resid(fit)$logs) %>%
  bind_cols() %>%
  setNames(unlist(names(fits)))

dat_logs <- bind_cols(dat, logs) 

png(here::here("figures",measure,"expanded_data","map_logs_6mods_mean.png"), height = 1000, width = 1600, res = 150)
dat_logs %>%
  pivot_longer(cols = base:BYM_geog) %>% 
  group_by(lad19cd, name) %>%
  summarise(value = -mean(value, na.rm = T)) %>%
  mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
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

for (s in seq_along(samples)){
  
  preds <- bind_cols(lapply(samples[[s]], get_preds))
  
  pdf(here::here("figures",measure,paste0("summ_post_", names(fits)[s],".pdf")), height = 8, width = 10)
  
  # par(mfrow = c(2,1))
  # hist(exp(as.matrix(preds))*dat$E_wk, breaks = 30, xlim = c(0,200), prob = T)
  # hist(dat$n, breaks = 30, xlim = c(0,200), prob = T)
  
  dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), preds) %>%
    pivot_longer(cols = -1:-8) %>%
    mutate(pred_n = exp(value)*E_wk)
  
  # png(here::here("figures",measure,"map_pred_quants.png"), height = 1000, width = 1600, res = 150)
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
    labs(title = "Predicted deaths per 100,000: Quantiles of 1000 posterior samples per local authority")
  )
# dev.off()

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

  la_samp_list <- list(random = sample(dat_tot$la,9),
                       range_rates = dat_tot$la[dat_tot$lad19nm %in% c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")],
                       positive_lag = dat_tot$la[dat_tot$lad19nm %in% c("Leeds", "Bradford") | dat$geography == "London Borough"])

  for (l in seq_along(la_samp_list)){
    plot_la_samp(dat_pred, la_samp_list[[l]])
  }
  
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
# SUMMARISE POSTERIOR SAMPLES - FINAL MODEL

nval <- nrow(dat)

preds <- bind_cols(lapply(samples_final, get_preds, dat))

dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), preds) %>%
  pivot_longer(cols = -1:-8) %>%
  mutate(pred = exp(value),
         pred_n = pred*E_wk)


dat_pred %>%
  # average over weeks
  group_by(lad19cd, name) %>%
  summarise(value = mean(pred)) %>%
  # summarise over samples
  group_by(lad19cd) %>%
  summarise(q50 = median(value),
            q01 = quantile(value, 0.01),
            q99 = quantile(value, 0.99)) -> quants_by_la

png(here::here("figures",measure,"expanded_data","map_pred_quants.png"), height = 1000, width = 1600, res = 150)
print(
  quants_by_la %>%
    pivot_longer(cols = contains("q")) %>% 
    left_join(regions) %>%
    basic_map(fill = "value", scale = F) +
    facet_wrap(~name) +
    labs(title = "Average predicted deaths relative to expected per LTLA",
         subtitle = "Quantiles of 1000 posterior samples")
)
dev.off()

png(here::here("figures",measure,"expanded_data","preds_geog.png"), height = 1000, width = 1200, res = 150)
print(
  dat_pred %>%
    group_by(week, name, geography) %>%
    summarise(pred = mean(pred),
              obs = sum(n)/sum(E_wk)) %>%
    ggplot() + 
    geom_line(aes(week, pred, group = name, col = geography), alpha = 0.1, col = "grey") +
    geom_point(aes(week, obs)) +
    geom_hline(yintercept = 1, col = "red",lty = "dashed") +
    facet_wrap(~geography, ) +
    labs(y = "Predicted versus expected deaths", x = "Week") +
    theme_minimal()
)
dev.off()

# LTLAs with low and high average relative risks
dat_tot %>% 
  mutate(crude_RR = n/E) %>%
  right_join(filter(quants_by_la, q50 > 2)) -> high_RR

dat_tot %>% 
  mutate(crude_RR = n/E) %>%
  right_join(filter(quants_by_la, q50 < 0.25)) -> low_RR

pdf(here::here("figures",measure,"expanded_data","high_low_RR.pdf"), width = 8, height = 6)
plot_la_samp(dat_pred, high_RR$la, pred = "pred_n", obs = "n") 
plot_la_samp(dat_pred, low_RR$la, pred = "pred_n", obs = "n") 
dev.off()

la_samp_list <- list(random = sample(dat_tot$la,9),
                     range_rates = dat_tot$la[dat_tot$lad19nm %in% c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")],
                     high_RR = high_RR$la)
  
  for (l in seq_along(la_samp_list)){
    png(here::here("figures",measure,paste0("preds_la",l,".png")), height = 1000, width = 1200, res = 150)
    plot_la_samp(dat_pred, la_samp_list[[1]])
    dev.off()
  }
  
  
  pdf(here::here("figures",measure,"parms_fixed.pdf"))
  lapply(names(fit_final_d$marginals.fixed), plot_parm, opt = 1)
  dev.off()
  pdf(here::here("figures",measure,"parms_hyper.pdf"))
  lapply(names(fit_final_c$marginals.hyperpar), plot_parm, opt = 2)
  dev.off()
  

  

## Fitted relative risks by LTLA

  dat$RR <- fit_final$summary.fitted.values[, "mean"]
  dat$LL <- fit_final$summary.fitted.values[, "0.025quant"]
  dat$UL <- fit_final$summary.fitted.values[, "0.975quant"]
  
  
  dat_la <- dat %>%
    group_by(lad19cd) %>%
    summarise(RR = mean(RR, na.rm = TRUE), LL = mean(LL, na.rm = TRUE), UL = mean(UL, na.rm = TRUE)) %>%
    pivot_longer(cols = c("LL","RR","UL"))
  
  
  png(here::here("figures",measure,"expanded_data","fitted_RR_CrI.png"), height = 1000, width = 1500, res = 150)
  regions %>%
    full_join(dat_la) %>%
    basic_map(fill = "value") +
    facet_wrap(~name) +
    scale_fill_gradient2(midpoint = 1) +
    labs(title = "Fitted relative risk with 95% credible interval")
  dev.off()
  

## Fitted temporal effects

  dat_pred %>%
    group_by(week, name) %>%
    summarise(pred = sum(pred_n),
              n = sum(n),
              pop = sum(la_pop)) %>%
    ggplot() + 
    geom_line(aes(week, pred*1e5/pop, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n*1e5/pop)) + 
    labs(x = "Calendar week", y = "Rate per 100,000", title = "Total fit over time, by calendar week", subtitle = "Observed rates shown in black, with 1000 posterior samples in grey") -> plot_fit_time
  
  dat_pred %>%
    group_by(week, name, geography) %>%
    summarise(pred = sum(pred_n),
              n = sum(n),
              pop = sum(la_pop)) %>%
    mutate(group = paste0(name, geography)) -> by_samp_geog
  by_samp_geog %>%
    group_by(week, geography) %>%
    summarise(pred = mean(pred),
              n = mean(n),
              pop = mean(pop)) -> by_geog
  

  ggplot(by_samp_geog) + 
    geom_line(aes(week, pred*1e5/pop, group = group, col = geography), alpha = 0.1) +
    geom_point(dat = by_geog, aes(week, n*1e5/pop, col = geography), pch = 21, fill = "white") + 
    labs(x = "Calendar week", y = "Rate per 100,000", title = "Total fit over time, by calendar week and geography", subtitle = "Observed rates shown in white, with 1000 posterior samples", col = "Geography") +
    theme(legend.position = c(0.2,0.7))  -> plot_fit_time_geog
  
  png(here::here("figures",measure,"expanded_data","temp_fit_wave2.png"), height = 1000, width = 1500, res = 150)
  plot_fit_time / plot_fit_time_geog
  dev.off()
  
  # length(fit_final$marginals.random$w)
  # length(fit_final$marginals.random$wk_since_first)
  
  dat_w <- data.frame(w = weekrange,
                      trend=fit_final$summary.random$w$mean,
                      Model = "Calendar week") %>%
    bind_rows(bind_cols(expand.grid(w = seq(min(dat$wk_since_first),max(dat$wk_since_first)), geography = sort(unique(dat$geography))),
                        data.frame(trend = fit_final$summary.random$wk_since_first$mean,
                                   Model = "Week since first death"))) %>%
    mutate(geography = factor(replace_na(as.character(geography), "All"),
                              levels = c("London Borough","Metropolitan District","Non-metropolitan District","Unitary Authority","All")))
  
  dat_w %>%
    ggplot(aes(w, trend, col = geography)) +
    geom_line() +
    geom_hline(yintercept = 0, col = "grey", lty = "dashed") +
    facet_grid(rows = vars(Model), scales = "free") + 
    labs(x = "", col = "Geography", y = "Trend") +
    theme_minimal() -> plot_rw
  
  png(here::here("figures",measure,"expanded_data","temp_re.png"), height = 1000, width = 1500, res = 150)
  plot_rw
  dev.off()
  
## Fitted spatial random effects

  sp_re <- data.frame(lad19cd = unique(dat$lad19cd), 
                      tau = fit_final$summary.hyperpar["Precision for la","mean"],
                      phi = fit_final$summary.hyperpar["Phi for la","mean"],
                      Total = fit_final$summary.random$la[1:nrow(regions),"mean"], 
                      Besag = fit_final$summary.random$la[(nrow(regions)+1):(2*nrow(regions)),"mean"]) %>% 
    mutate(IID = (Total*sqrt(tau) - Besag*sqrt(phi))/sqrt(1-phi)) %>%
    pivot_longer(c("Total","Besag","IID"))
  
  regions %>%
    full_join(sp_re) %>%
    basic_map(fill = "value") +
    geom_sf(data = border, aes(fill = NULL), alpha = 0, lwd = 0.5, col = "grey") +
    facet_wrap(~name) +
    labs(subtitle = "Decomposition of fitted spatial random effects", fill = "") +
    scale_fill_gradient2() -> map_sp_re
  
  png(here::here("figures",measure,"expanded_data","death_spatial_re.png"), height = 800, width = 1200, res = 150)
  map_sp_re
  dev.off()
  
  
## Highest IID effects
  sp_re %>%
    filter(name == "IID") %>%
    slice_max(order_by = abs(value), n = 9) %>%
    pull(lad19cd) %>%
    unique() -> hi_IID
  
  png(here::here("figures",measure,"expanded_data","ltla_hi_IID.png"), height = 800, width = 800)
  plot_la_samp(dat_pred, dat$la[dat$lad19cd %in% hi_IID])
  dev.off()
  
## Lowest IID effects
  sp_re %>%
    filter(name == "IID") %>%
    slice_min(order_by = abs(value), n = 9) %>%
    pull(lad19cd) %>%
    unique() -> low_IID
  
  png(here::here("figures",measure,"expanded_data","ltla_lo_IID.png"), height = 800, width = 800)
  plot_la_samp(dat_pred, dat$la[dat$lad19cd %in% low_IID])
  dev.off()
 
# ---------------------------------------------------------------------------- #
# OVER/UNDER-ESTIMATED PEAK
  
  
get_fit_peak <- function(dat, samples, offset){
    
    fitted <- bind_cols(lapply(samples, function(s) exp(s$latent[1:nrow(dat)]))) 
    
    fitted %>%
      mutate(across(everything(), function(x) x*offset)) %>%
      rowwise() %>%
      mutate(med = median(c_across(everything())),
             lo = quantile(c_across(everything()), p = 0.05),
             hi = quantile(c_across(everything()), p = 0.95)) %>%
      ungroup() -> fitted
    
    dat <- cbind(dat, fitted[,c("med","lo","hi")]) %>%
      mutate(obs_in_QR = between(n, lo, hi))
    
    dat %>%
      group_by(lad19cd) %>%
      summarise(fit_peak = week[which.max(med)],
                fit_peak_n = max(med),
                lo = lo[which.max(n)],
                hi = hi[which.max(n)],
                obs_peak = week[which.max(n)],
                obs_peak_n = max(n, na.rm = T),
                peak_n_err = (obs_peak_n - fit_peak_n)*100/obs_peak_n,
                peak_err = obs_peak - fit_peak,
                obs_in_QR = obs_in_QR[which.max(n)]) -> dat_peak
    
    return(list(dat_peak,dat))
  }
  
peak <- get_fit_peak(dat, samples_final, offset = dat$E_wk)

png(here::here("figures",measure,"peak_error_hist.png"), height = 1000, width = 1200, res = 150)
peak[[1]] %>%
  ggplot(aes(peak_n_err)) +
  geom_histogram(bins = 30, fill = "steelblue", col= "black") +
  labs(x = "Percentage error", title = "Percentage error in predicted peak deaths (median across 1000 posterior samples)", subtitle = paste0("Median error at +",round(median(peak[[1]]$peak_n_err),2))) + 
  geom_vline(xintercept = 0, lty = "dashed", lwd = 1)+
  geom_vline(aes(xintercept = median(peak_n_err)), col = "red", lwd = 1) 
dev.off()

png(here::here("figures",measure,"peak_error_map.png"), height = 1000, width = 1200, res = 150)
peak[[1]] %>%
  full_join(regions) %>%
  basic_map(fill = "peak_n_err") +
  labs(fill = "Percentage error", title = "Percentage error in predicted peak deaths (median across 1000 posterior samples)", subtitle = paste0("Median error at +",round(median(peak[[1]]$peak_n_err),2))) 
dev.off()


peak[[1]] %>%
  mutate(err_cat = case_when(peak_n_err < -50 ~ "< -50%",
                             (peak_n_err < -30 & peak_n_err >= -50) ~ "< -30%",
                             abs(peak_n_err) <= 30 ~ "within 30%",
                             (peak_n_err > 30 & peak_n_err <= 50) ~ "> 30%",
                             peak_n_err > 50 ~ "> 50%")) %>%
  mutate(err_cat = factor(err_cat, levels = c("< -50%","< -30%","within 30%","> 30%","> 50%"))) -> peak_err


png(here::here("figures",measure,"peak_error_cat_map.png"), height = 1000, width = 1200, res = 150)
regions %>%
  full_join(peak_err) %>% 
  basic_map(fill = "err_cat") +
  scale_fill_viridis_d() +
  labs(fill = "", title = "Percentage error in predicted peak deaths (median across 1000 posterior samples)")
dev.off()

# LTLAs with error > 50%

peak[[1]] %>%
  filter(peak_n_err < -30) %>%
  pull(lad19cd) -> overpredict

peak[[1]] %>%
  filter(peak_n_err > 30) %>%
  pull(lad19cd)-> underpredict

over <- unique(dat$la[dat$lad19cd %in% overpredict])
under <- unique(dat$la[dat$lad19cd %in% underpredict])


pdf(here::here("figures",measure,"overestimated_ltlas.pdf"), height = 8, width = 8)
plot_la_samp(dat_pred, over)
dev.off()

pdf(here::here("figures",measure,"underestimated_ltlas.pdf"), height = 15, width = 15)
plot_la_samp(dat_pred, under)
dev.off()

# ---------------------------------------------------------------------------- #
# Print sampled trajectories for each LTLA

preds <- bind_cols(lapply(samples[["BYM_geog"]], get_preds))

pdf(here::here("figures",measure,"summ_post_final_allLTLA.pdf"), height = 8, width = 10)

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

# IMD quintiles
post_fixed <- lapply(fit_final$marginals.fixed[[-6]], function(x) inla.tmarginal(exp, x)) 
lapply(post_fixed, inla.zmarginal)

# Prop minority (rescale to %)
inla.zmarginal(inla.tmarginal(function(x) exp(x/100), fit_final$marginals.fixed[[6]]))

png(here::here("figures",measure,"expanded_data","covariates_final.png"), height = 800, width = 1000, res = 150)
# tiff(filename = "./figures/final_altdata/map_resids.tif", height = 600, width = 1000)

fit_final$summary.fixed[-1,] %>% #rep_BYM
  rownames_to_column(var = "Covariate") %>%
  ggplot(aes(x = Covariate, y = mean,  ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_pointrange() +
  geom_hline(aes(yintercept = 0), lty = "dashed",col = "red") +
  labs(y = "Estimate") 
dev.off()


# ---------------------------------------------------------------------------- #
# Map fitted values over time

