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

list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

measure <- "deaths"
wave <- 1

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))
border <- st_union(regions)

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data",paste0(measure,".rds")))

dat <- dat_all[[wave]] 
period <- dat_all$breaks[[wave]]

linelist <- readRDS(paste0(datadir, sprintf("linelist_%s.rds",measure))) 

# covid_deaths_raw <- readRDS(paste0(datadir, "Deaths/covid_deaths_raw.rds"))
# covid_deaths_E <- readRDS(paste0(datadir, "Deaths/covid_deaths_E.rds"))

weekrange <- seq(min(dat$w), max(dat$w))

# Fitted models and posterior samples
fits <- readRDS(file = here::here("output/expanded_data",
                                  sprintf("fits_%s_%s.rds",measure, wave)))
fit_final <- fits[[6]]
samples <- readRDS(file = here::here("output/expanded_data",
                                     sprintf("samples_%s_%s.rds",measure, wave)))
samples_final <- samples[[6]]

# ---------------------------------------------------------------------------- #
# MODEL COMPARISON TABLE

model_comp <- data.frame(WAIC = sapply(fits, function(m){return(m$waic$waic)}),
                         logs = sapply(fits, function(m){return(-mean(log(m$cpo$cpo)))})) %>% 
  rownames_to_column(var = "Model")

model_comp %>%
  mutate(diff_WAIC = WAIC - min(WAIC),
         diff_logs = logs - min(logs)) %>%
  arrange(diff_WAIC) %>%
  mutate(across(-Model, function(x) round(x,2))) -> model_comp

model_comp

saveRDS(model_comp, here::here("output/expanded_data",paste0("model_comp_",measure,".rds")))

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
  facet_wrap(~name, ncol = 4) +
  labs(title = "Average pearson residual per local authority")
dev.off()

dat_resid %>%  pivot_longer(cols = base:BYM_geog) %>% group_by(name) %>% summarise(mean_resid = mean(value)) -> mean_resid
png(here::here("figures",measure,"map_resids_6mods.png"), height = 1000, width = 1600, res = 150)
# tiff(filename = "./figures/final_altdata/map_resids.tif", height = 600, width = 1000)
dat_resid %>%
  # dplyr::select(-BYM_geog_nocovs) %>%
  # pivot_longer(cols = base:BYM_geog) %>% 
  # mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
  pivot_longer(cols = base:BYM_geog) %>% 
  group_by(lad19cd, name) %>%
  summarise(value = mean(value)) %>%
  left_join(mean_resid) %>%
  mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
  group_by(name) %>%
  mutate(name2 = as.factor(paste0(name," (",round(mean_resid,4),")"))) %>% #View()
  left_join(regions) %>%
  basic_map(fill = "value") +
  geom_sf(data = border, aes(fill = NULL), alpha = 0, lwd = 0.5, col = "grey") +
  scale_fill_gradient2(midpoint = 0)+
  # scale_fill_viridis_c(trans = "log2") +
  facet_wrap(~name2) +
  labs(title = "Average pearson residual per local authority", subtitle = "Mean residual over all LTLAs given in brackets")
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
  facet_wrap(~name2, ncol = 4) +
  labs(title = "Average log score of prediction per local authority",
       subtitle = "Mean score over all LTLAs given in brackets")
dev.off()

png(here::here("figures",measure,"map_logs_6mods.png"), height = 1000, width = 1600, res = 150)
dat_logs %>%
  # dplyr::select(-BYM_geog_nocovs) %>%
  # pivot_longer(cols = base:BYM_geog) %>%
  pivot_longer(cols = base:BYM_geog) %>% 
  group_by(lad19cd, name) %>%
  summarise(value = -mean(value)) %>%
  # mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>% 
  mutate(name = factor(name, levels = c("base","base_geog","iid","iid_geog", "BYM","BYM_geog"), labels = c("Temporal only (A)", "Geography-specific temporal (B)","A + IID spatial", "B + IID spatial","A + BYM spatial","B + BYM spatial"))) %>%
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

  la_samp_list <- list(random = sample(dat$la,9),
                       range_rates = unique(dat$la[dat$lad19nm %in% c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")]),
                       positive_lag = unique(dat$la[dat$lad19nm %in% c("Leeds", "Bradford") | dat$geography == "London Borough"]))

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
  mutate(pred = exp(value))

png(here::here("figures",measure,"map_pred_quants.png"), height = 1000, width = 1600, res = 150)
print(
  dat_pred %>%
    group_by(lad19cd, name) %>%
    # average over weeks
    summarise(value = mean(pred)) %>%
    # average over samples
    group_by(lad19cd) %>%
    summarise(q50 = median(value),
              q01 = quantile(value, 0.01),
              q99 = quantile(value, 0.99)) %>% 
    pivot_longer(cols = contains("q")) %>% 
    left_join(regions) %>%
    basic_map(fill = "value", rate1e5 = TRUE) +
    facet_wrap(~name) +
    labs(title = "Average predicted deaths relative to expected per LTLA",
         subtitle = "Quantiles of 1000 posterior samples")
)
dev.off()

png(here::here("figures",measure,"preds_geog.png"), height = 1000, width = 1200, res = 150)
print(
  dat_pred %>%
    group_by(week, name, geography) %>%
    summarise(pred = mean(pred),
              obs = sum(n)/sum(E_wk)) %>%
    ggplot() + 
    geom_line(aes(week, pred, group = name, col = geography), alpha = 0.1, col = "grey") +
    geom_point(aes(week, obs)) + 
    facet_wrap(~geography, ) +
    labs(y = "Predicted versus expected deaths", x = "Week") +
    theme_minimal()
)
dev.off()


  
  la_samp_list <- list(random = sample(dat$la,9),
                       range_rates = unique(dat$la[dat$lad19nm %in% c("Liverpool", "Bromley","Bedford", "Allerdale","Wigan","Epping Forest")]),
                       positive_lag = unique(dat$la[dat$lad19nm %in% c("Leeds", "Bradford") | dat$geography == "London Borough"]))
  
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
  
  
  png(here::here("figures",measure,"fitted_RR_CrI.png"), height = 1000, width = 1500, res = 150)
  regions %>%
    full_join(dat_la) %>%
    basic_map(fill = "value") +
    facet_wrap(~name) +
    scale_fill_gradient2(midpoint = 1) +
    labs(title = "Fitted relative risk with 95% credible interval")
  dev.off()
  

## Fitted temporal random effects


  dat_pred %>%
    group_by(week, name) %>%
    summarise(pred = sum(pred_n),
              n = sum(n),
              pop = sum(la_pop)) %>% 
    ggplot() + 
    geom_line(aes(week, pred*1e5/pop, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n*1e5/pop)) + 
    scale_x_date(limits = as.Date(c('2020-06-28','2020-11-15'))) +
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
    scale_x_date(limits = as.Date(c('2020-06-28','2020-11-15'))) +
    labs(x = "Calendar week", y = "Rate per 100,000", title = "Total fit over time, by calendar week and geography", subtitle = "Observed rates shown in white, with 1000 posterior samples", col = "Geography") +
    theme(legend.position = c(0.2,0.7))  -> plot_fit_time_geog
  
  png(here::here("figures",measure,"temp_fit_wave2.png"), height = 1000, width = 1500, res = 150)
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
    # filter(Model == "Wk since first" & w>=0) %>%
    ggplot(aes(w, trend, col = geography)) +
    geom_line() +
    geom_hline(yintercept = 0, col = "grey", lty = "dashed") +
    facet_grid(rows = vars(Model), scales = "free") + 
    labs(x = "", col = "Geography", y = "Trend") +
    theme_minimal() -> plot_rw
  
  png(here::here("figures",measure,"temp_re.png"), height = 1000, width = 1500, res = 150)
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
  
  png(here::here("figures",measure,"death_spatial_re.png"), height = 800, width = 1200, res = 150)
  map_sp_re
  dev.off()
  
  
## Highest IID effects
  sp_re %>%
    filter(name == "IID") %>%
    slice_max(order_by = abs(value), n = 9) %>%
    pull(lad19cd) %>%
    unique() -> hi_IID
  
  png(here::here("figures",measure,"ltla_hi_IID.png"), height = 800, width = 800)
  plot_la_samp(dat_pred, dat$la[dat$lad19cd %in% hi_IID])
  dev.off()
  
## Lowest IID effects
  sp_re %>%
    filter(name == "IID") %>%
    slice_min(order_by = abs(value), n = 9) %>%
    pull(lad19cd) %>%
    unique() -> lo_IID
  
  png(here::here("figures",measure,"ltla_lo_IID.png"), height = 800, width = 800)
  plot_la_samp(dat_pred, dat$la[dat$lad19cd %in% lo_IID])
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

