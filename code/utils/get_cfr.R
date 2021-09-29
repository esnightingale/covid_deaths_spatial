get_cfr <- function(sims, cases, lag, plot = F){
  
  gc()
  sink(file = here::here(outdir,paste0("calcCFR_lag",lag,".txt")))
  
  nsims <- dplyr::n_distinct(sims$variable)
  
  # Join posterior predictions of deaths (averaging out covariates) with observed cases
  # Inner join to keep only overlapping weeks after lagging
  # Estimate CFR per week/per LTLA as cases/predicted deaths
  # Set CFR as NA if cases or pred deaths = 0
  sims %>%
    as.data.frame() %>%
    mutate(week = week-lag) %>%
    rename(sim = variable) %>%
    inner_join(dplyr::select(cases,lad19nm,geography,week,n), by = c("lad19nm","week","geography"), suffix = c("_d","_c")) %>% 
    mutate(CFR = (n_c/pred_n)*(pred_n > 0),
           period = case_when(week < ymd("2020-05-18") ~ 1,
                              week >= ymd("2020-05-18") ~ 2)) %>%
    mutate(CFR = na_if(CFR, 0)) %>%
    ungroup() -> ratio
  
  # Check output for one sim
  View(filter(ratio, sim == 1))
  
  # Summarise non-missing CFRs overall, by LA and by week
  ratio_nona <- filter(ratio, !is.na(CFR))
  ratio_nona %>% 
    summarise(med = median(CFR),
              q1 = quantile(CFR, p = 0.25),
              q3 = quantile(CFR, p = 0.75)) %>%
    print()
  
  ratio_nona %>% 
    group_by(period) %>%
    summarise(med = median(CFR),
              q1 = quantile(CFR, p = 0.25),
              q3 = quantile(CFR, p = 0.75)) %>%
    print()
  
  ratio_nona %>% 
    group_by(geography,period) %>%
    summarise(med = median(CFR),
              q1 = quantile(CFR, p = 0.25),
              q3 = quantile(CFR, p = 0.75)) %>%
    print()
  
  if (plot == T){
    
    # Ratio per LTLA
    ratio_nona %>%
      group_by(lad19nm, period, sim) %>%
      summarise(n_d = sum(n_d, na.rm = T),
                n_c = sum(n_c, na.rm = T),
                CFR = median(CFR)) %>%
      group_by(lad19nm, period) %>%
      summarise(n_d = median(n_d, na.rm = T),
                n_c = median(n_c, na.rm = T),
                CFR = median(CFR)) %>%
      ungroup() -> ratio_la
  
    print(summary(ratio_la$CFR[ratio_la$period == 1]))
    print(summary(ratio_la$CFR[ratio_la$period == 2]))
  
    regions %>%
      full_join(ratio_la) %>%
      mutate(period = factor(period,labels = c("Pre-pillar 2","Post-pillar 2"))) %>%
      basic_map(fill = "CFR", scale = F) +
      scale_fill_viridis_c(trans = "log2") +
      facet_wrap(~period, ncol = 2) +
      labs(title = "",
           fill = "Ratio",
           caption = paste("Ratios where no. deaths or no. cases = 0 are not calculated")) -> maps
  
    # Ratio per week
    ratio_nona %>%
      group_by(week, period, sim) %>%
      summarise(n_d = sum(n_d, na.rm = T),
                n_c = sum(n_c, na.rm = T),
                CFR = median(CFR)) %>%
      group_by(week, period) %>%
      summarise(med = median(CFR),
                l1 = quantile(CFR, 0.01),
                l2 = quantile(CFR, 0.25),
                h2 = quantile(CFR, 0.75),
                h1 = quantile(CFR, 0.99)) %>%
      ungroup() -> ratio_wk
    
    ratio_wk %>%
      ggplot(aes(x = week)) +
        geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
        geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
        geom_line(aes(y = med), col = "steelblue") +
        scale_y_continuous(trans = "log2") +
        geom_vline(xintercept = ymd("2020-05-18"), col = "indianred") +
        annotate("text", x = ymd("2020-05-19"), y = 0.25, label = "P2 available to all symptomatic cases", cex = 2, hjust = "left") +
        geom_vline(xintercept = ymd("2020-04-15"), col = "indianred") +
        annotate("text", x = ymd("2020-04-16"), y = 0.15, label = "P2 available to care home residents and staff", cex = 2, hjust = "left") +
        labs(x = "", y = "Observed CFR") +
        scale_x_date(limits = c(ymd("2020-02-05"),ymd("2020-06-17"))) -> time
    
    png(here::here(figdir,paste0("map_obs_cfr_",lag,".png")), height = 1000, width = 1500, res = 200)
    print(maps)
    dev.off()
  
    png(here::here(figdir,paste0("time_obs_cfr_",lag,".png")), height = 1000, width = 1500, res = 200)
    print(time)
    dev.off()
  }
 
  out <- list(ratio = ratio, lag = lag, nsims = nsims)
  saveRDS(out, here::here("output","reconstruct",paste0("cfr_lag",lag,".rds")))

  sink()
}


