get_cfr <- function(deaths, cases, lag, scale_quants){
  
  deaths %>%
    mutate(n = replace_na(n, 0),
           week = week-lag) %>%
    full_join(dplyr::select(cases,lad19nm,geography,week,n), by = c("lad19nm","week","geography"), suffix = c("_d","_c")) %>%
    mutate(CFR = n_c/n_d,
           wt = round(n_d*100000/sum(n_d, na.rm = T)),
           period = case_when(week < ymd("2020-05-18") ~ 1,
                              week >= ymd("2020-05-18") ~ 2)) %>%  
    filter(is.finite(CFR) & CFR > 0) %>%
    group_by(lad19nm, geography, week, period, CFR) %>% 
    tidyr::expand(wt = seq(1:wt)) %>%
    ungroup() -> ratio

  ratio %>% 
    summarise(med = median(CFR, na.rm = TRUE),
              q1 = quantile(CFR, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  ratio %>% 
    group_by(period) %>%
    summarise(med = median(CFR, na.rm = TRUE),
              q1 = quantile(CFR, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  
  ratio %>% 
    group_by(geography,period) %>%
    summarise(med = median(CFR, na.rm = TRUE),
              q1 = quantile(CFR, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  # Plot distribution of observed CFR
  CFR_wtd <- ratio$CFR[ratio$period == 2]
  lims <- c(0,quantile(CFR_wtd, 0.999, na.rm = T))
  x <- 0:max(CFR_wtd)
  png(here::here("figures","compare",paste0("cfr_dist_lag",lag,".png")), height = 600, width = 800, res = 150)
  hist(CFR_wtd, breaks = 200, prob = T, xlim = lims, main = "")
  lines(x, EnvStats::demp(x, CFR_wtd), col = "red")
  dev.off()
  
  quants <- quantile(ratio$CFR[ratio$period == 2], probs = scale_quants, na.rm = TRUE)
  quants
  
  return(list(quants = quants, ratio = ratio, lag = lag))
  
}
