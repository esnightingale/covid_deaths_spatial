get_cfr <- function(deaths, cases, lag, scale_quants, denom_cutoff = 1){
  
  deaths %>%
    mutate(n = replace_na(n, 0),
           week = week-lag) %>%
    full_join(dplyr::select(cases,lad19cd,geography,week,n), by = c("lad19cd","week","geography"), suffix = c("_d","_c")) %>%
    mutate(CFR_obs = n_c/n_d,
           period = case_when(week < ymd("2020-05-18") ~ 1,
                              week >= ymd("2020-05-18") ~ 2),
           ID = row_number()) -> ratio
  
  # Drop unstable ratios where denominator less than cutoff
  ratio$CFR_obs[ratio$n_d < denom_cutoff | ratio$n_c < denom_cutoff] <- NA
  
  ratio %>% 
    summarise(n_d = sum(n_d, na.rm = TRUE),
              n_c = sum(n_c, na.rm = TRUE),
              med = median(CFR_obs, na.rm = TRUE),
              q1 = quantile(CFR_obs, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR_obs, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  ratio %>% 
    group_by(period) %>%
    summarise(n_d = sum(n_d, na.rm = TRUE),
              n_c = sum(n_c, na.rm = TRUE),
              med = median(CFR_obs, na.rm = TRUE),
              q1 = quantile(CFR_obs, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR_obs, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  
  ratio %>% 
    group_by(geography,period) %>%
    summarise(n_d = sum(n_d, na.rm = TRUE),
              n_c = sum(n_c, na.rm = TRUE),
              med = median(CFR_obs, na.rm = TRUE),
              q1 = quantile(CFR_obs, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR_obs, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  # Plot distribution of observed CFR
  CFR_obs <- ratio$CFR_obs[!is.na(ratio$CFR_obs) & ratio$period == 2]
  x <- 0:max(CFR_obs)
  hist(CFR_obs, breaks = 100, prob = T)
  lines(x, EnvStats::demp(x, CFR_obs), col = "red")
  
  quants <- quantile(ratio$CFR_obs[ratio$period == 2], probs = scale_quants, na.rm = TRUE)
  quants
  
  return(list(quants = quants, ratio = ratio))
  
}
