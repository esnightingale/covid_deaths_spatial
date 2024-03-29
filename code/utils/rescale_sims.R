rescale_sims <- function(cfr_out, scale_quants = c(0.25, 0.75)){ 

  gc()
  
  # Average over weekly CFRs where cases and deaths > 0
  cfr_out$ratio %>%
    filter(period == 2) %>%
    group_by(lad19nm, sim) %>%
    summarise(q50 = quantile(CFR, probs = 0.5, na.rm = TRUE),
              qlow = quantile(CFR, probs = scale_quants[1], na.rm = TRUE),
              qhi = quantile(CFR, probs = scale_quants[2], na.rm = TRUE)) %>%
    pivot_longer(starts_with("q"), names_to = "scale_quant", values_to = "cfr") -> CFR_la
  
  # CFR_la consists of one CFR estimate per posterior sample for each LTLA
  # Want to apply these estimates to scale up posterior predicted deaths over whole period
  
  # Join back to full CFR dataset and multiply posterior predicted deaths by LA-averaged CFR
  cfr_out$ratio %>%
    full_join(CFR_la) %>% #slice(1:48940) %>% View()
    mutate(pred_c = pred_n*cfr) %>%
    group_by(lad19nm, week, scale_quant) %>%
    dplyr::mutate(sim = row_number(),
                  lag = paste(cfr_out$lag/7,"weeks")) %>% 
    dplyr::ungroup() -> rescaled
 
  saveRDS(rescaled, here::here("output","reconstruct",paste0("rescaled_lag",lag,".rds")))
  
  return(rescaled)
  
}
