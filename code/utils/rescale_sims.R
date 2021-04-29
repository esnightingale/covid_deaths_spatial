rescale_sims <- function(cfr_out, scale_quants = c(0.25, 0.75)){ #

  cfr_out$ratio %>%
    filter(period == 2) %>%
    group_by(lad19nm, sim) %>%
    summarise(q50 = quantile(CFR, probs = 0.5),
              qlow = quantile(CFR, probs = scale_quants[1]),
              qhi = quantile(CFR, probs = scale_quants[2])) %>%
    pivot_longer(starts_with("q"), names_to = "scale_quant", values_to = "cfr") -> CFR_la
  
  # CFR_la consists of one CFR estimate per posterior sample for each LTLA
  
  cfr_out$ratio %>%
    full_join(CFR_la) %>% #slice(1:48940) %>% View()
    mutate(pred_c = pred_n*cfr) %>%
    group_by(lad19nm, week, scale_quant) %>%
    dplyr::mutate(sim = row_number(),
                  lag = paste(cfr_out$lag/7,"weeks")) %>% 
    dplyr::ungroup() -> rescaled

  return(rescaled)
  
}
