
rescale_sims <- function(sims, ratio){
  
  # For each posterior sample, rescale by specified quantiles of the CFR distribution
  scaled_sims <- dplyr::bind_rows(lapply(ratio$quants, lag_rescale, lag = ratio$lag, sims = sims_long))
  
  scaled_sims$scale_quant <- group_indices(scaled_sims, scale)
  
  scaled_sims %>%
    # scaled_sims <- lag_rescale(lag = lag, ratiodist = CFR_obs, sims = sims) %>%
    dplyr::group_by(lad19nm, week) %>%
    dplyr::mutate(sim = row_number(),
                  lag = paste(ratio$lag/7,"weeks")) %>% 
    dplyr::ungroup() %>%
    dplyr::left_join(dplyr::select(cases, lad19nm, geography, week, n), by = c("lad19nm","geography","week")) -> scaled_sims_cases 
 
  return(scaled_sims_cases)
  
}
