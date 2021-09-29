agg_sims <- function(rescaled, type = c("la","geog","total")){
  
  gc()
  
  if (type == "la") {
    rescaled <- group_by(rescaled, lad19nm, week) %>%
      rename(obs = n_c)
  }
  if (type == "geog") {
    rescaled <- group_by(rescaled, geography, week, sim, scale_quant) %>%
      dplyr::summarise(pred_c = sum(pred_c),
                       obs = sum(n_c, na.rm = T),
                       lag = unique(lag)) %>%
      group_by(geography, week)
  }
  if (type == "total") {
    rescaled <- group_by(rescaled, week, sim, scale_quant) %>%
      dplyr::summarise(pred_c = sum(pred_c),
                       obs = sum(n_c, na.rm = T),
                       lag = unique(lag)) %>%
      group_by(week)
  }
  
  preds <- rescaled %>%
    dplyr::summarise(l1 = quantile(pred_c, plot_quants[1]),
                     l2 = quantile(pred_c, plot_quants[2]),
                     med = quantile(pred_c, 0.5),
                     h1 = quantile(pred_c, plot_quants[3]),
                     h2 = quantile(pred_c, plot_quants[4]),
                     obs = unique(obs),
                     lag = unique(lag)) %>%
    dplyr::ungroup() 
  
  return(list(preds = preds, type = type))
  
}
