agg_sims <- function(scaled_sims, type = c("la","geog","total")){
  
  if (type == "la") {
    scaled_sims <- group_by(scaled_sims, lad19nm, week) %>%
      rename(obs = n)
  }
  if (type == "geog") {
    scaled_sims <- group_by(scaled_sims, geography, week, sim, scale_quant) %>%
      dplyr::summarise(pred_n_scale = sum(pred_n_scale),
                       obs = sum(n, na.rm = T),
                       lag = unique(lag)) %>%
      group_by(geography, week)
  }
  if (type == "total") {
    scaled_sims <- group_by(scaled_sims, week, sim, scale_quant) %>%
      dplyr::summarise(pred_n_scale = sum(pred_n_scale),
                       obs = sum(n, na.rm = T),
                       lag = unique(lag)) %>%
      group_by(week)
  }
  
  scaled_quants <- scaled_sims %>%
    dplyr::summarise(l1 = quantile(pred_n_scale, plot_quants[1]),
                     l2 = quantile(pred_n_scale, plot_quants[2]),
                     med = quantile(pred_n_scale, 0.5),
                     h1 = quantile(pred_n_scale, plot_quants[3]),
                     h2 = quantile(pred_n_scale, plot_quants[4]),
                     obs = unique(obs),
                     lag = unique(lag)) %>%
    dplyr::ungroup() 
  
  return(list(preds = scaled_quants, type = type))
  
}
