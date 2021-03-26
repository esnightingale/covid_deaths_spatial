reconstruct <- function(deaths, cases, sims, lag, 
                        scale_quants = seq(0.1,0.9,0.1), 
                        plot_quants = c(0.01,0.25,0.75,0.99),
                        plot = TRUE,
                        suffix = ""){
  
  ratio <- get_cfr(deaths, cases, lag, scale_quants)
  rescaled <- rescale_sims(sims, ratio)
  
  # Aggregate simulations to different spatial scales
  agg_total <- agg_sims(rescaled, type = "total")
  agg_geog <- agg_sims(rescaled, type = "geog")
  agg_la <- agg_sims(rescaled, type = "la")
  
  # Plot reconstructed quantiles of cases alongside observed
  if (plot == TRUE){
    plot_reconst(agg_total, lag, h = 6, w = 10, suffix = suffix)
    plot_reconst(agg_geog, lag, h = 8, w = 10, suffix = suffix)
    plot_reconst(agg_la, lag, suffix = suffix)
  }
  
  return(list(total = agg_total, geog = agg_geog, la = agg_la, sims = rescaled))
  
}

