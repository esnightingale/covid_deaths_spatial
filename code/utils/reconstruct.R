reconstruct <- function(sims_long, cases, lag, 
                        plot_quants = c(0.01,0.25,0.75,0.99),
                        scale_quants = c(0.25,0.75),
                        plot = TRUE,
                        suffix = ""){
  
  cfr_out <- get_cfr(sims_long, cases, lag, plot = T)
  rescaled <- rescale_sims(cfr_out, scale_quants) 
  
  # Aggregate simulations to different spatial scales
  agg_total <- agg_sims(rescaled, type = "total")
  agg_geog <- agg_sims(rescaled, type = "geog")
  agg_la <- agg_sims(rescaled, type = "la")
  
  # Plot reconstructed quantiles of cases alongside observed
  if (plot == TRUE){
    plot_reconst(agg_total, lag, suffix = suffix)
    plot_reconst(agg_geog, lag, suffix = suffix)
    plot_reconst(agg_la, lag, h = 50, w = 50, format = "pdf", suffix = "wtd")
  }
  
  return(list(total = agg_total, geog = agg_geog, la = agg_la, sims = rescaled, lag = lag))
  
}

