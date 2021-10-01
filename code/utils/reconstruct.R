reconstruct <- function(cases, lag, 
                        plot_quants = c(0.01,0.25,0.75,0.99),
                        scale_quants = c(0.25,0.75),
                        plot = TRUE,
                        suffix = ""){
  
  gc() 
  
  file.cfr <- here::here("output","reconstruct",paste0("cfr_lag",lag,".rds"))
  if (!file.exists(file.cfr)){
    
    sims_long <- as.data.frame(readRDS(here::here("output","predict","sims_long_avgcov.rds")))

    get_cfr(sims_long, cases, lag, plot = T)
    
    rm(list = c("sims_long","regions","regions.df","boundary"))
  }
  
  cfr_out <- readRDS(file.cfr)
  
  rescaled <- rescale_sims(cfr_out, scale_quants) 
  
  # Aggregate simulations to different spatial scales
  agg_total <- agg_sims(rescaled, type = "total")
  agg_geog <- agg_sims(rescaled, type = "geog")
  agg_la <- agg_sims(rescaled, type = "la")
  
  # Plot reconstructed quantiles of cases alongside observed
  if (plot == TRUE){
    plot_reconst(agg_total, lag, suffix = suffix)
    plot_reconst(agg_geog, lag, suffix = suffix)
    plot_reconst(agg_la, lag, h = 50, w = 50, format = "pdf", suffix = suffix)
  }
  
  return(list(total = agg_total, geog = agg_geog, la = agg_la, sims = rescaled, lag = lag))
  
}

