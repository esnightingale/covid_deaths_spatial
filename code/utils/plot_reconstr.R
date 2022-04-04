plot_reconst <- function(agg_sims, lag, sample = NULL, order = T, title = T, caption = T,
                         save = T, figdir = "figures/reconstruct", suffix = "", h = 1400, w = 2400, res = 300, format = "png"){
  
  plot.data <- agg_sims$preds %>%
    mutate(coverage = factor(
      case_when(week < ymd("2020-05-17") ~ "Pre-P2 expansion",
                week >= ymd("2020-05-18") ~ "Post-P2 expansion"),
      levels = c("Pre-P2 expansion","Post-P2 expansion")))
  
  nsims <- agg_sims$nsims
  
  if (agg_sims$type == "la"){
    suffix <- paste0(suffix, "_la")
    plot.data$facet <- plot.data$lad19nm
    if (length(sample) > 0){
      plot.data <- dplyr::filter(plot.data, lad19nm %in% sample) 
      suffix <- paste0(suffix, "samp")
      if(order == F){
        plot.data <- plot.data %>%
          dplyr::mutate(facet = factor(facet, levels = sample))
      }
    }
  }
  
  else if (agg_sims$type == "geog"){
    suffix <- paste0(suffix, "_geog")
    plot.data$facet <- plot.data$geography
  }
  
  else if (agg_sims$type == "total"){
    suffix <- paste0(suffix, "_total")
  }
  
  
  p <- plot.data %>%
    as.data.frame() %>%
    ggplot(aes(week)) + 
    geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
    geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = med), col = "steelblue") +
    geom_point(aes(y = obs, pch = coverage)) +
    scale_shape_manual(values = c(1,19), na.translate = FALSE) +
    labs(x = "Calendar week", 
         y = "Rate per 100,000", 
         shape = "") +
    theme_minimal() 
  
  if (agg_sims$type != "total"){
    p <- p + 
      facet_wrap(~facet, scales = "free") +
      guides(shape = "none")
  }
  
  if (agg_sims$type == "total"){
    p <- p + 
      theme(legend.position = c(0.2,0.8))
  }
  
  if (title == T){
    p <- p + 
      labs(title = "Reconstruction of confirmed cases from COVID-19-related deaths")
  }
  
  if (caption == T){
    p <- p + 
      labs(caption = paste0("Median, ",
                            (plot_quants[3]-plot_quants[2])*100, " and ",
                            (plot_quants[4]-plot_quants[1])*100,
                            "% quantile intervals over ", nsims, " posterior simulations, scaled by ",
                            # samples_cfr, 
                            # "samples from",
                            min(scale_quants)*100,"% to ", max(scale_quants)*100,"% quantiles of",
                            " LTLA-specific CFR post-P2 expansion."))
  }
  
  if (save == T){
    if (format == "pdf"){  
      pdf(here::here(figdir,paste0("reconstr_lag",lag,suffix,".pdf")), height = h, width = w)
      print(p)
      dev.off()
    } else if(format == "png") {
      png(here::here(figdir,paste0("reconstr_lag",lag,suffix,".png")), height = h, width = w, res = res)
      print(p)
      dev.off()
    }
  } else {
    return(p)
  }
}
