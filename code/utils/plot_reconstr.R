plot_reconst <- function(agg_sims, lag, sample = NA, save = T, suffix = "", h = 50, w = 50, order = T, format = "pdf"){
  
  plot.data <- agg_sims$preds
  
  if (agg_sims$type == "la"){
    suffix <- paste0(suffix, "_la")
    plot.data$facet <- plot.data$lad19nm
    if (!is.na(sample)){
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
  
  else if (agg_sims$type == "_total"){
    suffix <- paste0(suffix, "total")
  }
  
  
  p <- plot.data %>%
    as.data.frame() %>%
    ggplot(aes(week)) + 
    geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
    geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = med), col = "steelblue") +
    geom_point(aes(y = obs)) +
    labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths",
         caption = paste0("Median, ",
                          (plot_quants[3]-plot_quants[2])*100, " and ",
                          (plot_quants[4]-plot_quants[1])*100,
                          "% quantile intervals over ", nsims, " posterior simulations, scaled by ",
                          # samples_cfr, 
                          # "samples from",
                          min(scale_quants)*100,"% to ", max(scale_quants)*100,"% quantiles of",
                          " the observed CFR distribution post-P2 expansion.")) +
    theme_minimal()
  
  if (agg_sims$type != "total"){
    p <- p + 
      facet_wrap(~facet, scales = "free") 
  }

  if (save == T){
    if (format == "pdf"){  
      pdf(here::here("figures","compare","expanded",paste0("reconstr_lag",lag,suffix,".pdf")), height = h, width = w)
      print(p)
      dev.off()
    } else if(format == "png") {
      png(here::here("figures","compare","expanded",paste0("reconstr_lag",lag,suffix,".png")), height = h, width = w, res = 200)
      print(p)
      dev.off()
    }
  } else {
    return(p)
  }
}
