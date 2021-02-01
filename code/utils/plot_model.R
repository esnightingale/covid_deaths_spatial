
plot.model <- function(dat, fit, samples, type = "time", la_samp = NULL){

  nval <- nrow(dat)

  preds <- bind_cols(lapply(samples, get_preds))

  dat_pred <- bind_cols(dplyr::select(dat, geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n), preds) %>%
    pivot_longer(cols = -1:-8) %>%
    mutate(pred_n = exp(value)*E_wk)
  
  map <- dat_pred %>%
      group_by(lad19cd, name) %>%
      # sum over weeks
      summarise(value = sum(pred_n)) %>%
      # average over samples
      group_by(lad19cd) %>%
      summarise(q50 = median(value),
                q01 = quantile(value, 0.01),
                q99 = quantile(value, 0.99)) %>% 
      pivot_longer(cols = contains("q")) %>% 
      left_join(regions) %>%
      basic_map(fill = "value", rate1e5 = TRUE) +
      facet_wrap(~name) +
      labs(title = "Predicted deaths per 100,000: Quantiles of 1000 posterior samples per local authority")

  time <- dat_pred %>%
      group_by(week, name) %>%
      summarise(pred_n = sum(pred_n, na.rm = T),
                n = sum(n, na.rm = T)) %>%
      ggplot() + 
      geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
      geom_point(aes(week, n)) + 
      theme_minimal()

  
  time_geog <- dat_pred %>%
      group_by(week, name, geography) %>%
      summarise(pred_n = sum(pred_n),
                n = sum(n)) %>%
      ggplot() + 
      geom_line(aes(week, pred_n, group = name, col = geography), alpha = 0.1, col = "grey") +
      geom_point(aes(week, n)) + 
      facet_wrap(~geography) +
      theme_minimal()
  
  
  time_la <- dat_pred %>%
        dplyr::filter(la %in% la_samp) %>%
        ggplot() +
        geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
        geom_point(aes(week, n)) +
        facet_wrap(~lad19nm) +
        theme_minimal()
  
  if (type == "map"){return(map)
    }else if (type == "time"){return(time)
    }else if (type == "time_geog"){return(time_geog)
        }else if (type == "time_la"){return(time_la)}

}

plot.model(dat, fit, samples, type = "time")
