################################################################################
# Description: Helper functions for spatio-temporal analysis of COVID-19 
# deaths in England
# 
# Author: Emily S Nightingale
# Date created: 30/09/2020
# 
################################################################################
################################################################################

## MAP PLOT ##
map_theme <- function () {
  theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.position = c(0.1,0.5),
          plot.margin = margin(0, 0, 0, 0, "cm"))
}

basic_map <- function(sf, fill, rate1e5 = F, plot.border = F, scale = T){
  
  if (rate1e5 == T){
    sf <- mutate(sf, fill = !!sym(fill)*1e5/la_pop)
  }else{sf <- mutate(sf, fill = !!sym(fill))}
  
  p <- ggplot() +
    # Base layer to avoid gaps with missing fill values when faceting
    geom_sf(data = regions, aes(geometry = geometry), fill = "grey", colour = NA) +
    geom_sf(data = sf, aes(geometry = geometry, fill = fill), colour = NA) +
    scale_fill_viridis_c(na.value = "grey") +
    labs(fill = "") +
    map_theme() 

  if (plot.border == T){ 
    p <- p + geom_sf(data = border, alpha = 0, lwd = 0.5, col = "grey") 
  }
  
  if (scale == T){ 
    p <- p + annotation_scale(location = "br") 
  }

  return(p)
}

inset_map <- function(sf, fill, fill.lab = "", legend.pos = c(0.2, 0.5), rate1e5 = F, plot.border = F, scale = T){
  
  if (rate1e5 == T){
    sf <- mutate(sf, fill = !!sym(fill)*1e5/la_pop)
  }else{sf <- mutate(sf, fill = !!sym(fill))}
  
  lims <- c(min(sf$fill, na.rm = T), max(sf$fill, na.rm = T))
  print(lims)
  
  p1 <- ggplot() +
    # Base layer to avoid gaps with missing fill values when faceting
    geom_sf(data = regions, aes(geometry = geometry), fill = "grey", colour = NA) +
    geom_sf(data = sf, aes(geometry = geometry, fill = fill), colour = NA) +
    scale_fill_viridis_c(na.value = "grey", limits = lims) +
    labs(fill = fill.lab) +
    map_theme() +
    theme(legend.position = legend.pos)
  
  sf.l <- filter(sf, geography == "London Borough")
  regions.l <- filter(regions, geography == "London Borough")
  
  p2 <- ggplot() +
    # Base layer to avoid gaps with missing fill values when faceting
    geom_sf(data = regions.l, aes(geometry = geometry), fill = "grey", colour = NA) +
    geom_sf(data = sf.l, aes(geometry = geometry, fill = fill), colour = NA) +
    scale_fill_viridis_c(na.value = "grey", limits = lims) +
    labs(fill = "", subtitle = "London Boroughs") +
    guides(fill = "none") +
    map_theme() +
    theme(plot.subtitle=element_text(size=8, hjust=-0.5, face="italic"))
  
  if (plot.border == T){ 
    p1 <- p1 + geom_sf(data = border, alpha = 0, lwd = 0.5, col = "grey") 
    p2 <- p2 + geom_sf(data = border, alpha = 0, lwd = 0.5, col = "grey") 
  }
  
  if (scale == T){ 
    p1 <- p1 + annotation_scale(location = "br") 
    p2 <- p2 + annotation_scale(location = "br") 
  }

  p <- ggdraw() +
    draw_plot(p1) +
    draw_plot(p2, x = 0.7, y = 0.65, width = 0.3, height = 0.3)
  
  return(p)
}


## INLA MODEL SPECIFICATION ##
fit_mod <- function(f, dat, expected = "E", family = "nbinomial"){
  
  if (expected == "E"){
  fit <- inla(f,
              family,
              data = dat,
              E = E_wk,
              control.compute=list(dic=TRUE, 
                                   waic=TRUE, 
                                   cpo = TRUE,
                                   config = TRUE),
              control.predictor = list(compute = TRUE, link = 1),
              control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
              verbose = T)
  }else{
    fit <- inla(f,
                family,
                data = dat,
                offset = log(la_pop),
                control.compute=list(dic=TRUE, 
                                     waic=TRUE, 
                                     cpo = TRUE,
                                     config = TRUE),
                control.predictor = list(compute = TRUE, link = 1),
                control.fixed=list(mean=0, prec=0.1, mean.intercept=0, prec.intercept=0.001),
                control.inla = list(cmin = 0),
                verbose = T)
  }
  return(fit)
}

## MODEL COMPARISON TABLE ##
model_comp <- function(fit.list){
  return(
  data.frame(DIC = sapply(fit.list, function(m){return(m$dic$dic)}),
             WAIC = sapply(fit.list, function(m){return(m$waic$waic)}),
             logs = sapply(fit.list, function(m){return(-mean(log(m$cpo$cpo), na.rm = T))})) %>% 
    rownames_to_column(var = "Model") %>%
    mutate(diff_WAIC = WAIC - min(WAIC),
           diff_DIC = DIC - min(DIC),
           diff_logs = logs - min(logs)) %>%
    arrange(diff_WAIC) %>%
    mutate(across(-Model, function(x) round(x,3))) 
  )
}

## PIT HISTOGRAM ##
pit_hist <- function(fit, bins = 30){
  print(summary(fit$cpo$failure))
  plotdata <- data.frame(pit = fit$cpo$pit)
  return(
    ggplot(plotdata, aes(pit, after_stat(density))) + 
      geom_histogram(fill = "white", col = "black", bins = bins) + 
      geom_hline(yintercept = 1, lty = "dashed", col = "red") + 
      ylim(0,1.5) +
      theme_minimal()+ 
      labs(x = "Probability Integral Transform (PIT)", y = "Density")
  )
}

print_coeffs <- function(m){
  exp(m$summary.fixed)
}

transform_coeff <- function(marg){
  emarg <- inla.tmarginal(exp, marg)
  zmarg <- inla.zmarginal(emarg, silent = T)
  summ.coeff <- setNames(t(zmarg), names(zmarg))
  return(summ.coeff)
}

get_coeffs <- function(fit){
  coeftab <- purrr::map_dfr(fit$marginals.fixed, transform_coeff) %>%
    mutate(effect = names(fit$marginals.fixed), .before = everything())
  return(coeftab)
}

get_preds <- function(sample,dat){
  pred <- sample$latent[1:nrow(dat)]
}


# For each simulated trajectory, scale by a sampled CFR from the observed empirical distribution post-P2
# lag_rescale <- function(lag, ratiodist, sims, nsamp = 1){
# 
#   sims_scale <- sims[1,]
#   for (i in 1:nsamp){
#   sims %>%
#     group_by(name) %>%
#     mutate(scale = EnvStats::remp(1, ratiodist),
#            week = week - lag,
#            pred_n_scale = pred_n*scale) -> temp
#   sims_scale <- bind_rows(sims_scale, temp)
#   }
#   return(sims_scale[-1,])
# }

lag_rescale <- function(scale, lag, sims){
  
  # setDT(sims)
  # sims_scale <- sims[,week := week - lag]
  # sims_scale <- sims_scale[, pred_n_scale := pred_n*scale]
  sims %>%
    mutate(week = week - lag,
           pred_n_scale = pred_n*scale,
           scale = scale) -> sims_scale
  return(sims_scale)
}
