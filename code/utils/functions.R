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
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.position = c(0.1,0.5)) 
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
  
  # annotation_north_arrow(location = "tl", which_north = "true", 
  #     # pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #     style = north_arrow_fancy_orienteering) 

  return(p)
}

map_w_trans <- function(sf, fill, alpha, rate1e5 = F){
  
  if (rate1e5 == T){
    sf <- mutate(sf, fill = !!sym(fill)*1e5/la_pop)
  }else{sf <- mutate(sf, fill = !!sym(fill))}
  
  p <- ggplot(sf, aes(geometry = geometry, fill = fill, alpha = !!sym(alpha))) +
    geom_sf() +
    scale_fill_viridis_c() +
    labs(fill = "") +
    annotation_scale(location = "br") +
    map_theme() 
  # annotation_north_arrow(location = "tl", which_north = "true", 
  #     # pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #     style = north_arrow_fancy_orienteering) 
  return(p)
}



map_wlondon <- function(sf, fill, rate1e5 = F){
  
  if (rate1e5 == T){
    sf <- mutate(sf, fill = !!sym(fill)*1e5/la_pop)
  }else{sf <- mutate(sf, fill = !!sym(fill))}
  
  london <- filter(sf, geography == "London borough")
  
  p <- ggplot(sf, aes(geometry = geometry, fill = fill)) +
    geom_sf() +
    scale_fill_viridis_c() +
    labs(fill = "") +
    annotation_scale(location = "br") +
    map_theme() 
  
  p_sub <- ggplot(london, aes(geometry = geometry, fill = fill)) +
    geom_sf() +
    scale_fill_viridis_c() +
    labs(fill = "") +
    annotation_scale(location = "br") +
    map_theme() 
  
  return(p + p_sub)
}


## INLA MODEL SPECIFICATION ##
fit_mod <- function(f, dat, expected = "E", family = "nbinomial"){
  
  if (expected == "E"){
  fit <- inla(f,
              family,
              data = dat,
              E = E_wk,
              # offset = log(la_age_pop),
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
                # E = E_wk,
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


get_preds <- function(sample,dat){
  pred <- sample$latent[1:nrow(dat)]
}

## PLOT MODEL SUMMARY ##
plot.model <- function(m){
  p <- autoplot(m)
  print(
    cowplot::plot_grid(plotlist = p)
  )
}

plot.resids <- function(m){
  p1 <- ggplot_inla_residuals(m, dat$SMR, binwidth = 0.1)
  p2 <- ggplot_inla_residuals2(m, dat$SMR, se = TRUE)
  return(list(p1, p2))
}

# 
# get_resid <- function(fit){
#   dat %>% 
#     mutate(mu = fit$summary.fitted.values$mean,
#            sigma2 = mu*(1 + mu/fit$summary.hyperpar[1,"mean"]),
#            resid = (SMR-mu)/sqrt(sigma2),
#            MSqE = (SMR-mu)^2,
#            logs = log(fit$cpo$cpo)) 
# }
# 
# resids <- lapply(fits, function(fit) get_resid(fit)$resid) %>%
#   bind_cols() %>%
#   setNames(unlist(names(fits)))
# 
# dat_resid <- bind_cols(dat, resids) 
# 
# png(here("figures","final_altdata","map_resids.png"), height = 800, width = 1200, res = 150)
# # tiff(filename = "./figures/final_altdata/map_resids.tif", height = 600, width = 1000)
# dat_resid %>%
#   pivot_longer(cols = base:BYM_geog) %>%
#   group_by(lad19cd, name) %>%
#   summarise(value = mean(value)) %>% 
#   left_join(regions) %>%
#   basic_map(fill = "value") +
#   # scale_fill_viridis_c(trans = "log10") +
#   facet_wrap(~name, ncol = 3, nrow = 2) +
#   labs(title = "Mean Squared Error per local authority")
# dev.off()
