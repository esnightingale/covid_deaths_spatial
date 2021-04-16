
plot_la_ts <- function(data, wave = 1, title1, title2){
  
  period <- data$breaks.first
  
  data[[wave]] %>%
    group_by(lad19cd) %>%
    mutate(peak = w[which.max(n)],
           rate_1e5 = n*1e5/la_pop) %>% 
    arrange(peak) %>% 
    ggplot(aes(w, rate_1e5, group = lad19cd, col = lad19cd)) +
    geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
    geom_line(alpha = 0.3) +
    labs(title = title1,
         subtitle = paste(period[1],"-",period[2]),
         x = "Calendar week",
         y = "Count per 100,000"
         # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
    ) +
    guides(col = F) +
    scale_colour_viridis_d(option = "cividis")  -> ts_by_wk
  
  data[[wave]] %>%
    group_by(lad19cd) %>%
    mutate(peak = wk_since_first[which.max(n)],
           rate_1e5 = n*1e5/la_pop) %>% 
    arrange(peak) %>%
    ggplot(aes(wk_since_first, rate_1e5, group = lad19cd, col = lad19cd)) +
    geom_vline(aes(xintercept = jitter(peak), col = lad19cd), alpha = 0.3) +
    geom_line(alpha = 0.3) +
    labs(title = title2,
         x = paste0("Weeks since first", measure),
         y = "Count per 100,000"
         # subtitle = paste0("Data up to ", max(lubridate::ymd(covid_deaths$dor))),
    ) +
    guides(col = F) +
    scale_colour_viridis_d(option = "cividis") -> ts_by_epiwk
  
  p <- ts_by_wk / ts_by_epiwk
  
  return(p)
  
}


plot_geog_ts <- function(data, wave = 1, title1, title2){
  
  period <- data$breaks[[wave]]
  
  data[[wave]] %>%
    group_by(w,week,geography) %>%
    summarise(n = sum(n, na.rm= T),
              geog_pop = sum(la_pop)) %>%
    ggplot(aes(week, n*1e5/geog_pop, group = geography, col = geography)) +
    geom_line() +
    labs(title = title1,
         subtitle = paste(period[1],"-",period[2]),
         x = "Calendar week",
         y = "Count per 100,000", 
         colour = "Geography type") +
    theme(legend.position = c(0.16,0.65), legend.text=element_text(size=8),  legend.title=element_text(size=8)) -> ts_geog_week
  
  
  data[[wave]] %>%
    group_by(wk_since_first,geography) %>%
    summarise(n = sum(n, na.rm= T),
              geog_pop = sum(la_pop)) %>%
    ggplot(aes(wk_since_first, n*1e5/geog_pop, group = geography, col = geography)) +
    geom_line() +
    labs(title = title2,
         x = "Weeks since first",
         y = "Count per 100,000",
         colour = "Geography type") +
    theme(legend.position = "none") -> ts_geog_epiwk
  
  p <- ts_geog_week / ts_geog_epiwk
  
  return(p)
  
}


plot_one_la <- function(laID){
  
  dat_pred_c %>%
    filter(la == laID) %>%
    ggplot() + 
    geom_line(aes(week, pred_n, group = name), alpha = 0.1, col = "grey") +
    geom_point(aes(week, n), col = "black") + 
    geom_line(data = filter(dat_pred_d,la == laID), aes(week-2, pred_n, group = name), alpha = 0.1, col = "steelblue3") +
    geom_point(data = filter(dat_pred_d,la == laID), aes(week-2, n), col = "steelblue4") +
    theme_minimal() -> p
  
  print(p)
}


plot_la_samp <- function(data, la_samp, pred = "pred_n", obs = "n"){
  
  data %>%
    filter(la %in% la_samp) %>%
    group_by(week, lad19nm) %>%
    summarise(med = median(!!sym(pred))) -> medians
  return(
    data %>%
      filter(la %in% la_samp) %>%
      ggplot() +
      geom_line(aes(week, !!sym(pred), group = name), alpha = 0.1, col = "grey") +
      geom_point(aes(week, !!sym(obs))) +
      facet_wrap(~lad19nm) +
      theme_minimal() +
      geom_line(data = medians, aes(week, med), col = "darkgrey") 
  )
}

plot_la_tot <- function(data, wave = 1, title){
  
  period <- data$breaks[[wave]]
  
  data[[wave]] %>%
    group_by(lad19cd) %>%
    summarise(n = sum(n, na.rm = TRUE)) %>% 
    full_join(regions) %>%
    basic_map(fill = "n", rate1e5 = TRUE, scale = F) +
    labs(title = title,
         subtitle = paste(period[1],"-",period[2])) + 
    map_theme() -> map_tot
  
  return(map_tot)
  
}

plot_la_mth <- function(data, title){
  
  data <- data %>%
    mutate(month = droplevels(month))
  
  combs <- data %>%
    tidyr::expand(lad19cd, month)
  agg <- data %>%
    group_by(lad19cd, month) %>%
    tally() %>%
    right_join(combs) 
  
  regions %>% 
    full_join(agg) %>%
    basic_map(fill = "n", rate1e5 = TRUE) +
    labs(title = title) + 
    map_theme() +
    facet_wrap(~month) -> map_mth
  
  return(map_mth)
  
}

plot_la_mth_animate <- function(data, title, file, path){

  data <- data %>%
    mutate(month = droplevels(month))
  
  combs <- data %>%
    tidyr::expand(lad19cd, month)
  agg <- data %>%
    group_by(lad19cd, month) %>%
    tally() %>%
    right_join(combs) 
  
  regions %>%
    full_join(agg) %>%
    # mutate(n = replace_na(n, 0)) %>%
    basic_map(fill = "n", rate1e5 = TRUE) +
    map_theme() +
    labs(title = title, 
         subtitle = '{current_frame}') +
    transition_manual(month) -> map_anim
  
  anim_save(
    file = file, 
    animate(
      map_anim,
      renderer = gifski_renderer(),
      height = 800, width = 800, res = 150
    ),
    path = path)
  
  
}

plot_epi_time <- function(data, wave = 1, measure){
  
  period <- data$breaks[[wave]]
  
  data[[wave]] %>%
    group_by(lad19cd) %>%
    summarise(first = unique(first)) %>% 
    full_join(regions) %>%
    basic_map(fill = "first") +
    scale_fill_viridis(option = "plasma", trans = "log2") +
    labs(title = paste0("Week of first ",measure),
         subtitle = paste(period[1],"-",period[2])) + 
    map_theme()  -> first_map
  
  data[[wave]] %>%
    group_by(lad19cd) %>%
    mutate(peak = w[which.max(n)]) %>% 
    arrange(peak) %>%
    full_join(regions) %>%
    basic_map(fill = "peak") +
    scale_fill_viridis(option = "plasma") +
    labs(title = paste0("Week of peak ",measure,"s")) +
    map_theme() -> peak_map
  
  p <- first_map + peak_map
  return(p)
  
}

plot_parm <- function(parm, fit, opt = 1){
  
  if (opt == 1){ d <- fit$marginals.fixed[[parm]]
  }else{ d <- fit$marginals.hyperpar[[parm]] }
  
  print(
    ggplot(data.frame(inla.smarginal(d)), aes(x, y)) +
      geom_line() +
      geom_vline(xintercept = 0, col = "red", lty = "dashed") +
      labs(title = parm) +
      theme_bw()
  )
}
