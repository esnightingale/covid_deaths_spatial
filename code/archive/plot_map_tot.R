plot_la_tot <- function(data, wave = 1, title){
  
  period <- data$breaks[[wave]]
  
  data[[wave]] %>%
    group_by(lad19cd) %>%
    summarise(n = sum(n)) %>% 
    full_join(regions) %>%
    basic_map(fill = "n", rate1e5 = TRUE) +
    labs(title = title,
         subtitle = paste(period[1],"-",period[2])) + 
    map_theme() -> map_tot
  
  return(map_tot)
  
  # dat2 %>%
  #   group_by(lad19cd) %>%
  #   summarise(n = sum(n)) %>% 
  #   full_join(regions) %>%
  #   basic_map(fill = "n", rate1e5 = TRUE) +
  #   labs(subtitle = paste(period2[1],"-",period2[2])) + 
  #   map_theme() -> map_tot2
  # 
  # dat %>%
  #   group_by(lad19cd) %>%
  #   summarise(n = unique(E)) %>% 
  #   full_join(regions) %>%
  #   basic_map(fill = "n", rate1e5 = TRUE) + 
  #   map_theme() +
  #   labs(title = "Expected deaths per 100,000", subtitle = "according to age distribution of population") -> map_E
  # 
  
}
