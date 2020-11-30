plot_epi_time <- function(data, wave = 1, measure){
  
  period <- data$breaks[[wave]]
  
data[[wave]] %>%
  group_by(lad19cd) %>%
  summarise(first = unique(first)) %>% 
  full_join(regions) %>%
  basic_map(fill = "first") +
  scale_fill_viridis(option = "plasma") +
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