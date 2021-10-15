################################################################################
# Description: Summarise reconstructed case counts
# 
# Author: Emily S Nightingale
# Date created: 30/04/2021
# 
################################################################################
################################################################################

outdir <- "output/reconstruct"
figdir <- "figures/reconstruct"

## Shapefiles
regions <- readRDS(here::here("data","LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
border <- sf::st_union(regions)

# Total population
tot_pop <- 56e6

plot_quants <- c(0.01,0.25,0.75,0.99)

# Load observed confirmed cases
cases <- readRDS(here::here("data","aggregated","cases.rds"))[[1]]

# Infection prevalence from ONS survey
infect_prev <- read.csv(here::here("data","infection_prevalence.csv"), header = T) %>%
  dplyr::select(-1) %>%
  dplyr::mutate(week = lubridate::dmy(week) + 2, # ONS week start monday but here week start wednesday, minus 1 week for onset-test delay
                across(estimate:hi, function(x) x/2),
                across(estimate:hi, function(x) (x/100)*tot_pop, .names = "{.col}_n"))

# Reconstructed cases
reconstruct7 <- readRDS(here::here(outdir,"reconstruct_lag7.rds"))
reconstruct14 <- readRDS(here::here(outdir,"reconstruct_lag14.rds"))
reconstruct21 <- readRDS(here::here(outdir,"reconstruct_lag21.rds"))

sink(here::here(outdir, "summary_log.txt"))

# ---------------------------------------------------------------------------- #
## Estimate % infections detected post-P2 ##

png(here::here(figdir,"reconstruct7_vs_infections.png"), height = 1500, width = 2000, res = 300)
reconstruct7$total$preds %>%
  inner_join(infect_prev) %>%
  ggplot(aes(week, estimate_n)) +
  geom_ribbon(aes(ymin = lo_n, ymax = hi_n), alpha = 0.2, fill = "indianred") +
  geom_line(col = "indianred") +
  geom_ribbon(aes(y = med, ymin = l1, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  labs(x = "", y = "Incidence", 
       caption = "Reconstructed P1+P2 cases in blue, estimated infections in red.")
dev.off()

# Sum predicted over LTLAs, join with estimated infection incidence and calculate
# proportion infections detected, then summarise over simulations
reconstruct7$sims %>%
  group_by(week, sim, scale_quant) %>%
  summarise(pred = sum(pred_c)) %>%
  inner_join(dplyr::select(infect_prev, week, estimate_n, lo_n, hi_n)) %>%
  mutate(p = pred*100/estimate_n,
         p_lo = pred*100/lo_n,
         p_hi = pred*100/hi_n) %>%
  pivot_longer(p:p_hi) %>%
  ungroup() -> p_detect

p_detect %>%
  group_by(week) %>%
  dplyr::summarise(lo = quantile(value, plot_quants[1]),
                   med = quantile(value, 0.5),
                   hi = quantile(value, plot_quants[4])) %>% 
  ungroup() -> p_detect_time

png(here::here(figdir,"perc_detect_time.png"), height = 1500, width = 2000, res = 300)
p_detect_time %>%
  ggplot(aes(week, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "steelblue") +
  geom_line(col = "steelblue") +
  labs(x = "", y = "% infections detected under Pillars 1 + 2")
dev.off()

# Overall
reconstruct7$sims %>%
  group_by(week, sim, scale_quant) %>%
  summarise(pred = sum(pred_c)) %>%
  inner_join(dplyr::select(infect_prev, week, estimate_n, lo_n, hi_n)) %>%
  filter(week > ymd("2020-05-18")) %>%
  group_by(sim, scale_quant) %>%
  summarise(across(pred:hi_n, sum)) %>% 
  mutate(p = pred*100/estimate_n,
         p_lo = pred*100/lo_n,
         p_hi = pred*100/hi_n) %>%
  pivot_longer(p:p_hi) %>%
  ungroup() %>%
  summarise(lo = quantile(value, plot_quants[1]),
            med = quantile(value, 0.5),
            hi = quantile(value, plot_quants[4])) %>%
  as.data.frame() -> p_detect_all

round(p_detect_all,2)
#   lo   med    hi
# 9.94 24.17 70.97

# ---------------------------------------------------------------------------- #
## Summarise percentage detected per week and plot ##

# Total observed
reconstruct7$total$preds$obs  -> obs

# Sum predicted over LTLAs, then summarise over posteriors/quantiles
reconstruct7$sims %>%
  dplyr::group_by(week, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c)) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(obs = reconstruct7$total$preds$obs) -> pred_time 

pred_time %>%
  mutate(across(l1:h2, function(x) x*100/p_detect_all$lo, .names = "{.col}_plo"),
         across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
         across(l1:h2, function(x) x*100/p_detect_all$hi,.names = "{.col}_phi"),
         across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_infect_time

png(here::here(figdir,"reconstruct7_vs_infections2.png"), height = 1500, width = 2000, res = 300)
infect_prev %>%
  right_join(pred_infect_time) %>%
  ggplot(aes(week)) +
  # # Infections - lower bound % detected
  # geom_ribbon(data = pred_infect_time,
  #             aes(ymin = l1_plo, ymax = h2_plo), alpha = 0.1, fill = "grey") +
  # geom_line(data = pred_infect_time,
  #           aes(y = med_plo), col = "lightgrey") +
  # # Infections - upper bound % detected
  # geom_ribbon(data = pred_infect_time,
  #             aes(ymin = l1_phi, ymax = h2_phi), alpha = 0.1, fill = "grey") +
  # geom_line(data = pred_infect_time,
  #           aes(y = med_phi), col = "lightgrey") +
  # Infections - estimate % detected
  geom_ribbon(aes(ymin = l1_est, ymax = h2_est), alpha = 0.3, fill = "grey") +
  geom_ribbon(aes(ymin = l2_est, ymax = h1_est), alpha = 0.3, fill = "grey") +
    geom_line(aes(y = med_est), col = "grey") +
  # ONS-estimated weekly incidence of test-positive infection
  geom_ribbon(aes(ymin = lo_n, ymax = hi_n), alpha = 0.2, fill = "indianred") +
  geom_line(aes(y = estimate_n), col = "indianred") +
  # Reconstructed cases detected under P1+2
  geom_ribbon(aes(y = med, ymin = l1, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(y = med, ymin = l2, ymax = h1), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  geom_point(aes(y = obs)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) +
  labs(x = "", y = "Incidence", 
       caption = "Reconstructed P1+2 cases in blue; Infections according to weekly ONS-estimated rates in red.\nGrey curve indicates potential infections under assumption of 23% detection under P1+2.")
dev.off()

pred_time %>%
  dplyr::mutate(obs = obs) %>%
  dplyr::mutate(across(l1:h2, function(pred) return(obs/pred))) -> obs_pred_time 

print("Summary observed versus predicted over time:")
summary(obs_pred_time)

# Plot relative differences in national total per week 
pred_infect_time %>%
  filter(obs > 0) %>%
  ggplot(aes(week)) +
  geom_ribbon(aes(ymin = l1_est_percdet, ymax = h2_est_percdet), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(ymin = l2_est_percdet, ymax = h1_est_percdet), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med_est_percdet), col = "steelblue") +
  scale_x_date(breaks = "month", date_labels = "%b") +
  # scale_y_continuous(trans = "log2",
  #                    labels = scales::number_format(accuracy = 0.01)) +
  geom_vline(xintercept = ymd("2020-03-13"), col = "red") +
  geom_vline(xintercept = ymd("2020-04-15"), col = "red") +
  geom_vline(xintercept = ymd("2020-05-18"), col = "red") +
  annotate("text", x = ymd("2020-05-19"), y = 15, 
           label = "P2 available to all symptomatic cases", 
           cex = 2, hjust = "left") +
  annotate("text", x = ymd("2020-04-16"), y = 10, 
           label = "P2 available to care home residents and staff", 
           cex = 2, hjust = "left") +
  annotate("text", x = ymd("2020-03-14"), y = 5, 
           label = "Testing pivoted from community to hospital need", 
           cex = 2, hjust = "left") +
  # geom_hline(yintercept = 1) +
  labs(x = "", y = "% infections detected") -> perc_detect_time #% total infections detected

png(here::here(figdir,"perc_detect_time.png"), height = 1200, width = 2000, res = 300)
perc_detect_time
dev.off()

# ---------------------------------------------------------------------------- #
# Compare P1 and P1+2

reconstruct7$total$preds %>%
  dplyr::mutate(period = case_when(week < ymd("2020-05-18") ~ "Pillar 1",
                                   week >= ymd("2020-05-18") ~ "Pillar 1+2")) %>% 
  dplyr::group_by(period) %>%
  dplyr::summarise(obs = sum(obs, na.rm = T)) %>%
  dplyr::ungroup() -> obs_period
  
reconstruct7$total$preds %>% 
  dplyr::mutate(period = case_when(week < ymd("2020-05-18") ~ "Pillar 1",
                                   week >= ymd("2020-05-18") ~ "Pillar 1+2")) %>%
  dplyr::group_by(period) %>%
  dplyr::summarise(across(l1:obs, sum)) -> pred_period

pred_period %>%
  mutate(across(l1:h2, function(x) x*100/p_detect_all$lo, .names = "{.col}_plo"),
         across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
         across(l1:h2, function(x) x*100/p_detect_all$hi,.names = "{.col}_phi"),
         across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_infect_period

png(here::here(figdir,"perc_detect_period.png"), height = 1000, width = 1200, res = 300)
pred_infect_period %>%
  ggplot(aes(period, med_est_percdet, ymin = l1_est_percdet, ymax = h2_est_percdet)) +
  geom_errorbar(width = 0.1) +
  geom_point() +
  labs(x = "", y = "Estimated % infections detected")
dev.off()

# March total
reconstruct7$total$preds %>%
  filter(week > ymd("2020-03-01") & week < ymd("2020-04-01")) %>%
  dplyr::summarise(obs = sum(obs)) %>% pull(obs) -> obs_march

reconstruct7$sims %>% 
  filter(week > ymd("2020-03-01") & week < ymd("2020-04-01")) %>%
  dplyr::group_by(sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c)) %>%
  dplyr::ungroup() %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>%
  ungroup() %>%
  mutate(obs = obs_march) %>% 
  dplyr::mutate(across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
                across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_march

print("Estimated total infections, March 2020:")
dplyr::select(pred_march, l1_est:h2_est)
# l1_est  l2_est med_est  h1_est  h2_est
# 239460. 253010. 306806. 378747. 398601.

print("% infections detected, March 2020:")
dplyr::select(pred_march, l1_est_percdet:h2_est_percdet)
# l1_est_percdet l2_est_percdet med_est_percdet h1_est_percdet h2_est_percdet
# 13.2           12.5            10.3           8.33           7.92


# ~10% infections detected by testing in March 2020

# ---------------------------------------------------------------------------- #
## Summarise percentage detected per LTLA ##

reconstruct7$la$preds %>%
  # dplyr::filter(week < ymd("2020-05-18")) %>%
  dplyr::group_by(lad19nm) %>%
  dplyr::summarise(obs = sum(obs, na.rm = T)) %>%
  ungroup() -> obs_la

# Sum predicted over LTLAs, then summarise over posteriors/quantiles
reconstruct7$sims %>%
  # dplyr::filter(week < ymd("2020-05-18")) %>%
  dplyr::group_by(lad19nm, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c)) %>%
  dplyr::group_by(lad19nm) %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>% 
  dplyr::ungroup() %>%
  dplyr::full_join(obs_la) -> pred_la

pred_la %>%
  mutate(across(l1:h2, function(x) x*100/p_detect_all$lo, .names = "{.col}_plo"),
         across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
         across(l1:h2, function(x) x*100/p_detect_all$hi,.names = "{.col}_phi"),
         across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_infect_la

summary(pred_infect_la$med_est_percdet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.964  13.675  20.944  24.447  30.434  86.606 

png(here::here(figdir,"perc_detect_space_hist_lag7.png"), height = 1500, width = 2000, res = 300)
pred_infect_la %>%
  ggplot() +
  geom_histogram(aes(x = med_est_percdet),alpha = 0.2) +
  geom_histogram(aes(x = l1_est_percdet),alpha = 0.2) +
  geom_histogram(aes(x = h2_est_percdet),alpha = 0.2) +
  labs(x = "% infections detected",
       y = "Frequency")
dev.off()

# Map out percentage infections detected in cumulative totals per LTLA
regions %>% 
  dplyr::full_join(pred_infect_la) %>%
  basic_map(fill = "med_est_percdet", plot.border = T) +
  # scale_fill_gradient2(midpoint = 0, trans = "log2") +
  labs(fill = "% infections\ndetected") +
  theme(legend.position = c(0.2, 0.5)) -> perc_detect_space

png(here::here(figdir,"perc_detect_space_lag7.png"), height = 1500, width = 2000, res = 300)
perc_detect_space
dev.off()

png(here::here(figdir,"perc_detect_time_space_lag7.png"), height = 1200, width = 2000, res = 300)
perc_detect_space + perc_detect_time + 
  plot_layout(widths = c(1,1.5), heights = c(1,0.8))
dev.off()

# ---------------------------------------------------------------------------- #
# Plot over time for sample with highest and lowest % detected

pred_infect_la %>%
  slice_max(order_by = med_est_percdet, n = 2) %>%
  pull(lad19nm) -> perc_hi
pred_infect_la %>%
  slice_min(order_by = med_est_percdet, n = 2) %>%
  pull(lad19nm) -> perc_lo
sample <- c(perc_hi,perc_lo)
sample_r <- sample(regions.df$lad19nm, 9)
  
# Predictions per week and per ltla
reconstruct7$la$preds %>%
  full_join(dplyr::select(regions.df, lad19nm, la_pop)) %>%
  left_join(infect_prev) %>%
  mutate(across(estimate:hi, function(x) (x/100)*la_pop, .names = "{.col}_n"),
         # across(l1:h2, function(x) x*100/p_detect_all$lo, .names = "{.col}_plo"),
         across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
         # across(l1:h2, function(x) x*100/p_detect_all$hi,.names = "{.col}_phi"),
         across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_propdet")) -> pred_infect

png(here::here(figdir,"reconstruct7_vs_infections_lasamp.png"), height = 2000, width = 2400, res = 300)
pred_infect %>%
  filter(lad19nm %in% sample_r) %>%
  ggplot(aes(week)) +
  geom_ribbon(aes(ymin = l1_est, ymax = h2_est), alpha = 0.3, fill = "grey") +
  geom_line(aes(y = med_est), col = "grey") +
  # Reconstructed cases detected under P1+2
  geom_ribbon(aes(y = med, ymin = l1, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  geom_point(aes(y = obs)) +
  # scale_y_continuous(trans = "log10") +
  labs(x = "", y = "Incidence", 
       caption = "Reconstructed P1+2 cases in blue. Grey curve indicates infection trajectory under the\nassumption of 23% detection overall under P1+2.") +
  facet_wrap(~lad19nm, scales = "free_y")
dev.off()

# ---------------------------------------------------------------------------- #
## Compare one/two/three week lags by national total time series ##

lagcomp <- dplyr::bind_rows(reconstruct7$total$preds, 
                            reconstruct14$total$preds, 
                            reconstruct21$total$preds) %>%
  rename(Lag = lag)

cases %>%
  group_by(week) %>%
  summarise(obs = sum(n)) -> obs_wk

png(here::here(figdir,"compare_lags.png"), height = 1400, width = 2200, res = 300)
ggplot(lagcomp, aes(week)) + 
  geom_ribbon(aes(ymin = l1, ymax = h1, fill = Lag), alpha = 0.2) +
  geom_ribbon(aes(ymin = l2, ymax = h2, fill = Lag), alpha = 0.2) +
  geom_line(aes(y = med, colour = Lag)) +
  geom_point(data = obs_wk, aes(y = obs)) +
  scale_fill_viridis_d(begin = 0.2, end = 0.8) +
  scale_colour_viridis_d(begin = 0.2, end = 0.8) +
  scale_x_date(breaks = "month", date_labels = "%b") +
  labs(x = "",y = "" 
       #, title = "Reconstruction of confirmed cases from COVID-19-related deaths across England",
       # subtitle = "Comparison of assumed lags between case confirmation and death"
  ) +
  theme_minimal() +
  theme(legend.position = c(0.8,0.8))
dev.off()

# ---------------------------------------------------------------------------- #
# Construct table of observed and predicted counts, nationally and by geog, with
# CrIs and percentage differences

make_table <- function(reconstructed) {
  
  lag <- reconstructed$lag
  
  ## By geography
  reconstructed$sims %>%
    dplyr::group_by(geography, sim, scale_quant) %>%
    dplyr::summarise(pred_c = sum(pred_c, na.rm = TRUE),
                     cases = sum(n_c, na.rm = TRUE)) %>%
    dplyr::ungroup() -> geog_totals
  
  reconstructed$sims %>%
    dplyr::group_by(sim, scale_quant) %>%
    dplyr::summarise(pred_c = sum(pred_c, na.rm = TRUE),
                     cases = sum(n_c, na.rm = TRUE),
                     geography = "England") %>% 
    dplyr::ungroup() %>%
    dplyr::bind_rows(geog_totals) %>%
    dplyr::group_by(geography) %>%
    dplyr::summarise(low = quantile(pred_c, plot_quants[2], na.rm = TRUE),
                     med = quantile(pred_c, 0.5, na.rm = TRUE),
                     high = quantile(pred_c, plot_quants[3], na.rm = TRUE),
                     observed = unique(cases, na.rm = TRUE)) %>%
    dplyr::ungroup() -> totals
  
  totals
  # geography                     low     med    high observed
  # 1 England                   276219. 336378. 420491.   231817
  # 2 London Borough             36124.  43969.  51631.    33399
  # 3 Metropolitan District      95779. 109869. 129519.    64007
  # 4 Non-metropolitan District  79038.  98271. 122657.    79441
  # 5 Unitary Authority          67947.  84017. 111498.    54970
  
  
  totals %>%
    dplyr::mutate(across(c(-geography,-observed), function(x) round((x - totals$observed)*100/totals$observed,1), .names = "diff_{.col}")) %>%
    dplyr::mutate(across(c("observed","med","low","high"), function(x) scales::comma(round(x, 0)))) %>%
    dplyr::mutate(pred = paste0(med, " [", low, " - ", high,"]"),
                  diff = paste0(diff_med, " [", diff_low, " - ", diff_high,"]")) %>%
    dplyr::select(geography, observed, pred, diff) %>%
    dplyr::rename_with(function(n) paste0(n,lag), .cols = !geography) -> tab
  
  tab
  # geography                 observed pred                        diff               
  # 1 England                     231817 336,378 [276,219 - 420,491] 45.1 [19.2 - 81.4] 
  # 2 London Borough               33399 43,969 [36,124 - 51,631]    31.6 [8.2 - 54.6]  
  # 3 Metropolitan District        64007 109,869 [95,779 - 129,519]  71.7 [49.6 - 102.4]
  # 4 Non-metropolitan District    79441 98,271 [79,038 - 122,657]   23.7 [-0.5 - 54.4] 
  # 5 Unitary Authority            54970 84,017 [67,947 - 111,498]   52.8 [23.6 - 102.8]
  
  return(tab)
}

tabs <- lapply(list(reconstruct7, reconstruct14, reconstruct21), make_table)

tab2 <- tabs[[1]] 
tabS2 <- full_join(tabs[[2]], tabs[[3]], by = "geography")

write.csv(tab2, here::here(outdir,"table2.csv"), row.names = FALSE)
write.csv(tabS2, here::here(outdir,"tableS2.csv"), row.names = FALSE)

################################################################################
sink()
################################################################################




# ---------------------------------------------------------------------------- #
# By geography
# 
# regions.df %>%
#   group_by(geography) %>%
#   summarise(geog_pop = sum(la_pop)) -> geog_pops
# 
# # Predictions per week by geography
# reconstruct7$geog$preds %>%
#   full_join(geog_pops) %>%
#   left_join(infect_prev) %>%
#   mutate(across(estimate:hi, function(x) (x/100)*geog_pop, .names = "{.col}_n"),
#          # across(l1:h2, function(x) x*100/p_detect_all$lo, .names = "{.col}_plo"),
#          across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
#          # across(l1:h2, function(x) x*100/p_detect_all$hi,.names = "{.col}_phi"),
#          across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_propdet")) -> pred_infect_geog
# 
# png(here::here(figdir,"reconstruct7_vs_infections_geog.png"), height = 1800, width = 2000, res = 300)
# pred_infect %>%
#   ggplot(aes(week)) +
#   # # Infections - lower bound % detected
#   # geom_ribbon(data = pred_infect,
#   #             aes(ymin = l1_plo, ymax = h2_plo), alpha = 0.1, fill = "grey") +
#   # geom_line(data = pred_infect,
#   #           aes(y = med_plo), col = "lightgrey") +
#   # # Infections - upper bound % detected
#   # geom_ribbon(data = pred_infect,
#   #             aes(ymin = l1_phi, ymax = h2_phi), alpha = 0.1, fill = "grey") +
#   # geom_line(data = pred_infect,
#   #           aes(y = med_phi), col = "lightgrey") +
#   # Infections - estimate % detected
# geom_ribbon(aes(ymin = l1_est, ymax = h2_est), alpha = 0.3, fill = "grey") +
#   geom_line(aes(y = med_est), col = "grey") +
#   # ONS-estimated weekly incidence of test-positive infection
#   geom_ribbon(aes(ymin = lo_n, ymax = hi_n), alpha = 0.2, fill = "indianred") +
#   geom_line(aes(y = estimate_n), col = "indianred") +
#   # Reconstructed cases detected under P1+2
#   geom_ribbon(aes(y = med, ymin = l1, ymax = h2), alpha = 0.2, fill = "steelblue") +
#   geom_line(aes(y = med), col = "steelblue") +
#   # scale_y_continuous(trans = "log10") +
#   labs(x = "", y = "Incidence", 
#        caption = "Reconstructed P1+2 cases in blue; infections based on national rates in red.\n
#        Grey curve indicates infection trajectory under the assumption of 23% detection overall under P1+2.") +
#   facet_wrap(~geography, scales = "free_y")
# dev.off()