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
regions.df <- sf::st_drop_geometry(regions)
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

sink(here::here(outdir, "summary_log.txt"))

# ---------------------------------------------------------------------------- #
## Estimate % infections detected post-P2 ##

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

reconstruct7$total$preds %>%
  inner_join(dplyr::select(infect_prev, week, estimate_n, lo_n, hi_n)) %>% 
  filter(week > ymd("2020-05-18")) %>%
  summarise(across(c(obs,estimate_n:hi_n), sum)) %>% 
  mutate(med = obs*100/estimate_n,
         hi = obs*100/lo_n,
         lo = obs*100/hi_n) %>%
  select(lo, med, hi) %>%
  as.data.frame() -> p_detect_all

round(p_detect_all,2)
#    lo   med    hi
# 12.89 25.32 58.01

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

pred_time %>%
  dplyr::select(-obs) %>%
  mutate(across(l1:h2, function(x) x*100/p_detect_all$med),
         variable = "Total infections") -> plot_infects

pred_time %>%
  dplyr::select(-obs) %>%
  mutate(variable = "Predicted-P1+P2 cases") -> plot_cases

pred_time %>%
  dplyr::select(week,obs) %>%
  mutate(variable = "Observed cases",
         coverage = factor(
           case_when(week < ymd("2020-05-17") ~ "Pre-P2 expansion",
                     week >= ymd("2020-05-18") ~ "Post-P2 expansion"),
           levels = c("Pre-P2 expansion","Post-P2 expansion"))) -> plot_obs

infect_prev %>%
  dplyr::select(week, lo_n, estimate_n, hi_n) %>%
  mutate(variable = "ONS Infection survey estimate") %>%
  bind_rows(plot_infects) %>%
  bind_rows(plot_cases) %>%
  bind_rows(plot_obs) -> plot_allvars
 
pal <- c("Total infections" = "grey",
         "Predicted-P1+P2 cases" = "steelblue",
         "Observed cases" = "black",
         "ONS Infection survey estimate" = "indianred")

png(here::here(figdir,"figure5.png"), height = 1500, width = 2400, res = 300)
tiff(here::here("figures","paper","fig5.tif"), height = 1500, width = 2400, res = 300)
ggplot(plot_allvars, aes(week, med, ymin = l1, ymax = h2)) +
  # Reconstructed cases and infections
  geom_ribbon(aes(fill = variable), alpha = 0.3) +
  geom_ribbon(aes(fill = variable, ymin = l2, ymax = h1), alpha = 0.3) +
  geom_line(aes(col = variable)) +
  # ONS-estimated weekly incidence of test-positive infection
  geom_errorbar(aes(ymin = lo_n, ymax = hi_n, col = variable), width = 1) +
  geom_point(aes(y = estimate_n, col = variable)) +
  # Observed cases (by surveillance period)
  geom_point(aes(y = obs, col = variable, pch = coverage)) +
  # Key time points
  geom_vline(xintercept = ymd("2020-03-13"), col = "grey30", lty = "dashed") +
  geom_vline(xintercept = ymd("2020-04-15"), col = "grey30", lty = "dashed") +
  geom_vline(xintercept = ymd("2020-05-18"), col = "grey30", lty = "dashed") +
  annotate("text", x = ymd("2020-05-19"), y = 28e4, 
           label = "P2 available to all symptomatic cases", col = "grey30",
           cex = 2, hjust = "left") +
  annotate("text", x = ymd("2020-04-16"), y = 26e4, 
           label = "P2 available to care home residents and staff", col = "grey30",
           cex = 2, hjust = "left") +
  annotate("text", x = ymd("2020-03-14"), y = 24e4, 
           label = "Testing pivoted from community to hospital need", col = "grey30",
           cex = 2, hjust = "left") +
  # Formatting
  scale_y_continuous(labels = scales::number_format(accuracy = 1000)) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_shape_manual(values = c(1,19), na.translate = FALSE) +
  theme(legend.position = c(0.2,0.7)) +
  labs(x = "", y = "Incidence", col = "", fill = "", shape = "")
dev.off()

# ---------------------------------------------------------------------------- #
# Infections detected per week 

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
  geom_vline(xintercept = ymd("2020-03-13"), col = "grey30", lty = "dashed") +
  geom_vline(xintercept = ymd("2020-04-15"), col = "grey30", lty = "dashed") +
  geom_vline(xintercept = ymd("2020-05-18"), col = "grey30", lty = "dashed") +
  annotate("text", x = ymd("2020-05-19"), y = 15, col = "grey30",
           label = "P2 available to all symptomatic cases", 
           cex = 2, hjust = "left") +
  annotate("text", x = ymd("2020-04-16"), y = 10, col = "grey30",
           label = "P2 available to care home residents and staff", 
           cex = 2, hjust = "left") +
  annotate("text", x = ymd("2020-03-14"), y = 5, col = "grey30",
           label = "Testing pivoted from community to hospital need", 
           cex = 2, hjust = "left") +
  # geom_hline(yintercept = 1) +
  labs(x = "", y = "% infections detected") -> perc_detect_time # % total infections detected per week

png(here::here(figdir,"perc_detect_time.png"), height = 1200, width = 2000, res = 300)
perc_detect_time
dev.off()

# ---------------------------------------------------------------------------- #
# Overall total estimated infections

reconstruct7$total$preds %>%
  dplyr::summarise(obs = sum(obs)) %>% pull(obs) -> obs_tot

reconstruct7$sims %>% 
  dplyr::group_by(sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c)) %>%
  dplyr::ungroup() %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>%
  ungroup() %>%
  mutate(obs = obs_tot) %>% 
  dplyr::mutate(across(l1:h2, function(x) obs*100/x, .names = "{.col}_obspred"),
                across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
                across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_tot

print("Estimated total infections:")
dplyr::select(pred_tot, l1_est:h2_est)
# l1_est   l2_est  med_est   h1_est   h2_est
# 1038213. 1088191. 1323622. 1654501. 1737564.

# ---------------------------------------------------------------------------- #
# Infections in March 2020

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
  dplyr::mutate(across(l1:h2, function(x) obs*100/x, .names = "{.col}_obspred"),
                across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
                across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_march

print("Estimated total infections, March 2020:")
dplyr::select(pred_march, l1_est:h2_est)
# l1_est  l2_est med_est  h1_est  h2_est
# 239460. 253010. 306806. 378747. 398601.

print("% infections detected, March 2020:")
dplyr::select(pred_march, l1_est_percdet:h2_est_percdet)
# l1_est_percdet l2_est_percdet med_est_percdet h1_est_percdet h2_est_percdet
# 13.9           13.1            10.8           8.72           8.21

# ~11% (9-13%) infections detected by testing in March 2020

print("% predicted cases detected, March 2020:")
dplyr::select(pred_march, l1_obspred:h2_obspred)
# l1_obspred l2_obspred med_obspred h1_obspred h2_obspred
# 55.1       51.6        42.5       34.5       32.4

# ---------------------------------------------------------------------------- #
# Total estimated infections by geography/region

reconstruct7$geog$preds %>%
  dplyr::group_by(geography) %>%
  dplyr::summarise(obs = sum(obs)) %>% pull(obs) -> obs_geog

reconstruct7$sims %>% 
  dplyr::group_by(geography, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c)) %>%
  dplyr::group_by(geography) %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>%
  ungroup() %>%
  mutate(obs = obs_geog) %>% 
  dplyr::mutate(across(l1:h2, function(x) obs*100/x, .names = "{.col}_obspred"),
                across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
                across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_geog

print("Estimated total infections by geography:")
dplyr::select(pred_geog, l1_est:h2_est)

reconstruct7$la$preds %>%
  dplyr::left_join(dplyr::select(regions, lad19nm, rgn19nm)) %>%
  dplyr::group_by(rgn19nm) %>%
  dplyr::summarise(obs = sum(obs)) %>% pull(obs) -> obs_rgn

reconstruct7$sims %>%
  dplyr::left_join(dplyr::select(regions, lad19nm, rgn19nm)) %>%
  dplyr::group_by(rgn19nm, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c)) %>%
  dplyr::group_by(rgn19nm) %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>%
  ungroup() %>%
  mutate(obs = obs_rgn) %>%
  dplyr::mutate(across(l1:h2, function(x) obs*100/x, .names = "{.col}_obspred"),
                across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
                across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_percdet")) -> pred_rgn

print("Estimated infections by region:")
dplyr::select(pred_rgn, rgn19nm, l1_est:h2_est)
# rgn19nm                   l1_est  l2_est med_est  h1_est  h2_est
# 1 East Midlands            131958. 139341. 176864. 252968. 271843.
# 2 East of England           88180.  93577. 116037. 141894. 151430.
# 3 London                   122976. 141735. 172478. 202789. 233570.
# 4 North East                36572.  39109.  48207.  60228.  64796.
# 5 North West               212280. 229268. 258741. 294696. 315380.
# 6 South East               106817. 113125. 138331. 175212. 184415.
# 7 South West                27010.  28917.  47713.  61298.  65821.
# 8 West Midlands            104508. 110805. 129636. 148193. 156690.
# 9 Yorkshire and The Humber 184770. 199236. 235422. 301367. 335078.

print("% infections detected, by region:")
dplyr::select(pred_rgn, rgn19nm, l1_est_percdet:h2_est_percdet)
# rgn19nm                  l1_est_percdet l2_est_percdet med_est_percdet h1_est_percdet h2_est_percdet
# 1 East Midlands                      15.2           14.4            11.3           7.93           7.38
# 2 East of England                    26.1           24.6            19.9          16.3           15.2 
# 3 London                             27.2           23.6            19.4          16.5           14.3 
# 4 North East                         41.0           38.3            31.1          24.9           23.1 
# 5 North West                         19.6           18.1            16.1          14.1           13.2 
# 6 South East                         31.1           29.4            24.0          19.0           18.0 
# 7 South West                         46.7           43.7            26.5          20.6           19.2 
# 8 West Midlands                      23.8           22.4            19.2          16.8           15.9 
# 9 Yorkshire and The Humber           15.1           14.0            11.9           9.28           8.35


# Commbine into supplementary table S4

format_tab <- function(tab){
  
  tab %>%
    dplyr::mutate(across(l1_est:h2_est, round),
                  across(l1_est_percdet:h2_est_percdet, round, digits = 1)) %>%
    dplyr::mutate(infect_CI98 = paste0(med_est, " [",l1_est," - ", h2_est,"]"),
                  infect_CI50 = paste0(med_est, " [",l2_est, " - ", h1_est,"]"),
                  pdet_CI98 = paste0(med_est_percdet," [",h2_est_percdet, " - ", l1_est_percdet,"]"),
                  pdet_CI50 = paste0(med_est_percdet, " [",h1_est_percdet, " - ", l2_est_percdet,"]")) %>%
    dplyr::select(obs, med_est, infect_CI98, infect_CI50, med_est_percdet, pdet_CI98, pdet_CI50) %>%
    return()
    
}

tabs <- list(pred_tot, pred_geog, pred_rgn)
tabS4 <- bind_rows(lapply(tabs, format_tab))

tabS4$Region <- c("England",pred_geog$geography, pred_rgn$rgn19nm)
tabS4 <- dplyr::select(tabS4, Region, everything())

write.csv(tabS4, here::here(outdir, "TableS4.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# Compare P1 and P1+2

reconstruct7$total$preds %>%
  dplyr::mutate(period = case_when(week < ymd("2020-05-18") ~ "Pillar 1",
                                   week >= ymd("2020-05-18") ~ "Pillar 1+2")) %>% 
  dplyr::group_by(period) %>%
  dplyr::summarise(obs = sum(obs, na.rm = T)) %>%
  dplyr::ungroup() -> obs_period

reconstruct7$total$preds %>% 
  dplyr::mutate(period = factor(case_when(week < ymd("2020-05-18") ~ "Pre-pillar 2 expansion",
                                          week >= ymd("2020-05-18") ~ "Post-pillar 2 expansion"),
                                levels = c("Pre-pillar 2 expansion","Post-pillar 2 expansion"))) %>%
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

# ---------------------------------------------------------------------------- #
## Infections detected per LTLA ##

covs <- readRDS(here::here("data","covs.rds"))

reconstruct7$la$preds %>%
  # dplyr::filter(week < ymd("2020-05-18")) %>%
  dplyr::group_by(lad19nm) %>%
  dplyr::summarise(obs = sum(obs, na.rm = T)) %>%
  ungroup() %>%
  dplyr::full_join(dplyr::select(regions.df, rgn19nm, geography, lad19nm, la_pop, med_age)) %>%
  dplyr::mutate(obs_inc = obs*1e4/la_pop) %>%
  dplyr::left_join(covs) -> obs_la

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

pred_infect_la %>%
  summarise(low = quantile(med_est_percdet, 0.01),
            med = median(med_est_percdet),
            high = quantile(med_est_percdet,0.99))
#  low   med  high
# 7.84  23.4  85.2

pred_infect_la %>%
  group_by(geography) %>%
  summarise(low = quantile(med_est_percdet, 0.01),
            med = median(med_est_percdet),
            high = quantile(med_est_percdet,0.99)) -> tab_geog
tab_geog
# geography                   low   med  high
# London Borough            11.0   19.3  43.4
# Metropolitan District      7.20  18.8  68.9
# Non-metropolitan District  8.40  26.2  93.1
# Unitary Authority          8.43  22.5  60.4

pred_infect_la %>%
  group_by(rgn19nm) %>%
  summarise(low = quantile(med_est_percdet, 0.01),
            med = median(med_est_percdet),
            high = quantile(med_est_percdet,0.99)) -> tab_rgn
tab_rgn
# rgn19nm                    low   med  high
# East Midlands             7.33  12.8  41.4
# East of England           9.90  29.0  67.5
# London                   11.0   19.3  43.4
# North East               15.6   29.6  73.1
# North West                8.24  17.4  74.6
# South East                9.10  28.7  75.9
# South West               10.1   38.9  96.5
# West Midlands            10.6   20.8  47.8
# Yorkshire and The Humber  7.29  16.2  44.6

# Absolute total infections by region
pred_infect_la %>%
  group_by(rgn19nm) %>%
  summarise(low = quantile(med_est, 0.01),
            med = median(med_est),
            high = quantile(med_est,0.99)) -> tab_rgn_abs
tab_rgn_abs


pred_infect_la$rgn19nm[pred_infect_la$rgn19nm == "Yorkshire and The Humber"] <- "Yorkshire and\nThe Humber"
ggplot(pred_infect_la, aes(x = rgn19nm, y = med_est_percdet, 
                           ymin = l1_est_percdet,
                           ymax = h2_est_percdet,
                           col = rgn19nm,
                           fill = rgn19nm)) +
  geom_jitter() +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  guides(col = "none", fill = "none") +
  labs(x = "", y = "% infections detected") +
  ylim(c(0,100)) -> perc_detect_region

png(here::here(figdir,"perc_detect_region.png"), height = 1500, width = 2500, res = 300)
perc_detect_region
dev.off()

png(here::here(figdir,"perc_detect_geog.png"), height = 1500, width = 2000, res = 300)
ggplot(pred_infect_la, aes(x = geography, y = med_est_percdet, 
                           ymin = l1_est_percdet,
                           ymax = h2_est_percdet,
                           col = geography,
                           fill = geography)) +
  geom_jitter() +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  guides(col = "none", fill = "none") +
  labs(x = "", y = "% infections detected") +
  ylim(c(0,100))
dev.off()

pred_infect_la %>%
  ggplot() +
  geom_histogram(aes(x = med_est_percdet, fill = geography), 
                 stat = "bin",alpha = 1, position = "stack") +
  # geom_histogram(aes(x = l1_est_percdet),alpha = 0.2) +
  # geom_histogram(aes(x = h2_est_percdet),alpha = 0.2) +
  labs(x = "% infections detected",
       y = "Frequency",
       fill = "") +
  guides(fill = "none") +
  facet_wrap(~geography, scales = "free_y") -> detect_hist

png(here::here(figdir,"perc_detect_space_hist_lag7.png"), height = 1500, width = 2000, res = 300)
detect_hist
dev.off()

# Compare average % detected by region, geography, LTLA covariates
pred_infect_la %>%
  ggplot(aes(y = med_est_percdet, x= obs_inc, col = rgn19nm, fill = rgn19nm)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", lwd = 0.5, alpha = 0.2) +
  labs(y = "% infections detected",
       x = "Observed incidence per 10,000",
       col = "") +
  guides(fill = "none",
         col = "none") +
  facet_wrap(~rgn19nm) -> perc_detect_inc_rgn

png(here::here(figdir,"perc_detect_la_inc_byrgn_lag7.png"), height = 2000, width = 3000, res = 300)
perc_detect_inc_rgn
dev.off()

pred_infect_la %>%
  ggplot(aes(y = med_est_percdet, x= obs_inc, col = geography, fill = geography)) +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm") +
  labs(y = "% infections detected",
       x = "Observed incidence per 10,000",
       col = "") +
  guides(fill = "none",
         col = "none") +
  facet_wrap(~geography) -> perc_detect_inc_geog
png(here::here(figdir,"perc_detect_la_inc_bygeog_lag7.png"), height = 2000, width = 3000, res = 300)
perc_detect_inc_geog
dev.off()

png(here::here(figdir,"perc_detect_la_imd_bygeog_lag7.png"), height = 2000, width = 3000, res = 300)
pred_infect_la %>%
  ggplot(aes(y = med_est_percdet, x= IMD, col = geography, fill = geography)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "% infections detected",
       x = "Index of multiple deprivation",
       col = "") +
  guides(fill = "none",
         col = "none") +
  facet_wrap(~geography)
dev.off()

png(here::here(figdir,"perc_detect_la_min_bygeog_lag7.png"), height = 2000, width = 3000, res = 300)
pred_infect_la %>%
  ggplot(aes(y = med_est_percdet, x= prop_minority, col = geography, fill = geography)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "% infections detected",
       x = "% minority ethnic population",
       col = "") +
  guides(fill = "none",
         col = "none") +
  facet_wrap(~geography)
dev.off()

png(here::here(figdir,"perc_detect_la_age_bygeog_lag7.png"), height = 2000, width = 3000, res = 300)
pred_infect_la %>%
  ggplot(aes(y = med_est_percdet, x= med_age, col = geography, fill = geography)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "% infections detected",
       x = "Median age",
       col = "") +
  guides(fill = "none",
         col = "none") +
  facet_wrap(~geography)
dev.off()


# Map out percentage infections detected in cumulative totals per LTLA
regions %>% 
  dplyr::full_join(dplyr::select(pred_infect_la, lad19nm, med_est_percdet)) %>% 
  inset_map(fill = "med_est_percdet", fill.lab = "% infections\ndetected", 
            plot.border = F, scale = F) -> perc_detect_space

png(here::here(figdir,"perc_detect_space_lag7.png"), height = 1500, width = 2000, res = 300)
perc_detect_space
dev.off()

library(grid)
library(gridExtra)
# png(here::here(figdir,"detect_space_rgn_time_layout.png"), height = 2600, width = 3000, res = 320)
tiff(here::here("figures","paper","fig6.tif"), height = 2600, width = 3000, res = 320)
gridExtra::grid.arrange(arrangeGrob(perc_detect_space,  left = textGrob("A",
                                                                              x = unit(0.5,"npc"),
                                                                              y = unit(0.9, "npc"),
                                                                              gp=gpar(fontsize=16))),
                        arrangeGrob(perc_detect_inc_rgn, left = textGrob("B",
                                                                       x = unit(0.5,"npc"),
                                                                       y = unit(0.9, "npc"),
                                                                       gp=gpar(fontsize=16))), 
                        arrangeGrob(perc_detect_time, left = textGrob("C",
                                                                            x = unit(0.5,"npc"),
                                                                            y = unit(0.9, "npc"),
                                                                            gp=gpar(fontsize=16))),
                        layout_matrix = rbind(c(1,1,2,2,2),c(3,3,3,3,3)))
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
sample <- sample(regions.df$lad19nm, 9)
  
# Predictions per week and per ltla
reconstruct7$la$preds %>%
  dplyr::mutate(period = factor(case_when(week < ymd("2020-05-18") ~ "Pre-pillar 2 expansion",
                                          week >= ymd("2020-05-18") ~ "Post-pillar 2 expansion"),
                                levels = c("Pre-pillar 2 expansion","Post-pillar 2 expansion"))) %>%
  full_join(dplyr::select(regions.df, lad19nm, la_pop)) %>%
  left_join(infect_prev) %>%
  mutate(across(estimate:hi, function(x) (x/100)*la_pop, .names = "{.col}_n"),
         # across(l1:h2, function(x) x*100/p_detect_all$lo, .names = "{.col}_plo"),
         across(l1:h2, function(x) x*100/p_detect_all$med, .names = "{.col}_est"),
         # across(l1:h2, function(x) x*100/p_detect_all$hi,.names = "{.col}_phi"),
         across(l1_est:h2_est, function(x) return(obs*100/x), .names = "{.col}_propdet")) -> pred_infect

png(here::here(figdir,"reconstruct7_vs_infections_lasamp.png"), height = 2000, width = 2600, res = 300)
# pdf(here::here(figdir,"reconstruct7_la_all.pdf"), height = 60, width = 60)
pred_infect %>%
  filter(lad19nm %in% sample) %>%
  ggplot(aes(week)) +
  geom_ribbon(aes(ymin = l1_est, ymax = h2_est), alpha = 0.3, fill = "grey") +
  geom_line(aes(y = med_est), col = "grey") +
  # Reconstructed cases detected under P1+2
  geom_ribbon(aes(y = med, ymin = l1, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  geom_point(aes(y = obs, shape = period)) +
  # scale_y_continuous(trans = "log10") +
  scale_shape_manual(values = c(1,19), na.translate = FALSE) +
  # theme(legend.position = c(0.1,0.9)) +
  labs(x = "", y = "Incidence", shape = "" ) +
  facet_wrap(~lad19nm, scales = "free_y")
dev.off()

# ---------------------------------------------------------------------------- #
## Compare one/two/three week lags by national total time series ##

reconstruct14 <- readRDS(here::here(outdir,"reconstruct_lag14.rds"))
reconstruct21 <- readRDS(here::here(outdir,"reconstruct_lag21.rds"))

lagcomp <- dplyr::bind_rows(reconstruct7$total$preds, 
                            reconstruct14$total$preds, 
                            reconstruct21$total$preds) %>%
  rename(Lag = lag)

cases %>%
  group_by(week) %>%
  summarise(obs = sum(n),
            coverage = factor(
              case_when(week < ymd("2020-05-18") ~ "Pre-P2 expansion",
                        week >= ymd("2020-05-18") ~ "Post-P2 expansion"),
              levels = c("Pre-P2 expansion","Post-P2 expansion"))) -> obs_wk

# png(here::here(figdir,"compare_lags.png"), height = 1400, width = 2200, res = 300)
tiff(here::here("figures","paper","fig4.tif"), height = 1400, width = 2200, res = 300)
ggplot(lagcomp, aes(week)) + 
  geom_ribbon(aes(ymin = l1, ymax = h1, fill = Lag), alpha = 0.2) +
  geom_ribbon(aes(ymin = l2, ymax = h2, fill = Lag), alpha = 0.2) +
  geom_line(aes(y = med, colour = Lag)) +
  geom_point(data = obs_wk, aes(y = obs, pch = coverage)) +
  scale_fill_viridis_d(begin = 0.2, end = 0.8) +
  scale_colour_viridis_d(begin = 0.2, end = 0.8) +
  scale_x_date(breaks = "month", date_labels = "%b") +
  scale_shape_manual(values = c(1,19), na.translate = FALSE) +
  labs(x = "",y = "Incidence", shape = "" 
       #, title = "Reconstruction of confirmed cases from COVID-19-related deaths across England",
       # subtitle = "Comparison of assumed lags between case confirmation and death"
  ) +
  theme_minimal() +
  theme(legend.position = c(0.8,0.7))
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
################################################################################
