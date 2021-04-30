################################################################################
# Description: Reconstruct case counts from model-predicted deaths
# 
# Author: Emily S Nightingale
# Date created: 18/03/2021
# 
################################################################################
################################################################################

library(tidyverse)

# Source functions
list.files(here::here("code","utils"), full.names = TRUE) %>% purrr::walk(source)

# Set directory for figures and outputs
figdir <- "figures/final"
outdir <- "output/final"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
border <- sf::st_union(regions)

# Specify quantiles for CFR scaling and plotting
scale_quants <- c(0.25, 0.75)
plot_quants <- c(0.01,0.25,0.75,0.99)

# Load posterior predictions from death model, given average values of 
# case-fatality risk factors
sims_long <- as.data.frame(readRDS(here::here("output","sims_long_avgcov.rds")))
nsims <- dplyr::n_distinct(sims_long$variable)

# Load observed confirmed cases
cases <- readRDS(here::here("data","cases.rds"))[[1]]

# # Check predictions at average covariates/relative to absolute population size 
# # rather than age-stratified
# sims_long %>% 
#   dplyr::filter(lad19nm %in% la_samp) %>% #(unique(deaths$lad19nm), 4)
#   dplyr::left_join(dplyr::select(deaths, lad19nm, week, n)) %>% 
#   dplyr::mutate(lad19nm = factor(lad19nm, levels = la_samp)) %>%
#   ggplot(aes(week, pred_n, group = variable)) +
#   geom_line(alpha = 0.2, col = "grey") +
#   geom_point(aes(y = n)) +
#   facet_wrap(~lad19nm, scales = "free")
# 
# sims_long %>%
#   dplyr::group_by(week, lad19nm, geography) %>%
#   dplyr::summarise(n = median(pred_n)) %>%
#   dplyr::ungroup() -> pred_deaths

# ---------------------------------------------------------------------------- #
## Reconstruct cases from predicted deaths ##

# Four steps performed by helper functions in /utils:
# + [get_cfr] Calculate post-P2 CFR between cases and predicted deaths per week/LTLA given 
# assumed lag
# + [rescale_sims] Scale each posterior sample of death time series by estimates post-pillar 2 
# CFR (according to specified scale_quants)
# + [agg_sims] Aggregate posteriors over country/geography/LTLA
# + [plot_reconstr] Plot quantiles alongside observed cases

reconstruct7 <- reconstruct(sims_long, cases, lag = 7, plot = F, suffix = "")
reconstruct14 <- reconstruct(sims_long, cases, lag = 14, plot = T, suffix = "")
reconstruct21 <- reconstruct(sims_long, cases, lag = 21, plot = T, suffix = "")

saveRDS(reconstruct7, here::here(outdir,"reconstruct_lag7.rds"))
saveRDS(reconstruct14, here::here(outdir,"reconstruct_lag14.rds"))
saveRDS(reconstruct21, here::here(outdir,"reconstruct_lag21.rds"))

# overall CFR quantiles for 1 week lag:
# period 1  4.11  2.62  6.75
# period 2  5.37  3.07  9.09

# Check plot for a sample of LTLAs
plot_reconst(reconstruct7$la, 7, sample = la_samp, save = F, h = 10, w = 12)

# ---------------------------------------------------------------------------- #
## Summarise ratio of observed:predicted per LTLA and plot ##

# For each LA/posterior/quantile, sum observed and predicted over weeks
reconstruct7$sims %>%
  dplyr::group_by(lad19nm, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c),
                   obs = sum(n_c, na.rm = T)) %>%
  dplyr::group_by(lad19nm) %>%
  dplyr::summarise(med = median(pred),
                   obs = median(obs)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(obs_med = obs/med) -> obs_pred_la

summary(obs_pred_la$obs_med)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2630  0.5960  0.9201  1.0722  1.3349  3.8052  

# Map out relative difference in cumulative totals per LTLA
png(here::here(figdir,"underascertainment_space_lag7.png"), height = 1500, width = 2000, res = 250)
regions %>% 
  dplyr::full_join(obs_pred_la) %>%
  basic_map(fill = "obs_med", plot.border = T) +
  scale_fill_gradient2(midpoint = 0, trans = "log2") +
  theme(legend.position = c(0.2, 0.5))
dev.off()

# ---------------------------------------------------------------------------- #
## Summarise observed/predicted per week and plot ##

# For each week/posterior/quantile, sum observed over LTLAs
reconstruct7$sims %>%
  dplyr::group_by(week, sim, scale_quant) %>%
  dplyr::summarise(obs = sum(n_c, na.rm = T)) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(obs = median(obs)) %>% pull(obs) -> obs
  
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
  dplyr::mutate(obs = obs) %>%
  dplyr::mutate(across(l1:h2, function(pred) return(obs/pred))) -> obs_pred_time 

summary(obs_pred_time)

# Plot relative differences in national total per week 
scaleFUN <- function(x) sprintf("%.1f", x)
png(here::here(figdir,"underascertainment_time_lag7.png"), height = 1200, width = 1800, res = 250)
obs_pred_time %>%
  filter(obs > 0) %>%
  ggplot(aes(week, med)) +
  geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  # scale_x_date(limits = c(ymd("2020-02-01",ymd("2020-06-17")))) +
  scale_y_continuous(trans = "log2", labels = scales::number_format(accuracy = 0.01)) +
  geom_vline(xintercept = ymd("2020-04-15"), col = "red") +
  annotate("text", x = ymd("2020-05-19"), y = 0.3, label = "P2 available to all symptomatic cases", cex = 2, hjust = "left") +
  geom_vline(xintercept = ymd("2020-03-13"), col = "red") +
  annotate("text", x = ymd("2020-04-16"), y = 0.2, label = "P2 available to care home residents and staff", cex = 2, hjust = "left") +
  geom_vline(xintercept = ymd("2020-05-18"), col = "red") +
  annotate("text", x = ymd("2020-03-14"), y = 0.1, label = "Testing pivoted from community to hospital need", cex = 2, hjust = "left") +
  geom_hline(yintercept = 1) +
  labs(x = "", y = "Observed/predicted cases") 
dev.off()

# ---------------------------------------------------------------------------- #
## Compare one/two/three week lags by national total time series ##

lagcomp <- dplyr::bind_rows(reconstruct7$total$preds, reconstruct14$total$preds, reconstruct21$total$preds) %>%
  rename(Lag = lag)

cases %>%
  group_by(week) %>%
  summarise(obs = sum(n)) -> obs_wk

png(here::here(figdir,"compare_lags.png"), height = 1400, width = 2200, res = 250)
ggplot(lagcomp, aes(week)) + 
  geom_ribbon(aes(ymin = l1, ymax = h1, fill = Lag), alpha = 0.2) +
  geom_ribbon(aes(ymin = l2, ymax = h2, fill = Lag), alpha = 0.2) +
  geom_line(aes(y = med, colour = Lag)) +
  geom_point(data = obs_wk, aes(y = obs)) +
  scale_fill_viridis_d(begin = 0.2, end = 0.8) +
  scale_colour_viridis_d(begin = 0.2, end = 0.8) +
  scale_x_date() +
  labs(x = "",y = "", title = "Reconstruction of confirmed cases from COVID-19-related deaths across England",
       subtitle = "Comparison of assumed lags between case confirmation and death") +
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

###############################################################################
###############################################################################