################################################################################
# Description: Reconstruct case counts from model-predicted deaths
# 
# Author: Emily S Nightingale
# Date created: 18/03/2021
# 
################################################################################
################################################################################

library(tidyverse)
list.files(here::here("code","utils"), full.names = TRUE) %>% purrr::walk(source)

figdir <- here::here("figures","compare")

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
border <- sf::st_union(regions)

set.seed(101)
scale_quants <- c(0.25, 0.75)
plot_quants <- c(0.01,0.25,0.75,0.99)


sims_long <- as.data.frame(readRDS(here::here("output","sims_long_avgcov.rds")))
nsims <- dplyr::n_distinct(sims_long$variable)

cases <- readRDS(here::here("data","cases.rds"))[[1]]

# # Check predictions at average covariates/relative to absolute population size rather than age-stratified
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

# + Calculate post-P2 CFR between cases and predicted deaths per week/LTLA with assumed lag
# + Rescale posterior samples of death time series according to CFR distribution
# + Aggregate over sims to country/geography/LTLA
# + Plot quantiles alongside observed cases

reconstruct7 <- reconstruct(sims_long, cases, lag = 7, plot = T, suffix = "nocov")
# 1  4.11  2.62  6.75
# 2  5.37  3.07  9.09
reconstruct14 <- reconstruct(sims_long, cases, lag = 14, plot = T, suffix = "")
reconstruct21 <- reconstruct(sims_long, cases, lag = 21, plot = T, suffix = "")

saveRDS(reconstruct7, here::here("output","reconstruct_lag7.rds"))
saveRDS(reconstruct14, here::here("output","reconstruct_lag14.rds"))
saveRDS(reconstruct21, here::here("output","reconstruct_lag21.rds"))

# Plot for a sample of LTLAs
# plot_reconst(reconstruct7$la, 7, sample = la_samp, save = F, h = 10, w = 12)

# ---------------------------------------------------------------------------- #
## Summarise observed/predicted per LTLA and map ##

reconstruct7$sims %>%
  dplyr::group_by(lad19nm, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c),
                   obs = sum(n_c, na.rm = T)) %>%
  dplyr::group_by(lad19nm) %>%
  dplyr::summarise(med = median(pred),
                   obs = median(obs)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(obs_med = obs/med) -> obs_pred_diffs

summary(obs_pred_diffs$obs_med)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2653  0.6046  0.9239  1.0890  1.3537  3.9174 
# 
png(here::here("figures","compare","underascertainment_lag_7.png"), height = 1500, width = 2000, res = 250)
regions %>% 
  dplyr::full_join(obs_pred_diffs) %>%
  basic_map(fill = "obs_med", plot.border = T) +
  scale_fill_gradient2(midpoint = 0, trans = "log2") +
  theme(legend.position = c(0.2, 0.5))
dev.off()

# ---------------------------------------------------------------------------- #
## Summarise observed/predicted per week and plot ##

reconstruct7$sims %>%
  dplyr::group_by(week, sim, scale_quant) %>%
  dplyr::summarise(obs = sum(n_c, na.rm = T)) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(obs = median(obs)) %>% pull(obs) -> obs
  # ggplot(aes(week, obs)) + geom_line() 

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
  dplyr::mutate(across(l1:h2, function(pred) return(obs/pred), .names = "{.col}_pred")) -> obs_pred_time

summary(obs_pred_time)

scaleFUN <- function(x) sprintf("%.1f", x)
png(here::here("figures","compare","underascertainment_time_lag7.png"), height = 1200, width = 1800, res = 250)
obs_pred_time %>%
  filter(med > 0) %>%
  ggplot(aes(week, med)) +
  geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  scale_y_continuous(trans = "log2", labels = scales::number_format(accuracy = 0.01)) +
  geom_vline(xintercept = ymd("2020-05-18"), col = "red") +
  annotate("text", x = ymd("2020-05-19"), y = 0.25, label = "P2 available to all symptomatic cases", cex = 2, hjust = "left") +
  geom_vline(xintercept = ymd("2020-04-15"), col = "red") +
  annotate("text", x = ymd("2020-04-16"), y = 0.15, label = "P2 available to care home residents and staff", cex = 2, hjust = "left") +
  geom_hline(yintercept = 1) +
  labs(x = "", y = "Observed/predicted cases") 
dev.off()

# ---------------------------------------------------------------------------- #
## Compare one/two week lags for total time series ##

# lagcomp <- dplyr::bind_rows(reconstruct7$total$preds, reconstruct14$total$preds)
# 
# png(here::here("figures","compare","reconstr_totals.png"), height = 1200, width = 2800, res = 250)
# lagcomp %>%
#   ggplot(aes(week)) + 
#   geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
#   geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
#   geom_line(aes(y = med), col = "steelblue") +
#   scale_x_date() +
#   geom_point(aes(y = obs)) +
#   facet_wrap(~lag) +
#   labs(x = "",y = "", title = "Reconstruction of confirmed cases from COVID-19-related deaths across England",
#        subtitle = "Comparison of assumed lags between case confirmation and death",
#        caption = paste0("Median, ",
#                         (plot_quants[3]-plot_quants[2])*100, " and ",
#                         (plot_quants[4]-plot_quants[1])*100,
#                         "% quantile intervals over ", nsims, " posterior simulations, scaled by ",
#                         # samples_cfr, 
#                         # "samples",
#                         min(scale_quants)*100,"% to ", max(scale_quants)*100, "% quantiles",
#                         " from the observed CFR distribution post-P2 expansion.")) +
#   theme_minimal() 
# dev.off()


lagcomp <- dplyr::bind_rows(reconstruct7$total$preds, reconstruct14$total$preds, reconstruct21$total$preds) %>%
  rename(Lag = lag)

cases %>%
  group_by(week) %>%
  summarise(obs = sum(n)) -> obs_wk

png(here::here("figures","compare","reconstr_totals.png"), height = 1400, width = 2200, res = 250)
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

## Total cases overall
reconstruc7$sims %>%
  dplyr::group_by(sim, scale_quant) %>%
  dplyr::summarise(pred_c = sum(pred_c, na.rm = TRUE),
                   cases = sum(n_c, na.rm = TRUE)) %>% #filter(sim == 1) %>% View
  dplyr::ungroup() %>%
  dplyr::summarise(low = quantile(pred_c, plot_quants[1], na.rm = TRUE),
                   med = quantile(pred_c, 0.5, na.rm = TRUE),
                   high = quantile(pred_c, plot_quants[4], na.rm = TRUE),
                   observed = median(cases, na.rm = TRUE)) -> totals
totals
# low     med    high observed
# <dbl>   <dbl>   <dbl>    <dbl>
# 262119. 333647. 437866.   231817

dplyr::select(totals, -observed) - totals$observed

totals %>%
  dplyr::mutate(across(-observed, function(x) (x - totals$observed)*100/totals$observed)) %>%
  as.data.frame()
# low      med     high observed
# 13.0715 43.9269 88.88451   231817

## By geography
reconstruct7$sims %>%
  dplyr::group_by(geography, sim, scale_quant) %>%
  dplyr::summarise(pred_c = sum(pred_c, na.rm = TRUE),
            cases = sum(n_c, na.rm = TRUE)) %>%
  dplyr::group_by(geography) %>%
  dplyr::summarise(low = quantile(pred_c, plot_quants[2], na.rm = TRUE),
                   med = quantile(pred_c, 0.5, na.rm = TRUE),
                   high = quantile(pred_c, plot_quants[3], na.rm = TRUE),
                   observed = unique(cases, na.rm = TRUE)) %>%
  dplyr::ungroup() -> geog_totals
geog_totals
# geography                    low    med    high observed
# 1 London Borough            35611.  43326.  50993.    33399
# 2 Metropolitan District     95379. 109349. 128949.    64007
# 3 Non-metropolitan District 78179.  97145. 121257.    79441
# 4 Unitary Authority         67491.  83499. 110771.    54970


geog_totals %>%
  dplyr::mutate(across(c(-geography,-observed), function(x) (x - geog_totals$observed)*100/geog_totals$observed)) %>%
  as.data.frame()
# geography       low       med      high observed
# 1            London Borough  6.622555 29.72227  52.67792    33399
# 2     Metropolitan District 49.013866 70.83977 101.46129    64007
# 3 Non-metropolitan District -1.588752 22.28592  52.63794    79441
# 4         Unitary Authority 22.778411 51.90000 101.51212    54970


###############################################################################
###############################################################################