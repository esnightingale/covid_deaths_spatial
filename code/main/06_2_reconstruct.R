################################################################################
# Description: Reconstruct case counts from model-predicted deaths
# 
# Author: Emily S Nightingale
# Date created: 18/03/2021
# 
################################################################################
################################################################################

library(tidyverse)
list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))
border <- st_union(regions)

set.seed(101)
scale_quants <- seq(0.1,0.9,0.1)
plot_quants <- c(0.01,0.25,0.75,0.99)

sims_long <- as.data.frame(readRDS(here::here("output","sims_long_avgcov.rds")))
nsims <- n_distinct(sims_long$variable)

deaths <- readRDS(here::here("data","expanded","deaths.rds"))[[1]]
cases <- readRDS(here::here("data","expanded","cases.rds"))[[1]]

# Check predictions at average covariates/relative to absolute population size rather than age-stratified
sims_long %>% 
  filter(lad19nm %in% sample) %>% #(unique(deaths$lad19nm), 4)
  left_join(dplyr::select(deaths, lad19nm, week, n)) %>% 
  mutate(lad19nm = factor(lad19nm, levels = sample)) %>%
  ggplot(aes(week, pred_n, group = variable)) +
  geom_line(alpha = 0.2, col = "grey") +
  geom_point(aes(y = n)) +
  facet_wrap(~lad19nm, scales = "free")

sims_long %>%
  group_by(week, lad19nm, geography) %>%
  summarise(n = median(pred_n)) %>%
  ungroup() -> pred_deaths

# ---------------------------------------------------------------------------- #
## Reconstruct cases from predicted deaths ##

# + Calculate post-P2 CFR between cases and predicted deaths per week/LTLA with assumed lag
# + Rescale posterior samples of death time series according to CFR distribution
# + Aggregate over sims to country/geography/LTLA
# + Plot quantiles alongside observed cases

reconstruct7 <- reconstruct(deaths, cases, sims_long, lag = 7, plot = T, suffix = "")
reconstruct14 <- reconstruct(deaths, cases, sims_long, lag = 14, plot = T, suffix = "")

saveRDS(reconstruct7, here::here("output","reconstruct_lag7.rds"))
saveRDS(reconstruct14, here::here("output","reconstruct_lag14.rds"))

# Plot for a sample of LTLAs
plot_reconst(reconstruct7$la, 7, sample = sample(unique(cases$lad19nm), 9), h = 10, w = 12)
plot_reconst(reconstruct7$total, 7, format = "png", h = 1200, w = 1800)

# ---------------------------------------------------------------------------- #
## Summarise observed/predicted per LTLA and map ##

dplyr::as_tibble(reconstruct7$sims) %>%
  group_by(lad19nm, sim, scale_quant) %>%
  summarise(pred = sum(pred_n_scale),
            obs = sum(n, na.rm = T)) %>%
  group_by(lad19nm) %>%
  summarise(med = median(pred),
            obs = median(obs)) %>% 
  ungroup() %>%
  mutate(obs_med = obs/med) -> obs_pred_diffs

summary(obs_pred_diffs$obs_med)

png(here::here("figures","compare","expanded","underascertainment_lag_7.png"), height = 1500, width = 2000, res = 250)
regions %>% 
  full_join(obs_pred_diffs) %>%
  basic_map(fill = "obs_med", plot.border = T) +
  scale_fill_gradient2(midpoint = 0, trans = "log2") +
  theme(legend.position = c(0.2, 0.5))
dev.off()

## By period ##
# dplyr::as_tibble(reconstruct7$sims) %>%
#   mutate(period = case_when(week < ymd("2020-05-18") ~ 1,
#                             week >= ymd("2020-05-18") ~ 2)) %>%
#   group_by(lad19nm, period, sim, scale_quant) %>%
#   summarise(pred = sum(pred_n_scale),
#             obs = sum(n, na.rm = T)) %>%
#   group_by(lad19nm, period) %>%
#   summarise(med = median(pred),
#             obs = median(obs)) %>% 
#   ungroup() %>%
#   mutate(obs_med = obs/med) -> obs_pred_diffs
# 
# obs_pred_diffs %>%
#   ggplot(aes(x = obs_med)) +
#   geom_histogram() +
#   geom_vline(aes(xintercept = 1)) +
#   geom_vline(aes(xintercept = median(obs_med)), col = "red",lty = "dashed") +
#   geom_vline(aes(xintercept = mean(obs_med)), col = "blue",lty = "dashed") +
#   scale_x_continuous(trans = "log2") +
#   facet_wrap(~period, scales = "free")
# 
# png(here::here("figures","compare","expanded","underascertainment_byper_lag_7.png"), height = 1500, width = 2000, res = 250)
# regions %>% 
#   full_join(obs_pred_diffs) %>%
#   basic_map(fill = "obs_med", plot.border = T) +
#   scale_fill_gradient2(midpoint = 0, trans = "log2") +
#   facet_wrap(~period)
# dev.off()

# dplyr::as_tibble(reconstruct7$la$preds) %>%
#   dplyr::mutate(obs = replace_na(obs, 0),
#                 obs_med_ratio = obs/med,
#                 obs_med_diff1 =  abs(obs_med_ratio - 1),
#                 se = (obs - med)^2,
#                 period = case_when(week < ymd("2020-05-18") ~ 1,
#                                    week >= ymd("2020-05-18") ~ 2)) -> obs_pred_diffs
# 
# obs_pred_diffs %>%
#   ggplot(aes(x = obs_med_ratio)) +
#   geom_histogram() +
#   geom_vline(aes(xintercept = 1)) +
#   geom_vline(aes(xintercept = median(obs_med_ratio)), col = "red",lty = "dashed") +
#   geom_vline(aes(xintercept = mean(obs_med_ratio)), col = "blue",lty = "dashed") +
#   scale_x_continuous(trans = "log2") +
#   facet_wrap(~period, scales = "free")
# 
# obs_pred_diffs %>%
#   group_by(lad19nm, period) %>%
#   summarise(avg_diff = median(obs_med_ratio*is.finite(obs_med_ratio), na.rm = T)) -> avg_diffs_period
# 
# png(here::here("figures","compare","expanded","underascertainment_byper_lag_7.png"), height = 1500, width = 2000, res = 250)
# regions %>% 
#   full_join(avg_diffs_period) %>%
#   basic_map(fill = "avg_diff", plot.border = T) +
#   scale_fill_gradient2(midpoint = 0, trans = "log2") +
#   facet_wrap(~period)
# dev.off()

# avg_diffs %>%
#   slice_min(avg_diff, n = 3) %>%
#   pull(lad19nm) -> low
# 
# avg_diffs %>%
#   mutate(diff_1 = abs(avg_diff - 1)) %>% #View()
#   slice_min(diff_1, n = 3) %>% #abs(mse-mean(mse)
#   pull(lad19nm) -> mid
# 
# avg_diffs %>%
#   slice_max(avg_diff, n = 3) %>%
#   pull(lad19nm) -> high
# 
# sample <- c(low, mid, high)
# 
# # Plot for a sample of high/low error LTLAs
# plot_reconst(reconstruct7$la, 7, sample = sample, suffix = "_lohi", h = 1800, w = 2000, order = F, format = "png")

# ---------------------------------------------------------------------------- #
## Compare one/two week lags for total time series ##

lagcomp <- bind_rows(reconstruct7$total$preds, reconstruct14$total$preds)

png(here::here("figures","compare","expanded","reconstr_totals.png"), height = 1200, width = 2800, res = 250)
lagcomp %>%
  ggplot(aes(week)) + 
  geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  scale_x_date() +
  geom_point(aes(y = obs)) +
  facet_wrap(~lag) +
  labs(x = "",y = "", title = "Reconstruction of confirmed cases from COVID-19-related deaths across England",
       subtitle = "Comparison of assumed lags between case confirmation and death",
       caption = paste0("Median, ",
                        (plot_quants[3]-plot_quants[2])*100, " and ",
                        (plot_quants[4]-plot_quants[1])*100,
                        "% quantile intervals over ", nsims, " posterior simulations, scaled by ",
                        # samples_cfr, 
                        # "samples",
                        min(scale_quants)*100,"% to ", max(scale_quants)*100, "% quantiles",
                        " from the observed CFR distribution post-P2 expansion.")) +
  theme_minimal() 
dev.off()


# ---------------------------------------------------------------------------- #

## Total cases overall
reconstruct14$sims %>%
  group_by(sim, scale_quant) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>% #filter(sim == 1) %>% View
  ungroup() %>%
  summarise(low = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[3], na.rm = TRUE),
            observed = median(cases, na.rm = TRUE)) -> totals
totals
# low     med    high observed
# <dbl>   <dbl>   <dbl>    <dbl>
#   1 160124. 263323. 410893.   234302

dplyr::select(totals, -observed) - totals$observed

totals %>%
  mutate(across(-observed, function(x) (x - totals$observed)*100/totals$observed)) %>%
  as.data.frame()
# low      med     high observed
# 1 -31.65919 12.38598 75.36876   234302

## By geography
reconstruct14$sims %>%
  group_by(geography, sim, scale_quant) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(geography) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[3], na.rm = TRUE),
            observed = unique(cases, na.rm = TRUE)) %>%
  ungroup() -> geog_totals
geog_totals
# geography                    low    med    high observed
# * <chr>                      <dbl>  <dbl>   <dbl>    <dbl>
#   1 London Borough            29956. 50710.  81360.    33582
# 2 Metropolitan District     38149. 62983.  98353.    64992
# 3 Non-metropolitan District 55499. 91161. 142132.    80077
# 4 Unitary Authority         35579. 58418.  90758.    55651


geog_totals %>%
  mutate(across(c(-geography,-observed), function(x) (x - geog_totals$observed)*100/geog_totals$observed)) %>%
  as.data.frame()
# geography       low       med      high observed
# 1            London Borough -10.79867 51.002137 142.27252    33582
# 2     Metropolitan District -41.30153 -3.091206  51.33037    64992
# 3 Non-metropolitan District -30.69301 13.841407  77.49418    80077
# 4         Unitary Authority -36.06744  4.972474  63.08426    55651

# ---------------------------------------------------------------------------- #


###############################################################################
###############################################################################