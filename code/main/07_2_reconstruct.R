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

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
border <- sf::st_union(regions)

set.seed(101)
scale_quants <- seq(0.1,0.9,0.1)
plot_quants <- c(0.01,0.25,0.75,0.99)

sims_long <- as.data.frame(readRDS(here::here("output","sims_long_avgcov.rds")))
nsims <- dplyr::n_distinct(sims_long$variable)

deaths <- readRDS(here::here("data","deaths.rds"))[[1]]
cases <- readRDS(here::here("data","cases.rds"))[[1]]

# Check predictions at average covariates/relative to absolute population size rather than age-stratified
sims_long %>% 
  dplyr::filter(lad19nm %in% la_samp) %>% #(unique(deaths$lad19nm), 4)
  dplyr::left_join(dplyr::select(deaths, lad19nm, week, n)) %>% 
  dplyr::mutate(lad19nm = factor(lad19nm, levels = la_samp)) %>%
  ggplot(aes(week, pred_n, group = variable)) +
  geom_line(alpha = 0.2, col = "grey") +
  geom_point(aes(y = n)) +
  facet_wrap(~lad19nm, scales = "free")

sims_long %>%
  dplyr::group_by(week, lad19nm, geography) %>%
  dplyr::summarise(n = median(pred_n)) %>%
  dplyr::ungroup() -> pred_deaths

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
# plot_reconst(reconstruct7$la, 7, sample = la_samp, save = F, h = 10, w = 12)

# ---------------------------------------------------------------------------- #
## Summarise observed/predicted per LTLA and map ##

reconstruct7$sims %>%
  dplyr::group_by(lad19nm, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_n_scale),
            obs = sum(n, na.rm = T)) %>%
  dplyr::group_by(lad19nm) %>%
  dplyr::summarise(med = median(pred),
            obs = median(obs)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(obs_med = obs/med) -> obs_pred_diffs

summary(obs_pred_diffs$obs_med)

png(here::here("figures","compare","expanded","underascertainment_lag_7.png"), height = 1500, width = 2000, res = 250)
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
  dplyr::summarise(obs = sum(n, na.rm = T)) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(obs = median(obs)) %>% 
  ggplot(aes(week, obs)) + geom_line() #pull(obs) -> obs

reconstruct7$sims %>%
  dplyr::group_by(week, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_n_scale)) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(across(l1:h2, function(pred) return(obs/pred))) -> obs_pred_time

summary(obs_pred_time)

scaleFUN <- function(x) sprintf("%.1f", x)
png(here::here("figures","compare","expanded","underascertainment_time_lag7.png"), height = 1000, width = 1400, res = 250)
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
  labs(x = "", y = "Ratio observed:predicted cases") 
dev.off()

# ---------------------------------------------------------------------------- #
## Compare one/two week lags for total time series ##

lagcomp <- dplyr::bind_rows(reconstruct7$total$preds, reconstruct14$total$preds)

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

# reconstruct7$sims %>% 
#   dplyr::filter(sim %in% 1:100) %>%
#   dplyr::group_by(week, sim, scale_quant) %>%
#   dplyr::summarise(pred_n_scale = sum(pred_n_scale),
#                    obs = sum(n, na.rm = T)) %>%
#   ggplot(aes(week, pred_n_scale, group = sim)) +
#   geom_line(alpha = 0.5, col = "grey") +
#   geom_point(aes(y = obs)) 


# ---------------------------------------------------------------------------- #

## Total cases overall
reconstruct7$sims %>%
  dplyr::group_by(sim, scale_quant) %>%
  dplyr::summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>% #filter(sim == 1) %>% View
  dplyr::ungroup() %>%
  dplyr::summarise(low = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
                   med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
                   high = quantile(pred_n_scale, plot_quants[3], na.rm = TRUE),
                   observed = median(cases, na.rm = TRUE)) -> totals
totals
# low     med    high observed
# <dbl>   <dbl>   <dbl>    <dbl>
#   1 160124. 263323. 410893.   234302

dplyr::select(totals, -observed) - totals$observed

totals %>%
  dplyr::mutate(across(-observed, function(x) (x - totals$observed)*100/totals$observed)) %>%
  as.data.frame()
# low      med     high observed
# 1 -31.65919 12.38598 75.36876   234302

## By geography
reconstruct7$sims %>%
  dplyr::group_by(geography, sim, scale_quant) %>%
  dplyr::summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  dplyr::group_by(geography) %>%
  dplyr::summarise(low = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
                   med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
                   high = quantile(pred_n_scale, plot_quants[3], na.rm = TRUE),
                   observed = unique(cases, na.rm = TRUE)) %>%
  dplyr::ungroup() -> geog_totals
geog_totals
# geography                    low    med    high observed
# * <chr>                      <dbl>  <dbl>   <dbl>    <dbl>
#   1 London Borough            29956. 50710.  81360.    33582
# 2 Metropolitan District     38149. 62983.  98353.    64992
# 3 Non-metropolitan District 55499. 91161. 142132.    80077
# 4 Unitary Authority         35579. 58418.  90758.    55651


geog_totals %>%
  dplyr::mutate(across(c(-geography,-observed), function(x) (x - geog_totals$observed)*100/geog_totals$observed)) %>%
  as.data.frame()
# geography       low       med      high observed
# 1            London Borough -10.79867 51.002137 142.27252    33582
# 2     Metropolitan District -41.30153 -3.091206  51.33037    64992
# 3 Non-metropolitan District -30.69301 13.841407  77.49418    80077
# 4         Unitary Authority -36.06744  4.972474  63.08426    55651


###############################################################################
###############################################################################