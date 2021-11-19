################################################################################
# Description: Summarise reconstructed case counts
# 
# Author: Emily S Nightingale
# Date created: 30/04/2021
# 
################################################################################
################################################################################

library(EpiNow2)

outdir <- "output/reconstruct"
figdir <- "figures/reconstruct"

## Shapefiles
regions <- readRDS(here::here("data","LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
border <- sf::st_union(regions)

plot_quants <- c(0.01,0.25,0.75,0.99)

# Load observed confirmed cases
cases <- readRDS(here::here("data","aggregated","cases.rds"))[[1]]

reconstruct7 <- readRDS(here::here(outdir,"reconstruct_lag7.rds"))
reconstruct14 <- readRDS(here::here(outdir,"reconstruct_lag14.rds"))
reconstruct21 <- readRDS(here::here(outdir,"reconstruct_lag21.rds"))

# ---------------------------------------------------------------------------- #
## Calculate quantiles of reconstructed cases

# Sum predicted over LTLA/weeks, then summarise over posteriors/quantiles
reconstruct7$sims %>%
  dplyr::group_by(week, sim, scale_quant) %>%
  dplyr::summarise(pred = sum(pred_c)) %>%
  dplyr::group_by(week) %>%
  dplyr::summarise(l1 = quantile(pred, plot_quants[1]),
                   l2 = quantile(pred, plot_quants[2]),
                   med = quantile(pred, 0.5),
                   h1 = quantile(pred, plot_quants[3]),
                   h2 = quantile(pred, plot_quants[4])) %>% 
  dplyr::ungroup() -> summ_sims

reconstructed_cases <- data.frame(
  date = summ_sims$week,
  confirm = round(summ_sims$med)) 

# set up example generation time
generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
# set delays between infection and case report 
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")
reporting_delay <- list(mean = convert_to_logmean(3, 1), mean_sd = 0.1,
                        sd = convert_to_logsd(3, 1), sd_sd = 0.1, max = 10)

reconstruct_infects <- estimate_infections(
  reported_cases = reconstructed_cases,
  generation_time = generation_time,
  delays = delay_opts(incubation_period, reporting_delay),
  rt = rt_opts(prior = list(mean = 2, sd = 0.1)),
  stan = stan_opts(control = list(adapt_delta = 0.95))
)
saveRDS(reconstruct_infects, here::here(outdir, "reconstruct_infects7.rds"))

backcalc <- estimate_infections(reported_cases = reconstructed_cases, 
                                generation_time = generation_time,
                                delays = delay_opts(incubation_period),
                                rt = NULL, backcalc = backcalc_opts(),
                                obs = obs_opts(week_effect = FALSE),
                                horizon = 0)

# real time estimates
summary(backcalc)

# summary plot
plot(backcalc)





reconstruct7$sims %>%
  dplyr::group_by(lad19nm, week) %>%
  dplyr::summarise(l1 = quantile(pred_c, plot_quants[1]),
                   l2 = quantile(pred_c, plot_quants[2]),
                   med = median(pred_c),
                   h1 = quantile(pred_c, plot_quants[3]),
                   h2 = quantile(pred_c, plot_quants[4]),
                   obs = median(n_c)) %>% 
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(regions.df, lad19nm, lad19cd, la_pop)) %>%
  dplyr::rename(region_code = lad19cd,
                date = week) -> summ_sims

pred_underrep_join <- summ_sims %>%
  dplyr::mutate(date = date - 2) %>%
  dplyr::full_join(under_reporting_fits) %>%
  dplyr::mutate(pred_infect = med/estimate,
                pred_lag = lead(med, 1),
                pred_infect_crude = pred_lag/0.25) %>%
  dplyr::filter(date > ymd("2020-01-14"))

pred_underrep_join %>%
  group_by(date) %>%
  summarise(obs = sum(obs),
            pred = sum(med),
            pred_lag = sum(pred_lag),
            pred_infect_crude = sum(pred_infect_crude),
            pred_infect = sum(pred_infect)) %>%
  ggplot(aes(date, pred_infect_crude)) +
  geom_line(col = "indianred") +
  geom_line(aes(y = pred), col = "steelblue") +
  geom_point(aes(y = obs))

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