library(tidyverse)
library(lubridate)
library(dtplyr)

sims_long <- readRDS(here::here("output","sims_long_avgcov.rds"))

deaths <- readRDS(here::here("data","expanded","deaths.rds"))
deaths <- deaths[[1]]
cases <- readRDS(here::here("data","expanded","cases.rds"))
cases <- cases[[1]]

# ---------------------------------------------------------------------------- #
## Scale predictions by death:case ratio ##

# lag <- 14
set.seed(101)

scale_quants <- seq(0.1,0.9,0.1)
plot_quants <- c(0.05,0.25,0.75,0.95)

get_cfr <- function(deaths, cases, lag, scale_quants, denom_cutoff = 1){
  
  deaths %>%
    mutate(n = replace_na(n, 0),
           week = week-lag) %>%
    full_join(dplyr::select(cases,lad19cd,geography,week,n), by = c("lad19cd","week","geography"), suffix = c("_d","_c")) %>%
    mutate(CFR_obs = n_c/n_d,
           period = case_when(week < ymd("2020-05-18") ~ 1,
                              week >= ymd("2020-05-18") ~ 2),
           ID = row_number()) -> ratio
  
  # Drop unstable ratios where denominator less than cutoff
  ratio$CFR_obs[ratio$n_d < denom_cutoff | ratio$n_c < denom_cutoff] <- NA
  
  ratio %>% 
    summarise(n_d = sum(n_d, na.rm = TRUE),
              n_c = sum(n_c, na.rm = TRUE),
              med = median(CFR_obs, na.rm = TRUE),
              q1 = quantile(CFR_obs, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR_obs, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  ratio %>% 
    group_by(period) %>%
    summarise(n_d = sum(n_d, na.rm = TRUE),
              n_c = sum(n_c, na.rm = TRUE),
              med = median(CFR_obs, na.rm = TRUE),
              q1 = quantile(CFR_obs, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR_obs, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  
  ratio %>% 
    group_by(geography,period) %>%
    summarise(n_d = sum(n_d, na.rm = TRUE),
              n_c = sum(n_c, na.rm = TRUE),
              med = median(CFR_obs, na.rm = TRUE),
              q1 = quantile(CFR_obs, p = 0.25, na.rm = TRUE),
              q3 = quantile(CFR_obs, p = 0.75, na.rm = TRUE)) %>%
    print()
  
  # Plot distribution of observed CFR
  CFR_obs <- ratio$CFR_obs[!is.na(ratio$CFR_obs) & ratio$period == 2]
  x <- 0:max(CFR_obs)
  hist(CFR_obs, breaks = 100, prob = T)
  lines(x, EnvStats::demp(x, CFR_obs), col = "red")
  
  quants <- quantile(ratio$CFR_obs[ratio$period == 2], probs = scale_quants, na.rm = TRUE)
  quants
  
  return(list(quants = quants, ratio = ratio))
  
}
  
ratio7 <- get_cfr(deaths, cases, lag = 7, scale_quants)
ratio14 <- get_cfr(deaths, cases, lag = 14, scale_quants)



reconstruct <- function(sims, quants){
  
  # For each posterior sample, rescale by specified quantiles of the CFR distribution
  scaled_sims <- bind_rows(lapply(quants, lag_rescale, lag = lag, sims = sims)) %>%
    lazy_dt() %>%
    # scaled_sims <- lag_rescale(lag = lag, ratiodist = CFR_obs, sims = sims) %>%
    group_by(lad19nm, week) %>%
    mutate(sim = row_number(),
           lag = paste(lag/7,"weeks")) %>%
    ungroup() %>%
    rename(n_d = n) %>%
    left_join(dplyr::select(cases[["first"]], lad19nm, geography, week, n))
  
  scaled_quants <- scaled_sims %>%
    group_by(lad19nm,la, week) %>%
    summarise(l1 = quantile(pred_n_scale, plot_quants[1]),
              l2 = quantile(pred_n_scale, plot_quants[2]),
              med = quantile(pred_n_scale, 0.5),
              h1 = quantile(pred_n_scale, plot_quants[3]),
              h2 = quantile(pred_n_scale, plot_quants[4]),
              obs = unique(n))
  
  return(list(scaled_sims = scaled_sims, scaled_quants = scaled_quants))
  
}

reconstruct7 <- reconstruct(sims_long, ratio7$quants)

plot_reconst <- function(scaled_sims, plot_quants, suffix = ""){
  
  la_samp <- sample(scaled_quants$lad19nm, size =  4)
  png(here::here("figures","compare","expanded",paste0("reconstr_lasamp_lag_",lag,"_",suffix,".png")), height = 1000, width = 1500, res = 150)
  print(
    scaled_quants %>%
      filter(lad19nm %in% la_samp) %>% #View()
      as.data.frame() %>%
      ggplot(aes(week)) + 
      geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
      geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
      geom_line(aes(y = med), col = "steelblue") +
      scale_x_date(limits = range(data$week)) +
      geom_point(aes(y = obs)) +
      facet_wrap(~lad19nm, scales = "free") +
      labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths",
           subtitle = "Four sampled LTLAs",
           caption = paste0("Median, ",
                            (plot_quants[2]-plot_quants[1])*100, " and ",
                            (plot_quants[2]-plot_quants[1])*100,
                            "% quantile intervals over ", nsims, " posterior simulations, scaled by ",
                            # samples_cfr, 
                            # "samples from",
                            min(scale_quants)*100,"% to ", max(scale_quants)*100,"% quantiles of",
                            " the observed CFR distribution post-P2 expansion."
           )
      ) +
      theme_minimal()
  )
  dev.off()
  
  scaled_quants_geog <- scaled_sims %>%
    group_by(geography, week, sim) %>%
    summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
              cases = sum(n, na.rm = TRUE)) %>%
    group_by(geography, week) %>%
    summarise(l1 = quantile(pred_n_scale, plot_quants[1]),
              l2 = quantile(pred_n_scale, plot_quants[2]),
              med = quantile(pred_n_scale, 0.5),
              h1 = quantile(pred_n_scale, plot_quants[3]),
              h2 = quantile(pred_n_scale, plot_quants[4]),
              obs = unique(cases, na.rm = TRUE)) 
  
  png(here::here("figures","compare","expanded",paste0("reconstr_geog_lag_",lag,"_",suffix,".png")), height = 1000, width = 1500, res = 150)
  print(
    scaled_quants_geog %>%
      as.data.frame() %>%
      ggplot(aes(week)) + 
      geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
      geom_line(aes(y = med), col = "steelblue") +
      scale_x_date(limits = range(data$week)) +
      geom_point(aes(y = obs)) +
      facet_wrap(~geography, scales = "free_y") +
      labs(x = "",y = "Confirmed case count", title = "Reconstruction of confirmed cases from COVID-19-related deaths",
           subtitle = "by geography type",
           caption = paste0("Median, ",
                            (plot_quants[2]-plot_quants[1])*100, " and ",
                            (plot_quants[2]-plot_quants[1])*100,
                            "% quantile intervals over ", nsims, " posterior simulations, scaled by ",
                            # samples_cfr, 
                            # "samples from",
                            min(scale_quants)*100,"% to ", max(scale_quants)*100, "% quantiles of",
                            " the observed CFR distribution post-P2 expansion.")) +
      theme_minimal()
  )
  dev.off()
}

plot_quants
scaled_sims0 <- reconstruct(sims = dat_sims, data = dat, lag = 0, scale_quants = scale_quants, plot_quants = plot_quants, suffix = "avgcov")
scaled_sims7 <- reconstruct(sims = dat_sims, data = dat, lag = 7, scale_quants = scale_quants, plot_quants = plot_quants, suffix = "avgcov")
scaled_sims14 <- reconstruct(sims = dat_sims, data = dat, lag = 14, scale_quants = scale_quants, plot_quants = plot_quants, suffix = "avgcov")
scaled_sims21 <- reconstruct(sims = dat_sims, data = dat, lag = 21, scale_quants = scale_quants, plot_quants = plot_quants, suffix = "avgcov")

lagcomp <- bind_rows(scaled_sims7, scaled_sims14)

png(here::here("figures","compare","expanded","reconstr_totals.png"), height = 1200, width = 2800, res = 250)
lagcomp %>%
  group_by(lag, week, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(lag,week) %>%
  summarise(l1 = quantile(pred_n_scale, plot_quants[1]),
            l2 = quantile(pred_n_scale, plot_quants[2]),
            med = quantile(pred_n_scale, 0.5),
            h1 = quantile(pred_n_scale, plot_quants[3]),
            h2 = quantile(pred_n_scale, plot_quants[4]),
            cases = unique(cases, na.rm = TRUE)) %>%
  ggplot(aes(week, colour = lag)) + 
  geom_ribbon(aes(ymin = l1, ymax = h1), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(ymin = l2, ymax = h2), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = med), col = "steelblue") +
  scale_x_date() +
  geom_point(aes(y = cases)) +
  facet_wrap(~lag) +
  labs(x = "",y = "", title = "Reconstruction of confirmed cases from COVID-19-related deaths across England",
       subtitle = "Comparison of assumed lags between case confirmation and death",
       caption = paste0("Median, ",
                        (plot_quants[2]-plot_quants[1])*100, " and ",
                        (plot_quants[2]-plot_quants[1])*100,
                        "% quantile intervals over ", nsims, " posterior simulations, scaled by ",
                        # samples_cfr, 
                        "samples",
                        # min(scale_quants)*100,"% to ", max(scale_quants)*100, "% quantiles",
                        " from the observed CFR distribution post-P2 expansion.")) +
  theme_minimal() 
dev.off()

## Total cases overall
scaled_sims7 %>%
  group_by(sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            observed = unique(cases, na.rm = TRUE)) -> totals
totals
#     low     med    high observed
# 107783. 334155. 884277.   229758

dplyr::select(totals, -observed) - totals$observed
#       low      med     high
# -121974.7 104396.6 654519.5

totals %>%
  mutate(across(-observed, function(x) (x - totals$observed)*100/totals$observed)) %>%
  as.data.frame()

## By geography
scaled_sims7 %>%
  group_by(geography, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(geography) %>%
  summarise(low = quantile(pred_n_scale, plot_quants[1], na.rm = TRUE),
            med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            high = quantile(pred_n_scale, plot_quants[2], na.rm = TRUE),
            observed = unique(cases, na.rm = TRUE)) %>%
  ungroup() -> geog_totals
geog_totals

geog_totals %>%
  mutate(across(c(-geography,-observed), function(x) (x - geog_totals$observed)*100/geog_totals$observed)) %>%
  as.data.frame()

# ---------------------------------------------------------------------------- #

# Map out the median differences
scaled_sims7 %>%
  group_by(lad19cd, sim) %>%
  summarise(pred_n_scale = sum(pred_n_scale, na.rm = TRUE),
            cases = sum(n, na.rm = TRUE)) %>%
  group_by(lad19cd) %>%
  summarise(med = quantile(pred_n_scale, 0.5, na.rm = TRUE),
            cases = unique(cases, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(ratio = (med/cases)) -> scale_by_la

summary(scale_by_la$ratio)

scale_by_la %>%
  summarise(low = quantile(ratio, plot_quants[1], na.rm = TRUE),
            med = quantile(ratio, 0.5, na.rm = TRUE),
            high = quantile(ratio, plot_quants[2], na.rm = TRUE)) %>%
  as.data.frame()

png(here::here("figures","compare","expanded","underascertainment_lag_7.png"), height = 1500, width = 2000, res = 250)
regions %>%
  full_join(scale_by_la) %>%
  basic_map(fill = "ratio", plot.border = T) + 
  scale_fill_gradient2(midpoint = 0, trans = "log2") + 
  labs(fill = "Ratio", title = "Relative difference between median inferred total cases\nand observed test positives",
       subtitle = paste0("Predicted deaths back-dated by one week"))
dev.off()

###############################################################################
###############################################################################