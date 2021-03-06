################################################################################
# Description: Summarise and visualise model output for selected model
# 
# Author: Emily S Nightingale
# Date created: 08/04/2021
# 
################################################################################
################################################################################

################################################################################
# SETUP
################################################################################

# LTLA-week-aggregated observed deaths, expected deaths and LTLA covariates
# (first and second waves)
dat_all <- readRDS(here::here("data","deaths.rds"))

dat <- dat_all[[wave]] 
dat$n[dat$wk_since_first < 0] <- NA
period <- dat_all$breaks[[wave]]

dat_tot <- dat %>%
  group_by(la, lad19cd, lad19nm, geography) %>%
  summarise(n = sum(n, na.rm = T),
            E = unique(E))

weekrange <- seq(min(dat$w), max(dat$w))
weekseq <- seq(min(dat$week), max(dat$week), by = "week")

# Fitted models and posterior samples
fit_final <- readRDS(file = here::here("output",
                                  sprintf("fits_%s_%s.rds","deaths", wave)))[[6]]
samples_final <- readRDS(file = here::here("output",
                                     sprintf("samples_%s_%s.rds","deaths", wave)))[[6]]

nsims <- length(samples_final)

# ---------------------------------------------------------------------------- #
# All parameter posteriors

pdf(here::here(figdir,"parms_fixed.pdf"))
lapply(names(fit_final$marginals.fixed), plot_parm, fit = fit_final, opt = 1)
dev.off()
pdf(here::here(figdir,"parms_hyper.pdf"))
lapply(names(fit_final$marginals.hyperpar), fit = fit_final, plot_parm, opt = 2)
dev.off()

# ---------------------------------------------------------------------------- #
# Fixed effects

# IMD quintiles
post_fixed <- lapply(fit_final$marginals.fixed, function(x) inla.tmarginal(exp, x)) 
lapply(post_fixed, inla.zmarginal)

# Prop minority (rescale to %)
inla.zmarginal(inla.tmarginal(function(x) exp(x/100), fit_final$marginals.fixed[[6]]))

get_coeffs(fit_final)

png(here::here(figdir,"covariates_final.png"), height = 800, width = 1000, res = 150)
fit_final %>% 
  get_coeffs() %>% 
  slice(-1) %>%
  ggplot(aes(x = effect, y = mean,  ymin = `quant0.025`, ymax = `quant0.975`)) +
  geom_pointrange() +
  geom_hline(aes(yintercept = 1), lty = "dashed",col = "red") +
  labs(y = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, vjust =0.6))
dev.off()

# ---------------------------------------------------------------------------- #
# Random effects

## Temporal

dat_w <- data.frame(w = weekrange,
                    trend=fit_final$summary.random$w$`0.5quant`,
                    Model = "Calendar week") %>%
  bind_rows(bind_cols(expand.grid(w = seq(min(dat$wk_since_first),max(dat$wk_since_first)), 
                                  geography = sort(unique(dat$geography))),
                      data.frame(trend = fit_final$summary.random$wk_since_first$`0.5quant`,
                                 Model = "Week since first death"))) %>%
  mutate(geography = factor(replace_na(as.character(geography), "All"),
                            levels = c("London Borough","Metropolitan District","Non-metropolitan District","Unitary Authority","All")))

dat_w %>%
  ggplot(aes(w, trend, col = geography)) +
  geom_line() +
  geom_hline(yintercept = 0, col = "grey", lty = "dashed") +
  facet_grid(rows = vars(Model), scales = "free") + 
  labs(x = "", col = "Geography", y = "Trend") +
  theme_minimal() -> plot_rw

png(here::here(figdir,"temp_re.png"), height = 1000, width = 1800, res = 200)
plot_rw
dev.off()


## Spatial

sp_re <- data.frame(lad19cd = unique(dat$lad19cd), 
                    tau = fit_final$summary.hyperpar["Precision for la","0.5quant"],
                    phi = fit_final$summary.hyperpar["Phi for la","0.5quant"],
                    Total = fit_final$summary.random$la[1:nrow(regions),"0.5quant"], 
                    Structured = fit_final$summary.random$la[(nrow(regions)+1):(2*nrow(regions)),"0.5quant"]) %>% 
  mutate(Unstructured = (Total*sqrt(tau) - Structured*sqrt(phi))/sqrt(1-phi)) %>%
  pivot_longer(c("Total","Structured","Unstructured")) %>%
  mutate(name = factor(name, levels = c("Structured","Unstructured","Total")))

regions %>%
  full_join(sp_re) %>%
  basic_map(fill = "value", scale = FALSE, plot.border = T) + 
  facet_wrap(~name) +
  labs(subtitle = "Decomposition of fitted spatial random effects", fill = "") +
  scale_fill_gradient2() +
  theme(legend.position = c(0.05,0.5)) -> map_sp_re

png(here::here(figdir,"death_spatial_re.png"), height = 1200, width = 2400, res = 300)
map_sp_re
dev.off()


## Highest IID effects
sp_re %>%
  filter(name == "IID") %>%
  slice_max(order_by = abs(value), n = 9) %>%
  pull(lad19cd) %>%
  unique() -> hi_IID

plot_la_fit <- function(id){
  dat %>%
    bind_cols(fit_final$summary.fitted.values) %>%
    filter(lad19cd %in% id) %>%
    ggplot(aes(week, `0.5quant`)) +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.2, fill = "steelblue") +
    geom_line(col = "steelblue") +
    geom_point(aes(y = n/E_wk), cex = 0.5) + 
    facet_wrap(~lad19nm, scales = "free") +
    labs(x = "Calendar week", 
         y = "Observed/Expected", 
         title = "Model fit over time, by calendar week", 
         subtitle = "Observed rates shown in black, with 50-95% posterior quantiles") %>%
    return()
}

png(here::here(figdir,"ltla_hi_IID.png"), height = 1200, width = 1200, res = 200)
plot_la_fit(hi_IID)
dev.off()

## Lowest IID effects
sp_re %>%
  filter(name == "IID") %>%
  slice_min(order_by = abs(value), n = 9) %>%
  pull(lad19cd) %>%
  unique() -> low_IID

png(here::here(figdir,"ltla_lo_IID.png"), height = 1200, width = 1200, res = 200)
plot_la_fit(low_IID)
dev.off()

# ---------------------------------------------------------------------------- #
# Posterior predictions

sims <- as.data.table(bind_cols(lapply(samples_final, get_preds, dat)))

setDT(dat)
dat_sims <-  dplyr::bind_cols(dat[,.(geography, lad19cd, lad19nm, la, la_pop, week, E_wk, n)], 
                              sims)

dat_sims_long <- reshape2::melt(
  dat_sims, 
  id.vars = 1:9)

dat_sims_long <- dat_sims_long[, pred := exp(value)]
dat_sims_long <- dat_sims_long[, pred_n := exp(value)*E_wk]

dat_plot_tot <- dat_sims_long[,.(pred_n = sum(pred_n), 
                                 obs = sum(n, na.rm = T), 
                                 pop = sum(la_pop)), by = .(week, variable)] 

dat_plot_geog <- dat_sims_long[,.(pred_n = sum(pred_n), 
                                  obs = sum(n, na.rm = T), 
                                  pop = sum(la_pop)), by = .(week, variable, geography)] 
dat_plot_geog <- dat_plot_geog[, group := paste(variable, geography)]


## Overall and by geography ##

ggplot(as.data.frame(dat_plot_tot)) + 
  geom_line(aes(week, pred_n*1e5/pop, group = variable), alpha = 0.1, col = "grey") +
  geom_point(aes(week, obs*1e5/pop)) + 
  scale_x_date(breaks = "month", date_labels = "%b") +
  labs(
       # title = "Total fit over time, by calendar week", 
       # subtitle = paste0("Observed rates shown in black, with ", nsims, " posterior samples in grey")
       x = "", y = "Rate per 100,000") -> plot_fit_time

ggplot(as.data.frame(dat_plot_geog)) + 
  geom_line(aes(week, pred_n*1e5/pop, group = group, col = geography), alpha = 0.1) +
  geom_point(aes(week, obs*1e5/pop, col = geography), pch = 21, fill = "white") +
  scale_x_date(breaks = "month", date_labels = "%b") +
  labs(
    # title = "Total fit over time, by calendar week and geography", 
    # subtitle = paste0("Observed rates shown in white, with ", nsims, " posterior samples"),
    x = "Calendar week", y = "Rate per 100,000", col = "Geography") +
  theme(legend.position = c(0.2,0.7))  -> plot_fit_time_geog

png(here::here(figdir,"fit_total_geog.png"), height = 1800, width = 2400, res = 300)
plot_fit_time / plot_fit_time_geog
dev.off()


## By LTLA ##

agg_sims <- dat_sims_long[,.(q01 = quantile(pred_n, 0.01),
                             q25 = quantile(pred_n, 0.25),
                             q50 = quantile(pred_n, 0.5),
                             q75 = quantile(pred_n, 0.75),
                             q99 = quantile(pred_n, 0.99),
                             obs = mean(n)),
                          by = .(la, lad19cd, lad19nm, la_pop, geography, week)]

calc_inc <- function(x) x*1e5/agg_sims$la_pop
agg_sims_plot <- mutate(agg_sims, across(q05:obs, calc_inc))

plot_fit_la <- function(plotdata){
  ggplot(plotdata, aes(x = week)) + 
    geom_ribbon(aes(ymin = q01, ymax = q98), alpha = 0.2, fill = "steelblue") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = q50), col = "steelblue") +
    geom_point(aes(y = obs), cex = 0.5) + 
    facet_wrap(~lad19nm, scales = "free") +
    labs(x = "Calendar week", 
         y = "Rate per 100,000", 
         title = "Model fit over time, by calendar week", 
         subtitle = paste0("Observed rates shown in black, with 50-98% quantiles over ", nsims, " posterior samples")) %>%
    return()
}

la_samp <- sample(unique(dat$lad19nm),9)
png(here::here(figdir,"fit_lasamp.png"), height = 1000, width = 1500, res = 200)
plot_fit_la(as.data.frame(agg_sims[lad19nm %in% la_samp]))
dev.off()

pdf(here::here(figdir,"fit_all_ltlas.pdf"), height = 25, width = 25)
plot_fit_la(as.data.frame(agg_sims))
dev.off()


# LTLAs with low and high average relative risks

dat$RR <- fit_final$summary.fitted.values[, "0.5quant"]
dat$LL <- fit_final$summary.fitted.values[, "0.025quant"]
dat$UL <- fit_final$summary.fitted.values[, "0.975quant"]

dat_la <- dat %>%
  group_by(lad19cd) %>%
  summarise(Median = mean(RR, na.rm = TRUE), 
            Lower = mean(LL, na.rm = TRUE), 
            Upper = mean(UL, na.rm = TRUE)) %>%
  as.data.frame()

png(here::here(figdir,"fitted_RR_CrI.png"), height = 1000, width = 1800, res = 200)
regions %>%
  full_join(pivot_longer(dat_la, cols = c("Lower","Median","Upper"))) %>%
  basic_map(fill = "value", plot.border = T, scale = F) +
  facet_wrap(~name) +
  scale_fill_gradient2(midpoint = 1) +
  labs(title = "Median fitted relative risk with 95% credible interval")
dev.off()


## LTLAs with highest/lowest fitted relative risk ##
dat_la %>%
  slice_max(n = 6, Median) -> high_RR

dat_la %>%
  slice_min(n = 6, Median) -> low_RR

png(here::here(figdir,"high_RR.png"), width = 1800, height = 1000, res = 200)
plot_fit_la(as.data.frame(agg_sims[lad19cd %in% high_RR$lad19cd]))
dev.off()

png(here::here(figdir,"low_RR.png"), width = 1800, height = 1000, res = 200)
plot_fit_la(as.data.frame(agg_sims[lad19cd %in% low_RR$lad19cd]))
dev.off()

################################################################################
################################################################################
