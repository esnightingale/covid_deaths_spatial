################################################################################
# Description: Summarise and visualise model output
# 
# Author: Emily S Nightingale
# Date created: 30/09/2020
# 
################################################################################
################################################################################

################################################################################
# SETUP
################################################################################

library(data.table)
library(tidyverse)
library(dtplyr)

dat <- readRDS(here::here("data","expanded","deaths.rds"))
dat <- dat[[wave]] 
# dat$n[dat$wk_since_first < 0] <- NA

# Setup adjusted prediction data
setDF(dat) %>%
  mutate(prop_minority = median(prop_minority),
         IMD_quint = "(21.3,29.2]",
         E_wk = E_wk_unstrat,
         n = NA)  -> dat_avgcov

# Make pred data with average covariate values
dat_pred <- bind_rows(dat, dat_avgcov) %>%
  mutate(IMD_quint = factor(IMD_quint, levels = levels(dat$IMD_quint)))

pred_avgcov <- readRDS(here::here("output","fit_samples_expanded_avgcov.rds"))
fit_pred <- pred_avgcov$fit
samples_pred <- pred_avgcov$samples

sims <- as.data.table(bind_cols(lapply(samples_pred, get_preds, dat_pred)))
saveRDS(sims, file = here::here("output","sims_expanded_avgcov.rds"))

sims <- readRDS(file = here::here("output","sims_expanded_avgcov.rds"))
data.table::setDT(sims)

cases <- readRDS(here::here("data","cases.rds"))[[1]] %>%
  rename(cases = n) %>%
  dplyr::select(week, lad19nm, cases)

# ---------------------------------------------------------------------------- #
## Extract relevant samples, reshape and summarise

nsims <- ncol(sims)
idx <- -1:-(nrow(sims)/2)

sims <- sims[idx]
sims$week <- dat$week
sims$lad19nm <- dat$lad19nm
sims$geography <- dat$geography
sims$E_wk <- dat$E_wk_unstrat
sims_long <- reshape2::melt(sims, id.vars = nsims+1:4)

sims_long <- sims_long[, pred := exp(value)]
sims_long <- sims_long[, pred_n := exp(value)*E_wk]

saveRDS(sims_long, file = here::here("output","sims_long_avgcov.rds"))

agg_sims <- sims_long[,.(q05 = quantile(pred_n, 0.05),
                             q25 = quantile(pred_n, 0.25),
                             q50 = quantile(pred_n, 0.5),
                             q75 = quantile(pred_n, 0.75),
                             q95 = quantile(pred_n, 0.95)),
                          by = .(lad19nm, geography, week)]

write.csv(agg_sims, file = here::here("output","pred_quants_avgcov.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
## Check fit by geography and for a samples of LTLAs

plot_geog <- sims_long[,.(pred_n = sum(pred_n)), by = .(week, variable, geography)] 

dat %>%
  dplyr::group_by(week, geography) %>%
  dplyr::summarise(n = sum(n, na.rm = T),
                   pop = sum(la_pop)) -> dat_geog

plot_geog <- merge(dat_geog, plot_geog)

png(here::here("figures","deaths","pred_avgcov_geog.png"), height = 1000, width = 1200, res = 150)
print(
  data.table::setDF(plot_geog) %>%
    ggplot() + 
      geom_line(aes(week, pred_n*1e5/pop, group = variable, col = geography), alpha = 0.05) +
      geom_point(aes(week, n*1e5/pop)) +
      facet_wrap(~geography) + #
      labs(y = "Rate per 100,000", x = "Week", title = "Predicted COVID-19-related deaths at median covariate values") +
      guides(col = FALSE) +
      theme_minimal()
)
dev.off()

la_samp <- sample(unique(dat$lad19nm),9)
agg_sims$lad19nm <- dat$lad19nm
agg_sims$n <- dat$n
png(here::here("figures","deaths","pred_avgcov_lasamp.png"), height = 1000, width = 1200, res = 150)
ggplot(as.data.frame(agg_sims[lad19nm %in% la_samp]),
       aes(week, q50)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2) +
  geom_line() +
  geom_point(aes(y = n)) +
  facet_wrap(~lad19nm) +
  theme_minimal()
dev.off()

lag <- 7
agg_sims_lag <- agg_sims
agg_sims_lag$week <- agg_sims_lag$week - lag
agg_sims_lag <- merge(agg_sims_lag, cases[,c("week","lad19nm","cases")], by = c("week","lad19nm"))
pdf(here::here("figures","deaths","pred_avgcov_wcases.pdf"), height = 50, width = 50)
ggplot(as.data.frame(agg_sims_lag),
       aes(week, q50*5)) +
  geom_ribbon(aes(ymin = q05*5, ymax = q95*5, fill = geography), alpha = 0.2) +
  geom_ribbon(aes(ymin = q25*5, ymax = q75*5, fill = geography), alpha = 0.2) +
  geom_line(aes(col = geography)) +
  geom_line(aes(y = cases), lty = "dashed", col = "darkgrey", lwd = 0.5) +
  geom_point(aes(y = n*5), cex = 0.5, col= "darkgrey",pch = 20) +
  facet_wrap(~lad19nm, scales = "free") +
  labs(title = "Covariate-adjusted predicted deaths lagged by one week and rescaled by the median CFR post-P2 expansion (5 cases per death)",
       subtitle = "Confirmed cases indicated by the dashed line. Observed deaths indicated by points. Predictions are coloured by geography type.",
       y = "Rescaled weekly deaths",
       x = "") +
  theme_minimal()
dev.off()


###############################################################################
###############################################################################
