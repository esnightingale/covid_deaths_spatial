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

outdir <- "output/predict"
figdir <- "figures/predict"

dat <- readRDS(here::here("data","aggregated","deaths.rds"))
dat <- dat[["first"]] 

pred_avgcov <- readRDS(here::here(outdir,"fit_samples_avgcov.rds"))
fit_pred <- pred_avgcov$fit
samples_pred <- pred_avgcov$samples
dat_pred <- pred_avgcov$dat

sims <- as.data.table(dplyr::bind_cols(lapply(samples_pred, get_preds, dat_pred)))
saveRDS(sims, file = here::here(outdir,"sims_avgcov.rds"))

# Predictions at average covariates
# sims <- readRDS(file = here::here(outdir,"sims_avgcov.rds"))
data.table::setDT(sims)

cases <- readRDS(here::here("data","aggregated","cases.rds"))[[1]] %>%
  dplyr::rename(cases = n) %>%
  dplyr::select(week, lad19nm, cases)

# ---------------------------------------------------------------------------- #
## Extract relevant samples, reshape and summarise

nsims <- ncol(sims)
idx <- -1:-(nrow(sims)/2)

names(sims)[1:nsims] <- 1:nsims

sims <- sims[idx]
sims$week <- dat$week
sims$lad19nm <- dat$lad19nm
sims$geography <- dat$geography
sims$E_wk <- dat$E_wk_unstrat
sims$n <- dat$n
sims_long <- setDT(reshape2::melt(sims, id.vars = nsims+1:5))

sims_long <- sims_long[, pred := exp(value)]
sims_long <- sims_long[, pred_n := exp(value)*E_wk]

saveRDS(sims_long, file = here::here(outdir,"sims_long_avgcov.rds"))

agg_sims <- sims_long[,.(q01 = quantile(pred_n, 0.01),
                             q25 = quantile(pred_n, 0.25),
                             q50 = quantile(pred_n, 0.5),
                             q75 = quantile(pred_n, 0.75),
                             q99 = quantile(pred_n, 0.99)),
                          by = .(lad19nm, geography, week)]

write.csv(agg_sims, file = here::here(outdir,"pred_quants_avgcov.csv"), row.names = F)


################################################################################
################################################################################

## Check fit by geography and for a samples of LTLAs

plot_geog <- sims_long[,.(pred_n = sum(pred_n)), by = .(week, variable, geography)] 

dat %>%
  dplyr::group_by(week, geography) %>%
  dplyr::summarise(n = sum(n, na.rm = T),
                   pop = sum(la_pop)) -> dat_geog

plot_geog <- merge(dat_geog, plot_geog)

png(here::here(figdir,"pred_avgcov_geog.png"), height = 1000, width = 1200, res = 150)
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
png(here::here(figdir,"pred_avgcov_lasamp.png"), height = 1000, width = 1200, res = 150)
ggplot(as.data.frame(agg_sims[lad19nm %in% la_samp]),
       aes(week, q50)) +
  geom_ribbon(aes(ymin = q01, ymax = q99), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "steelblue") +
  geom_line(col = "steelblue") +
  geom_point(aes(y = n)) +
  facet_wrap(~lad19nm, scales = "free") +
  theme_minimal()
dev.off()


###############################################################################
###############################################################################
