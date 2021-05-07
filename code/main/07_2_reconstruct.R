################################################################################
# Description: Reconstruct case counts from model-predicted deaths
# 
# Author: Emily S Nightingale
# Date created: 18/03/2021
# 
################################################################################
################################################################################

## Shapefiles
regions <- readRDS(here::here("data","LA_shp_wpops.rds")) %>%
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

# plot_reconst(reconstruct7$total, 7, save = T, title = F, caption = F)
# plot_reconst(reconstruct7$geog, 7, save = T, title = F, caption = F)

###############################################################################
###############################################################################