################################################################################
# Description: Reconstruct case counts from model-predicted deaths
# 
# Author: Emily S Nightingale
# Date created: 18/03/2021
# 
################################################################################
################################################################################

outdir <- "output/reconstruct"
figdir <- "figures/reconstruct"

# Specify quantiles for CFR scaling and plotting
scale_quants <- c(0.25, 0.75)
plot_quants <- c(0.01,0.25,0.75,0.99)

# Load posterior predictions from death model, given average values of 
# case-fatality risk factors
# sims_long <- as.data.frame(readRDS(here::here("output","predict","sims_long_avgcov.rds")))
# nsims <- dplyr::n_distinct(sims_long$variable)

# Load observed confirmed cases
cases <- readRDS(here::here("data","aggregated","cases.rds"))[[1]]

# Sample of LAs to check reconstruction
la_samp <- sample(unique(cases$lad19nm),9)
la_samp
# [1] "Knowsley"           "Winchester"         "Copeland"           "Southwark"          "Lewes"              "Bromsgrove"         "Forest of Dean"    
# [8] "Broxbourne"         "North Warwickshire"

# ---------------------------------------------------------------------------- #
## Reconstruct cases from predicted deaths ##

# Four steps performed by helper functions in /utils:
# + [get_cfr] Calculate post-P2 CFR between cases and predicted deaths per week/LTLA given 
# assumed lag
# + [rescale_sims] Scale each posterior sample of death time series by estimates post-pillar 2 
# CFR (according to specified scale_quants)
# + [agg_sims] Aggregate posteriors over country/geography/LTLA
# + [plot_reconstr] Plot quantiles alongside observed cases

for (lag in c(7,14,21)){
  result <- reconstruct(cases, lag = lag, plot = T, suffix = "")
  saveRDS(result, here::here(outdir,paste0("reconstruct_lag",lag,".rds")))
  rm(result) 
}

# overall CFR quantiles for 1 week lag:
# period 1  4.11  2.62  6.75
# period 2  5.37  3.07  9.09


# ---------------------------------------------------------------------------- #
## Plot for lag 7 ##

reconstruct <- readRDS(here::here(outdir, "reconstruct_lag7.rds"))
plot_reconst(reconstruct$la, 7, sample = la_samp, save = T, figdir = figdir, h = 10, w = 12)
# plot_reconst(reconstruct$total, 7, save = T, figdir = figdir, title = F, caption = F)
# plot_reconst(reconstruct$geog, 7, save = T, figdir = figdir, title = F, caption = F)

###############################################################################
###############################################################################