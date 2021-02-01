################################################################################
# Description: Back-calculate onset of symptoms for each death, from estimated
# distribution.
# 
# Author: Emily S Nightingale
# Date created: 13/01/2021
# 
################################################################################
################################################################################

################################################################################
# SETUP
################################################################################

library(tidyverse)
library(lubridate)
library(linelist)
library(INLA)
library(ggplot2)
library(patchwork)
library(viridis)
library(here)
library(data.table)
# devtools::install_github("reconhub/distcrete")
library(distcrete)
library(fitdistrplus)
library(DiscreteWeibull)
library(mixdist)

# Local data directory
datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

theme_set(theme_minimal())

wave <- 1

# source(here::here("code","main","functions.R"))
list.files(here::here("code","utils"), full.names = TRUE) %>% walk(source)

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

## Shapefiles
regions <- readRDS(paste0(datadir,"maps/LA_shp_wpops.rds")) %>%
  filter(grepl("E", lad19cd))

# LTLA-week-aggregated observed deaths
deaths <- readRDS(here::here("data","deaths.rds")) 


# Fitted models and posterior samples
fit <- readRDS(file = here::here("output",
                                    sprintf("fits_deaths_%s.rds", wave)))[["BYM_geog"]]

samples <- readRDS(file = here::here("output",
                                       sprintf("samples_deaths_%s.rds", wave)))[["BYM_geog"]]

# outputdir <- "deaths"
dat <- deaths[[wave]]

################################################################################
# Delay distribution
################################################################################

# From CO-CIN data on symptom onset to death during first wave
# https://www.gov.uk/government/publications/co-cin-covid-19-time-from-symptom-onset-until-death-in-uk-hospitalised-patients-7-october-2020

med <- 14
IQ <- 14


weib_med <- function(scale, shape){
  return(scale * log(2)^(1/shape))
}

weib_med(1,1)

med_iqr <- function(theta){
  quants <- qweibull(p = c(0.25,0.5,0.75), shape = theta[1], scale = theta[2])
  return(c(quants[2],quants[3]-quants[1]))
}

med_iqr(c(1.5,18))
# [1] 14.09796 14.53493

# Delay distribution:
sample_lag <- rweibull(1,shape = 1.5,scale = 18)











######### FIT IQR'S ########

# # calculate the difference between quantiles and weibell output for theta
min_quantiles <- function(theta, dist_x, iqr = c(0.25,0.5,0.75)){
  #calculate quantiles with test parameters
  test_quantiles <- pweibull(q = dist_x, shape = theta[1], scale = theta[2])
  actual_quantiles <- iqr
  # calculate difference between actualy quantiles and tested quantiles
  diff_quantiles <- c(test_quantiles[2] - actual_quantiles[2],
                      (test_quantiles[3]-test_quantiles[1]) - actual_quantiles[3]-actual_quantiles[1])
  # make it absolute
  optim_value <- sum(diff_quantiles^2)
  return(optim_value)
}

# # optimising wrapper over the min_quantiles function - weibull
optimise_quantiles <- function(init_values, dist_x, iqr = c(0.25,0.5,0.75)){
  #optimise the min_quantiles function.
  #suppressWarnings so don't get message about trying values that don't work
  par_out <- suppressWarnings(optim(par = c(init_values[1], init_values[2]),
                                    fn = min_quantiles,
                                    dist_x = c(dist_x[c("median", "IQR")]),
                                    iqr=iqr))

  return(par_out)
}
