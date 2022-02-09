################################################################################
# DATA SETUP - ALL
# 
# Take cleaned linelists of deaths and cases, calculate E over specified time
# period, aggregate by week and LTLA and merge all
# 
################################################################################

## LTLA covariates
covs <- readRDS(here::here("data","covs.rds"))

## Shapefiles
regions <- readRDS(here::here("data","LA_shp_wpops.rds")) %>%
  dplyr::filter(grepl("E", lad19cd))
regions.df <- sf::st_drop_geometry(regions)

# Define time period by week start
waves <- list(first = c(lubridate::ymd("2020-01-01"), lubridate::ymd("2020-06-30")),
              second = c(lubridate::ymd("2020-07-01"), lubridate::ymd("2020-12-02")))

## CASES ##
source(here::here("code","main","00_3_setup_cases.R"))

cases <- lapply(waves,
                function(wave) setup_analysis_data(alldata, 
                                                   measure = "cases",
                                                   start = wave[1], 
                                                   end = wave[2]))

cases <- rlist::list.append(cases, breaks = waves)

saveRDS(cases, here::here("data","aggregated","cases.rds"))

## DEATHS ##
source(here::here("code","main","00_4_setup_deaths.R"))

deaths <- lapply(waves,
                 function(wave) setup_analysis_data(alldata, 
                                                    measure = "deaths",
                                                    start = wave[1], 
                                                    end = wave[2]))

deaths <- rlist::list.append(deaths, breaks = waves)

saveRDS(deaths, here::here("data","aggregated","deaths.rds"))

# range(alldata_sub$date) 
# "2020-01-05" "2020-06-30"

# summary(as.factor(alldata_sub$death_type2))
# Lab confirmed Not lab confirmed
#         39332             13228

################################################################################
################################################################################
