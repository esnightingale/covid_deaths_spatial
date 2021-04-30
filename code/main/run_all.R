################################################################################
# Description: Run whole analysis pipeline
# 
# Author: Emily S Nightingale
# Date created: 30/04/2021
# 
################################################################################
################################################################################

source(here::here("code","setup_env.R"))
list.files(here::here("code","utils"), full.names = TRUE) %>% purrr::walk(source)

# Set directory for figures and outputs
figdir <- "figures"
outdir <- "output"

## SETUP ##
# Tidy LTLA population data and add to shapefile 
source(here::here("code","main","00_1_setup_pop_la.R"))
# Tidy and merge covariate data for each LTLA
source(here::here("code","main","00_2_setup_covariates.R"))
# Tidy case and death linelists, aggregate and add population/covariate data
# Setup variables for analysis
source(here::here("code","main","01_setup_analysis_data.R"))

## DESCRIPTIVE ##
# Summaries and figures
source(here::here("code","main","02_descriptive.R"))

## MODEL FITTING ##
# Setup, fit and compare models
source(here::here("code","main","03_run_models.R"))
source(here::here("code","main","04_compare_models.R"))
# Summarise fit for final model
source(here::here("code","main","05_plot_final_model.R"))

## PREDICTION ##
# Predict from final model at averaged population case-fatality risk factors
source(here::here("code","main","06_predict.R"))
# Reshape dataset of posterior predictions using data.table
source(here::here("code","main","07_1_reshape_sims.R"))
# Calculate CFR and rescale posterior predictions of deaths to represent cases
source(here::here("code","main","07_2_reconstruct.R"))

## SUMMARISE ##
# Generate figures and tables to summarise final case predictions 
source(here::here("code","main","08_summarise.R"))

################################################################################
################################################################################

