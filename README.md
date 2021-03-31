# Reconstructing the symptomatic-case epidemic from COVID-19-related deaths in England

This repository contains the code used to produce estimated trajectories of detectable, symptomatic cases of COVID-19 during the first epidemic wave in England, based on reported COVID-19-related deaths per week and per lower tier local authority (LTLA). This is acheived by estimating a confirmed-case-fatality ratio according to cases detected under the expanded pillar 1 + pillar 2 testing system in place from 18th May 2020, and applying this to rescale smoothed trajectories of deaths from January to June, accounting for case-fatality risk factors of the local population. 

The analysis consists of the following steps:

0. **Import/clean data**
   Import and process population estimates (by age) and covariate values (deprivation index and proportion of ethnic minority population) for each LTLA in England.
   Import and clean linelist of confirmed cases (incl. pillars 1&2) and COVID-19-related deaths (incl. as underlying cause or mentioned on death certificate).
   *Scripts: 00_pop_by_la.R, 00_setup_covariates.R, 00_setup_cases.R, 00_setup_deaths.R*
   
1. **Analysis data setup**
   Aggregate death linelist to week and LTLA, join with covariates/population estimates and calculate age-stratified expected deaths.
   Aggregate case linelist to week and LTLA.
   *Script: 01_setup_analysis_data_all.R*
   
2. **Descriptive** Generate descriptive summaries and figures
   *Script: 02_descriptive.R*
   
3. **Model fitting** Specify priors and model formulae, then fit and generate posterior samples.
   *Script: 03_run_models.R*
   
4. **Model comparison** Summarise and compare fit of candidate models
   *Script: 04_output_figs_tabs.R*
   
5. **Predict** Generate predictions from selected model at average covariate values, by posterior sampling.
   *Script: 05_predict.R*
   
6. **Reconstruct cases** Lag and rescale samples to represent cases detectable under the expanded testing system, throughout the whole time period.
   *Scripts: 06_1_reshape_sims.R, 06_2_reconstruct.R*

Helper functions are included in /utils. 
