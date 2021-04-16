################################################################################
# DATA SETUP - Covariates
# 
# Read covariate datasets, clean and match on LAD19CD.
# 
################################################################################

library(tidyverse)
library(linelist)
library(INLA)
library(spdep)
library(sf)

datadir <- "C:/Users/phpuenig/Documents/COVID-19/Data/"

################################################################################

health <- readxl::read_xlsx(paste0(datadir,"covariates/PHE_localhealth_indicators.xlsx"), sheet = "Data", skip = 3) %>% 
  as_tibble() %>%
  rename(lad19nm = Label, lad19cd = Code)
  
health$lad19nm[health$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
health$lad19nm[health$lad19nm == "Shepway"] <- "Folkestone and Hythe"
health$lad19nm[health$lad19nm %in% c("Bournemouth","Christchurch", "Poole")] <- "Bournemouth, Christchurch and Poole"
health$lad19nm[grepl("Dorset", health$lad19nm) | health$lad19nm %in% c("Weymouth and Portland","Purbeck")] <- "Dorset"
health$lad19nm[health$lad19nm %in% c("Forest Heath","St Edmundsbury")] <- "West Suffolk"
health$lad19nm[health$lad19nm %in% c("Waveney","Suffolk Coastal")] <- "East Suffolk"
health$lad19nm[health$lad19nm %in% c("West Somerset","Taunton Deane")] <- "Somerset West and Taunton"
health$lad19nm[health$lad19nm == "City of London"] <- "Westminster"

health <- health %>% 
  group_by(lad19nm) %>% 
  summarise(across(c(`Total population`,`Older People in Deprivation, Number`), 
                   sum),
            across(starts_with(c("Emergency","Limiting","Incidence","Population",
                                 "Index","Deaths","Black","Older","Healthy","Low")), 
                               mean)) %>%
  rename(older_depr_n = `Older People in Deprivation, Number`,
         adm_chd = `Emergency hospital admissions for CHD`,
         adm_copd = `Emergency hospital admissions for Chronic Obstructive Pulmonary Disease (COPD)`,
         LTI_disability = `Limiting long-term illness or disability`, 
         lung_cancer = `Incidence of lung cancer`,
         pop_nonwhite = `Population whose ethnicity is not 'White UK'`,
         death_all = `Deaths from all causes, all ages`, 
         life_exp_m = `Healthy life expectancy for males, 2009-2013`,
         life_exp_f = `Healthy life expectancy for females, 2009-2013`)

saveRDS(health, paste0(datadir,"health_indicators.rds"))

## Health workers (census 2011)
kw <- readxl::read_xlsx(paste0(datadir,"covariates/keyworkersreferencetableupdated.xlsx"), sheet = "Table 19", range = "A5:C383") %>% #View()
  as_tibble() %>%
  rename(lad19nm = `...1`,
         pop_kw = population,
         prop_kw = percentage) %>%
  mutate(lad19nm = gsub(" UA","",lad19nm),
         prop_kw = round(prop_kw/100,2),
         n_kw = pop_kw*prop_kw) 

# Remove note rows
# kw <- kw[-which(is.na(kw$lad19nm))[1]:-nrow(kw),]

kw$lad19nm[kw$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
kw$lad19nm[kw$lad19nm == "Westminster and the City of London"] <- "Westminster"
kw$lad19nm[kw$lad19nm == "Cornwall and the Isles of Scilly"] <- "Cornwall"
kw$lad19nm[kw$lad19nm == "Shepway"] <- "Folkestone and Hythe"
kw$lad19nm[kw$lad19nm == "St. Albans"] <- "St Albans"
kw$lad19nm[kw$lad19nm %in% c("Bournemouth","Christchurch", "Poole")] <- "Bournemouth, Christchurch and Poole"
kw$lad19nm[grepl("Dorset", kw$lad19nm) | kw$lad19nm %in% c("Weymouth and Portland","Purbeck")] <- "Dorset"
kw$lad19nm[kw$lad19nm %in% c("Forest Heath","St. Edmundsbury")] <- "West Suffolk"
kw$lad19nm[kw$lad19nm %in% c("Waveney","Suffolk Coastal")] <- "East Suffolk"
kw$lad19nm[kw$lad19nm %in% c("West Somerset","Taunton Deane")] <- "Somerset West and Taunton"
kw$lad19nm[kw$lad19nm == "Rhondda, Cynon, Taf"] <- "Rhondda Cynon Taf"

kw <- 
  kw %>%
  group_by(lad19nm) %>%
  summarise(pop_kw = sum(pop_kw),
            prop_kw = mean(prop_kw)) 

# Cornwall and Westminster missing - impute Cornwall with median of neighbouring regions and leave Westminster as 0 since very small population
kw$lad19nm[kw$prop_kw == 0]
missing_kw <- which(regions$lad19nm %in% c("Cornwall","Westminster"))

# neighbour IDs
kw$prop_kw[kw$lad19nm == "Cornwall"] <- median(kw$prop_kw[kw$lad19nm %in% regions$lad19nm[unlist(g$nbs[missing_kw[1]])]])
kw$prop_kw[kw$lad19nm == "Westminster"] <- median(kw$prop_kw[kw$lad19nm %in% regions$lad19nm[unlist(g$nbs[missing_kw[2]])]])

# View(kw)


## Index of multiple deprivation
imd <- read.csv(paste0(datadir,"covariates/Indices_of_Multiple_Deprivation_(IMD)_2019.csv"), header = T) %>%
  mutate(lad19nm = as.character(LADnm))
imd$lad19nm[imd$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"

imd <- imd %>% 
  group_by(lad19nm) %>% 
  summarise(IMD = median(IMDScore))


## Ethnicity 

ethn <- readxl::read_xlsx(paste0(datadir,"covariates/ethnicity.xlsx"), skip = 9) %>%
  separate(Area, into = c("code","lad19nm"), sep = ":") %>% 
  filter(code == "ualad09") %>%
  dplyr::select(-code)

ethn$lad19nm[ethn$lad19nm == "City of London"] <- "Westminster"
ethn$lad19nm[ethn$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
ethn$lad19nm[ethn$lad19nm %in% c("Bournemouth","Christchurch", "Poole")] <- "Bournemouth, Christchurch and Poole"
ethn$lad19nm[grepl("Dorset", ethn$lad19nm) | ethn$lad19nm %in% c("Weymouth and Portland","Purbeck")] <- "Dorset"
ethn$lad19nm[ethn$lad19nm %in% c("Forest Heath","St Edmundsbury")] <- "West Suffolk"
ethn$lad19nm[ethn$lad19nm == "Rhondda Cynon Taff"] <- "Rhondda Cynon Taf"
ethn$lad19nm[ethn$lad19nm %in% c("Waveney","Suffolk Coastal")] <- "East Suffolk"
ethn$lad19nm[ethn$lad19nm %in% c("West Somerset","Taunton Deane")] <- "Somerset West and Taunton"

ethn <- 
  ethn %>%
  group_by(lad19nm) %>%
  summarise_all(sum) 
ethn$tot_pop <- rowSums(ethn[,-1])

ethn <- 
  ethn %>%
  mutate(prop_white = `White: Total`/tot_pop,
         prop_minority = 1-prop_white,
         prop_black = `Black/African/Caribbean/Black British: Total`/tot_pop,
         prop_asian = `Asian/Asian British: Total`/tot_pop,
         prop_mixed = `Mixed/multiple ethnic group: Total`/tot_pop,
         prop_oth = `Other ethnic group: Total`/tot_pop) %>%
  dplyr::select(lad19nm, prop_white:prop_oth) 


## Rural/urban proportion
# rural <- read.csv(paste0(datadir,"covariates/rural_urban.csv")) %>% View()
#   mutate(lad19nm = as.character(Name),
#          tot_pop = as.numeric(gsub("[^0-9.-]","",trimws(`Total.Population1`))),
#          rural_pop = as.numeric(gsub("[^0-9.-]","",trimws(`Total.Rural.Population..including.Large.Market.Town.population.2`)))) %>% 
#   dplyr::select(lad19nm, tot_pop, rural_pop)
# 
# rural$lad19nm[rural$lad19nm == "City of London"] <- "Westminster"
# rural$lad19nm[rural$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
# rural$lad19nm[rural$lad19nm %in% c("Bournemouth","Christchurch", "Poole")] <- "Bournemouth, Christchurch and Poole"
# rural$lad19nm[grepl("Dorset", rural$lad19nm) | rural$lad19nm %in% c("Weymouth and Portland","Purbeck")] <- "Dorset"
# rural$lad19nm[rural$lad19nm %in% c("Forest Heath","St. Edmundsbury")] <- "West Suffolk"
# rural$lad19nm[rural$lad19nm %in% c("Waveney","Suffolk Coastal")] <- "East Suffolk"
# rural$lad19nm[rural$lad19nm %in% c("West Somerset","Taunton Deane")] <- "Somerset West and Taunton"
# rural$lad19nm[rural$lad19nm == "Durham"] <- "County Durham"
# rural$lad19nm[rural$lad19nm == "Shepway"] <- "Folkestone and Hythe"
# rural$lad19nm <- gsub(" City of", ", City of", rural$lad19nm)
# rural$lad19nm <- gsub(" County of", ", County of", rural$lad19nm)
# 
# rural <- 
#   rural %>%
#   mutate(rural_pop = replace_na(rural_pop, 0)) %>%
#   group_by(lad19nm) %>%
#   summarise_all(sum) %>%
#   mutate(prop_rural = rural_pop/tot_pop)


## Land use
landuse <- read.csv(paste0(datadir,"covariates/land_use.csv")) %>% 
  mutate(lad19nm = as.character(Local.Authority.Name)) %>%
  filter(lad19nm != "") %>%
  mutate_at(vars(-ONS.Code, -lad19nm, -Local.Authority.Name), function(x) as.numeric(gsub(",","",gsub("-","",x)))) %>% 
  dplyr::select(lad19nm, Total, Developed.Use, Non.developed.Use, Total.residential)

landuse$lad19nm[landuse$lad19nm == "City of London"] <- "Westminster"
landuse$lad19nm[landuse$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
landuse$lad19nm[landuse$lad19nm %in% c("Bournemouth","Christchurch", "Poole")] <- "Bournemouth, Christchurch and Poole"
landuse$lad19nm[grepl("Dorset", landuse$lad19nm) | landuse$lad19nm %in% c("Weymouth and Portland","Purbeck")] <- "Dorset"
landuse$lad19nm[landuse$lad19nm %in% c("Forest Heath","St. Edmundsbury")] <- "West Suffolk"
landuse$lad19nm[landuse$lad19nm %in% c("Waveney","Suffolk Coastal")] <- "East Suffolk"
landuse$lad19nm[landuse$lad19nm %in% c("West Somerset","Taunton Deane")] <- "Somerset West and Taunton"
landuse$lad19nm[landuse$lad19nm == "Durham"] <- "County Durham"
landuse$lad19nm[landuse$lad19nm == "Shepway"] <- "Folkestone and Hythe"
landuse <- 
  landuse %>%
  group_by(lad19nm) %>%
  summarise_all(sum) %>%
  mutate(prop_dev_land = Developed.Use/Total,
         total_area_km2 = Total/100)

################################################################################

# Merge all

covs <- kw %>%
  inner_join(imd) %>%
  inner_join(ethn) %>%
  inner_join(landuse) 


saveRDS(covs, paste0(datadir,"covs.rds"))
