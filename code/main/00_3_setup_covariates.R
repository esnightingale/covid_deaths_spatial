################################################################################
# DATA SETUP - Covariates
# 
# Read covariate datasets, clean and match on LAD19CD.
# 
################################################################################
################################################################################

# ---------------------------------------------------------------------------- #
# Index of multiple deprivation

imd <- read.csv(here::here("data","raw","IMD.csv"), header = T) %>%
  dplyr::mutate(lad19nm = as.character(LADnm))
imd$lad19nm[imd$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"

imd <- imd %>% 
  dplyr::group_by(lad19nm) %>% 
  dplyr::summarise(IMD = median(IMDScore))

# ---------------------------------------------------------------------------- #
# Ethnicity 

ethn <- readxl::read_xlsx(here::here("data","raw","ethnicity.xlsx"), skip = 9) %>%
  dplyr::separate(Area, into = c("code","lad19nm"), sep = ":") %>% 
  dplyr::filter(code == "ualad09") %>%
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
  dplyr::group_by(lad19nm) %>%
  dplyr::summarise_all(sum) 
ethn$tot_pop <- rowSums(ethn[,-1])

ethn <- 
  ethn %>%
  dplyr::mutate(prop_white = `White: Total`/tot_pop,
         prop_minority = 1-prop_white,
         prop_black = `Black/African/Caribbean/Black British: Total`/tot_pop,
         prop_asian = `Asian/Asian British: Total`/tot_pop,
         prop_mixed = `Mixed/multiple ethnic group: Total`/tot_pop,
         prop_oth = `Other ethnic group: Total`/tot_pop) %>%
  dplyr::select(lad19nm, prop_white:prop_oth) 


# ---------------------------------------------------------------------------- #
# Merge all

covs <- kw %>%
  inner_join(imd) %>%
  inner_join(ethn) 


saveRDS(covs, here::here("data","covs.rds"))

################################################################################
################################################################################
