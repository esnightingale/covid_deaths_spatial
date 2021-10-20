################################################################################
# DATA SETUP - CASES
# 
# Call cases linelist from data pipeline, clean and aggregate to LTLA/week.
# 
################################################################################

## Load/tidy linelist

path_to_factory <- "C:/Users/emily/Documents/COVID-19/covid19_automation"
file_path <- file.path(path_to_factory, "data", "rds", "linelist_latest.rds")
key <- cyphr::data_key(file.path(path_to_factory, "data"))
linelist_raw <- cyphr::decrypt(readRDS(file_path), key)

# Filter to deaths in England with non-missing LTLA
linelist <- filter(linelist_raw, grepl("e",ltla_code)) %>% # filter to England. 
  dplyr::mutate(lad19cd = toupper(ltla_code),
                ltla_name = gsub(" And "," and ",
                             gsub(" Of "," of ",
                               stringr::str_to_title(gsub("_"," ",ltla_name)))),
                lad19nm = ltla_name)


# Aggregate Aylsbury Vale, Chiltern, South Bucks and Wycombe into 
# Buckinghamshire unitary authority (created April 2020)
linelist$lad19nm[linelist$lad19nm %in% 
                   c("Aylesbury Vale","Chiltern", 
                     "South Bucks", "Wycombe")] <- "Buckinghamshire"
linelist$lad19cd[linelist$ltla_name %in% 
                   c("Aylesbury Vale","Chiltern", 
                     "South Bucks", "Wycombe")] <- "E06000060"

# Join city of London into westminster
linelist$lad19cd[linelist$lad19cd == "E09000001"] <- "E09000033"
linelist$lad19nm[linelist$lad19nm == "City of London"] <- "Westminster"

# Exclude Isles of Scilly
linelist <- dplyr::filter(linelist, lad19nm != "Isles of Scilly")

# Tidy variables and define age group
linelist %>%
  dplyr::rename(date_report = date_lab_report,
                date = date_specimen) %>%
  dplyr::mutate_at(vars(date_report, date), lubridate::ymd) %>%
  dplyr::filter(!is.na(age), !is.na(date)) %>% 
  dplyr::mutate(result_delay = as.integer(difftime(date_report, 
                                                   date, 
                                                   units = "days")),
                week = lubridate::floor_date(date, 
                                             unit = "week", 
                                             week_start = 3),
                age_group = as.character(cut(age, 
                                             breaks = c(0,seq(10,90,10),120), 
                                             right = FALSE))) %>%
  dplyr::group_by(lad19cd) %>%
  dplyr::mutate(ltla_first = min(date, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(age) %>%
  dplyr::select(lad19cd, lad19nm, week, date, date_report, ltla_first, 
                result_delay, age_group, pillar) -> alldata

alldata$age_group[alldata$age_group == "[90,120)"] <- "[90,NA)"
alldata$age_group <- as.factor(alldata$age_group)
# summary(alldata$age_group)

################################################################################

saveRDS(alldata, here::here("data","linelist","linelist_cases.rds"))

################################################################################
################################################################################


