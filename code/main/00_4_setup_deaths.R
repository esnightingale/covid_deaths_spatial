################################################################################
# DATA SETUP - PIPELINE
# 
# Call deaths linelist from data pipeline and clean - defining 10-yr age groups
# and LTLA names to match with shapefiles and pops.
# 
################################################################################

## Load/tidy linelist

# READ ONS LINELIST
path_to_factory <- "C:/Users/emily/Documents/COVID-19/covid19_automation"
file_path <- file.path(path_to_factory, "data", "rds", "deaths_eng_latest.rds")
key <- cyphr::data_key(file.path(path_to_factory, "data"))
covid_deaths_raw <- cyphr::decrypt(readRDS(file_path), key)

# Filter to deaths in England with non-missing LTLA
# death_type28 == 1 indicates official covid-related death type
covid_deaths <- filter(covid_deaths_raw, grepl("e",ltla_code)) %>% 
  mutate(lad19cd = toupper(ltla_code),
         ltla_name = gsub(" And "," and ",
                        gsub(" Of "," of ",
                        stringr::str_to_title(gsub("_"," ",ltla_name)))),
         lad19nm = ltla_name) 

covid_deaths$date_swab[lubridate::ymd(covid_deaths$date_swab) < "2020-01-01"] <- NA

# Aggregate Aylsbury Vale, Chiltern, South Bucks and Wycombe into 
# Buckinghamshire unitary authority (created April 2020)
covid_deaths$lad19nm[covid_deaths$lad19nm %in% 
                       c("Aylesbury Vale","Chiltern", 
                         "South Bucks", "Wycombe")] <- "Buckinghamshire"
covid_deaths$lad19cd[covid_deaths$ltla_name %in% 
                       c("Aylesbury Vale","Chiltern", 
                         "South Bucks", "Wycombe")] <- "E06000060"

# Join city of London into westminster
covid_deaths$lad19cd[covid_deaths$lad19cd == "E09000001"] <- "E09000033"
covid_deaths$lad19nm[covid_deaths$lad19cd == "E09000033"] <- "Westminster"

# Tidy variables and define age group
covid_deaths %>%
  dplyr::rename(date = date_death,
                date_report = ons_date_report) %>%
  dplyr::mutate_at(vars(date_swab, date, date_report), lubridate::ymd) %>%
  dplyr::filter(!is.na(age) & !is.na(date)) %>% 
  dplyr::mutate(swab_death = as.integer(difftime(date, 
                                                 date_swab, 
                                                 units = "days")),
                report_delay = as.integer(difftime(date_report, 
                                                   date, 
                                                   units = "days")),
                death_type2 = factor(case_when(
                  death_type == "lab_confirmed" ~ "Lab confirmed",
                  grepl("not_lab",death_type) ~ "Not lab confirmed")),
                week = lubridate::floor_date(date, 
                                             unit = "week", 
                                             week_start = 3),
                age_group = as.character(cut(age, 
                                             breaks = c(0,seq(10,90,10),120), #
                                             right = FALSE))) %>%
  dplyr::group_by(lad19cd) %>%
  dplyr::mutate(ltla_first = min(date, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(age) %>%
  dplyr::select(lad19cd, ltla_name, ltla_first, week, date, date_report, pillars, 
                date_swab, swab_death, report_delay, age_group, death_type, 
                death_type2) -> alldata

alldata$age_group[alldata$age_group == "[90,120)"] <- "[90,NA)"
alldata$age_group <- as.factor(alldata$age_group)
# summary(alldata$age_group)

################################################################################

saveRDS(alldata, here::here("data","linelist","linelist_deaths.rds"))

################################################################################
################################################################################


