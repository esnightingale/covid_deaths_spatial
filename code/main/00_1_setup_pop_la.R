################################################################################
# Set up population data and add to shapefile
################################################################################
################################################################################

# Function for LTLA pops: join City of London into Westminster, remove note rows and sum population > 65yrs
join_city <- function(d){
  
  rmv_rows <- min(which(rowSums(!is.na(d)) == 0)):nrow(d) # rows with exactly 0 non-missing entries
  d <- d[-rmv_rows,]  
  
  d$Name[d$Name == "City of London"] <- "Westminster"
  d$Code[d$Name == "Westminster"] <- "E09000033"
  
  d <- d %>%
    dplyr::group_by(Code, Name, Geography1) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::ungroup()
  
  d$gt65 <- rowSums(dplyr::select(d,`65`:`90`))
  
  return(d)
}

# ---------------------------------------------------------------------------- #
# Load shapefiles

regions_raw <- rgdal::readOGR(dsn = here::here("data","raw","Local_Authority_Districts_(April_2019)_Boundaries_UK_BFC-shp"),
                              layer = "Local_Authority_Districts_(April_2019)_Boundaries_UK_BFC-shp")

# Region lookup

nhsregion <- read_csv(here::here("data","raw","Local_Authority_District_to_Region__December_2019__Lookup_in_England.csv")) %>%
  dplyr::select(-FID)

# ---------------------------------------------------------------------------- #
# Raw population data

pops_raw <- readxl::read_xls(here::here("data","raw","ukmidyearestimates20192020ladcodes.xls"), 
                             sheet = "MYE2 - Persons",
                             skip = 4) 

male <- readxl::read_xls(here::here("data","raw","ukmidyearestimates20192020ladcodes.xls"), 
                         sheet = "MYE2 - Males",
                         skip = 4)  %>% 
  join_city()

female <- readxl::read_xls(here::here("data","raw","ukmidyearestimates20192020ladcodes.xls"), 
                           sheet = "MYE2 - Females",
                           skip = 4) %>% 
  join_city()

# ---------------------------------------------------------------------------- #
# Tidy data

rmv_rows <- min(which(rowSums(!is.na(pops_raw)) == 0)):nrow(pops_raw) # Note rows with exactly 0 non-missing entries
pops <- pops_raw[-rmv_rows,]  

# Calculate average population age
agevars <- 6:ncol(pops)
pops$mean_age <- rowSums(pops[,agevars]*as.integer(names(pops)[agevars]))/pops$`All ages`

# Join City of London into Westminster
pops$Name[pops$Name == "City of London"] <- "Westminster"
pops$Code[pops$Name == "Westminster"] <- "E09000033"

pops.long <- pops %>%
  dplyr::rename(geography = Geography1,
         lad19cd = Code,
         lad19nm = Name,
         la_pop = `All ages`) %>%
  dplyr::select(lad19cd:la_pop, mean_age, everything()) %>%
  tidyr::pivot_longer(cols = -lad19cd:-mean_age, names_to = "age", values_to = "pop") %>% # View()
  dplyr::mutate(age = as.integer(age)) %>%
  dplyr::group_by(lad19cd,lad19nm,geography, age) %>%
  dplyr::summarise(la_pop = sum(unique(la_pop)),
            pop = sum(pop),
            mean_age = mean(mean_age)) %>%
  dplyr::group_by(lad19cd,lad19nm,la_pop, mean_age) %>% 
  dplyr::mutate(age_group = cut(age, breaks = c(seq(0,90,5),110), right = F, include.lowest = TRUE, ordered_result = TRUE),
                age_group10 = cut(age, breaks = c(seq(0,90,10),110), right = F, include.lowest = TRUE, ordered_result = TRUE),
                cum.prop = cumsum(pop)/la_pop,
                cum.prop = case_when(cum.prop <= 0.5 ~ cum.prop,
                                     cum.prop > 0.5 ~ 0),
                q50 = which.max(cum.prop),
                med_age = age[q50],
                med_agegrp = age_group[q50]) %>% 
  dplyr::group_by(lad19cd,lad19nm, geography, la_pop, age_group, age_group10) %>%
  dplyr::summarise(la_age_pop = sum(pop),
                   mean_age = unique(mean_age),
                   med_age = unique(med_age)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate_at(vars("age_group","age_group10"),as.character) 


pops.long$age_group[pops.long$age_group == "[90,110]"] <- "[90,NA)"
pops.long$age_group10[pops.long$age_group10 == "[90,110]"] <- "[90,NA)"

pops.long <-
  pops.long %>%
  mutate_at(vars("age_group","age_group10"),as.factor) 

pops.wide <- pops.long %>%
  dplyr::select(-age_group10) %>% #-age_quintile
  tidyr::pivot_wider(id_cols = c(-age_group,-la_age_pop), names_from = age_group, values_from = la_age_pop)

# ---------------------------------------------------------------------------- #
# Add tidied population data to shapefile

regions <- st_as_sf(regions_raw) %>%
  dplyr::full_join(nhsregion) %>%
  dplyr::rename_all(tolower) %>%
  dplyr::filter(grepl("E",lad19cd)) %>%  # remove NI, Wales and Scotland
  dplyr::mutate_if(is.factor, as.character) # make factors into character to edit levels for aggregation

# Aggregate Aylsbury Vale, Chiltern, South Bucks and Wycombe into Buckinghamshire unitary authority (created April 2020)
regions$lad19nm[regions$lad19nm %in% c("Aylesbury Vale","Chiltern", "South Bucks", "Wycombe")] <- "Buckinghamshire"
regions$lad19cd[regions$lad19nm == "Buckinghamshire"] <- "E06000060"

# Join City of London into Westminster
regions$lad19nm[regions$lad19nm == "City of London"] <- "Westminster"
regions$lad19cd[regions$lad19nm == "Westminster"] <- "E09000033"

# Join Isles of Scilly into Cornwall
regions$lad19nm[regions$lad19nm == "Isles of Scilly"] <- "Cornwall"
regions$lad19cd[regions$lad19nm == "Cornwall"] <- "E06000052"

regions <- regions %>%
  dplyr::group_by(rgn19cd, rgn19nm, lad19cd, lad19nm) %>%
  dplyr::summarise() %>%
  dplyr::ungroup()

# Merge population data with shapefile and add projection
regions <- dplyr::left_join(regions, pops.wide) %>%
  sf::st_transform("+init=epsg:27700 +units=km +no_defs")


# ---------------------------------------------------------------------------- #
# Calculate sex ratio  

prop_male <- male %>%
  dplyr::rename(lad19cd = Code,
         lad19nm = Name,
         prop_male_all = `All ages`,
         prop_male_gt65 = gt65)
prop_male[,-1:-3] <- prop_male[,-1:-3]/(prop_male[,-1:-3]+female[,-1:-3])

## Add to regions df:
regions <- 
  regions %>%
  dplyr::left_join(dplyr::select(prop_male, lad19cd, lad19nm, prop_male_all, prop_male_gt65))
regions$area_km2 <- sf::st_area(regions)


# ---------------------------------------------------------------------------- #
# Save outputs

saveRDS(pops.long, here::here("data","pops.long.rds"))

saveRDS(prop_male,  here::here("data","prop_male.rds"))

saveRDS(regions,  here::here("data","LA_shp_wpops.rds"))

################################################################################
################################################################################