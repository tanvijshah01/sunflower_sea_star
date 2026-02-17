################## ESM 211: Sunflower Seastar Homeranges #################################
### Danielle Turner
### February 2026

##########################################################################################
##########################################################################################
############################# EXPLORATORY DATA ###########################################

# get Sunflower Seastar data from iNaturalist
library(rinat)

obs_seastar <- get_inat_obs(taxon_id = 47673, 
                            quality = "research", maxresults = 10000)

# cleaning data dates
obs_seastar_clean <- obs_seastar |>
  mutate(
    datetime = as.numeric(substr(observed_on, 1, 4))
  )


# split data frame into 3 new data frames: 
# before (1973-2013), during (2013-2017), after (2017-present)
library(dplyr)
library(lubridate)

obs_seastar_before <- obs_seastar |>
  filter(year(datetime) < 2013)

obs_seastar_during <- obs_seastar |>
  filter(year(datetime) >= 2013 & year(datetime) <= 2017)

obs_seastar_after <- obs_seastar |>
  filter(year(datetime) > 2017)

