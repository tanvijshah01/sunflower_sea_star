################## ESM 211: Sunflower Seastar Homeranges #################################
### Danielle Turner
### February 2026

##########################################################################################
##########################################################################################
############################# SEPARATE DATA BY YEARS ####################################

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
  filter(year(datetime) >= 2013 & year(datetime) <= 2017)|>
  filter(!is.na(latitude) & !is.na(longitude))

obs_seastar_after <- obs_seastar |>
  filter(year(datetime) > 2017)|>
  filter(!is.na(latitude) & !is.na(longitude))


################################################################################
######################### BEFORE SWD HOMERANGE #################################
# Packages
library(here) #localized file paths
library(tidyverse) #compendium of useful packages
library(janitor) #cleans data
library(rinat) #data interface for iNaturalist
library(sf) # for converting and using simple features in maps
library(rnaturalearth) #map data and functions
library(rnaturalearthdata) #map data
library(patchwork) # for stacking figures
library(sp)#For dealing with spatial data
library(adehabitatHR)# for home range estimation

# Convert to simple feature object.
seastar_before_sf<-st_as_sf(obs_seastar_before, coords = c("longitude", "latitude"),crs = 4326) #Default CRS for many applications including R

# Convert to spatial points
seastar_before_proj<-st_transform(seastar_before_sf, 5070) # CRS ideal for HR estimation in the United States.  Area estimate is accurate. 
seastar_before_sp_points<-SpatialPoints(coords=st_coordinates(seastar_before_proj), 
                                proj4string = CRS(st_crs(seastar_before_proj)$proj4string))

# Calculate the MCPs
seastar_before_mcp_90 <- mcp(
  seastar_before_sp_points,
  percent = 90)

# Convert spatial points to simple feature for plotting
seastar_before_mcp_sf <- st_as_sf(seastar_before_mcp_95) |>
  st_transform(4326)

# Get state boundaries
states <- ne_states(country = "United States of America", returnclass = "sf")

# Define boundaries of interest
bbox <- list(
  swlat = 32.5,     # Southern California
  swlng = -180,     # Far western Aleutians
  nelat = 71.5,     # Northern Alaska coast
  nelng = -117      # Inland edge of California coast
)

# Create figure
before_combined_mcp <- ggplot() +
  geom_sf(data = states, fill = "white", color = "gray10", size=0.5) +
  geom_sf(data = seastar_before_sf, color="steelblue", alpha = 0.6,size = 0.5)  +
  geom_sf(data = seastar_before_mcp_sf, color="steelblue",size = 1)  +
  coord_sf(xlim = c(bbox$swlng-1, bbox$nelng+1),
           ylim = c(bbox$swlat-1, bbox$nelat+1)) +
  labs(title = "Sunflower Seastar MCP",
       subtitle = "Before SWD: 1978-2013") +
  theme_void(base_size = 14)

before_combined_mcp


################################################################################
######################### DURING SWD HOMERANGE #################################

# Convert to simple feature object.
seastar_during_sf<-st_as_sf(obs_seastar_during, coords = c("longitude", "latitude"),crs = 4326)

# Convert to spatial points
seastar_during_proj<-st_transform(seastar_during_sf, 5070) # CRS ideal for HR estimation in the United States.  Area estimate is accurate. 
seastar_during_sp_points<-SpatialPoints(coords=st_coordinates(seastar_during_proj), 
                                        proj4string = CRS(st_crs(seastar_during_proj)$proj4string))

# Calculate the MCPs
seastar_during_mcp_90 <- mcp(
  seastar_during_sp_points,
  percent = 90)

# Convert spatial points to simple feature for plotting
seastar_during_mcp_sf <- st_as_sf(seastar_during_mcp_95) |>
  st_transform(4326)

# # Get state boundaries
# states <- ne_states(country = "United States of America", returnclass = "sf")
# 
# # Define boundaries of interest
# bbox <- list(
#   swlat = 32.5,     # Southern California
#   swlng = -180,     # Far western Aleutians
#   nelat = 71.5,     # Northern Alaska coast
#   nelng = -117      # Inland edge of California coast
# )

# Create figure
during_combined_mcp <- ggplot() +
  geom_sf(data = states, fill = "white", color = "gray10", size=0.5) +
  geom_sf(data = seastar_during_sf, color="darkorange", alpha = 0.6,size = 0.5)  +
  geom_sf(data = seastar_during_mcp_sf, color="darkorange",size = 1)  +
  coord_sf(xlim = c(bbox$swlng-1, bbox$nelng+1),
           ylim = c(bbox$swlat-1, bbox$nelat+1)) +
  labs(title = "Sunflower Seastar MCP",
       subtitle = "During SWD: 2014-2017") +
  theme_void(base_size = 14)

during_combined_mcp

################################################################################
######################### AFTER SWD HOMERANGE #################################

# Convert to simple feature object.
seastar_after_sf<-st_as_sf(obs_seastar_after, coords = c("longitude", "latitude"),crs = 4326)

# Convert to spatial points
seastar_after_proj<-st_transform(seastar_after_sf, 5070) # CRS ideal for HR estimation in the United States.  Area estimate is accurate. 
seastar_after_sp_points<-SpatialPoints(coords=st_coordinates(seastar_after_proj), 
                                        proj4string = CRS(st_crs(seastar_after_proj)$proj4string))

# Calculate the MCPs
seastar_after_mcp_90 <- mcp(
  seastar_after_sp_points,
  percent = 90)

# Convert spatial points to simple feature for plotting
seastar_after_mcp_sf <- st_as_sf(seastar_after_mcp_95) |>
  st_transform(4326)

# # Get state boundaries
# states <- ne_states(country = "United States of America", returnclass = "sf")
# 
# # Define boundaries of interest
# bbox <- list(
#   swlat = 32.5,     # Southern California
#   swlng = -180,     # Far western Aleutians
#   nelat = 71.5,     # Northern Alaska coast
#   nelng = -117      # Inland edge of California coast
# )

# Create figure
after_combined_mcp <- ggplot() +
  geom_sf(data = states, fill = "white", color = "gray10", size=0.5) +
  geom_sf(data = seastar_after_sf, color="darkblue", alpha = 0.6,size = 0.5)  +
  geom_sf(data = seastar_after_mcp_sf, color="darkblue",size = 1)  +
  coord_sf(xlim = c(bbox$swlng-1, bbox$nelng+1),
           ylim = c(bbox$swlat-1, bbox$nelat+1)) +
  labs(title = "Sunflower Seastar MCP",
       subtitle = "After SWD: 2017-2026") +
  theme_void(base_size = 14)

after_combined_mcp
