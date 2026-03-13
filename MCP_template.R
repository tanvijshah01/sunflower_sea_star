# get the data from iNat for coyotes and roadrunners (MCP only)
# work through the KDE first. This is shorten code without as many plots.
# Jerde
# 28 January 2026
# #########################################################

#clear the environment
rm(list = ls()) #best practice when working with many scripts

#packages
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

# look up taxon_id on the iNaturalist website
# We use the Greater Roadrunner (beep, beep) for this example
# taxon id is the number on the URL for the speices' web page
coyote_id<-42051 #https://www.inaturalist.org/taxa/42051-Canis-latrans
roadrunner_id<-1986 #https://www.inaturalist.org/taxa/1986-Geococcyx-californianus



#define the area you are interested in studying with a bounding box
bbox <- list(
  swlat = 34,
  swlng = -116,
  nelat = 38.0,
  nelng = -108.0
)

# get the coyote data
obs_coyote_2025 <- get_inat_obs(
  taxon_id = coyote_id,
  quality  = "research",    #research grade observations
  year=2025,
  bounds = bbox,
  maxresults = 10000)

# get the roadrunner data
obs_roadrunner_2025 <- get_inat_obs(
  taxon_id = roadrunner_id,
  quality  = "research",  
  year=2025,
  bounds = bbox,
  maxresults = 10000)

#clean the data (note one of the packages loaded interferes with select. Must specify package)
obs_c_2025_clean<- obs_coyote_2025 |> dplyr::select(latitude, longitude) |>drop_na()
obs_rr_2025_clean <- obs_roadrunner_2025 |> dplyr::select(latitude, longitude) |>drop_na()

#convert to simple feature object.
c_2025_sf<-st_as_sf(obs_c_2025_clean, coords = c("longitude", "latitude"),crs = 4326) #Default CRS for many applications including R
rr_2025_sf<-st_as_sf(obs_rr_2025_clean, coords = c("longitude", "latitude"),crs = 4326)

#convert to spatial points
c_2025_proj<-st_transform(c_2025_sf, 5070) # CRS ideal for HR estimation in the United States.  Area estimate is accurate. 
c_2025_sp_points<-SpatialPoints(coords=st_coordinates(c_2025_proj), 
                                proj4string = CRS(st_crs(c_2025_proj)$proj4string))
rr_2025_proj<-st_transform(rr_2025_sf, 5070) 
rr_2025_sp_points<-SpatialPoints(coords=st_coordinates(rr_2025_proj), 
                                proj4string = CRS(st_crs(rr_2025_proj)$proj4string))

#calculate the MCPs
c_2025_mcp_95 <- mcp(
  c_2025_sp_points,
  percent = 95)

rr_2025_mcp_95 <- mcp(
  rr_2025_sp_points,
  percent = 95)

#convert spatial points to simple feature for plotting
c_2025_mcp_sf <- st_as_sf(c_2025_mcp_95) |>
  st_transform(4326)
rr_2025_mcp_sf <- st_as_sf(rr_2025_mcp_95) |>
  st_transform(4326)

# Get state boundaries
states <- ne_states(country = "United States of America", returnclass = "sf")

#overlapping figure
######################################################
p_combined_mcp <- ggplot() +
  geom_sf(data = states, fill = "white", color = "gray10", size=1) +
  geom_sf(data = c_2025_sf, color="brown", alpha = 0.6,size = 1)  +
  geom_sf(data = c_2025_mcp_sf, color="brown",size = 1)  +
  geom_sf(data = rr_2025_sf, color="blue", alpha = 0.6,size = 1)  +
  geom_sf(data = rr_2025_mcp_sf, color="blue",size = 1)  +
  coord_sf(xlim = c(bbox$swlng-1, bbox$nelng+1),
           ylim = c(bbox$swlat-1, bbox$nelat+1)) +
  labs(title = "Coyote (brown) & Roadrunner (blue) MCP",
       subtitle = "Year:2025") +
  theme_void(base_size = 14)
ggsave(here("figures","coyote_roadrunner_MCP.jpg"), p_combined_mcp, width=6, height=4, unit="in",dpi=300)


