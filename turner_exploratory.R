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

# plot the data
plot(obs_seastar$longitude, obs_seastar$latitude,
     xlab = "Longitude",
     ylab = "Latitude",
     pch = 20,
     col = "steelblue",
     main = "iNaturalist Sea Star Observations")

# plot for visualizing data on map
library(ggplot2)
library(maps)

world <- map_data("world")

ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray90",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = obs_seastar,
    aes(x = longitude, y = latitude),
    size = 0.6,
    alpha = 0.6,
    color = "steelblue"
  ) +
  coord_fixed(
    xlim = c(-170, -100),   # Alaska to Mexico (longitude)
    ylim = c(10, 75)        # Mexico to Arctic Alaska (latitude)
  ) +
  theme_minimal() +
  labs(
    title = "iNaturalist Sea Star Observations",
    x = "Longitude",
    y = "Latitude"
  )

##########################################################################################
##########################################################################################
############################# DISEASE IMPACTS ###########################################

# get pre- and post-seastar wasting disease data
library(dplyr)
library(hexbin)

obs_seastar_clean <- obs_seastar |>
  mutate(
    datetime = as.numeric(substr(observed_on, 1, 4))
  )

pre_2014 <- obs_seastar_clean %>%
  filter(datetime <= 2014)

post_2014 <- obs_seastar_clean %>%
  filter(datetime > 2014)

# plot the pre-disease data
ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = pre_2014,
    aes(x = longitude, y = latitude),
    size = 0.6,
    alpha = 0.6,
    color = "darkblue"
  ) +
  coord_fixed(
    xlim = c(-170, -100),
    ylim = c(10, 75)
  ) +
  theme_minimal() +
  labs(
    title = "Sunflower Sea Star Observations (≤ 2013)",
    x = "Longitude",
    y = "Latitude"
  )

# plot the post-disease data
ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = post_2014,
    aes(x = longitude, y = latitude),
    size = 0.6,
    alpha = 0.6,
    color = "red"
  ) +
  coord_fixed(
    xlim = c(-170, -100),
    ylim = c(10, 75)
  ) +
  theme_minimal() +
  labs(
    title = "Sunflower Sea Star Observations (> 2013)",
    x = "Longitude",
    y = "Latitude"
  )

# plot side by side
obs_seastar_clean_v2 <- obs_seastar %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  mutate(
    year = as.numeric(substr(observed_on, 1, 4)),
    period = factor(
      ifelse(datetime <= 2014, "≤ 2014", "> 2014"),
      levels = c("≤ 2014", "> 2014")
    )
  )


seastar_pre_post_2014 <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = obs_seastar_clean_v2,
    aes(x = longitude, y = latitude),
    size = 0.5,
    alpha = 0.5,
    color = "steelblue"
  ) +
  facet_wrap(~ period, nrow = 1) +
  coord_fixed(
    xlim = c(-155, -115),   # tighter longitudinal coastal band
    ylim = c(28, 62)        # SoCal → Gulf of Alaska
  ) +
  theme_minimal() +
  labs(
    title = "Sunflower Sea Star Observations Before/During vs After 2014",
    x = "Longitude",
    y = "Latitude"
  )

# save plot
ggsave(
  filename = "seastar_pre_post_2014.png",
  plot = seastar_pre_post_2014,
  width = 10,
  height = 5,
  dpi = 300
)

