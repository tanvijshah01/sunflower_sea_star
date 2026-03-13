################## ESM 211: Sunflower Seastar Homeranges #################################
### Danielle Turner
### Updated to match current plotting theme and date windows
##########################################################################################

# -----------------------------
# Packages
# -----------------------------
library(here)
library(tidyverse)
library(janitor)
library(rinat)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(sp)
library(adehabitatHR)
library(lubridate)
library(showtext)
library(sysfonts)

# -----------------------------
# Theme setup
# -----------------------------
font_add_google("EB Garamond", "eb_garamond")

showtext_auto(TRUE)
showtext_opts(dpi = 300)

col_title <- "#003660"
col_text  <- "black"

col_panel_bg <- "#FFFFFF"
col_plot_bg  <- "#FFFFFF"

col_axis     <- "black"
col_grid_maj <- "#D9D9D9"

theme_bolb <-
  theme_minimal(base_size = 12, base_family = "eb_garamond") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      family = "eb_garamond",
      face = "bold",
      size = 14,
      color = col_title,
      hjust = 0,
      margin = margin(b = 8)
    ),
    plot.subtitle.position = "plot",
    plot.subtitle = element_text(
      family = "eb_garamond",
      face = "plain",
      size = 12,
      color = col_text,
      hjust = 0,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(
      family = "eb_garamond",
      face = "bold",
      size = 12,
      color = col_text,
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      family = "eb_garamond",
      face = "bold",
      size = 12,
      color = col_text,
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      family = "eb_garamond",
      size = 12,
      color = col_text,
      margin = margin(t = 4)
    ),
    axis.text.y = element_text(
      family = "eb_garamond",
      size = 12,
      color = col_text,
      margin = margin(r = 4)
    ),
    legend.title = element_text(
      family = "eb_garamond",
      face = "plain",
      size = 12,
      color = col_text
    ),
    legend.text = element_text(
      family = "eb_garamond",
      face = "plain",
      size = 12,
      color = col_text
    ),
    strip.text = element_text(
      family = "eb_garamond",
      face = "plain",
      size = 12,
      color = col_text
    ),
    plot.background  = element_rect(fill = col_plot_bg, color = NA),
    panel.background = element_rect(fill = col_panel_bg, color = NA),
    axis.line        = element_line(color = col_axis, linewidth = 0.5),
    axis.ticks       = element_line(color = col_axis, linewidth = 0.5),
    axis.ticks.length = unit(3, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = col_grid_maj, linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

theme_bolb_map <- theme_bolb +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

theme_set(theme_bolb)

# -----------------------------
# Time windows
# -----------------------------
yr_before_lo <- 2008
yr_before_hi <- 2012

yr_during_lo <- 2013
yr_during_hi <- 2017

yr_after_lo  <- 2018
yr_after_hi  <- 2022

# -----------------------------
# Study extent
# -----------------------------
bbox <- list(
  swlat = 30,
  swlng = -180,
  nelat = 60,
  nelng = -110
)

# -----------------------------
# Get Sunflower Seastar data from iNaturalist
# -----------------------------
obs_seastar <- get_inat_obs(
  taxon_id = 47673,
  quality = "research",
  maxresults = 10000
) |>
  janitor::clean_names()

# -----------------------------
# Clean dates + coordinates
# -----------------------------
obs_seastar_clean <- obs_seastar |>
  mutate(
    observed_on = as.Date(observed_on),
    year = lubridate::year(observed_on)
  ) |>
  filter(
    !is.na(year),
    !is.na(latitude),
    !is.na(longitude)
  ) |>
  filter(year >= yr_before_lo, year <= yr_after_hi)

# -----------------------------
# Split data by time period
# -----------------------------
obs_seastar_before <- obs_seastar_clean |>
  filter(year >= yr_before_lo, year <= yr_before_hi)

obs_seastar_during <- obs_seastar_clean |>
  filter(year >= yr_during_lo, year <= yr_during_hi)

obs_seastar_after <- obs_seastar_clean |>
  filter(year >= yr_after_lo, year <= yr_after_hi)

# -----------------------------
# Basemap
# -----------------------------
countries_3 <- ne_countries(scale = "medium", returnclass = "sf") |>
  filter(admin %in% c("United States of America", "Canada", "Mexico"))

states_3 <- ne_states(returnclass = "sf") |>
  filter(geonunit %in% c("United States of America", "Canada", "Mexico"))

# -----------------------------
# Output folder + save helper
# -----------------------------
fig_dir <- here("figures")
dir_create(fig_dir)

save_plot <- function(p, filename, w = 10, h = 7, dpi = 400) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = p,
    width = w,
    height = h,
    units = "in",
    dpi = dpi
  )
}

# -----------------------------
# Helper for MCP
# -----------------------------
make_mcp_objects <- function(df, percent = 95, proj_epsg = 5070) {
  if (nrow(df) < 3) return(list(points_sf = NULL, mcp_sf = NULL))
  
  pts_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  pts_proj <- st_transform(pts_sf, proj_epsg)
  
  pts_sp <- SpatialPoints(
    coords = st_coordinates(pts_proj),
    proj4string = CRS(st_crs(pts_proj)$proj4string)
  )
  
  mcp_poly <- mcp(pts_sp, percent = percent)
  
  mcp_sf <- st_as_sf(mcp_poly) |>
    st_transform(4326)
  
  list(points_sf = pts_sf, mcp_sf = mcp_sf)
}

plot_mcp_map <- function(df, subtitle_txt) {
  res <- make_mcp_objects(df, percent = 95)
  
  if (is.null(res$mcp_sf)) {
    return(
      ggplot() +
        labs(
          title = "Sunflower sea star MCP",
          subtitle = paste0(subtitle_txt, " (not enough points for MCP)")
        ) +
        theme_bolb
    )
  }
  
  ggplot() +
    geom_sf(data = countries_3, fill = "white", color = "#434343") +
    geom_sf(data = states_3, fill = NA, color = "#434343", linewidth = 0.3) +
    geom_sf(data = res$mcp_sf, fill = "#faedcd", alpha = 0.5, color = NA) +
    geom_sf(data = res$mcp_sf, fill = NA, color = "#d4a373", linewidth = 1) +
    geom_sf(data = res$points_sf, color = "#003660", alpha = 0.6, size = 1.2) +
    coord_sf(
      xlim = c(bbox$swlng - 1, bbox$nelng + 1),
      ylim = c(bbox$swlat - 1, bbox$nelat + 1),
      expand = FALSE
    ) +
    labs(
      title = "Sunflower sea star MCP",
      subtitle = subtitle_txt
    ) +
    theme_bolb_map
}

# -----------------------------
# MCP maps
# -----------------------------
before_combined_mcp <- plot_mcp_map(obs_seastar_before, "Before SWD: 2008–2012")
during_combined_mcp <- plot_mcp_map(obs_seastar_during, "During SWD: 2013–2017")
after_combined_mcp  <- plot_mcp_map(obs_seastar_after,  "After SWD: 2018–2022")

before_combined_mcp
during_combined_mcp
after_combined_mcp

# -----------------------------
# Lollipop plot of annual latitude range
# -----------------------------
lat_summary <- obs_seastar_clean |>
  group_by(year) |>
  summarise(
    min_lat = min(latitude, na.rm = TRUE),
    median_lat = median(latitude, na.rm = TRUE),
    max_lat = max(latitude, na.rm = TRUE),
    .groups = "drop"
  )

lat_lollipop <- ggplot(lat_summary, aes(x = year)) +
  annotate(
    "rect",
    xmin = yr_during_lo,
    xmax = yr_during_hi,
    ymin = -Inf,
    ymax = Inf,
    fill = "#faedcd",
    alpha = 0.5
  ) +
  geom_text(
    aes(x = 2015, y = Inf, label = "SWD"),
    vjust = 1.5,
    family = "eb_garamond",
    face = "bold",
    color = "#d4a373",
    size = 4
  ) +
  geom_segment(
    aes(xend = year, y = min_lat, yend = max_lat),
    linewidth = 0.6,
    color = "#434343",
    alpha = 0.8
  ) +
  geom_point(aes(y = min_lat), size = 2, color = "#003660") +
  geom_point(aes(y = median_lat), size = 2.2, color = "#03859c") +
  geom_point(aes(y = max_lat), size = 2, color = "#09a99a") +
  scale_x_continuous(breaks = seq(min(lat_summary$year), max(lat_summary$year), by = 2)) +
  labs(
    x = NULL,
    y = "Latitude (°N)"
  ) +
  theme_bolb +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

lat_lollipop

# -----------------------------
# Plot of total latitudinal range per year
# -----------------------------
range_value <- lat_summary |>
  mutate(
    year = as.numeric(year),
    lat_range = max_lat - min_lat
  )

lat_range_plot <- ggplot(range_value, aes(x = year, y = lat_range)) +
  annotate(
    "rect",
    xmin = yr_during_lo,
    xmax = yr_during_hi,
    ymin = -Inf,
    ymax = Inf,
    fill = "#faedcd",
    alpha = 0.5
  ) +
  geom_text(
    aes(x = 2015, y = Inf, label = "SWD"),
    vjust = 1.5,
    family = "eb_garamond",
    face = "bold",
    color = "#d4a373",
    size = 4
  ) +
  geom_line(color = "#434343", linewidth = 0.6, alpha = 0.8) +
  geom_point(color = "#03859c", size = 2.2) +
  scale_x_continuous(breaks = seq(min(range_value$year), max(range_value$year), by = 2)) +
  labs(
    x = NULL,
    y = "Latitudinal range (max °N - min °N)"
  ) +
  theme_bolb +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

lat_range_plot

# -----------------------------
# Save plots
# -----------------------------
save_plot(before_combined_mcp, "before_combined_mcp.png", w = 6.5, h = 3.25)
save_plot(during_combined_mcp, "during_combined_mcp.png", w = 6.5, h = 3.25)
save_plot(after_combined_mcp,  "after_combined_mcp.png",  w = 6.5, h = 3.25)

save_plot(lat_lollipop, "lat_lollipop.png", w = 6.5, h = 3.25)
save_plot(lat_range_plot, "lat_range_plot.png", w = 6.5, h = 3.25)

save_plot(before_combined_mcp, "before_combined_mcp.png", w = 6.5, h = 3.25)
save_plot(during_combined_mcp, "during_combined_mcp.png", w = 6.5, h = 3.25)
save_plot(after_combined_mcp,  "after_combined_mcp.png",  w = 6.5, h = 3.25)

save_plot(lat_lollipop, "lat_lollipop.png", w = 6.5, h = 3.25)
save_plot(lat_range_plot, "lat_range_plot.png", w = 6.5, h = 3.25)