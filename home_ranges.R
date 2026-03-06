# ESM 211: Sunflower Seastar Homeranges
# Authors: Danielle Turner & Lilia Mourier
# Last updated on: Mar 5, 2026


# -----------------------------
# Load Libraries
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
library(dplyr)
library(lubridate)
library(fs)

# Theme & design libraries
library(showtext)
library(sysfonts)

# -----------------------------
# Theme setup (your BOLB theme)
# -----------------------------
font_add_google("Merriweather", "merriweather")
font_add_google("Source Sans 3", "source_sans")

showtext_auto(TRUE)
showtext_opts(dpi = 300)

col_title <- "#003660"
col_text  <- "#434343"

col_panel_bg <- "#FFFFFF"
col_plot_bg  <- "#FFFFFF"

col_axis     <- "#434343"
col_grid_maj <- "#D9D9D9"

theme_bolb <-
  theme_minimal(base_size = 12, base_family = "source_sans") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      family = "merriweather",
      face = "plain",
      size = 18,
      color = col_title,
      hjust = 0,
      margin = margin(b = 8)
    ),
    plot.subtitle.position = "plot",
    plot.subtitle = element_text(
      family = "source_sans",
      face = "plain",
      size = 16,
      color = col_text,
      hjust = 0,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(
      family = "source_sans",
      face = "bold",
      size = 16,
      color = col_text,
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      family = "source_sans",
      face = "bold",
      size = 16,
      color = col_text,
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      family = "source_sans",
      size = 14,
      color = col_text,
      margin = margin(t = 4)
    ),
    axis.text.y = element_text(
      family = "source_sans",
      size = 14,
      color = col_text,
      margin = margin(r = 4)
    ),
    legend.title = element_text(
      family = "source_sans",
      face = "plain",
      size = 16,
      color = col_text
    ),
    legend.text = element_text(
      family = "source_sans",
      face = "plain",
      size = 16,
      color = col_text
    ),
    strip.text = element_text(
      family = "source_sans",
      face = "plain",
      size = 14,
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
    plot.margin = margin(t = 10, r = 0, b = 0, l = 0)
  )

# A map-friendly variant: keeps title/subtitle styling, removes axes (like theme_void)
theme_bolb_map <- theme_bolb +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Default theme (used by histogram plots)
theme_set(theme_bolb)

# -----------------------------
# Tunable parameters
# -----------------------------
acc_thresh_m     <- 5000     # positional_accuracy threshold in meters
fixed_longitude  <- -128.5   # chosen longitude for "simple" workflow
mcp_percent      <- 95

# Time windows
yr_min_before <- 2004
yr_before_end <- 2013   # before: 2004–2012
yr_during_lo  <- 2013   # during: 2013–2017
yr_during_hi  <- 2017
yr_after_lo   <- 2018   # after: 2018–present

# Map extent (for clean home range maps)
bbox <- list(
  swlat = 30,
  swlng = -180,
  nelat = 60,
  nelng = -110
)

# Histogram controls (for simple latitude plots)
bin_width_deg  <- 0.5
use_lat_window <- TRUE
lat_min <- 30
lat_max <- 60

# Latitude boundary lines (lower bounds for each region + upper Alaska bound)
regions_bounds <- tibble::tribble(
  ~label,          ~lat,
  "Baja Mexico",    22.876,
  "California",     32.533,
  "Oregon",         42.000,
  "Washington",     45.5436,
  "Canada",         48.30,
  "Alaska",         51.262,
  "Alaska max",     71.390
)

if (use_lat_window) {
  regions_bounds <- regions_bounds |>
    dplyr::filter(lat >= lat_min, lat <= lat_max)
}

# -----------------------------
# Basemaps (US + Canada + Mexico)
# -----------------------------
countries_3 <- ne_countries(scale = "medium", returnclass = "sf") |>
  filter(admin %in% c("United States of America", "Canada", "Mexico"))

states_3 <- ne_states(returnclass = "sf") |>
  filter(geonunit %in% c("United States of America", "Canada", "Mexico"))

# -----------------------------
# Load data
# -----------------------------
obs_seastar_raw <- read.csv(here("obs_seastar.csv"))

# -----------------------------
# Tidy + filter (quality + obscured + positional accuracy) then create CLEAN + SIMPLE
# -----------------------------
obs_seastar_filt <- obs_seastar_raw |>
  filter(!is.na(observed_on), observed_on != "") |>
  mutate(
    year = suppressWarnings(as.integer(substr(observed_on, 1, 4))),
    positional_accuracy = suppressWarnings(as.numeric(positional_accuracy)),
    coordinates_obscured = tolower(as.character(coordinates_obscured)) == "true"
  ) |>
  filter(
    !is.na(latitude), !is.na(longitude),
    !is.na(positional_accuracy),
    !is.na(year),
    positional_accuracy < acc_thresh_m,
    quality_grade == "research",
    coordinates_obscured == FALSE
  )

# Clean: keep true lon/lat
obs_seastar_clean <- obs_seastar_filt

# Simple: overwrite longitude with chosen constant
obs_seastar_simple <- obs_seastar_filt |>
  mutate(longitude = fixed_longitude)

# -----------------------------
# Split into periods (CLEAN + SIMPLE)
# -----------------------------
obs_before_clean <- obs_seastar_clean |>
  filter(year >= yr_min_before, year < yr_before_end)

obs_during_clean <- obs_seastar_clean |>
  filter(year >= yr_during_lo, year <= yr_during_hi)

obs_after_clean <- obs_seastar_clean |>
  filter(year >= yr_after_lo)

obs_before_simple <- obs_seastar_simple |>
  filter(year >= yr_min_before, year < yr_before_end)

obs_during_simple <- obs_seastar_simple |>
  filter(year >= yr_during_lo, year <= yr_during_hi)

obs_after_simple <- obs_seastar_simple |>
  filter(year >= yr_after_lo)

# -----------------------------
# Helpers
# -----------------------------
make_mcp_sf <- function(df, percent = 99, proj_epsg = 5070) {
  df <- df |>
    filter(!is.na(latitude), !is.na(longitude))
  
  if (nrow(df) < 3) {
    return(list(points_sf = NULL, mcp_sf = NULL))
  }
  
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

plot_clean_map <- function(df_clean, subtitle_txt) {
  res <- make_mcp_sf(df_clean, percent = mcp_percent)
  
  if (is.null(res$mcp_sf)) {
    return(
      ggplot() +
        labs(
          title = paste0(subtitle_txt, " (not enough points for MCP)")
          #subtitle = paste0(subtitle_txt, " (not enough points for MCP)")
        ) +
        theme_bolb
    )
  }
  
  ggplot() +
    geom_sf(data = countries_3, fill = "white", color = "#434343") +
    geom_sf(data = states_3, fill = NA, color = "#434343", linewidth = 0.3) +
    geom_sf(data = res$mcp_sf, fill = "#faedcd", alpha = 0.5, color = NA) +
    geom_sf(data = res$mcp_sf, fill = NA, color = "#d4a373", linewidth = 1) +
    geom_sf(data = res$points_sf, color = "#003660", fill = "#003660", alpha = 0.6, size = 2) +
    coord_sf(
      xlim = c(bbox$swlng - 1, bbox$nelng + 1),
      ylim = c(bbox$swlat - 1, bbox$nelat + 1),
      expand = FALSE
    ) +
    labs(
      title = subtitle_txt,
      #subtitle = subtitle_txt
    ) + 
    theme_bolb_map +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank()
    )
}

plot_simple_lat_hist <- function(df_simple, subtitle_txt) {
  # MCP first (on simple points), then lat extent from the MCP polygon bbox
  res <- make_mcp_sf(df_simple, percent = mcp_percent)
  
  plot_df <- df_simple |>
    dplyr::filter(!is.na(latitude))
  
  if (use_lat_window) {
    plot_df <- plot_df |>
      dplyr::filter(latitude >= lat_min, latitude <= lat_max)
  }
  
  x_limits <- if (use_lat_window) c(lat_min, lat_max) else range(plot_df$latitude, na.rm = TRUE)
  
  # MCP latitude range
  lat_line <- tibble::tibble(lat_lo = NA_real_, lat_hi = NA_real_)
  if (!is.null(res$mcp_sf)) {
    bb <- sf::st_bbox(res$mcp_sf)
    lat_line <- tibble::tibble(
      lat_lo = as.numeric(bb["ymin"]),
      lat_hi = as.numeric(bb["ymax"])
    )
  }
  
  # Keep only boundary labels/lines that are in (or near) the plotting window
  regions_bounds_plot <- regions_bounds |>
    dplyr::filter(lat >= (x_limits[1] - 1), lat <= (x_limits[2] + 1))
  
  regions_bounds_plot <- regions_bounds_plot |>
    dplyr::mutate(
      vjust_use = dplyr::case_when(
        label == "Washington" ~ 1.0,  # slightly higher
        label == "Canada"     ~ 2.5,  # slightly lower
        TRUE                  ~ 1.75   # default for everyone else
      )
    )
  
  # ----- Add padding so labels aren't clipped -----
  x_pad <- 0.8   # degrees of latitude padding on left/right (tweak)
  y_pad <- 0.18  # percent headroom on top for labels (tweak)
  
  # Compute max histogram bin count to set a y-limit with headroom
  y_max <- plot_df |>
    dplyr::mutate(lat_bin = floor(latitude / bin_width_deg) * bin_width_deg) |>
    dplyr::count(lat_bin, name = "n") |>
    dplyr::summarise(mx = max(n, na.rm = TRUE)) |>
    dplyr::pull(mx)
  
  x_limits_pad <- c(x_limits[1] - x_pad, x_limits[2] + x_pad)
  y_limits_pad <- c(0, y_max * (1 + y_pad))
  
  # Slightly larger margins for whitespace + label room
  theme_bolb_more_space <- theme_bolb +
    theme(
      plot.margin = margin(t = 40, r = 70, b = 30, l = 55)
    )
  
  ggplot(plot_df, aes(x = latitude)) +
    # MCP shaded band (background)
    geom_rect(
      data = lat_line,
      aes(xmin = lat_lo, xmax = lat_hi, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#faedcd",
      alpha = 0.5,
      color = NA
    ) +
    # Region boundary lines
    geom_vline(
      data = regions_bounds_plot,
      aes(xintercept = lat),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.5,
      color = "#434343"
    ) +
    geom_histogram(
      binwidth = bin_width_deg, boundary = 0, closed = "left",
      fill = "#003660", color = "#003660", alpha = 0.25
    ) +
    geom_density(
      aes(y = after_stat(count) * bin_width_deg),
      linewidth = 1, alpha = 0.25,
      color = "#003660", fill = "#003660"
    ) +
    geom_vline(
      xintercept = lat_line$lat_lo,
      linetype = "dashed", linewidth = 0.6,
      color = "#d4a373",
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = lat_line$lat_hi,
      linetype = "dashed", linewidth = 0.6,
      color = "#d4a373",
      na.rm = TRUE
    ) +
    # Labels just to the right of each boundary line
    geom_text(
      data = regions_bounds_plot,
      aes(x = lat + 0.15, y = Inf, label = label, vjust = vjust_use),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4,
      family = "source_sans",
      color = "#434343"
    ) +
    # Label for MCP band (near top, inside panel; clipping off allows margin space too)
    annotate(
      "text",
      x = mean(c(lat_line$lat_lo, lat_line$lat_hi), na.rm = TRUE),
      y = y_limits_pad[2],
      label = "MCP home range",
      vjust = 3.5,
      hjust = -.8,
      family = "source_sans",
      fontface = "bold",
      size = 4,
      color = "#d4a373"
    ) +
    coord_cartesian(
      xlim = x_limits_pad,
      ylim = y_limits_pad,
      clip = "off"
    ) +
    labs(
      title = subtitle_txt,
      x = "Latitude (°N)",
      y = "Number of observations"
    ) +
    theme_bolb_more_space
}
# -----------------------------
# Build ALL plots (clean map + simple histogram) for each period
# -----------------------------
# BEFORE
shr_before_clean_map <- plot_clean_map(obs_before_clean, "Before SWD: 2004–2012")
shr_before_simple_hist <- plot_simple_lat_hist(obs_before_simple, "Before SWD: 2004–2012")

# DURING
shr_during_clean_map <- plot_clean_map(obs_during_clean, "During SWD: 2013–2017")
shr_during_simple_hist <- plot_simple_lat_hist(obs_during_simple, "During SWD: 2013–2017")

# AFTER
shr_after_clean_map <- plot_clean_map(obs_after_clean, "After SWD: 2018–present")
shr_after_simple_hist <- plot_simple_lat_hist(obs_after_simple, "After SWD: 2018–present")

# -----------------------------
# Print plots (individually)
# -----------------------------
shr_before_clean_map
shr_before_simple_hist

shr_during_clean_map
shr_during_simple_hist

shr_after_clean_map
shr_after_simple_hist

# -----------------------------
# Optional: stack each period as "map over histogram"
# -----------------------------
(homerange <- shr_before_clean_map / shr_during_clean_map / shr_after_clean_map)
(lat_gradient <- shr_before_simple_hist / shr_during_simple_hist / shr_after_simple_hist)

# Optional: stack all periods
# (before_combo) / (during_combo) / (after_combo)

# -----------------------------
# Save plots to /figures in project folder
# -----------------------------

fig_dir <- here("figures")
dir_create(fig_dir)  # creates folder if it doesn't exist

# helper to save with consistent settings
save_plot <- function(p, filename, w = 10, h = 7, dpi = 300) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = p,
    width = w, height = h, units = "in",
    dpi = dpi
  )
}

# --- individual plots ---
save_plot(shr_before_clean_map,   "map_before_clean_mcp.png",   w = 10, h = 7)
save_plot(shr_during_clean_map,   "map_during_clean_mcp.png",   w = 10, h = 7)
save_plot(shr_after_clean_map,    "map_after_clean_mcp.png",    w = 10, h = 7)

save_plot(shr_before_simple_hist, "hist_before_latitude.png",   w = 10, h = 7)
save_plot(shr_during_simple_hist, "hist_during_latitude.png",   w = 10, h = 7)
save_plot(shr_after_simple_hist,  "hist_after_latitude.png",    w = 10, h = 7)

# --- stacked plots (if you created them) ---
save_plot(homerange,    "stack_homerange_maps.png",   w = 10, h = 14)
save_plot(lat_gradient, "stack_latitude_gradients.png", w = 10, h = 14)

# Optional: also save PDFs (nice for papers/slides)
save_plot(homerange,    "stack_homerange_maps.pdf",   w = 10, h = 14)
save_plot(lat_gradient, "stack_latitude_gradients.pdf", w = 10, h = 14)
