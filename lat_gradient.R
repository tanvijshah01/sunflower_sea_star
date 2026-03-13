# ESM 211: Sunflower Seastar Home Ranges
# Authors: Danielle Turner & Lilia Mourier
# Rebuilt: Mar 12, 2026

# -----------------------------
# Load libraries
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
library(fs)
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
# Output folder + save helper
# -----------------------------
fig_dir <- here("figures")
dir_create(fig_dir)

save_plot <- function(p, filename, w = 10, h = 7, dpi = 300) {
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
# Tunable parameters
# -----------------------------
acc_thresh_m     <- 5000
fixed_longitude  <- -128.5
mcp_percent      <- 95

# Equal 5-year periods
yr_before_lo <- 2008
yr_before_hi <- 2012

yr_during_lo <- 2013
yr_during_hi <- 2017

yr_after_lo  <- 2018
yr_after_hi  <- 2022

# Study extent
bbox <- list(
  swlat = 30,
  swlng = -180,
  nelat = 60,
  nelng = -110
)

# Histogram controls
bin_width_deg  <- 0.5
use_lat_window <- TRUE
lat_min <- 30
lat_max <- 60

# Region labels
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
# Basemaps
# -----------------------------
countries_3 <- ne_countries(scale = "medium", returnclass = "sf") |>
  filter(admin %in% c("United States of America", "Canada", "Mexico"))

states_3 <- ne_states(returnclass = "sf") |>
  filter(geonunit %in% c("United States of America", "Canada", "Mexico"))

# -----------------------------
# Load and clean sunflower sea star data
# -----------------------------
obs_seastar_raw <- read.csv(here("obs_seastar.csv")) |>
  janitor::clean_names()

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
  ) |>
  filter(year >= yr_before_lo, year <= yr_after_hi)

# -----------------------------
# Plot all filtered sea star observations together
# -----------------------------
all_seastar_points_sf <- st_as_sf(
  obs_seastar_filt,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

all_points_map <- ggplot() +
  geom_sf(data = countries_3, fill = "white", color = "#434343") +
  geom_sf(data = states_3, fill = NA, color = "#434343", linewidth = 0.3) +
  geom_sf(
    data = all_seastar_points_sf,
    color = "#003660",
    fill = "#003660",
    alpha = 0.35,
    size = 1.2
  ) +
  coord_sf(
    xlim = c(bbox$swlng - 1, bbox$nelng + 1),
    ylim = c(bbox$swlat - 1, bbox$nelat + 1),
    expand = FALSE
  )+
  theme_bolb_map+
  theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 20))
  # theme(
  #   axis.title = element_blank(),
  #   axis.text  = element_blank(),
  #   axis.text.x = element_blank(),
  #   axis.text.y = element_blank(),
  #   axis.ticks = element_blank(),
  #   axis.line  = element_blank()
  # )

# -----------------------------
# Pull and clean anemone data for effort standardization
# -----------------------------
# Aggregating anemone
anemone_id <- 52564

inat_bounds <- c(
  bbox$swlat,
  bbox$swlng,
  bbox$nelat,
  bbox$nelng
)

obs_anemone_raw <- get_inat_obs(
  taxon_id = anemone_id,
  quality = "research",
  bounds = inat_bounds,
  maxresults = 10000
) |>
  janitor::clean_names()

obs_anemone_filt <- obs_anemone_raw |>
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
  ) |>
  filter(year >= yr_before_lo, year <= yr_after_hi)

# -----------------------------
# Annual effort standardization
# -----------------------------
yearly_sunflower <- obs_seastar_filt |>
  count(year, name = "sunflower_obs")

yearly_anemone <- obs_anemone_filt |>
  count(year, name = "anemone_obs")

median_anemone_effort <- median(yearly_anemone$anemone_obs, na.rm = TRUE)

yearly_effort_index <- yearly_sunflower |>
  full_join(yearly_anemone, by = "year") |>
  arrange(year) |>
  mutate(
    sunflower_obs = coalesce(sunflower_obs, 0L),
    anemone_obs   = coalesce(anemone_obs, 0L),
    effort_index = if_else(anemone_obs > 0, sunflower_obs / anemone_obs, NA_real_),
    norm_weight = if_else(anemone_obs > 0, median_anemone_effort / anemone_obs, NA_real_),
    period = case_when(
      year >= yr_before_lo & year <= yr_before_hi ~ "Before SWD",
      year >= yr_during_lo & year <= yr_during_hi ~ "During SWD",
      year >= yr_after_lo & year <= yr_after_hi ~ "After SWD",
      TRUE ~ NA_character_
    )
  )

readr::write_csv(
  yearly_effort_index,
  file.path(fig_dir, "yearly_effort_standardized_index.csv")
)

# -----------------------------
# Raw and normalized sea star plotting datasets
# -----------------------------
obs_seastar_raw_plot <- obs_seastar_filt

obs_seastar_norm_plot <- obs_seastar_filt |>
  left_join(
    yearly_effort_index |>
      select(year, norm_weight),
    by = "year"
  ) |>
  filter(!is.na(norm_weight))

# -----------------------------
# Split raw and normalized data into equal periods
# -----------------------------
obs_before_clean_raw <- obs_seastar_raw_plot |>
  filter(year >= yr_before_lo, year <= yr_before_hi)

obs_during_clean_raw <- obs_seastar_raw_plot |>
  filter(year >= yr_during_lo, year <= yr_during_hi)

obs_after_clean_raw <- obs_seastar_raw_plot |>
  filter(year >= yr_after_lo, year <= yr_after_hi)

obs_before_clean_norm <- obs_seastar_norm_plot |>
  filter(year >= yr_before_lo, year <= yr_before_hi)

obs_during_clean_norm <- obs_seastar_norm_plot |>
  filter(year >= yr_during_lo, year <= yr_during_hi)

obs_after_clean_norm <- obs_seastar_norm_plot |>
  filter(year >= yr_after_lo, year <= yr_after_hi)

obs_before_simple_raw <- obs_before_clean_raw |>
  mutate(longitude = fixed_longitude)

obs_during_simple_raw <- obs_during_clean_raw |>
  mutate(longitude = fixed_longitude)

obs_after_simple_raw <- obs_after_clean_raw |>
  mutate(longitude = fixed_longitude)

obs_before_simple_norm <- obs_before_clean_norm |>
  mutate(longitude = fixed_longitude)

obs_during_simple_norm <- obs_during_clean_norm |>
  mutate(longitude = fixed_longitude)

obs_after_simple_norm <- obs_after_clean_norm |>
  mutate(longitude = fixed_longitude)

# -----------------------------
# Helpers
# -----------------------------
make_mcp_sf <- function(df, percent = 95, proj_epsg = 5070, use_weights = FALSE) {
  df <- df |>
    filter(!is.na(latitude), !is.na(longitude))
  
  if (use_weights && "norm_weight" %in% names(df)) {
    min_w <- suppressWarnings(min(df$norm_weight[df$norm_weight > 0], na.rm = TRUE))
    
    if (is.finite(min_w)) {
      df <- df |>
        mutate(weight_n = pmax(1, round(norm_weight / min_w))) |>
        tidyr::uncount(weights = weight_n)
    }
  }
  
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

plot_clean_map <- function(df_clean, title_txt, use_weights = FALSE) {
  res <- make_mcp_sf(df_clean, percent = mcp_percent, use_weights = use_weights)
  
  if (is.null(res$mcp_sf)) {
    return(
      ggplot() +
        labs(title = paste0(title_txt, " (not enough points for MCP)")) +
        theme_bolb
    )
  }
  
  ggplot() +
    geom_sf(data = countries_3, fill = "white", color = "#434343") +
    geom_sf(data = states_3, fill = NA, color = "#434343", linewidth = 0.3) +
    geom_sf(data = res$mcp_sf, fill = "#faedcd", alpha = 0.5, color = NA) +
    geom_sf(data = res$mcp_sf, fill = NA, color = "#d4a373", linewidth = 1) +
    geom_sf(data = res$points_sf, color = "#003660", fill = "#003660", alpha = 0.5, size = 1.5) +
    coord_sf(
      xlim = c(bbox$swlng - 1, bbox$nelng + 1),
      ylim = c(bbox$swlat - 1, bbox$nelat + 1),
      expand = FALSE
    ) +
    labs(title = title_txt) +
    theme_bolb_map+
    theme(plot.margin = margin(t = 5, r = 10, b = 10, l = 5))
    # theme(
    #   axis.title = element_blank(),
    #   axis.text  = element_blank(),
    #   axis.text.x = element_blank(),
    #   axis.text.y = element_blank(),
    #   axis.ticks = element_blank(),
    #   axis.line  = element_blank()
    # )
}

plot_simple_lat_hist <- function(df_simple, title_txt, use_weights = FALSE) {
  res <- make_mcp_sf(df_simple, percent = mcp_percent, use_weights = use_weights)
  
  plot_df <- df_simple |>
    filter(!is.na(latitude))
  
  if (use_lat_window) {
    plot_df <- plot_df |>
      filter(latitude >= lat_min, latitude <= lat_max)
  }
  
  if (use_weights && "norm_weight" %in% names(plot_df)) {
    plot_df <- plot_df |>
      mutate(plot_weight = norm_weight)
  } else {
    plot_df <- plot_df |>
      mutate(plot_weight = 1)
  }
  
  x_limits <- if (use_lat_window) c(lat_min, lat_max) else range(plot_df$latitude, na.rm = TRUE)
  
  lat_line <- tibble::tibble(lat_lo = NA_real_, lat_hi = NA_real_)
  if (!is.null(res$mcp_sf)) {
    bb <- sf::st_bbox(res$mcp_sf)
    lat_line <- tibble::tibble(
      lat_lo = as.numeric(bb["ymin"]),
      lat_hi = as.numeric(bb["ymax"])
    )
  }
  
  regions_bounds_plot <- regions_bounds |>
    filter(lat >= (x_limits[1] - 1), lat <= (x_limits[2] + 1)) |>
    mutate(
      vjust_use = case_when(
        label == "Washington" ~ 2.0,
        label == "Canada" ~ 4,
        TRUE ~ 3
      )
    )
  
  x_pad <- 0.8
  y_pad <- 0.18
  
  hist_counts <- ggplot_build(
    ggplot(plot_df, aes(x = latitude, weight = plot_weight)) +
      geom_histogram(binwidth = bin_width_deg, boundary = 0, closed = "left")
  )$data[[1]]$count
  
  y_max <- if (length(hist_counts) == 0) 1 else max(hist_counts, na.rm = TRUE)
  x_limits_pad <- c(x_limits[1] - x_pad, x_limits[2] + x_pad)
  if (use_weights) {
    y_limits_pad <- c(0, 125)
    label_y <- 95
  } else {
    y_limits_pad <- c(0, 125)
    label_y <- 475
  }
  
  theme_bolb_more_space <- theme_bolb +
    theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 20))
  
  ggplot(plot_df, aes(x = latitude, weight = plot_weight)) +
    geom_rect(
      data = lat_line,
      aes(xmin = lat_lo, xmax = lat_hi, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#faedcd",
      alpha = 0.5,
      color = NA
    ) +
    geom_vline(
      data = regions_bounds_plot,
      aes(xintercept = lat),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.5,
      color = "#434343"
    ) +
    geom_histogram(
      binwidth = bin_width_deg,
      boundary = 0,
      closed = "left",
      fill = "#003660",
      color = "#003660",
      alpha = 0.25
    ) +
    # geom_density(
    #   aes(y = after_stat(count) * bin_width_deg),
    #   linewidth = 1,
    #   alpha = 0.25,
    #   color = "#003660",
    #   fill = "#003660"
    # ) +
    geom_vline(
      xintercept = lat_line$lat_lo,
      linetype = "dashed",
      linewidth = 0.6,
      color = "#d4a373",
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = lat_line$lat_hi,
      linetype = "dashed",
      linewidth = 0.6,
      color = "#d4a373",
      na.rm = TRUE
    ) +
    geom_text(
      data = regions_bounds_plot,
      aes(x = lat + 0.15, y = Inf, label = label, vjust = vjust_use),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4,
      family = "eb_garamond",
      color = "#434343"
    ) +
    # annotate(
    #   "text",
    #   x = mean(c(lat_line$lat_lo, lat_line$lat_hi), na.rm = TRUE),
    #   y = label_y,
    #   label = "MCP home range",
    #   vjust = 3.5,
    #   hjust = -0.8,
    #   family = "eb_garamond",
    #   fontface = "bold",
    #   size = 4,
    #   color = "#d4a373"
    # ) +
    coord_cartesian(
      xlim = x_limits_pad,
      ylim = y_limits_pad,
      clip = "on"
    ) +
    labs(
      title = title_txt,
      x = "Latitude (°N)",
      y = ifelse(use_weights, "Weighted\n observations", "Number of observations")
    ) +
    theme_bolb_more_space
}

# -----------------------------
# Standardization figure
# -----------------------------
p_sunflower_raw <- ggplot(yearly_effort_index, aes(x = year, y = sunflower_obs)) +
  geom_col(fill = "#003660", alpha = 0.85) +
  geom_vline(xintercept = c(2013, 2018), linetype = "dashed", color = "#434343") +
  scale_x_continuous(breaks = seq(min(yearly_effort_index$year), max(yearly_effort_index$year), by = 2)) +
  labs(
    title = "A. Sunflower seastar observations",
    x = "Year",
    y = "Observations"
  ) +
  theme_bolb

p_anemone_raw <- ggplot(yearly_effort_index, aes(x = year, y = anemone_obs)) +
  geom_col(fill = "#d4a373", alpha = 0.85) +
  geom_vline(xintercept = c(2013, 2018), linetype = "dashed", color = "#434343") +
  scale_x_continuous(breaks = seq(min(yearly_effort_index$year), max(yearly_effort_index$year), by = 2)) +
  labs(
    title = "B. Anemone observations",
    x = "Year",
    y = "Observations"
  ) +
  theme_bolb

p_standardized <- ggplot(yearly_effort_index, aes(x = year, y = effort_index)) +
  geom_line(color = "#003660", linewidth = 1.2) +
  geom_point(color = "#003660", size = 2.5) +
  geom_vline(xintercept = c(2013, 2018), linetype = "dashed", color = "#434343") +
  scale_x_continuous(breaks = seq(min(yearly_effort_index$year), max(yearly_effort_index$year), by = 2)) +
  labs(
    title = "C. Standardized sunflower seastar occurrence",
    x = "Year",
    y = "Relative observations\n (seastar / anemone)"
  ) +
  theme_bolb

standardization_figure <- p_sunflower_raw / p_anemone_raw / p_standardized

save_plot(standardization_figure, "standardization_workflow_anemone.png", w = 6.5, h = 8)
save_plot(standardization_figure, "standardization_workflow_anemone.pdf", w = 6.5, h = 8)

# -----------------------------
# Build RAW plots
# -----------------------------
raw_before_map <- plot_clean_map(obs_before_clean_raw, "Before", use_weights = FALSE)
raw_during_map <- plot_clean_map(obs_during_clean_raw, "During", use_weights = FALSE)
raw_after_map  <- plot_clean_map(obs_after_clean_raw,  "After", use_weights = FALSE)

raw_before_hist <- plot_simple_lat_hist(obs_before_simple_raw, "Before", use_weights = FALSE)
raw_during_hist <- plot_simple_lat_hist(obs_during_simple_raw, "During", use_weights = FALSE)
raw_after_hist  <- plot_simple_lat_hist(obs_after_simple_raw,  "After", use_weights = FALSE)

# -----------------------------
# Build NORMALIZED plots
# -----------------------------
norm_before_map <- plot_clean_map(obs_before_clean_norm, "Before", use_weights = TRUE)
norm_during_map <- plot_clean_map(obs_during_clean_norm, "During", use_weights = TRUE)
norm_after_map  <- plot_clean_map(obs_after_clean_norm,  "After", use_weights = TRUE)

norm_before_hist <- plot_simple_lat_hist(obs_before_simple_norm, "Before", use_weights = TRUE)
norm_during_hist <- plot_simple_lat_hist(obs_during_simple_norm, "During", use_weights = TRUE)
norm_after_hist  <- plot_simple_lat_hist(obs_after_simple_norm,  "After", use_weights = TRUE)

# -----------------------------
# Stack plots
# -----------------------------
homerange_raw <- raw_before_map + raw_during_map + raw_after_map +
  plot_layout(axes = "collect_y", axis_titles = "collect_y")

lat_gradient_raw <- raw_before_hist / raw_during_hist / raw_after_hist +
  plot_layout(axes = "collect_x", axis_titles = "collect_x")

homerange_norm <- norm_before_map + norm_during_map + norm_after_map +
  plot_layout(axes = "collect_y", axis_titles = "collect_y")

lat_gradient_norm <- norm_before_hist / norm_during_hist / norm_after_hist +
  plot_layout(axes = "collect_x", axis_titles = "collect_x")

# Optional side-by-side comparison layouts
homerange_compare <- homerange_raw | homerange_norm
lat_gradient_compare <- lat_gradient_raw | lat_gradient_norm

# -----------------------------
# Print plots
# -----------------------------
all_points_map

standardization_figure

homerange_raw
lat_gradient_raw

homerange_norm
lat_gradient_norm

homerange_compare
lat_gradient_compare

# -----------------------------
# Save plots
# -----------------------------

save_plot(all_points_map, "all_seastar_points_map.png", w = 6.5, h = 4.5)
save_plot(all_points_map, "all_seastar_points_map.pdf", w = 6.5, h = 4.5)

save_plot(raw_before_map,  "raw_map_before_clean_mcp.png",  w = 6.5, h = 3.25)
save_plot(raw_during_map,  "raw_map_during_clean_mcp.png",  w = 6.5, h = 3.25)
save_plot(raw_after_map,   "raw_map_after_clean_mcp.png",   w = 6.5, h = 3.25)

save_plot(raw_before_hist, "raw_hist_before_latitude.png",  w = 6.5, h = 3.25)
save_plot(raw_during_hist, "raw_hist_during_latitude.png",  w = 6.5, h = 3.25)
save_plot(raw_after_hist,  "raw_hist_after_latitude.png",   w = 6.5, h = 3.25)

save_plot(norm_before_map,  "norm_map_before_clean_mcp.png", w = 6.5, h = 3.25)
save_plot(norm_during_map,  "norm_map_during_clean_mcp.png", w = 6.5, h = 3.25)
save_plot(norm_after_map,   "norm_map_after_clean_mcp.png",  w = 6.5, h = 3.25)

save_plot(norm_before_hist, "norm_hist_before_latitude.png", w = 6.5, h = 3.25)
save_plot(norm_during_hist, "norm_hist_during_latitude.png", w = 6.5, h = 3.25)
save_plot(norm_after_hist,  "norm_hist_after_latitude.png",  w = 6.5, h = 3.25)

save_plot(homerange_raw,      "stack_homerange_maps_raw.png",         w = 6.5, h = 8)
save_plot(lat_gradient_raw,   "stack_latitude_gradients_raw.png",     w = 6.5, h = 8)
save_plot(homerange_norm,     "stack_homerange_maps_normalized.png",  w = 6.5, h = 8)
save_plot(lat_gradient_norm,  "stack_latitude_gradients_normalized.png", w = 6.5, h = 8)
save_plot(standardization_figure,  "standardization_figure.png", w = 6.5, h = 8)

save_plot(homerange_compare,   "compare_homerange_maps_raw_vs_norm.png",   w = 9, h = 6.5)
save_plot(lat_gradient_compare, "compare_latitude_raw_vs_norm.png",        w = 9, h = 6.5)

save_plot(homerange_raw,      "stack_homerange_maps_raw.pdf",         w = 6.5, h = 8)
save_plot(lat_gradient_raw,   "stack_latitude_gradients_raw.pdf",     w = 6.5, h = 8)
save_plot(homerange_norm,     "stack_homerange_maps_normalized.pdf",  w = 6.5, h = 8)
save_plot(lat_gradient_norm,  "stack_latitude_gradients_normalized.pdf", w = 6.5, h = 8)

save_plot(homerange_compare,   "compare_homerange_maps_raw_vs_norm.pdf",   w = 8, h = 6.5)
save_plot(lat_gradient_compare, "compare_latitude_raw_vs_norm.pdf",        w = 8, h = 6.5)
