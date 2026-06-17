library(sf)
library(tidyverse)
library(terra)
library(ncdf4)
library(yaml)

cfg       <- yaml.load_file('config.yml')
data_root <- cfg$environments[[cfg$environment]]$data_root
d         <- function(subpath) file.path(data_root, subpath)

OUT_DIR <- cfg$results$dir
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BUFFER_M           <- cfg$params$buffer_m
SDD_MIN            <- cfg$params$sdd_window_min
SDD_MAX            <- cfg$params$sdd_window_max
LIGHTNING_WINDOW   <- cfg$params$lightning_window_days
LIGHTNING_BUFFER_M <- cfg$params$lightning_buffer_m


#  Lightning exclusion 
# Returns a logical vector (length = nrow(hotspots_sf)):
#   TRUE  → lightning strike within buffer_m and window_days → exclude
#   FALSE → keep
#
# NetCDF cg_flashes has dimensions (time, lat, lon) in C order.
# ncdf4 reverses to Fortran order, so R array is [lon, lat, time].

check_lightning <- function(hotspots_sf, nc_path, buffer_m, window_days) {
  nc <- nc_open(nc_path)
  on.exit(nc_close(nc))

  lat_vals   <- ncvar_get(nc, "lat")
  lon_vals   <- ncvar_get(nc, "lon")
  time_raw   <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value
  n_time     <- length(time_raw)

  origin_str <- sub("^.+ since ", "", time_units)
  origin     <- as.POSIXct(origin_str, tz = "UTC")
  time_dates <- as.Date(origin + as.numeric(time_raw) * 3600)

  # Clip to hotspot bounding box + 0.3° pad (>> 2 km lightning buffer at this latitude).
  # Reads only the study-area subset instead of the full continental grid.
  bbox     <- st_bbox(st_transform(hotspots_sf, 4326))
  pad      <- 0.3
  lon_idx  <- which(lon_vals >= bbox["xmin"] - pad & lon_vals <= bbox["xmax"] + pad)
  lat_idx  <- which(lat_vals >= bbox["ymin"] - pad & lat_vals <= bbox["ymax"] + pad)

  if (length(lon_idx) == 0L || length(lat_idx) == 0L) return(rep(FALSE, nrow(hotspots_sf)))

  lon_start <- min(lon_idx);  lon_count <- max(lon_idx) - lon_start + 1
  lat_start <- min(lat_idx);  lat_count <- max(lat_idx) - lat_start + 1
  lon_sub   <- lon_vals[lon_start:(lon_start + lon_count - 1)]
  lat_sub   <- lat_vals[lat_start:(lat_start + lat_count - 1)]

  # ncdf4 start/count follow R (Fortran) dimension order: [lon, lat, time]
  flashes <- ncvar_get(nc, "cg_flashes",
                       start = c(lon_start, lat_start, 1),
                       count = c(lon_count, lat_count, n_time))

  # Reshape to [n_cells, n_time]; lon varies fastest (matches expand.grid row order)
  flash_mat <- matrix(as.numeric(flashes), nrow = lon_count * lat_count, ncol = n_time)
  rm(flashes)

  first_idx <- apply(flash_mat > 0, 1, function(row) {
    idx <- which(row)
    if (length(idx) == 0L) NA_integer_ else idx[1L]
  })
  rm(flash_mat)

  grid <- expand.grid(lon = lon_sub, lat = lat_sub)
  grid$first_date <- as.Date(
    ifelse(is.na(first_idx), NA_real_, as.numeric(time_dates[first_idx])),
    origin = "1970-01-01"
  )
  grid <- grid[!is.na(grid$first_date), ]

  if (nrow(grid) == 0L) return(rep(FALSE, nrow(hotspots_sf)))

  grid_sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326) |>
    st_transform(3005)

  hs_3005       <- st_transform(hotspots_sf, 3005)
  hotspot_dates <- as.Date(hotspots_sf$acq_date)

  hs_buf        <- st_buffer(hs_3005, dist = buffer_m)
  hs_buf$hs_idx <- seq_len(nrow(hs_buf))

  flagged_idx <- st_join(grid_sf, hs_buf["hs_idx"], join = st_within) |>
    filter(!is.na(hs_idx)) |>
    st_drop_geometry() |>
    mutate(
      hotspot_date = hotspot_dates[hs_idx],
      days_diff    = abs(as.numeric(first_date - hotspot_date))
    ) |>
    filter(days_diff <= window_days) |>
    pull(hs_idx)

  result          <- rep(FALSE, nrow(hotspots_sf))
  result[flagged_idx] <- TRUE
  result
}


#  Load qualifying 2024 spring hotspots (from 06a) 
spring_hs_2024 <- st_read(cfg$intermediate$spring_hotspots_2024_in_2023, quiet = TRUE)

#  Buffer 2024 hotspots 
hs2024_buffered <- st_buffer(spring_hs_2024, dist = BUFFER_M)

#  Load 2025 hotspots 
hotspots_2025 <- st_read(d(cfg$data$fires$all_hotspots), quiet = TRUE) |>
  filter(year == 2025) |>
  st_transform(3005) |>
  mutate(
    acq_date = as.Date(acq_date),
    doy      = as.integer(format(acq_date, '%j'))
  ) |>
  select(acq_date, doy, geometry)

#  Extract mean SDD 2025 for 2023 perimeters 
clean_names <- function(df, year) {
  names(df) <- tolower(names(df))
  df <- df |>
    mutate(
      fire_id  = paste0(fire_year, '_', fire_no),
      geo_cat  = paste0('perim_', substr(as.character(year), 3, 4))
    )
  if (!'fire_cause' %in% names(df)) {
    df <- df |> mutate(fire_cause = NA_character_)
  }
  df |> select(fire_id, fire_year, fire_cause, size_ha, geo_cat, area_sqm, geometry)
}

perims_2023 <- st_read(d(cfg$data$perimeters$perims_2023), quiet = TRUE) |>
  rename_with(tolower) |> clean_names(2023) |>
  select(fire_id, geometry) |>
  filter(fire_id %in% unique(spring_hs_2024$fire_id_2023))

sdd_rast_2025 <- rast(d(cfg$data$snow_dd$sdd_raster_2025))
perims_v      <- project(vect(perims_2023), crs(sdd_rast_2025))
sdd_vals      <- terra::extract(sdd_rast_2025, perims_v, fun = mean, na.rm = TRUE)

sdd_2025_lookup <- tibble(
  fire_id_2023 = perims_2023$fire_id,
  sdd_2025     = sdd_vals[, 2]
)

#  Spatial join: 2025 hotspots within buffered 2024 hotspot areas 
hs2024_join <- hs2024_buffered |>
  select(hotspot_id_2024, fire_id_2023, acq_date_2024 = acq_date, geometry)

hs2025_in_buf <- st_join(
  hotspots_2025,
  hs2024_join,
  join = st_within,
  left = FALSE
)


#  SDD filter for 2025: -5 to 60 days after snow melt 
hs2025_filtered <- hs2025_in_buf |>
  left_join(sdd_2025_lookup, by = 'fire_id_2023') |>
  mutate(days_after_sdd_2025 = doy - sdd_2025) |>
  filter(days_after_sdd_2025 >= SDD_MIN, days_after_sdd_2025 <= SDD_MAX) |>
  mutate(hotspot_id_2025 = paste0('hs25_', row_number()))


#  Lightning exclusion for 2025 hotspots 
if (nrow(hs2025_filtered) > 0) {
  lightning_flag <- check_lightning(
    hotspots_sf = hs2025_filtered,
    nc_path     = d(cfg$data$lightning$cg_flashes_2025),
    buffer_m    = LIGHTNING_BUFFER_M,
    window_days = LIGHTNING_WINDOW
  )
  hs2025_filtered <- hs2025_filtered[!lightning_flag, ]
}


#  Save 2025 hotspot locations for mapping 
hs2025_spatial <- hs2025_filtered |>
  rename(acq_date_2025 = acq_date) |>
  select(fire_id_2023, hotspot_id_2024, hotspot_id_2025,
         acq_date_2025, days_after_sdd_2025, geometry)

st_write(hs2025_spatial, cfg$results$spring_hotspots_2025_in_chains, delete_dsn = TRUE, quiet = TRUE)

#  Build and save tabular output 
chains <- hs2025_filtered |>
  st_drop_geometry() |>
  rename(acq_date_2025 = acq_date) |>
  left_join(
    spring_hs_2024 |> st_drop_geometry() |>
      select(hotspot_id_2024, days_after_sdd_2024 = days_after_sdd),
    by = 'hotspot_id_2024'
  ) |>
  mutate(
    days_2024_to_2025 = as.integer(as.Date(acq_date_2025) - as.Date(acq_date_2024))
  ) |>
  select(
    fire_id_2023,
    hotspot_id_2024, acq_date_2024, days_after_sdd_2024,
    hotspot_id_2025, acq_date_2025, days_after_sdd_2025,
    days_2024_to_2025
  )

write_csv(chains, cfg$results$overwintering_perimeter_chains)

# Per-fire summary 
# Multiple hotspots on the same day within a perimeter are evidence of activity
# extent, not independent events. Summarise to one row per fire_id_2023:
#   - first detection date and days-after-SDD for each year (timing)
#   - n_hotspot_days (unique detection dates, a measure of persistence)
#   - n_hotspots (total detections, a measure of intensity/spatial extent)
#   - is_multiyear (TRUE when qualifying hotspots exist in both 2024 and 2025)

summary_2024 <- spring_hs_2024 |>
  st_drop_geometry() |>
  group_by(fire_id_2023) |>
  summarise(
    sdd_2024                  = first(sdd_2024),
    first_acq_date_2024       = min(acq_date),
    days_after_sdd_first_2024 = min(days_after_sdd),
    n_hotspot_days_2024       = n_distinct(acq_date),
    n_hotspots_2024           = n(),
    .groups = 'drop'
  )

summary_2025 <- hs2025_filtered |>
  st_drop_geometry() |>
  group_by(fire_id_2023) |>
  summarise(
    sdd_2025                  = first(sdd_2025),
    first_acq_date_2025       = min(acq_date),
    days_after_sdd_first_2025 = min(days_after_sdd_2025),
    n_hotspot_days_2025       = n_distinct(acq_date),
    n_hotspots_2025           = n(),
    .groups = 'drop'
  )

chains_summary <- summary_2024 |>
  left_join(summary_2025, by = 'fire_id_2023') |>
  mutate(
    is_multiyear            = !is.na(first_acq_date_2025),
    days_first_2024_to_2025 = as.integer(as.Date(first_acq_date_2025) - as.Date(first_acq_date_2024))
  ) |>
  select(
    fire_id_2023,
    sdd_2024, first_acq_date_2024, days_after_sdd_first_2024,
    n_hotspot_days_2024, n_hotspots_2024,
    sdd_2025, first_acq_date_2025, days_after_sdd_first_2025,
    n_hotspot_days_2025, n_hotspots_2025,
    is_multiyear, days_first_2024_to_2025
  )

write_csv(chains_summary, cfg$results$overwintering_perimeter_chains_summary)

