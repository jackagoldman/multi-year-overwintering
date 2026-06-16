library(sf)
library(tidyverse)
library(terra)
library(ncdf4)
library(yaml)

cfg       <- yaml.load_file('config.yml')
data_root <- cfg$environments[[cfg$environment]]$data_root
d         <- function(subpath) file.path(data_root, subpath)

INT_DIR <- cfg$intermediate$dir
dir.create(INT_DIR, showWarnings = FALSE, recursive = TRUE)

SDD_MIN            <- cfg$params$sdd_window_min
SDD_MAX            <- cfg$params$sdd_window_max
LIGHTNING_WINDOW   <- cfg$params$lightning_window_days
LIGHTNING_BUFFER_M <- cfg$params$lightning_buffer_m


# ── Lightning exclusion helper ─────────────────────────────────────────────────
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

  # First time step with flashes > 0 per cell
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

  # Buffer all hotspots at once; join grid cells into buffers
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


# ── Load 2023 perimeters ───────────────────────────────────────────────────────
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
  select(fire_id, geometry)

# ── Extract mean SDD 2024 per 2023 perimeter ──────────────────────────────────
sdd_rast  <- rast(d(cfg$data$snow_dd$sdd_raster_2024))
perims_v  <- project(vect(perims_2023), crs(sdd_rast))
sdd_vals  <- terra::extract(sdd_rast, perims_v, fun = mean, na.rm = TRUE)
perims_2023$sdd_2024 <- sdd_vals[, 2]

# ── Load 2024 hotspots ─────────────────────────────────────────────────────────
hotspots_2024 <- st_read(d(cfg$data$fires$all_hotspots), quiet = TRUE) |>
  filter(year == 2024) |>
  st_transform(3005) |>
  mutate(
    acq_date = as.Date(acq_date),
    doy      = as.integer(format(acq_date, '%j'))
  ) |>
  select(acq_date, doy, geometry)

# ── Spatial join: 2024 hotspots within 2023 perimeters ────────────────────────
hs_in_perims <- st_join(
  hotspots_2024,
  perims_2023 |> select(fire_id_2023 = fire_id, sdd_2024, geometry),
  join = st_within,
  left = FALSE
)

cat("2024 hotspots inside 2023 perimeters:", nrow(hs_in_perims), "\n")

# ── SDD filter: -5 to 60 days after snow melt ─────────────────────────────────
hs_spring <- hs_in_perims |>
  mutate(days_after_sdd = doy - sdd_2024) |>
  filter(days_after_sdd >= SDD_MIN, days_after_sdd <= SDD_MAX) |>
  mutate(hotspot_id_2024 = paste0('hs24_', row_number()))

cat("After SDD filter (days_after_sdd", SDD_MIN, "to", SDD_MAX, "):", nrow(hs_spring), "\n")

# ── Lightning exclusion ────────────────────────────────────────────────────────
if (nrow(hs_spring) > 0) {
  cat("Running lightning check...\n")
  lightning_flag <- check_lightning(
    hotspots_sf = hs_spring,
    nc_path     = d(cfg$data$lightning$cg_flashes_2024),
    buffer_m    = LIGHTNING_BUFFER_M,
    window_days = LIGHTNING_WINDOW
  )
  cat("Lightning-flagged (excluded):", sum(lightning_flag), "\n")
  hs_spring <- hs_spring[!lightning_flag, ]
}

cat("Qualifying spring 2024 hotspots:", nrow(hs_spring), "\n")

# ── Save ───────────────────────────────────────────────────────────────────────
hs_out <- hs_spring |>
  select(hotspot_id_2024, fire_id_2023, acq_date, doy, sdd_2024, days_after_sdd, geometry)

st_write(hs_out, cfg$intermediate$spring_hotspots_2024_in_2023, delete_dsn = TRUE, quiet = TRUE)
cat("Saved:", cfg$intermediate$spring_hotspots_2024_in_2023, "\n")
