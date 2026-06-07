library(sf)
library(tidyverse)

OUT_DIR <- 'workflow/results/overwintering_2'
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ──────────────────────────────────────────────────────────────────
sdd_2024 <- read_csv('workflow/results/intermediate_products/sdd_2024.csv', show_col_types = FALSE)
sdd_2025 <- read_csv('workflow/results/intermediate_products/sdd_2025.csv', show_col_types = FALSE)

perims_2023 <- st_read('data/processed_fire_perimeters/perims_2023_clipped_3005.shp', quiet = TRUE)
perims_2024 <- st_read('data/processed_fire_perimeters/perims_2024_clipped_3005.shp', quiet = TRUE)
perims_2025 <- st_read('data/processed_fire_perimeters/perims_2025_clipped_3005.shp', quiet = TRUE)

hotspots_all <- st_read('data/raw/all_hotspots_2023_2025.geojson', quiet = TRUE) |>
  st_transform(3005) |>
  mutate(acq_date = as.Date(acq_date)) # non-conservative hotspots.... would be data/analysis/all_hotspots_2023_2025_night_nh.geojson better? but that would require re-running the fall hotspot candidate analysis to get the clipped version. For now, just load all and filter to fall later.

hotspots_2023 <- hotspots_all |> filter(year == 2023)
hotspots_2024 <- hotspots_all |> filter(year == 2024)
hotspots_2025 <- hotspots_all |> filter(year == 2025)

# ── Clean perimeter names ──────────────────────────────────────────────────────
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

perims_2023 <- clean_names(perims_2023, 2023)
perims_2024 <- clean_names(perims_2024, 2024)
perims_2025 <- clean_names(perims_2025, 2025)

# join SDD to perimeters
perims_2024 <- perims_2024 |> left_join(sdd_2024, by = 'fire_id')
perims_2025 <- perims_2025 |> left_join(sdd_2025, by = 'fire_id')

# ── Build fall 2023/2024 last hotspots ──────────────────────────────────────────────
# Last detected fall hotspot per 2023 fire (DOY >= 213, Aug 1 onward).
# This is where the fire was confirmed smoldering before winter dormancy.
# Used as the spatial reference for finding overwintering 2024 ignitions.

get_fall_confirmed <- function(hotspots, perims, min_hotspots = 2) {

  # all fall hotspots (DOY >= 213) inside each perimeter
  hs_fall <- st_join(
    hotspots |>
      mutate(doy = as.integer(format(acq_date, '%j'))) |>
      filter(doy >= 213) |>
      select(acq_date, doy, geometry),
    perims |> select(fire_id, geometry),
    join = st_within,
    left = FALSE
  )

  # count per fire — used for both the filter and the output
  counts <- hs_fall |>
    st_drop_geometry() |>
    count(fire_id, name = 'n_fall_hotspots')

  # confirmed fire ids: >= min_hotspots detections
  confirmed_ids <- counts |>
    filter(n_fall_hotspots >= min_hotspots) |>
    pull(fire_id)

  # last-day hotspots for confirmed fires only
  last_hotspots <- hs_fall |>
    filter(fire_id %in% confirmed_ids) |>
    group_by(fire_id) |>
    filter(acq_date == max(acq_date)) |>
    ungroup() |>
    left_join(counts, by = 'fire_id')

  # confirmed perimeters
  confirmed_perims <- perims |>
    filter(fire_id %in% confirmed_ids)

  return(list(
    last_hotspots     = last_hotspots,
    confirmed_perims  = confirmed_perims,
    counts            = counts
  ))
}

# run for 2023 and 2024
fall_2023 <- get_fall_confirmed(hotspots_2023, perims_2023)
fall_2024 <- get_fall_confirmed(hotspots_2024, perims_2024)

# extract what you need
fall_2023_last       <- fall_2023$last_hotspots
perims_2023_confirmed <- fall_2023$confirmed_perims

fall_2024_last       <- fall_2024$last_hotspots
perims_2024_confirmed <- fall_2024$confirmed_perims


# ── Ignition point function ────────────────────────────────────────────────────
#
#' Get ignition points for fires in a given year.
#'
#' For each fire perimeter, finds all hotspots on the earliest detected date.
#' Among those first-day hotspots, selects the single point closest to the
#' last confirmed fall hotspot from the previous year. This anchors the
#' ignition point spatially to where overwintering smoldering was last seen,
#' rather than taking a centroid that could fall anywhere in a large fire.
#'
#' @param perims     sf — current-year fire perimeters with sdd column joined
#' @param hotspots   sf — current-year VIIRS hotspots (already filtered to year)
#' @param fall_last  sf — previous-year last fall hotspot per fire
#' @param sdd_col    column name (unquoted) of the SDD DOY column in perims
#'
#' @return sf with one row per fire:
#'   fire_id, fire_year, acq_date, ignition_doy, sdd_doy, days_after_sdd,
#'   n_first_day_hotspots, dist_to_fall_m, nearest_fall_fire_id, geometry


get_ignition_points <- function(perims, hotspots, prev_perims, sdd_col) {

  sdd_lookup <- perims |>
    st_drop_geometry() |>
    select(fire_id, sdd_doy = {{ sdd_col }})

  # Step 1: hotspots inside each current-year perimeter
  hs_in_perims <- st_join(
    hotspots |> select(acq_date, geometry),
    perims   |> select(fire_id, fire_year, geometry),
    join = st_within,
    left = FALSE
  ) |>
    mutate(acq_date = as.Date(acq_date))

  # Step 2: first-day hotspots only
  earliest <- hs_in_perims |>
    st_drop_geometry() |>
    group_by(fire_id) |>
    summarise(earliest_date = min(acq_date), .groups = 'drop')

  first_day_hs <- hs_in_perims |>
    inner_join(earliest, by = 'fire_id') |>
    filter(acq_date == earliest_date)

  # Step 3: for each first-day hotspot, find distance to nearest
  # previous-year perimeter boundary
  prev_geoms <- prev_perims |>
    select(prev_fire_id = fire_id, geometry)

  ignition_pts <- first_day_hs |>
    rowwise() |>
    mutate(
      dists_to_prev_perim  = list(as.numeric(st_distance(geometry, prev_geoms$geometry))),
      nearest_prev_idx     = which.min(dists_to_prev_perim),
      dist_to_prev_perim_m = dists_to_prev_perim[[nearest_prev_idx]],
      nearest_prev_fire_id = prev_geoms$prev_fire_id[nearest_prev_idx]
    ) |>
    ungroup() |>
    # Step 4: per fire, keep the first-day hotspot closest to any previous-year perimeter
    group_by(fire_id) |>
    slice_min(dist_to_prev_perim_m, n = 1, with_ties = FALSE) |>
    ungroup() |>
    left_join(sdd_lookup, by = 'fire_id') |>
    left_join(
      first_day_hs |>
        st_drop_geometry() |>
        count(fire_id, name = 'n_first_day_hotspots'),
      by = 'fire_id'
    ) |>
    mutate(
      ignition_doy   = as.integer(format(acq_date, '%j')),
      days_after_sdd = ignition_doy - sdd_doy
    ) |>
    select(fire_id, fire_year, acq_date, ignition_doy, sdd_doy,
           days_after_sdd, n_first_day_hotspots,
           dist_to_prev_perim_m, nearest_prev_fire_id,
           geometry) |>
    st_as_sf(crs = 3005)

  return(ignition_pts)
}

# now pass confirmed perimeters as prev_perims
ignitions_2024 <- get_ignition_points(
  perims      = perims_2024,
  hotspots    = hotspots_2024,
  prev_perims = perims_2023_confirmed,   # only confirmed smoldering 2023 fires
  sdd_col     = sdd_2024_mean
) |>
  filter(days_after_sdd >= 0, days_after_sdd <= 50)

ignitions_2025 <- get_ignition_points(
  perims      = perims_2025,
  hotspots    = hotspots_2025,
  prev_perims = perims_2024_confirmed,   # only confirmed smoldering 2024 fires
  sdd_col     = sdd_2025_mean
) |>
  filter(days_after_sdd >= 0, days_after_sdd <= 50)

# ── Save ───────────────────────────────────────────────────────────────────────
st_write(ignitions_2024, file.path(OUT_DIR, 'ignitions_2024.geojson'), delete_dsn = TRUE, quiet = TRUE)
st_write(ignitions_2025, file.path(OUT_DIR, 'ignitions_2025.geojson'), delete_dsn = TRUE, quiet = TRUE)
