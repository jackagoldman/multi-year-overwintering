library(sf)
library(dplyr)
library(readr)

INT_DIR <- 'results/'
OUT_DIR <- 'results/'
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Load ignition points ───────────────────────────────────────────────────────
ignitions_2024 <- st_read(file.path(INT_DIR, 'ignitions_2024.geojson'), quiet = TRUE)
ignitions_2025 <- st_read(file.path(INT_DIR, 'ignitions_2025.geojson'), quiet = TRUE)


# ── Spatial overlap: 2023->2024 ────────────────────────────────────────────────
# A 2024 fire is a candidate overwinter if its ignition point
# is within threshold distance of a confirmed 2023 perimeter

DIST_THRESH_M <- 1000

candidates_23_24 <- ignitions_2024 |>
  filter(dist_to_prev_perim_m <= DIST_THRESH_M)

# ── Spatial overlap: 2024->2025 ────────────────────────────────────────────────
candidates_24_25 <- ignitions_2025 |>
  filter(dist_to_prev_perim_m <= DIST_THRESH_M)

# ── Multi-year chains: 2023->2024->2025 ───────────────────────────────────────
# A chain exists when the 2024 fire that a 2023 fire reactivated into
# also appears as a confirmed 2024->2025 candidate.
# The link is: candidates_23_24$fire_id == candidates_24_25$nearest_prev_fire_id
# the 2024 fire must itself be a confirmed smoldering fire
# i.e. it must appear in perims_2024_confirmed
multiyear <- candidates_23_24 |>
  st_drop_geometry() |>
  rename(
    fire_id_2023         = nearest_prev_fire_id,
    fire_id_2024         = fire_id,
    acq_date_2024        = acq_date,
    ignition_doy_2024    = ignition_doy,
    sdd_doy_2024         = sdd_doy,
    days_after_sdd_2324  = days_after_sdd,
    n_hotspots_2024      = n_first_day_hotspots,
    dist_to_2023_perim_m = dist_to_prev_perim_m
  ) |>
  # 2024 fire must have confirmed fall smoldering
  filter(fire_id_2024 %in% perims_2024_confirmed$fire_id) |>
  inner_join(
    candidates_24_25 |>
      st_drop_geometry() |>
      rename(
        fire_id_2024         = nearest_prev_fire_id,
        fire_id_2025         = fire_id,
        acq_date_2025        = acq_date,
        ignition_doy_2025    = ignition_doy,
        sdd_doy_2025         = sdd_doy,
        days_after_sdd_2425  = days_after_sdd,
        n_hotspots_2025      = n_first_day_hotspots,
        dist_to_2024_perim_m = dist_to_prev_perim_m
      ),
    by = 'fire_id_2024'
  ) |> select(
    fire_id_2023, fire_id_2024,fire_id_2025, 
    acq_date_2024, acq_date_2025, 
    ignition_doy_2024,ignition_doy_2025,
    sdd_doy_2024, sdd_doy_2025,
    days_after_sdd_2324, days_after_sdd_2425,
    dist_to_2023_perim_m,
    dist_to_2024_perim_m
  ) 
# ── Single-year ────────────────────────────────────────────────────────────────
multiyear_ids <- unique(multiyear$fire_id_2024)

single_23_24 <- candidates_23_24 |>
  filter(!fire_id %in% multiyear_ids)

single_24_25 <- candidates_24_25 |>
  filter(!nearest_prev_fire_id %in% multiyear_ids)

# ── Save ───────────────────────────────────────────────────────────────────────
write_csv(multiyear,    file.path(OUT_DIR, 'multiyear_overwinter.csv'))
st_write(single_23_24, file.path(OUT_DIR, 'single_23_24_overwinter.geojson'), delete_dsn = TRUE)
st_write(single_24_25, file.path(OUT_DIR, 'single_24_25_overwinter.geojson'), delete_dsn = TRUE)
