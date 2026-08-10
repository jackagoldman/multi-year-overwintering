library(readr)
library(dplyr)
library(yaml)

cfg <- yaml.load_file('config.yml')

summary <- read_csv(cfg$results$overwintering_perimeter_chains_summary, show_col_types = FALSE)

table_out <- summary |>
  filter(is_multiyear == TRUE) |>
  transmute(
    `Fire ID (2023)`               = fire_id_2023,
    `First detection 2024`         = as.character(first_acq_date_2024),
    `Days after SDD (2024)`        = round(days_after_sdd_first_2024, 1),
    `Detection days 2024`          = n_hotspot_days_2024,
    `Hotspot count 2024`           = n_hotspots_2024,
    `First detection 2025`         = as.character(first_acq_date_2025),
    `Days after SDD (2025)`        = round(days_after_sdd_first_2025, 1),
    `Detection days 2025`          = n_hotspot_days_2025,
    `Hotspot count 2025`           = n_hotspots_2025,
    `Days between first detections` = days_first_2024_to_2025
  )

out_path <- 'supplementary/inplace_multiyear_fires_table.csv'
write_csv(table_out, out_path)
