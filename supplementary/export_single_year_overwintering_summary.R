library(sf)
library(readr)
library(dplyr)
library(yaml)

cfg <- yaml.load_file('config.yml')

# New fire ID: 05_overwinter.R already excludes 3-year-chain fires from these.
n_new_23_24 <- nrow(st_read(cfg$results$single_23_24_overwinter, quiet = TRUE))
n_new_24_25 <- nrow(st_read(cfg$results$single_24_25_overwinter, quiet = TRUE))

# In-place 2023-24: all qualifying reignitions, minus the ones that continue
# into a 2023->2024->2025 three-year chain.
inplace_2024 <- st_read(cfg$intermediate$spring_hotspots_2024_in_2023, quiet = TRUE)
chains_summary <- read_csv(cfg$results$overwintering_perimeter_chains_summary, show_col_types = FALSE)

n_inplace_2024_total <- n_distinct(inplace_2024$fire_id_2023)
n_inplace_multiyear  <- sum(chains_summary$is_multiyear, na.rm = TRUE)
n_inplace_23_24      <- n_inplace_2024_total - n_inplace_multiyear

# In-place 2024-25: independent census (see detect_inplace_2025_in_2024.R).

inplace_2025_path <- cfg$intermediate$spring_hotspots_2025_in_2024
if (!file.exists(inplace_2025_path)) {
  stop(
    "Missing ", inplace_2025_path,
    " — run supplementary/detect_inplace_2025_in_2024.R first (requires the mounted data drive)."
  )
}
inplace_2025          <- st_read(inplace_2025_path, quiet = TRUE)
multiyear_overwinter  <- read_csv(cfg$results$multiyear_overwinter, show_col_types = FALSE)
inplace_2025_ids      <- unique(inplace_2025$fire_id_2024)
already_chained_ids   <- intersect(inplace_2025_ids, multiyear_overwinter$fire_id_2024)
n_inplace_24_25       <- length(setdiff(inplace_2025_ids, already_chained_ids))

if (length(already_chained_ids) > 0) {
  message(
    "Excluded from in-place 2024-25 (already in the new-fire-ID 3-year chain): ",
    paste(already_chained_ids, collapse = ", ")
  )
}

summary_table <- tibble(
  Transition = c('2023-24', '2023-24', '2024-25', '2024-25'),
  Category   = c('New fire ID', 'In-place', 'New fire ID', 'In-place'),
  `N single-year overwintering fires` = c(n_new_23_24, n_inplace_23_24, n_new_24_25, n_inplace_24_25)
)

out_path <- 'supplementary/single_year_overwintering_summary.csv'
write_csv(summary_table, out_path)

print(summary_table)
