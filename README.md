# Multi-year overwintering fires search algorithm and analyses

This repository contains the code and documentation for the multi-year overwintering fire detection algorithm and analyses presented in Wilkinson et al.

- Algorithm description: [`detection_algorithm/description.md`](detection_algorithm/description.md)
- Datasets used: [`supplementary/datasets_used.md`](supplementary/datasets_used.md)
- Overwintering detection criteria: [`supplementary/overwintering_criteria.md`](supplementary/overwintering_criteria.md)

## Repository structure

| Folder                   | Contents                                                                                    |
| ------------------------ | ------------------------------------------------------------------------------------------- |
| `prep/`                | Data preparation run before the algorithm (snow disappearance date rasters computed in GEE) |
| `detection_algorithm/` | The detection algorithm, steps 01–06, run in order                                         |
| `results/`             | Algorithm outputs; `results/intermediate_products/` holds files passed between steps      |
| `supplementary/`       | Summary tables, sensitivity results, dataset and criteria docs                              |
| `visualizations/`      | Manuscript maps (`maps.qmd`)                                                              |
| `air_quality/`         | Seasonal PM2.5 analysis for the Fort Nelson region                                          |

## Setup

### 1. Software

Python packages are managed with [uv](https://docs.astral.sh/uv/) (`pyproject.toml`) and R packages with uvr (`uvr.toml`). The R version is pinned in `.r-version`.

```bash
uv sync     # Python environment in .venv
uvr sync    # R library in .uvr/library
```

The SDD step also needs a Google Earth Engine account (`earthengine-api`).

### 2. Configure paths in `config.yml`

All scripts read file paths and parameters from `config.yml`, so no paths are hard-coded in the scripts.

- Set `environment` to `home` (data in `./data`) or `remote` (data on the external drive), and check that `data_root` points to your data folder.
- Data paths under `data:` are relative to `data_root`. Do not start them with `/`, or `data_root` is ignored.
- Algorithm parameters (buffer distances, SDD window, lightning window) are under `params:`. Change them here for sensitivity analyses.

### 3. Input data

Download the datasets listed in `supplementary/datasets_used.md` and place them under `data_root` at the paths given in `config.yml`:

| Dataset                                                              | `config.yml` key                                      |
| -------------------------------------------------------------------- | ------------------------------------------------------- |
| Study area boundary (fire zones)                                     | `data.plotting.study_zones`                           |
| BC fire perimeters 2023–2025, clipped to the study zones, EPSG:3005 | `data.perimeters`, `data.fires.all_fires_processed` |
| VIIRS nighttime hotspots 2023–2025                                  | `data.fires.all_hotspots`                             |
| CLDN cloud-to-ground lightning 2024–2025 (NetCDF)                   | `data.lightning`                                      |
| Road network clipped to the study zones, EPSG:3005                   | `data.roads.roads_study_zones`                        |

## Running the workflow

Run every script from the project root, since scripts open `config.yml` relative to the working directory:

```bash
uv run python detection_algorithm/01_sdd.py
uvr run detection_algorithm/02_hotspotsb.R
```

### Prep

| Step                    | Script                 | Requires                                                 | Produces                                                                                                |
| ----------------------- | ---------------------- | -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Snow disappearance date | `prep/sdd_acq.ipynb` | Earth Engine account and project ID; study area boundary | Annual MODIS SDD rasters, exported to Google Drive (`SDDperYear`). Download them to `data.snow_dd`. |
|                         |                        |                                                          |                                                                                                         |

### Detection algorithm

Steps 01–05 detect overwintering fires from new ignitions near the previous year's perimeter (adapted from Scholten et al. 2021). Steps 06a–06b detect in-place reactivation within the original perimeter.

| Step | Script                         | Requires                                                   | Produces                                                                                                                                                                           |
| ---- | ------------------------------ | ---------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 01   | `01_sdd.py`                  | Fire perimeters; SDD rasters 2024–2025                    | `intermediate_products/sdd_2024.csv`, `sdd_2025.csv` (mean SDD per fire)                                                                                                       |
| 02   | `02_hotspotsb.R`             | Perimeters 2023–2025; hotspots; output of 01              | `ignitions_2024.geojson`, `ignitions_2025.geojson`; `intermediate_products/perims_2023_confirmed.geojson`, `perims_2024_confirmed.geojson` (perimeters with fall activity) |
| 03   | `03_lightning.py`            | Output of 02; lightning NetCDF                             | `ignitions_2024_lightning.csv`, `ignitions_2025_lightning.csv` (lightning flag per ignition)                                                                                   |
| 04   | `04_roads.py`                | `ignitions_2025.geojson`; road network                   | `ignitions_2025_roads.csv` (road distance as a proxy for human cause, 2025 only)                                                                                                 |
| 05   | `05_overwinter.R`            | Output of 02                                               | `multiyear_overwinter.csv`, `single_23_24_overwinter.geojson`, `single_24_25_overwinter.geojson`                                                                             |
| 06a  | `06a_spring_hotspots.R`      | 2023 perimeters; hotspots; SDD raster 2024; lightning 2024 | `intermediate_products/spring_hotspots_2024_in_2023_perims.geojson`                                                                                                              |
| 06b  | `06b_overwintering_chains.R` | Output of 06a; hotspots; SDD raster 2025; lightning 2025   | `inplace_overwintering_perimeter_chains.csv`, `..._summary.csv`, `inplace_spring_hotspots_2025.geojson`                                                                      |

All outputs are written to `results/` unless marked `intermediate_products/`.

### Supplementary and figures

After the algorithm has run:

- `supplementary/` scripts build the summary tables and the KML export from the files in `results/`.
- `visualizations/maps.qmd` makes the manuscript maps. It is rendered from inside `visualizations/` and reads `../config.yml`.
- `air_quality/seasonal_pm25_2019_2025.R` runs the PM2.5 analysis, separate from the algorithm.

## Notes

- **Sensitivity analysis:** change `params.buffer_m` and rerun from step 02. Previous runs at 500 m and 2000 m are in `supplementary/sensitivity_results/`.
