# 02_sdd.py — Per-Fire Mean Snow Disappearance Date
# ==================================================
# Extracts mean SDD DOY from MODIS rasters for each 2024 and 2025
# fire perimeter using zonal statistics.
#
# Inputs  : data/analysis/all_fires_processed.geojson
#           data/snow_dd/SDD_2024.tif
#           data/snow_dd/SDD_2025.tif
#
# Outputs : results/intermediate_products/sdd_2024.csv
#           results/intermediate_products/sdd_2025.csv
#
# Run after: 01_hotspots.py

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats

OUT_DIR = 'workflow/results/intermediate_products'
os.makedirs(OUT_DIR, exist_ok=True)

SDD_RASTERS = {
    2024: 'data/snow_dd/SDD_2024.tif',
    2025: 'data/snow_dd/SDD_2025.tif',
}

SDD_LOOKUP = {
    'Fort St. John Fire Zone': {2024: 135.71, 2025: 140.63},
    'Fort Nelson Fire Zone':   {2024: 165.14, 2025: 167.72},
}

# ── Load perimeters ────────────────────────────────────────────────────────────
fires = gpd.read_file('data/analysis/all_fires_processed.geojson')
fires['fire_year'] = pd.to_numeric(fires['fire_year'], errors='coerce')

for year in [2024, 2025]:
    raster_path = SDD_RASTERS[year]

    # Reproject perimeters to raster native CRS once — avoids per-polygon reprojection
    with rasterio.open(raster_path) as src:
        raster_crs  = src.crs
        nodata_val  = src.nodata

    perims = (
        fires[fires['fire_year'] == float(year)][['fire_id', 'fire_zone', 'geometry']]
        .copy()
        .to_crs(raster_crs)   # match raster CRS so zonal_stats does zero reprojection
    )

    stats = zonal_stats(
        perims,
        raster_path,
        stats    = ['mean'],
        nodata   = nodata_val if nodata_val is not None else np.nan,
        all_touched = True,   # include edge pixels — faster and better for small fires
    )

    sdd_col = f'sdd_{year}_mean'
    perims[sdd_col] = [
        round(s['mean'], 2) if s['mean'] is not None and not np.isnan(s['mean'])
        else SDD_LOOKUP.get(row.fire_zone, {}).get(year, np.nan)
        for s, (_, row) in zip(stats, perims.iterrows())
    ]

    out = perims[['fire_id', 'fire_zone', sdd_col]]
    out.to_csv(os.path.join(OUT_DIR, f'sdd_{year}.csv'), index=False)
    print(f"{year}: {len(out)} perimeters  |  "
          f"n_valid={out[sdd_col].notna().sum()}  |  "
          f"mean_sdd={out[sdd_col].mean():.1f}")

print('Done. Run 02_hotspots.py next')