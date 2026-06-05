#  OVERWINTERING FIRE DETECTION ANALYSIS: Spatial and Temporal Proximity Criteria
# Author: Jack A. Goldman
#
# Detects fires that overwintered in northeastern BC boreal
# plains between 2023-2025 using VIIRS FIRMS hotspots and BC fire perimeters.
#
# Algorithm based on Scholten et al. with the following adaptations:
#   - VIIRS FIRMS hotspots used as proxy for ignition points (no large fire
#     database ignition points available for this region/period)
#   - Fall hotspot persistence check added (>= 2 detections within perimeter
#     after Aug 1) to confirm active smoldering before winter dormancy
#   - Per-fire SDD extracted from MODIS raster (vs. regional mean in Scholten)
#   - Fire cause filter replaces Scholten road proximity filter: spring-year
#     perimeters with fire_cause == 'Person' are excluded as human ignitions
#   - Lightning filter: spring hotspots within 7 days of a CG lightning strike
#     at that location are excluded (Scholten uses 6 days + 1 day uncertainty)
#
# Criterion (perimeter-distance version):
#   1. >= 2 fall VIIRS hotspots within previous-year perimeter (DOY >= 213)
#   2. Spring VIIRS hotspot within 1000m buffer of previous-year perimeter,
#      between SDD and SDD + 50 days
#   3. That spring hotspot is within 1000m of a current-year fire perimeter
#   4. Spring-year perimeter fire_cause != 'Person'
#   5. No CG lightning strike within 7 days before spring hotspot detection
#
# Study region: Fort St. John and Fort Nelson Fire Zones, northeastern BC
# 

# --- imports 
import os
import datetime
import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
from shapely.geometry import Point
from pyproj import Transformer
from rasterstats import zonal_stats
import rasterio
import logging

# 
# OUTPUT DIRECTORY
# 

RUN_DATE = datetime.date.today().strftime('%Y%m%d')
OUT_DIR  = "results/results_stp"
os.makedirs(OUT_DIR, exist_ok=True)
print(f"Outputs will be saved to: {OUT_DIR}/")
LOG_PATH = os.path.join(OUT_DIR, "last_run.log")

logging.basicConfig(
    filename=LOG_PATH,
    filemode="w",  # 'w' overwrites the file so it only shows the latest run info
    level=logging.INFO,
    format="%(message)s",  # Keeps the format clean to just show your message
)

# 3. Write the run date into the file
logging.info(f"run_date: {RUN_DATE}")

# 
# CONSTANTS
# 

TRANSFORMER    = Transformer.from_crs('EPSG:4326', 'EPSG:3005', always_xy=True)
BUFFER_DISTANCES   = [500, 1000, 2000]
FALL_DOY_START     = 213          # Aug 1
MAX_PERIM_DIST_M   = 1000         # max distance from spring hotspot to current-year perimeter
SDD_MAX_DAYS       = 50           # spring window: SDD to SDD + 50 days
LIGHTNING_MAX_DAYS = 7            # max lag between lightning strike and spring hotspot

SDD_LOOKUP = {
    'Fort St. John Fire Zone': {2024: 135.71, 2025: 140.63},
    'Fort Nelson Fire Zone':   {2024: 165.14, 2025: 167.72}
}

SDD_RASTERS = {
    2024: '../data/sdd/sdd_2024.tif',
    2025: '../data/sdd/sdd_2025.tif'
}

# 
# STEP 1 — LOAD DATA
# 

# --- fire perimeters (2023-2025) ---
fires = gpd.read_file('data/analysis/all_fires_processed.geojson')
fires['fire_year'] = pd.to_numeric(fires['fire_year'], errors='coerce')

print(f"Fires loaded: {len(fires)}")
print(f"Years: {sorted(fires['fire_year'].unique())}")
print(f"Zones: {fires['fire_zone'].unique()}")
print(f"CRS: {fires.crs}")
print(f"fire_cause values: {fires['fire_cause'].unique()}")

# --- buffered fire perimeters with per-buffer SDD ---
buffered_dfs = {
    500:  gpd.read_file('data/analysis/fires_buffer_500m_with_sdd.geojson'),
    1000: gpd.read_file('data/analysis/fires_buffer_1000m_with_sdd.geojson'),
    2000: gpd.read_file('data/analysis/fires_buffer_2000m_with_sdd.geojson')
}

for dist in BUFFER_DISTANCES:
    sdd_cols = [c for c in buffered_dfs[dist].columns if 'sdd' in c]
    print(f"Buffer {dist}m: {len(buffered_dfs[dist])} fires | SDD cols: {sdd_cols}")

# --- VIIRS FIRMS hotspots ---
# pre-filtered: nighttime only (daynight == 'N'),
#               confidence != low, scan/track <= 0.5
hotspots = gpd.read_file('data/analysis/all_hotspots_2023_2025_night_nh.geojson')
hotspots['datetime'] = pd.to_datetime(hotspots['acq_date'])
hotspots['doy']      = hotspots['datetime'].dt.dayofyear
hotspots['year']     = hotspots['datetime'].dt.year
hotspots = hotspots.to_crs('EPSG:3005')

print(f"\nHotspots loaded: {len(hotspots)}")
print(f"Years: {sorted(hotspots['year'].unique())}")
print(f"CRS: {hotspots.crs}")

# --- fire cause lookup (fire_id -> fire_cause) ---
cause_lookup = fires.set_index('fire_id')['fire_cause'].to_dict()

# 
# STEP 2 — SPLIT HOTSPOTS BY YEAR AND SEASON
# 

def split_hotspots_fall_spring(hotspots_year, year, sdd_lookup,
                                fall_doy_start=FALL_DOY_START):
    """
    Split annual hotspots into fall and spring subsets.

    Fall:   DOY >= fall_doy_start (Aug 1)
    Spring: DOY >= min SDD across zones for that year
            Upper bound (SDD + 50 days) applied per-fire in detect_overwintering
    """
    fall = hotspots_year[hotspots_year['doy'] >= fall_doy_start].copy()

    sdd_values = [
        sdd_lookup[zone][year]
        for zone in sdd_lookup
        if year in sdd_lookup[zone]
    ]

    if sdd_values:
        min_sdd = min(sdd_values)
        spring  = hotspots_year[hotspots_year['doy'] >= min_sdd].copy()
    else:
        print(f"  No SDD found for year {year}, spring subset will be empty")
        spring  = gpd.GeoDataFrame(columns=hotspots_year.columns, crs=hotspots_year.crs)
        min_sdd = 'N/A'

    print(f"Year {year} | fall: {len(fall)} hotspots (DOY >= {fall_doy_start}) "
          f"| spring: {len(spring)} hotspots (DOY >= {min_sdd})")

    return fall, spring


hotspots_2023 = hotspots[hotspots['year'] == 2023].copy()
hotspots_2024 = hotspots[hotspots['year'] == 2024].copy()
hotspots_2025 = hotspots[hotspots['year'] == 2025].copy()

fall_2023, _           = split_hotspots_fall_spring(hotspots_2023, 2023, SDD_LOOKUP)
fall_2024, spring_2024 = split_hotspots_fall_spring(hotspots_2024, 2024, SDD_LOOKUP)
_,         spring_2025 = split_hotspots_fall_spring(hotspots_2025, 2025, SDD_LOOKUP)

# 
# STEP 3 — DETECT OVERWINTERING FIRES (perimeter-distance criterion)
# 
# For each fire perimeter:
#   1. Find fall hotspots within perimeter (>= 2 detections)
#   2. Find spring hotspots within buffer between SDD and SDD + 50 days
#   3. Check spring hotspot is within MAX_PERIM_DIST_M of a current-year perimeter
#   4. Record linked current-year perimeter ID for cause filtering downstream

def detect_overwintering(fires, buffered_dfs, fall_hotspots, spring_hotspots,
                         buffer_dist, fire_year, spring_year,
                         max_perim_dist_m=MAX_PERIM_DIST_M):
    """
    Detect overwintering fires using perimeter proximity criterion.

    Parameters
    ----------
    fires : GeoDataFrame
        Fire perimeters with fire_id, fire_year, fire_zone, fire_cause columns
    buffered_dfs : dict
        Buffered perimeters keyed by buffer distance, with SDD columns
    fall_hotspots : GeoDataFrame
        Hotspots from fall of fire_year (DOY >= 213)
    spring_hotspots : GeoDataFrame
        Hotspots from spring of spring_year (DOY >= min SDD)
    buffer_dist : int
        Buffer distance in metres
    fire_year : int
        Year fire originated (e.g. 2023)
    spring_year : int
        Year of spring reactivation (e.g. 2024)
    max_perim_dist_m : int
        Max distance from spring hotspot to current-year perimeter edge

    Returns
    -------
    list of dicts, one record per fire
    """
    perims         = fires[fires['fire_year'] == float(fire_year)].copy()
    buffers        = buffered_dfs[buffer_dist][
        buffered_dfs[buffer_dist]['fire_year'] == float(fire_year)
    ].copy()
    current_perims = fires[fires['fire_year'] == float(spring_year)].copy()

    results = []

    for _, fire in perims.iterrows():
        fid       = fire['fire_id']
        fire_zone = fire['fire_zone']
        fire_geom = fire['geometry']

        buf_row = buffers[buffers['fire_id'] == fid]
        if buf_row.empty:
            continue
        buf_row  = buf_row.iloc[0]
        buf_geom = buf_row['geometry']

        sdd_col = f'sdd_{spring_year}_mean'
        if sdd_col not in buf_row.index or pd.isna(buf_row[sdd_col]):
            continue
        sdd_doy = buf_row[sdd_col]

        # --- fall hotspots within perimeter (>= 2 required) ---
        fall_in_perim = fall_hotspots[fall_hotspots.within(fire_geom)]
        has_fall      = len(fall_in_perim) >= 2

        # last fall hotspot
        if has_fall:
            last_fall      = fall_in_perim.loc[fall_in_perim['datetime'].idxmax()]
            last_fall_date = last_fall['datetime']
            last_fall_lat  = last_fall['latitude']
            last_fall_lon  = last_fall['longitude']
        else:
            last_fall_date = last_fall_lat = last_fall_lon = None

        # --- spring hotspots within buffer in SDD window ---
        spring_in_window = spring_hotspots[
            (spring_hotspots['doy'] >= sdd_doy) &
            (spring_hotspots['doy'] <= sdd_doy + SDD_MAX_DAYS)
        ]
        spring_in_buffer = spring_in_window[spring_in_window.within(buf_geom)]
        has_spring       = len(spring_in_buffer) > 0

        # --- spring hotspot proximity to current-year perimeter ---
        if has_spring:
            spring_in_buffer = spring_in_buffer.copy()
            spring_in_buffer['dist_to_current_perim'] = spring_in_buffer['geometry'].apply(
                lambda pt: current_perims['geometry'].distance(pt).min()
                if not current_perims.empty else np.inf
            )
            spring_in_buffer['nearest_current_perim_id'] = spring_in_buffer['geometry'].apply(
                lambda pt: current_perims.loc[
                    current_perims['geometry'].distance(pt).idxmin(), 'fire_id'
                ] if not current_perims.empty else None
            )

            spring_near_perim     = spring_in_buffer[
                spring_in_buffer['dist_to_current_perim'] <= max_perim_dist_m
            ]
            has_spring_near_perim = len(spring_near_perim) > 0

            if has_spring_near_perim:
                best_spring            = spring_near_perim.loc[
                    spring_near_perim['dist_to_current_perim'].idxmin()
                ]
                first_spring_date          = best_spring['datetime']
                first_spring_lat           = best_spring['latitude']
                first_spring_lon           = best_spring['longitude']
                first_spring_dist_to_perim = best_spring['dist_to_current_perim']
                linked_perim_id            = best_spring['nearest_current_perim_id']
                dormancy_days              = (
                    pd.to_datetime(first_spring_date) - pd.to_datetime(last_fall_date)
                ).days if last_fall_date is not None else None
                days_after_sdd             = (
                    pd.to_datetime(first_spring_date).dayofyear - sdd_doy
                )
            else:
                first_spring_date = first_spring_lat = first_spring_lon = None
                first_spring_dist_to_perim = linked_perim_id = None
                dormancy_days = days_after_sdd = None

        else:
            has_spring_near_perim      = False
            spring_near_perim          = gpd.GeoDataFrame()
            first_spring_date          = None
            first_spring_lat           = None
            first_spring_lon           = None
            first_spring_dist_to_perim = None
            linked_perim_id            = None
            dormancy_days              = None
            days_after_sdd             = None

        overwintered = has_fall and has_spring_near_perim

        results.append({
            'fire_id':                    fid,
            'fire_year':                  fire_year,
            'fire_zone':                  fire_zone,
            'buffer_dist':                buffer_dist,
            'spring_year':                spring_year,
            'sdd_doy':                    sdd_doy,
            'sdd_max_doy':                sdd_doy + SDD_MAX_DAYS,
            'n_fall_hotspots':            len(fall_in_perim),
            'n_spring_hotspots':          len(spring_in_buffer) if has_spring else 0,
            'n_spring_near_perim':        len(spring_near_perim) if has_spring else 0,
            'has_fall_hotspots':          has_fall,
            'has_spring_hotspots':        has_spring,
            'has_spring_near_perim':      has_spring_near_perim,
            'overwintered':               overwintered,
            'linked_perim_id':            linked_perim_id,
            'first_spring_dist_to_perim': first_spring_dist_to_perim,
            'last_fall_date':             last_fall_date,
            'last_fall_lat':              last_fall_lat,
            'last_fall_lon':              last_fall_lon,
            'first_spring_date':          first_spring_date,
            'first_spring_lat':           first_spring_lat,
            'first_spring_lon':           first_spring_lon,
            'dormancy_days':              dormancy_days,
            'days_after_sdd':             days_after_sdd,
            'exclusion_reason':           None,
            'flag_near_lightning':        None,
            'flag_near_infrastructure':   None,
        })

    return results


# run detection for all buffer distances
all_results = {}

for dist in BUFFER_DISTANCES:
    results_23_24 = detect_overwintering(
        fires, buffered_dfs, fall_2023, spring_2024,
        dist, fire_year=2023, spring_year=2024
    )
    results_24_25 = detect_overwintering(
        fires, buffered_dfs, fall_2024, spring_2025,
        dist, fire_year=2024, spring_year=2025
    )

    df_23_24 = pd.DataFrame(results_23_24)
    df_24_25 = pd.DataFrame(results_24_25)

    all_results[dist] = pd.concat([df_23_24, df_24_25], ignore_index=True)

    print(f"\nBuffer {dist}m (perimeter distance criterion):")
    print(f"  2023->2024 overwintering: "
          f"{df_23_24['overwintered'].sum()} of {len(df_23_24)} fires")
    print(f"  2024->2025 overwintering: "
          f"{df_24_25['overwintered'].sum()} of {len(df_24_25)} fires")

# 
# STEP 4 — LIGHTNING FILTER
# 
# Exclude spring hotspots where a CG lightning strike occurred at or near
# the detection location within LIGHTNING_MAX_DAYS before the hotspot date.

def get_lightning_date_at_point(lat, lon, lightning_dates, lats, lons):
    """
    Return the first lightning strike date at the nearest grid cell to lat/lon.
    Returns None if no lightning recorded at that location.
    """
    lat_idx = np.argmin(np.abs(lats - lat))
    lon_idx = np.argmin(np.abs(lons - lon))
    val     = lightning_dates[lat_idx, lon_idx]
    return None if pd.isnull(val) else pd.Timestamp(val)


def add_lightning_filter(all_results, lightning_by_year, max_days=LIGHTNING_MAX_DAYS):
    """
    Flag overwintering candidates where a lightning strike preceded the
    spring hotspot detection by <= max_days at the same location.

    Parameters
    ----------
    all_results : dict
        Overwintering results keyed by buffer distance
    lightning_by_year : dict
        Keys: year (int). Values: dict with 'dates', 'lats', 'lons' arrays
    max_days : int
        Maximum lag in days between lightning and hotspot (default 7)

    Returns
    -------
    all_results with flag_near_lightning, lightning_date,
    lightning_days_before columns populated
    """
    for dist in BUFFER_DISTANCES:
        df = all_results[dist].copy()

        lightning_flag = []
        lightning_date = []
        lightning_days = []

        for _, row in df.iterrows():
            if not row['overwintered']:
                lightning_flag.append(None)
                lightning_date.append(None)
                lightning_days.append(None)
                continue

            spring_lat  = row['first_spring_lat']
            spring_lon  = row['first_spring_lon']
            spring_date = pd.to_datetime(row['first_spring_date'])
            spring_year = int(row['spring_year'])

            if pd.isna(spring_lat) or pd.isna(spring_lon):
                lightning_flag.append(None)
                lightning_date.append(None)
                lightning_days.append(None)
                continue

            if spring_year not in lightning_by_year:
                lightning_flag.append(None)
                lightning_date.append(None)
                lightning_days.append(None)
                continue

            ld      = lightning_by_year[spring_year]
            lt_date = get_lightning_date_at_point(
                spring_lat, spring_lon, ld['dates'], ld['lats'], ld['lons']
            )

            if lt_date is None:
                lightning_flag.append(False)
                lightning_date.append(None)
                lightning_days.append(None)
                continue

            days_diff = (spring_date - lt_date).days
            flagged   = 0 <= days_diff <= max_days

            lightning_flag.append(flagged)
            lightning_date.append(lt_date)
            lightning_days.append(days_diff)

        df['flag_near_lightning']   = lightning_flag
        df['lightning_date']        = lightning_date
        df['lightning_days_before'] = lightning_days
        all_results[dist]           = df

        ow        = df[df['overwintered'] == True]
        n_flagged = ow['flag_near_lightning'].sum()
        print(f"\nBuffer {dist}m:")
        print(f"  Overwintering candidates: {len(ow)}")
        print(f"  Flagged near lightning:   {n_flagged}")

    return all_results


# load lightning first-occurrence grids
lightning_by_year = {}
for lightning_year in [2024, 2025]:
    ds              = xr.open_dataset(
        f'data/processed_lightning/cg_first_occurrence_{lightning_year}.nc'
    )
    lightning_dates = ds['time'].values
    lats            = ds['lat'].values
    lons            = ds['lon'].values

    print(f"\nLightning {lightning_year}:")
    print(f"  Grid shape: {lightning_dates.shape}")
    print(f"  Lat range:  {lats.min():.2f} to {lats.max():.2f}")
    print(f"  Lon range:  {lons.min():.2f} to {lons.max():.2f}")

    lightning_by_year[lightning_year] = {
        'dates': lightning_dates,
        'lats':  lats,
        'lons':  lons,
    }

all_results = add_lightning_filter(all_results, lightning_by_year)

# apply lightning exclusion
for dist in BUFFER_DISTANCES:
    df           = all_results[dist].copy()
    lightning_mask = df['flag_near_lightning'] == True
    df.loc[lightning_mask, 'overwintered']    = False
    df.loc[lightning_mask, 'exclusion_reason'] = 'lightning_strike'
    all_results[dist] = df

    ow = df[df['overwintered'] == True]
    print(f"\nBuffer {dist}m — after lightning filter: {len(ow)} overwintering fires")

# 
# STEP 5 — FIRE CAUSE FILTER
# 
# Exclude confirmed OW fires where the linked spring-year perimeter has
# fire_cause == 'Person'. This replaces Scholten's road proximity filter
# with direct cause attribution from the BC Wildfire Service database.

for dist in BUFFER_DISTANCES:
    df = all_results[dist].copy()

    df['linked_perim_cause'] = df['linked_perim_id'].map(cause_lookup)

    person_mask = df['linked_perim_cause'] == 'Person'

    df.loc[person_mask, 'overwintered']     = False
    df.loc[person_mask, 'exclusion_reason'] = 'human_caused_perimeter'

    all_results[dist] = df

    ow     = df[df['overwintered'] == True]
    n_excl = person_mask.sum()
    print(f"\nBuffer {dist}m — after fire cause filter:")
    print(f"  Human-caused perimeters excluded: {n_excl}")
    print(f"  Overwintering fires remaining:    {len(ow)}")
    print(ow[['fire_id', 'fire_year', 'spring_year',
              'dormancy_days', 'days_after_sdd',
              'linked_perim_cause']].to_string())

# 
# STEP 6 — BUILD MULTI-YEAR CHAINS (2023 -> 2024 -> 2025)
# 
# A confirmed multi-year chain requires:
#   1. 2023 fire overwintered into 2024 (confirmed in Step 3-5)
#   2. The linked 2024 perimeter also overwintered into 2025 (confirmed in Step 3-5)

for dist in BUFFER_DISTANCES:
    df = all_results[dist]

    ow_23_24 = df[
        (df['fire_year'] == 2023) &
        (df['overwintered'] == True) &
        (df['spring_year'] == 2024)
    ][['fire_id', 'fire_zone', 'linked_perim_id', 'dormancy_days',
       'days_after_sdd', 'first_spring_dist_to_perim']].rename(columns={
        'fire_id':                    'fire_id_2023',
        'linked_perim_id':            'fire_id_2024',
        'dormancy_days':              'dormancy_days_2324',
        'days_after_sdd':             'days_after_sdd_2324',
        'first_spring_dist_to_perim': 'dist_to_perim_2324'
    })

    ow_24_25 = df[
        (df['fire_year'] == 2024) &
        (df['overwintered'] == True) &
        (df['spring_year'] == 2025)
    ][['fire_id', 'linked_perim_id', 'dormancy_days',
       'days_after_sdd', 'first_spring_dist_to_perim']].rename(columns={
        'fire_id':                    'fire_id_2024',
        'linked_perim_id':            'fire_id_2025',
        'dormancy_days':              'dormancy_days_2425',
        'days_after_sdd':             'days_after_sdd_2425',
        'first_spring_dist_to_perim': 'dist_to_perim_2425'
    })

    multiyear = ow_23_24.merge(ow_24_25, on='fire_id_2024', how='inner')

    print(f"\nBuffer {dist}m — multi-year chains after all filters: {len(multiyear)}")
    if not multiyear.empty:
        print(multiyear[[
            'fire_id_2023', 'fire_id_2024', 'fire_id_2025',
            'dormancy_days_2324', 'dormancy_days_2425',
            'days_after_sdd_2324', 'days_after_sdd_2425',
            'dist_to_perim_2324',  'dist_to_perim_2425'
        ]].to_string())

    # save multi-year chains for 1000m buffer
    if dist == 1000 and not multiyear.empty:
        multiyear.to_csv(
            os.path.join(OUT_DIR, 'multiyear_chains_1000m.csv'), index=False
        )
        print(f"Saved multiyear_chains_1000m.csv")

# 
# STEP 7 — SINGLE-YEAR CLASSIFICATION
# 
# Fires that overwintered once only — not part of a multi-year chain.

multiyear_1000 = pd.read_csv(os.path.join(OUT_DIR, 'multiyear_chains_1000m.csv')) \
    if os.path.exists(os.path.join(OUT_DIR, 'multiyear_chains_1000m.csv')) \
    else pd.DataFrame(columns=['fire_id_2023', 'fire_id_2024'])

multiyear_ids_2023 = set(multiyear_1000['fire_id_2023']) if not multiyear_1000.empty else set()
multiyear_ids_2024 = set(multiyear_1000['fire_id_2024']) if not multiyear_1000.empty else set()

for dist in BUFFER_DISTANCES:
    df = all_results[dist]

    ow_23_24 = df[(df['fire_year'] == 2023) & (df['overwintered'] == True) & (df['spring_year'] == 2024)]
    ow_24_25 = df[(df['fire_year'] == 2024) & (df['overwintered'] == True) & (df['spring_year'] == 2025)]

    single_23_24 = ow_23_24[~ow_23_24['fire_id'].isin(multiyear_ids_2023)]
    single_24_25 = ow_24_25[~ow_24_25['fire_id'].isin(multiyear_ids_2024)]

    print(f"\nBuffer {dist}m:")
    print(f"  2023->2024 total: {len(ow_23_24)} | single-year only: {len(single_23_24)}")
    print(f"  2024->2025 total: {len(ow_24_25)} | single-year only: {len(single_24_25)}")

# 
# STEP 8 — SAVE ALL RESULTS
# 

for dist in BUFFER_DISTANCES:
    all_results[dist].to_csv(
        os.path.join(OUT_DIR, f"overwintering_results_{dist}m.csv"), index=False
    )
    print(f"Saved overwintering_results_{dist}m.csv")

print("\nAll results saved.")
