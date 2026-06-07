# 03_lightning.py — adapted for ignition points
# change SPRING_FILES to read ignition points
import os
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import Point


IGNITION_FILES = {
    2024: 'workflow/results/overwintering_2/ignitions_2024.geojson',
    2025: 'workflow/results/overwintering_2/ignitions_2025.geojson',
}

LIGHTNING_FILES = {
    2024: 'data/raw/lightning/cg_flashes_3hr_0.1-deg_2024.nc',
    2025: 'data/raw/lightning/cg_flashes_3hr_0.1-deg_2025.nc',
}

LIGHTNING_MAX_DAYS = 7   # days before OR after ignition

LIGHTNING_BUFFER_M = 2000   # metres — projected buffer radius

def load_lightning_grid(year):
    """
    Load raw CG flash data and compute first-occurrence date per grid cell.

    Returns
    -------
    first_occ : xr.DataArray (lat, lon) — datetime64, NaT where no strike
    lats      : np.ndarray 1-D
    lons      : np.ndarray 1-D
    """
    ds  = xr.open_dataset(LIGHTNING_FILES[year], engine='netcdf4')
    da  = ds['cg_flashes']   # (time, lat, lon)

    # First time step with cg_flashes > 0 per cell
    mask     = da > 0
    has_any  = mask.any(dim='time')
    idx      = mask.argmax(dim='time')           # index of first True per cell
    first_occ = da['time'].isel(time=idx).where(has_any)   # NaT where no strike

    return first_occ, da['lat'].values, da['lon'].values


def first_lightning_in_buffer(lat_hs, lon_hs, first_occ, lats, lons):
    """
    Find the earliest CG lightning first-occurrence date within a 2 km
    buffer of the hotspot location.

    Steps
    -----
    1. Build 2 km buffer around hotspot (EPSG:3005 → WGS84 for grid matching)
    2. Identify grid cell centres that fall inside the buffer
    3. Return the minimum (earliest) first-occurrence date among those cells

    Parameters
    ----------
    lat_hs, lon_hs : float — hotspot coords in WGS84
    first_occ      : xr.DataArray (lat, lon) — first strike datetime per cell
    lats, lons     : np.ndarray — 1-D grid coordinate arrays

    Returns
    -------
    pd.Timestamp | None
    """
    # Build 2 km buffer: project point to EPSG:3005, buffer, reproject to WGS84
    pt_proj   = gpd.GeoSeries([Point(lon_hs, lat_hs)], crs='EPSG:4326').to_crs('EPSG:3005')
    buf_proj  = pt_proj.buffer(LIGHTNING_BUFFER_M)
    buf_wgs84 = buf_proj.to_crs('EPSG:4326').iloc[0]

    # Find grid cells whose centres fall within the buffer
    # Use bounding box first for speed, then exact within check
    minx, miny, maxx, maxy = buf_wgs84.bounds
    lon_mask = (lons >= minx) & (lons <= maxx)
    lat_mask = (lats >= miny) & (lats <= maxy)

    if not lon_mask.any() or not lat_mask.any():
        return None

    # Subset first_occ to bounding box
    subset = first_occ.sel(lat=lats[lat_mask], lon=lons[lon_mask])

    # Check each candidate cell centre against the actual buffer polygon
    candidate_dates = []
    for lat_val in lats[lat_mask]:
        for lon_val in lons[lon_mask]:
            if not buf_wgs84.contains(Point(lon_val, lat_val)):
                continue
            val = subset.sel(lat=lat_val, lon=lon_val).values
            if not pd.isnull(val):
                candidate_dates.append(pd.Timestamp(val))

    if not candidate_dates:
        return None

    return min(candidate_dates)

def add_lightning_flag(df, year):
    first_occ, lats, lons = load_lightning_grid(year)

    lt_dates, lt_days, lt_flags = [], [], []

    for _, row in df.iterrows():
        # use ignition point coords instead of first_spring_lat/lon
        lat = row.get('ignition_lat')
        lon = row.get('ignition_lon')
        doy = row.get('ignition_doy')

        if pd.isna(lat) or pd.isna(lon) or pd.isna(doy):
            lt_dates.append(None)
            lt_days.append(None)
            lt_flags.append(None)
            continue

        lt_date = first_lightning_in_buffer(float(lat), float(lon), first_occ, lats, lons)

        if lt_date is None:
            lt_dates.append(None)
            lt_days.append(None)
            lt_flags.append(False)
            continue

        lt_doy    = int(lt_date.timetuple().tm_yday)
        days_diff = int(doy) - lt_doy   # positive = lightning before ignition
                                         # negative = lightning after ignition
        # flag if within 7 days before OR after
        flagged   = abs(days_diff) <= LIGHTNING_MAX_DAYS

        lt_dates.append(lt_date.date())
        lt_days.append(days_diff)
        lt_flags.append(flagged)

    df = df.copy()
    df['first_lightning_date_2km']  = lt_dates
    df['days_lightning_to_hotspot'] = lt_days   # + = strike before ignition
                                                 # - = strike after ignition
    df['lightning_flag']            = lt_flags
    return df

for year in [2024, 2025]:
    df = gpd.read_file(f'workflow/results/overwintering_2/ignitions_{year}.geojson')
    
    # extract lat/lon from geometry for the lightning buffer lookup
    df['ignition_lon'] = df.geometry.to_crs('EPSG:4326').x
    df['ignition_lat'] = df.geometry.to_crs('EPSG:4326').y
    
    df = add_lightning_flag(df, year)
    df.drop(columns='geometry').to_csv(
        f'workflow/results/overwintering_2/ignitions_{year}_lightning.csv',
        index=False
    )
