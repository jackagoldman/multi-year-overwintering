# 04_roads.py — Road Proximity for 2025 Ignition Points
# =======================================================
# For 2024 ignitions: fire_cause is already populated — no road check needed.
# For 2025 ignitions: fire_cause is blank, so road distance is used as a
# proxy for human-caused ignitions.
#
# For each 2025 ignition point, computes distance to nearest road and
# records the road class. Flags ignitions within 1000m as potentially
# human-caused.
#
# Inputs
#   workflow/results/intermediate_products/ignitions_2025.geojson
#   data/roads_study_zones_3005.geojson
#
# Outputs
#   workflow/results/intermediate_products/ignitions_2025_roads.csv
#
# New columns
#   road_dist_m        : distance in metres from ignition point to nearest road
#   nearest_road_class : RD_CLASS value of the nearest road segment
#   road_flag          : True if road_dist_m <= 1000, else False

import os
import numpy as np
import geopandas as gpd

ROAD_DIST_M    = 1000
IGNITION_PATH  = 'results/ignitions_2025.geojson'
ROADS_PATH     = 'data/roads_study_zones_3005.geojson'
OUT_PATH       = 'results/ignitions_2025_roads.csv'

ignitions = gpd.read_file(IGNITION_PATH).to_crs('EPSG:3005')
roads     = gpd.read_file(ROADS_PATH).to_crs('EPSG:3005')

road_dists   = []
road_classes = []
road_flags   = []

for _, row in ignitions.iterrows():
    pt = row.geometry

    if pt is None or pt.is_empty:
        road_dists.append(None)
        road_classes.append(None)
        road_flags.append(None)
        continue

    dists       = roads.geometry.distance(pt)
    nearest_idx = dists.idxmin()
    dist_m      = dists[nearest_idx]

    road_dists.append(round(dist_m, 1))
    road_classes.append(roads.loc[nearest_idx, 'RD_CLASS'])
    road_flags.append(dist_m <= ROAD_DIST_M)

ignitions['road_dist_m']        = road_dists
ignitions['nearest_road_class'] = road_classes
ignitions['road_flag']          = road_flags

ignitions.drop(columns=['geometry']).to_csv(OUT_PATH, index=False)
