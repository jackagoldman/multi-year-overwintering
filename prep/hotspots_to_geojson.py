# Convert the VIIRS hotspot GeoParquet to the GeoJSON read by the detection algorithm.

import os
import yaml
import geopandas as gpd

cfg = yaml.safe_load(open('config.yml'))
data_root = cfg['environments'][cfg['environment']]['data_root']

src = '' #= path to your file
dst = os.path.join(data_root, cfg['data']['fires']['all_hotspots'])

os.makedirs(os.path.dirname(dst), exist_ok=True)
hotspots = gpd.read_parquet(src)
hotspots.to_file(dst, driver='GeoJSON')

print(f'{len(hotspots)} hotspots written to {dst}')
