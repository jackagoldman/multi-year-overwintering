# clip nbac files to the study zones 
# read in the study zones and then dissolve them into one big zone, 
# clip the study nbac .shp to all the study zones, and then save the clipped files to a new folder.
# THIS CODE ASSUMES YOU HAVE A INDEPENDENT SHAPEFILE FOR EACH NBAC YEAR - IN THIS CASE 2023-2025 
import geopandas as gpd
import os
import geopandas as pd
import yaml
from pathlib import Path

# Find the project root by walking up from the notebook's folder until config.yml turns up
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / 'config.yml').exists())

with open(ROOT / 'config.yml') as f:
    cfg = yaml.safe_load(f)

data_root = cfg['environments'][cfg['environment']]['data_root']

def d(subpath):
    return os.path.join(data_root, subpath)
# read in the study zones
study_zones = gpd.read_file(d(cfg['data']['plotting']['study_zones']))

# study zones crs
print(study_zones.crs)
# dissolve the study zones into one big zone
study_zone = study_zones.dissolve()
# read in the nbac .shp files
nbac_folder = d(cfg['data']['fires']['nbac'])
nbac_files = [f for f in os.listdir(nbac_folder) if f.endswith('.shp')]
# create a new folder to save the clipped files
clipped_folder = d(cfg['data']['fires']['nbac_clipped'])
if not os.path.exists(clipped_folder):
    os.makedirs(clipped_folder)

# nbac crs
for nbac_file in nbac_files:
    nbac_path = os.path.join(nbac_folder, nbac_file)
    nbac_gdf = gpd.read_file(nbac_path)
    print(f"{nbac_file} crs: {nbac_gdf.crs}")

# loop through the nbac files and clip them to the study zone
for nbac_file in nbac_files:
    nbac_path = os.path.join(nbac_folder, nbac_file)
    nbac_gdf = gpd.read_file(nbac_path)
    nbac_gdf = nbac_gdf.to_crs(study_zone.crs)
    clipped_gdf = gpd.clip(nbac_gdf, study_zone)
    clipped_path = os.path.join(clipped_folder, nbac_file)
    clipped_gdf.to_file(clipped_path)



# read in a single nbac and get columns
nbac_gdf = gpd.read_file(os.path.join(nbac_folder, nbac_files[0]))
print(nbac_gdf.columns)


#read in files from clipped folder and get columns
clipped_files = [f for f in os.listdir(clipped_folder) if f.endswith('.shp')]
# remove if starts with ._
clipped_files = [f for f in clipped_files if not f.startswith('._')]
for clipped_file in clipped_files:
    clipped_path = os.path.join(clipped_folder, clipped_file)
    clipped_gdf = gpd.read_file(clipped_path)
    print(f"{clipped_file} columns: {clipped_gdf.columns}")


# check columns for BC perims
old_fire = d(cfg['data']['perimeters']['perims_folder'])
perims_files = [f for f in os.listdir(old_fire) if f.endswith('.shp')]
perims_files = [f for f in perims_files if not f.startswith('._')]
for perims_file in perims_files:
    perims_path = os.path.join(old_fire, perims_file)
    perims_gdf = gpd.read_file(perims_path)
    print(f"{perims_file} columns: {perims_gdf.columns}")

# rename FIRECAUS to fire_cause, YEAR to fire_year, GID to fire_id, POLY_HA to siza_ha, build area_sqm
for clipped_file in clipped_files:
    clipped_path = os.path.join(clipped_folder, clipped_file)
    clipped_gdf = gpd.read_file(clipped_path)
    clipped_gdf = clipped_gdf.rename(columns={
        "FIRECAUS": "fire_cause",
        "YEAR": "fire_year",
        "GID": "fire_id",
        "POLY_HA": "size_ha"
    })
    clipped_gdf["area_sqm"] = clipped_gdf["size_ha"] * 10000
    print(f"{clipped_file} columns after renaming: {clipped_gdf.columns}")
    #save to clipped_folder with same name but with _renamed at the end
    renamed_path = os.path.join(clipped_folder, clipped_file.replace('.shp', '_renamed.shp'))
    clipped_gdf.to_file(renamed_path)

# combine all the fires into a single dataframe and write as .geojson
clipped_files = [f for f in os.listdir(clipped_folder) if f.endswith('_renamed.shp')]
clipped_files = [f for f in clipped_files if not f.startswith('._')]
all_fires_list = []
for clipped_file in clipped_files:
    clipped_path = os.path.join(clipped_folder, clipped_file)
    clipped_gdf = gpd.read_file(clipped_path)
    all_fires_list.append(clipped_gdf)

all_fires = gpd.GeoDataFrame(pd.concat(all_fires_list, ignore_index=True), geometry="geometry")

# if fire intersects with study zones, add a column called fire_zone and assign the value of the fire zone from the study zones columns MFFRZNNM
all_fires['fire_zone'] = None
for index, fire in all_fires.iterrows():
    for _, zone in study_zones.iterrows():
        if fire.geometry.intersects(zone.geometry):
            all_fires.at[index, 'fire_zone'] = zone['MFFRZNNM']
            break

# write as geojson in clipped_folder
all_fires.to_file(os.path.join(clipped_folder, "all_fires_processed.geojson"), driver='GeoJSON')


