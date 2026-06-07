# Multi-year overwintering detection
This algorithm is used to detect overwintering and specifically multi-year overwintering fires

The algorithm is based on previous work described by Scholten et al. 2021 but with key metholodological differences.

Step 1: Quantified Snow Disappearance Date for each fire perimeter (adapted from Scholten et al. 2021)

We extract the mean SDD from the MODIS raster over each 2024 and 2025 fire perimeter using zonal statistics. This gives a per-fire snowmelt date rather than Scholten's regional mean. Hotspots before SDD are excluded from the spring pool, the temporal window is SDD to SDD + 50 days, matching Scholten's 50-day threshold.


Step 2: Confirm fall smouldering (not in Scholten et al. 2021)

Because fall smouldering (or burning into the fall) is a key indicator of potential overwintering, it is essential to confirm its presence before proceeding with the detection of multi-year overwintering fires. This is because if a fire perimeter in the preceeding year had no activity detected in spring, it may be unlikely that the fire could be an overwintering fire. Therefore for 2023 and 2024 fire perimeters, which represent the fires that precede overwintering identification, we extract all VIIRS hotspots detected after August 1 (DOY ≥ 213) that fall within the perimeter boundary. This produces a set of confirmed smouldering or fall active fire perimeter and only these perimeters are used as a spatial references in the ignition point search. In scholten, any previous-year fire periemter is a valid reference in their framework regardless of whether it was still burning before winter. 

Step 3: Derive ignition points (adapted from Scholten et al. 2021)

For each 2024 and 2025 fire perimeter, we find all VIIRS hotspots detected on the earliest date inside the fire perimeter. Among those fire-dat hotspots, we select the single point closest to any confirmed previous-year perimeter boundary. This is the ignition point. 

Step 4: spatial overlap with previous year burn scars (adapted from Scholten et al. 2021)

We check whetehr the ignition point falls within a distance threshold of a confirmed previous-year perimeter. In scholten this threshold is 1000m, which they derived from MODIS picel uncertainty (926m) plus a small buffer. We compute a column called `dist_to_prev_perim_m` which is the key output of this step, it is the Euclidean distance from the chosen ignition point to the nearest confirmed previous-year perimeter boundary. This allows us the run sensitivity analyses to see how our detection algorithm differs based on different distances. 

Step 5: temporal filter (matches scholten et al. 2021)

We keep only ignition points where `days_after_sdd` falls between 0 and 50. This ensures that only fires occurring within the defined temporal window are considered for multi-year overwintering detection.

Step 6: Lightning exclusion (adapted from Scholten et al. 2021)

For each canadiate ignition point, we find the earliest CG strike within a 2km buffer using the raw 3-hourly lightning data from CLDN for 2024 and 2025 respectively. We flag ignition if a strike occured within 7 days before or after detection. We use a geometrically precise 2km buffer projected in EPSG:3005 and search all grid cells whose centres fall within the buffer, returning the earliest strike date among them.

Step 7: Human ignition exclusion (departure from Scholten et al. 2021)

For 2024 ignitions, we use the BC Wildfire Service `fire_cause` attribute directly, whereby fires attributed to a human igition `person` are exlude. For 2025 ignitions, `fire_cause` is not available and in this case we use road proximity as a proxy: ignitions within 100m of a road are flagged (analogous to scholten) as potentially human-caused. 

Step 8: Multi-year chain detection 

Our mutli-year chains required three confirmed links:
1. A 2023 fire with confirmed fall smouldering/hotspots
2. A 2024 fire whose ignition point is within XXXXm if that 2023 perimeter and which itself has confirmed fall smouldering
3. A 2025 fire whose ignition is within XXXXm of that same 2024 perimeter

The chain is then identified by joining the 2023->2024 and 2024->2025 candidate tables on the 2024 `fire_id`. Fires that overwintered exactly once and are not part of a chain are classified as single-year overwinters. 