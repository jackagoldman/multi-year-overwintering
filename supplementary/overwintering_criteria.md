# Multi-year Overwintering Fire Detection Algorithm Criteria

| **Criterion** | **Detection Algorithm** |
|---|---|
| Minimum fall VIIRS hotspot detections within previous-year perimeter (DOY ≥ 213) | ≥ 2 detections |
| Spring VIIRS hotspot within buffer of previous-year perimeter | ≤ 2000 m |
| Spring detection window relative to snow disappearance date (SDD) | SDD − 5 days to SDD + 60 days |
| Ignition point selection | First-day hotspot closest to confirmed previous-year perimeter boundary |
| Lightning exclusion (CG strike within 2 km of ignition point) | Excluded if strike within ± 7 days of ignition |
| Human ignition exclusion — 2024 | FIRE\_CAUSE = Person excluded |
| Human ignition exclusion — 2025 | Road proximity ≤ 1000 m excluded |
| Snow disappearance date (SDD) source | MODIS fractional snow cover, per-fire zonal mean |
| VIIRS hotspot pre-filtering | Nighttime only, ≥ medium confidence, scan/track ≤ 0.5 |
| Multi-year chain detection | 2024 fire must have confirmed fall smoldering (≥ 2 hotspots DOY ≥ 213) and a 2025 ignition point within 2000 m of its perimeter |