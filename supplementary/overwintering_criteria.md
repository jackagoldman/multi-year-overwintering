# Overwintering Fire Detection Criteria

| **Criterion** | **Candidate (Spatial and Temporal Proximity)** | **Confirmed (Spatial Continuity)** |
|---|---|---|
| Minimum fall VIIRS hotspot detections within previous-year perimeter (DOY ≥ 213) | ≥ 2 | ≥ 2 |
| Spring VIIRS hotspot within buffer of previous-year perimeter | ≤ 1000 m | ≤ 1000 m |
| Spring detection window relative to snow disappearance date (SDD) | SDD to SDD + 50 days | SDD to SDD + 50 days |
| Spring hotspot distance to current-year fire perimeter | ≤ 1000 m | ≤ 1000 m |
| Fall-to-spring spatial continuity (distance from last fall hotspot to spring hotspot) | Not required | ≤ 1000 m |
| Lightning exclusion (CG strike at spring hotspot location prior to detection) | ≤ 7 days excluded | ≤ 7 days excluded |
| Human ignition exclusion (current-year perimeter fire cause) | FIRE\_CAUSE = Person excluded | FIRE\_CAUSE = Person excluded |
| Snow disappearance date (SDD) source | MODIS fractional snow cover, per-fire raster mean | MODIS fractional snow cover, per-fire raster mean |
| VIIRS hotspot pre-filtering | Nighttime only, ≥ medium confidence, scan/track ≤ 0.5 | Nighttime only, ≥ medium confidence, scan/track ≤ 0.5 |
| Spring hotspot selection method | Hotspot closest to current-year perimeter boundary | Hotspot closest to last fall hotspot location |
