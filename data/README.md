# Data folder

## The Mendeley Data workbook (required)

Download `Shenyang_metro_entrance_acoustic_environment_data.xlsx` from Mendeley Data —
**DOI [10.17632/7rwt35tchz](https://doi.org/10.17632/7rwt35tchz)** (CC BY 4.0) — and place it in **this folder**:

```
data/Shenyang_metro_entrance_acoustic_environment_data.xlsx
```

The workbook is not tracked in this repository — it is archived at Mendeley Data, which is its citable home.

Sheets the analysis reads (`code/load_data.py` extracts them to `intermediate/clean/`):

| Sheet | What it holds |
|---|---|
| `station_master` | one row per station site (111): station attributes, rounded coordinates, the dominant sound event, the four-group sound-source shares, and the land-use, road-density, point-of-interest and NDVI indicators at each buffer distance |
| `station_hour_panel` | one row per station and hour (2,664): hourly LAeq, primary and secondary sound-event labels, day type and the population heat values at 20 m and 50 m |
| `noise_metrics_station` | the station-level acoustic indicators derived from the hourly LAeq series (24 h, daytime, night-time, peak and inter-peak levels, variability and exceedance shares) |
| `landuse_name_map` | the map from the original land-use indicator names to the cleaned column names, used to group the indicators by category and scale |

The workbook also carries a README sheet, a sheet summary, the raw-form hourly tables (`laeq_hourly`, `sound_event_hourly`, `people_heat_hourly_*`), the `sites` and `landuse` tables from which `station_master` was assembled, and the dictionary sheets that define every column and code. None of those is read here.

## Shipped with the repository

- `osm_roads_shenyang.csv` — the major-road polylines (motorway, trunk and primary) drawn as the grey context layer of Figures 1 and 5. Fetched from the OpenStreetMap Overpass API by `code/fetch_osm.py` on 2026-09-12. Data © OpenStreetMap contributors, available under the [Open Database License](https://www.openstreetmap.org/copyright). Rerun the fetch script only to refresh the cache; the figures draw nothing if the file is absent.

## Optional inputs

`validation/` holds two small tables that are not part of the public deposit — see `validation/README.md`. Everything else runs without them.
