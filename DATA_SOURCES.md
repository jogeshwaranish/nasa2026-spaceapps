# Data sources and links

Everything the project downloads or loads, with where it is used. Links were checked on 2026-11-14 and returned HTTP 200 (the FPA FOD file URL was downloaded; the NASA POWER API URL needs query parameters, so it answers 422 without them).

## Fire detections (NASA FIRMS)

| What | Link | Used by |
|---|---|---|
| FIRMS area API (docs) | https://firms.modaps.eosdis.nasa.gov/api/ | `fetch.py` downloads three 5-day snippets per month (days 1-5, 11-15, 21-25, about half of every month) |
| FIRMS map key (free, needed for the API) | https://firms.modaps.eosdis.nasa.gov/api/map_key/ | `.env` -> `FIRMS_MAP_KEY` |
| FIRMS archive download (full CSVs by country) | https://firms.modaps.eosdis.nasa.gov/download/ | input for `viirs_from_modis/fire_spike_matcher.py` in daily mode (`--firms`) |

Feeds requested in `fetch.py`:

| Sensor | Feeds | Coverage |
|---|---|---|
| MODIS (Terra + Aqua, 1 km), Collection 6.1 | `MODIS_SP` (standard), `MODIS_NRT` (near-real-time) | Nov 2000 - present |
| VIIRS S-NPP (375 m) | `VIIRS_SNPP_SP`, `VIIRS_SNPP_NRT` | Jan 2012 - present |

## Known fire incidents (used to check the predictions)

| What | Link | Used by |
|---|---|---|
| FPA FOD, "Spatial wildfire occurrence data for the United States, 1992-2020" (Short, Forest Service RDS-2013-0009.6) | https://www.fs.usda.gov/rds/archive/catalog/RDS-2013-0009.6 | catalog page, citation and metadata |
| FPA FOD SQLite download (about 224 MB zipped, 958 MB unzipped; table `Fires`) | https://www.fs.usda.gov/rds/archive/products/RDS-2013-0009.6/RDS-2013-0009.6_Data_Format4_SQLITE.zip | `fire_spike_matcher.py --fod ...`; not committed to this repo |
| Other FPA FOD formats on the same catalog page | `..._Data_Format1_ACCDB.zip`, `..._Data_Format2_GDB.zip`, `..._Data_Format3_GPKG.zip` in the same `products/RDS-2013-0009.6/` folder | alternatives |

## Weather (experiment)

| What | Link | Used by |
|---|---|---|
| NASA POWER daily point API (community AG) | https://power.larc.nasa.gov/api/temporal/daily/point | `weather.py` (optional; no longer used by the model) |
| NASA POWER docs | https://power.larc.nasa.gov/docs/ | reference |

## Optional lookups

| What | Link | Used by |
|---|---|---|
| Wikipedia search API (one request per matched fire) | https://en.wikipedia.org/w/api.php | `fire_spike_matcher.py --wiki`; `known_fires.py` (fires of 5,000+ acres, cached in `data/wiki_cache.json`) |

## Front-end libraries and map tiles (loaded by the HTML pages)

| What | Link |
|---|---|
| Leaflet 1.9.4 (JS) | https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.js |
| Leaflet 1.9.4 (CSS) | https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.css |
| Esri World Light Gray base tiles | https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer |
| Esri World Light Gray reference (labels) | https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Reference/MapServer |

## Product documentation (background)

| What | Link |
|---|---|
| FIRMS: MODIS active fire, attribute descriptions | https://www.earthdata.nasa.gov/data/tools/firms/faq |
| FIRMS home / data catalog | https://firms.modaps.eosdis.nasa.gov/ |

## Files in this repo derived from the sources above

| File | Derived from |
|---|---|
| `monthly.csv`, `fires.html` | FIRMS (sampled days 1-5, 11-15 and 21-25 of each month) |
| `viirs_from_modis/out/*.json`, `*.csv` | FIRMS data, via `train.py` |
| `viirs_from_modis/out/validation/` | FIRMS-derived predictions compared with FPA FOD incidents |
| `app/data/fires.json` | FPA FOD fires matched to observed spikes (`fire_spike_matcher.py` output) + Wikipedia links, via `known_fires.py` |
| `app/data/hotspots.csv.gz`, `app/data/meta.json` | FIRMS hotspots, VIIRS consolidated to 1 km and MODIS weighted by the back-cast into VIIRS units, via `hotspots.py` (read by the Fire Activity Explorer in `app/`) |

Terms of use: FIRMS data is free to use with NASA acknowledgement; FPA FOD is a US Forest Service dataset with its own citation requirement (see its catalog page).
