# MODIS vs VIIRS fire detections (NASA FIRMS)

These scripts compare MODIS (since Nov 2000) and VIIRS S-NPP (since Jan 2012) active-fire detections over California, using the [NASA FIRMS API](https://firms.modaps.eosdis.nasa.gov/api/).

Open `fires.html` in a browser. The page shows:
- a monthly/yearly time series of detections per day,
- two synced detection maps by year,
- **VIIRS in MODIS units** (`consolidate.py`), drawn as a model line you can turn off. MODIS is the reference, so 2000–2011 is shown as measured.
  - **Merge:** VIIRS detections in the same 1 km cell within one overpass are merged into one MODIS-pixel-sized fire. That cuts VIIRS from 4.4× to 2.3× MODIS.
  - **Calibrate:** a fit handles what's left, using how clustered the fires are (VIIRS pixels per merged cell). It is scored leave-one-year-out: typical monthly error ×/÷ 1.20, and the 90% band covers 90% of held-out months.
  - **What the larger sample showed:** a simple fit on raw VIIRS totals is about as accurate (×/÷ 1.19). Merging plus clustering is only slightly better in quiet months (×/÷ 1.24 vs 1.25). Season and NASA POWER weather (`weather.py`: cloud, rain, humidity, temperature, wind) don't improve it enough to be used.
- raw FIRMS rows and `df.describe()` tables for months where both sensors have data.

## Fire Activity Explorer (`app/`)

A web app for the harmonized hotspots: pick a region (preset, typed coordinates, or a box drawn on the map) and a date range, and get:
- four headline tiles, a hotspot density map, a fire-activity chart, a fire-season chart, and the recorded fires in the selection;
- downloads: an HTML report (charts, map, fires and method), hotspots as CSV or GeoJSON, a monthly summary CSV, a recorded-fires CSV, analytics JSON, and each chart as SVG or PNG.

More options let you pick sensors and a minimum confidence.

The selection is kept in the URL, so a view can be shared as a link.

Data: `app/data/hotspots.csv.gz` + `app/data/meta.json`, written by `python3 hotspots.py` (after `build.py`). One row per hotspot: MODIS detections as-is, and VIIRS pixels merged to 1 km per overpass with a `modis_equiv` weight from the calibration, so summing `modis_equiv` gives MODIS-equivalent counts. Any CSV with the same columns can be loaded with "Load hotspot file…".

Two scales (chart header):
- **MODIS units** (default): the line is observed MODIS. VIIRS, merged and calibrated to MODIS units, is drawn dashed as the model's prediction with its 90% range. It would also fill any month without MODIS (e.g. after MODIS is retired), shaded as estimated.
- **Unconverted**: each sensor as detected.

Recorded fires: `app/data/fires.json`, written by `python3 known_fires.py`. It takes the US Forest Service fires (FPA FOD, to 2020) that `viirs_from_modis/fire_spike_matcher.py` matched to spikes in the observed MODIS/VIIRS detections, and links large ones to Wikipedia using the matcher's `wiki_lookup` plus a stricter title check. The app lists them for the selection, marks them on the map, and labels the largest on the activity chart.

Run it: `python3 -m http.server` in the repo, then open http://localhost:8000/app/. (Opened straight from disk, the browser blocks the data fetch; the app then asks you to pick `app/data/hotspots.csv.gz`.)

## `viirs_from_modis/`

Vaangmaya's folder. Its `fire_spike_matcher.py` and validation matches are what `known_fires.py` uses for the named fires. Its MODIS→VIIRS back-cast is no longer used by the main page or the app, and its outputs were built from the earlier 5-day-per-month sample.

## Rebuild

```sh
echo "FIRMS_MAP_KEY=<your key>" > .env
python3 fetch.py   # downloads days 1-5, 11-15 and 21-25 of every month per feed into data/raw (~35 min)
python3 weather.py # optional: NASA POWER weather for the model comparison into data/weather (~1 min, no key)
python3 build.py   # writes fires.html and monthly.csv
python3 hotspots.py    # writes app/data/hotspots.csv.gz and meta.json
python3 known_fires.py # writes app/data/fires.json (named fires + Wikipedia)
```

The main page needs only the Python standard library. The downloaded data is sampled (about half of every month), so the page reports rates as detections per sampled day.
