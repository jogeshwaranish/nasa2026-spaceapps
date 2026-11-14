# MODIS vs VIIRS fire detections (NASA FIRMS)

These scripts compare MODIS (since Nov 2000) and VIIRS S-NPP (since Jan 2012) active-fire detections over California, using the [NASA FIRMS API](https://firms.modaps.eosdis.nasa.gov/api/).

Open `fires.html` in a browser. The page shows:
- a monthly/yearly time series of detections per day,
- two synced detection maps by year,
- a model line you can switch between:
  - **VIIRS in MODIS units** (`consolidate.py`): VIIRS detections in the same 1 km cell within one overpass are merged into one MODIS-pixel-sized detection, which cuts VIIRS from 4.4× to 2.3× MODIS. A calibration fit then handles what's left, using how clustered the fires are (VIIRS pixels per merged cell) and the season. Weather on the sampled days from NASA POWER (`weather.py`: cloud, rain, humidity, temperature, wind) was tested and doesn't improve it. It is scored leave-one-year-out, and its 90% band covers 90% of held-out months.
  - **MODIS in VIIRS units** (`harmonize.py`): the back-cast prototype, which estimates what VIIRS would have seen in 2000–2011.
- raw FIRMS rows and `df.describe()` tables for months where both sensors have data.

## Fire Activity Explorer (`app/`)

A web app for the harmonized hotspots: pick a region (preset, typed coordinates, or a box drawn on the map) and a date range, choose sensors, units (MODIS-equivalent or raw) and minimum confidence, and get:
- summary tiles, a hotspot density map, monthly and seasonal charts, a fire-intensity (FRP) distribution and the most active places;
- downloads: an HTML report (all charts, map and tables), hotspots as CSV or GeoJSON, a monthly summary CSV, analytics JSON, and each chart as SVG or PNG.

The selection is kept in the URL, so a view can be shared as a link.

Data: `app/data/hotspots.csv.gz` + `app/data/meta.json`, written by `python3 hotspots.py` (after `build.py`). One row per hotspot: MODIS detections as-is, and VIIRS pixels merged to 1 km per overpass with a `modis_equiv` weight from the calibration, so summing `modis_equiv` gives MODIS-equivalent counts. Any CSV with the same columns can be loaded with "Load hotspot file…".

Two harmonized scales (chart header):
- **MODIS units** (default): MODIS is the reference. VIIRS hotspots are merged to 1 km and calibrated to MODIS (`modis_equiv`). The line is observed MODIS; converted VIIRS fills any month without MODIS (e.g. after MODIS is retired) and is drawn dashed as the model prediction, with its 90% range.
- **VIIRS units**: VIIRS is the reference. A model trained on 2012+ (`viirs_from_modis/train.py`, leave-one-year-out) predicts VIIRS from MODIS. The line is observed VIIRS from 2012 and the prediction before that, shaded as estimated; the prediction is also drawn dashed with its 90% range across all years, so it can be compared with observed VIIRS from 2012. `hotspots.py` stores the per-month factors in `meta.json` (`backcast`), so each MODIS hotspot carries a VIIRS-equivalent weight for any region or date selection.

Recorded fires: `app/data/fires.json`, written by `python3 known_fires.py`. It takes the US Forest Service fires (FPA FOD, to 2020) that `viirs_from_modis/fire_spike_matcher.py` matched to spikes in the observed MODIS/VIIRS detections, and links large ones to Wikipedia using the matcher's `wiki_lookup` plus a stricter title check. The app lists them for the selection, marks them on the map, and labels the largest on the activity chart.

Run it: `python3 -m http.server` in the repo, then open http://localhost:8000/app/. (Opened straight from disk, the browser blocks the data fetch; the app then asks you to pick `app/data/hotspots.csv.gz`.)

## Back-cast models

`viirs_from_modis/` holds separate back-cast models (ridge and gradient boosting, plus a per-cell spatial model). They need numpy, scipy and scikit-learn.

## Rebuild

```sh
echo "FIRMS_MAP_KEY=<your key>" > .env
python3 fetch.py   # downloads one 5-day snippet per month per feed into data/raw (~10 min)
python3 weather.py # optional: NASA POWER weather for the model comparison into data/weather (~1 min, no key)
python3 build.py   # writes fires.html and monthly.csv
```

The main page needs only the Python standard library. The downloaded data is sampled (5 days per month), so the page reports rates as detections per sampled day.
