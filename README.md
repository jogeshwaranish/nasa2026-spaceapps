# MODIS vs VIIRS fire detections (NASA FIRMS)

These scripts compare MODIS (since Nov 2000) and VIIRS S-NPP (since Jan 2012) active-fire detections over California, using the [NASA FIRMS API](https://firms.modaps.eosdis.nasa.gov/api/).

Open `fires.html` in a browser. The page shows:
- a monthly/yearly time series of detections per day,
- two synced detection maps by year,
- **Harmonization in VIIRS units** (`consolidate.py`), drawn as a dashed model line with a 90% band:
  - **Consolidate:** VIIRS detections in the same 1 km cell within one overpass count as one fire. That cuts VIIRS from 4.4× to 2.3× MODIS, and this consolidated count is the true VIIRS record (2012 on).
  - **Back-cast:** on the 179 months both satellites flew, a model learns consolidated VIIRS from MODIS (+ neighbouring months + season, R² 0.94). Run on every MODIS month, it predicts the VIIRS that would probably have been seen in 2000–2011. From 2012 it can be checked against true VIIRS.
  - **How good it is:** scored leave-one-year-out, the typical monthly error is ×/÷ 1.17, and the 90% band covers 90% of held-out months. A plain log-log fit is nearly as good (×/÷ 1.17); a constant ratio is much worse (×/÷ 1.44).
- raw FIRMS rows and `df.describe()` tables for months where both sensors have data.

## Fire Activity Explorer (`app/`)

A web app for the harmonized hotspots. From top to bottom:
- **Map and selection:** pick a region (preset, or a box drawn on the map) and dates. Exact coordinates and a confidence filter are tucked under the map.
- **Four tiles:** fire activity in VIIRS units, the **model prediction** from MODIS (compared with true VIIRS where both exist), the peak month, and the largest recorded fire.
- **Fire activity over time:** true VIIRS from 2012 (solid) and the model prediction from MODIS back to 2000 (dashed) with its 90% range. Months before 2012 are shaded: there, the prediction is the only record. Labelled with the largest recorded fires.
- **Recorded fires and fire season.**
- **How the data was harmonized**, for the current selection:
  - three step charts on one scale: raw VIIRS (≈4.4× MODIS), consolidated VIIRS (≈2.3×), and the model prediction vs true VIIRS (≈0.9–1.0×);
  - a scatter of model prediction vs true VIIRS per month;
  - a table of how close the prediction is, alongside the held-out test.
- **Downloads:** an HTML report (charts, map, fires, harmonization check, method), hotspots as CSV or GeoJSON, a monthly CSV, a recorded-fires CSV, analytics JSON, and each chart as SVG or PNG.

The selection is kept in the URL, so a view can be shared as a link.

Data: `app/data/hotspots.csv.gz` + `app/data/meta.json`, written by `python3 hotspots.py` (after `build.py`). One row per hotspot: VIIRS pixels consolidated to 1 km per overpass (`viirs_equiv` = 1), and MODIS detections with a `viirs_equiv` weight from the back-cast, so summing `viirs_equiv` gives counts in consolidated-VIIRS units. Any CSV with the same columns can be loaded with "Load hotspot file…".

Recorded fires: `app/data/fires.json`, written by `python3 known_fires.py`. It takes the US Forest Service fires (FPA FOD, to 2020) that `viirs_from_modis/fire_spike_matcher.py` matched to spikes in the observed MODIS/VIIRS detections, cleans up the record names (e.g. `MCNALLY` → McNally Fire, `DOE` → August Complex (Doe Fire), out-of-state fires get their state), and links large ones to Wikipedia: an exact article title whose intro names the fire's year and state, or the matcher's `wiki_lookup` with a stricter title check. The app lists them for the selection, marks them on the map, and labels the largest on the activity chart.

Run it: `python3 -m http.server` in the repo, then open http://localhost:8000/app/. (Opened straight from disk, the browser blocks the data fetch; the app then asks you to pick `app/data/hotspots.csv.gz`.)

## `viirs_from_modis/`

Vaangmaya's folder. Its `fire_spike_matcher.py` and validation matches are what `known_fires.py` uses for the named fires. Its own MODIS→VIIRS back-cast (on raw VIIRS) is not used; the main page and app use `consolidate.py`'s back-cast on consolidated VIIRS. Its outputs were built from the earlier 5-day-per-month sample.

## Rebuild

```sh
echo "FIRMS_MAP_KEY=<your key>" > .env
python3 fetch.py   # downloads days 1-5, 11-15 and 21-25 of every month per feed into data/raw (~35 min)
python3 build.py   # writes fires.html and monthly.csv
python3 hotspots.py    # writes app/data/hotspots.csv.gz and meta.json
python3 known_fires.py # writes app/data/fires.json (named fires + Wikipedia)
```

The main page needs only the Python standard library. The downloaded data is sampled (about half of every month), so the page reports rates as detections per sampled day.
