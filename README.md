# MODIS vs VIIRS fire detections (NASA FIRMS)

These scripts compare MODIS (since Nov 2000) and VIIRS S-NPP (since Jan 2012) active-fire detections over California, using the [NASA FIRMS API](https://firms.modaps.eosdis.nasa.gov/api/).

Open `fires.html` in a browser. The page shows:
- a monthly/yearly time series of detections per day,
- two synced detection maps by year,
- a model line you can switch between:
  - **VIIRS in MODIS units** (`consolidate.py`): VIIRS detections in the same 1 km cell within one overpass are merged into one MODIS-pixel-sized detection, which cuts VIIRS from 4.4× to 2.3× MODIS. A calibration fit then handles what's left, using how clustered the fires are (VIIRS pixels per merged cell) and the season. Weather on the sampled days from NASA POWER (`weather.py`: cloud, rain, humidity, temperature, wind) was tested and doesn't improve it. It is scored leave-one-year-out, and its 90% band covers 90% of held-out months.
  - **MODIS in VIIRS units** (`harmonize.py`): the back-cast prototype, which estimates what VIIRS would have seen in 2000–2011.
- raw FIRMS rows and `df.describe()` tables for months where both sensors have data.

`viirs_from_modis/` holds separate back-cast models (ridge and gradient boosting, plus a per-cell spatial model). They need numpy, scipy and scikit-learn.

## Rebuild

```sh
echo "FIRMS_MAP_KEY=<your key>" > .env
python3 fetch.py   # downloads one 5-day snippet per month per feed into data/raw (~10 min)
python3 weather.py # optional: NASA POWER weather for the model comparison into data/weather (~1 min, no key)
python3 build.py   # writes fires.html and monthly.csv
```

The main page needs only the Python standard library. The downloaded data is sampled (5 days per month), so the page reports rates as detections per sampled day.
