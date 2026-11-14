# MODIS vs VIIRS fire detections (NASA FIRMS)

These scripts compare MODIS (since Nov 2000) and VIIRS S-NPP (since Jan 2012) active-fire detections over California, using the [NASA FIRMS API](https://firms.modaps.eosdis.nasa.gov/api/).

Open `fires.html` in a browser. The page shows:
- a monthly/yearly time series of detections per day,
- two synced detection maps by year,
- raw FIRMS rows and `df.describe()` tables for months where both sensors have data.

## Rebuild

```sh
echo "FIRMS_MAP_KEY=<your key>" > .env
python3 fetch.py   # downloads one 5-day snippet per month per feed into data/raw (~10 min)
python3 build.py   # writes fires.html and monthly.csv
```

Python standard library only. The downloaded data is sampled (5 days per month), so the page reports rates as detections per sampled day.
