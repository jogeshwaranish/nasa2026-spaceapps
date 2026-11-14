#!/usr/bin/env python3
"""Daily weather from NASA POWER for the fire-prone parts of the region, averaged
over the same sampled days as the FIRMS snippets (13th-17th of each month).

Used by consolidate.py to test whether weather explains what's left of the
MODIS/VIIRS gap after merging pixels. No key needed.

Run `python3 weather.py` to download (cached in data/weather/, ~1 min); build.py
uses the cache if it's there and skips the weather model if it isn't.
"""
import json, statistics, time, urllib.request
from datetime import date
from pathlib import Path

import fetch

ROOT = Path(__file__).parent
CACHE = ROOT / "data" / "weather"
START = "20120101"
PARAMS = {  # POWER name -> short name used in the model
    "CLOUD_AMT": "cloud",       # cloud amount, %
    "PRECTOTCORR": "precip",    # precipitation, mm/day
    "RH2M": "humidity",         # relative humidity at 2 m, %
    "T2M_MAX": "tmax",          # daily max temperature at 2 m, °C
    "WS10M": "wind",            # wind speed at 10 m, m/s
}
# Fire-prone parts of California, equally weighted (the bbox itself is mostly ocean and desert).
POINTS = {
    "Klamath / North Coast": (41.0, -123.2), "Northeast plateau": (41.0, -121.0),
    "Sacramento Valley": (39.5, -122.0), "Northern Sierra": (38.8, -120.5),
    "San Joaquin Valley": (36.5, -119.8), "Central Coast": (35.6, -120.8),
    "Southern California": (34.3, -118.0), "San Diego / Inland": (33.0, -116.8),
}
URL = ("https://power.larc.nasa.gov/api/temporal/daily/point?parameters={p}&community=AG"
       "&longitude={lon}&latitude={lat}&start={s}&end={e}&format=JSON")


def download():
    CACHE.mkdir(parents=True, exist_ok=True)
    end = date.today().strftime("%Y%m%d")
    for name, (lat, lon) in POINTS.items():
        out = CACHE / f"{lat}_{lon}.json"
        url = URL.format(p=",".join(PARAMS), lat=lat, lon=lon, s=START, e=end)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=180) as r:
                    data = json.load(r)["properties"]["parameter"]
                break
            except Exception as ex:
                print(f"  {name}: {ex}, retrying"); time.sleep(5)
        else:
            raise SystemExit(f"could not download {name}")
        out.write_text(json.dumps(data))
        print(f"{name}: {len(data['RH2M'])} days")


def monthly():
    """{"YYYY-MM": {short name: mean over sampled days and points}}; None if not downloaded."""
    files = [CACHE / f"{lat}_{lon}.json" for lat, lon in POINTS.values()]
    if not all(f.exists() for f in files):
        return None
    acc = {}
    for f in files:
        data = json.loads(f.read_text())
        for p, short in PARAMS.items():
            for d, v in data[p].items():
                if v is None or v <= -999 or not (fetch.SAMPLE_DAY <= int(d[6:]) < fetch.SAMPLE_DAY + fetch.WINDOW):
                    continue
                acc.setdefault(f"{d[:4]}-{d[4:6]}", {}).setdefault(short, []).append(v)
    return {m: {k: statistics.fmean(v) for k, v in vals.items()} for m, vals in sorted(acc.items())}


if __name__ == "__main__":
    download()
    w = monthly()
    print(f"{len(w)} months; latest {max(w)}: {w[max(w)]}")
