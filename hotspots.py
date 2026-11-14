#!/usr/bin/env python3
"""Write the harmonized hotspot dataset the web app (app/) reads.

One row per hotspot, both sensors on one scale (consolidated VIIRS fires):
  VIIRS  pixels in the same 1 km cell within one overpass are consolidated into one fire
         (consolidate.merged_groups); lat/lon = mean of the pixels, frp = sum, pixels =
         how many were merged. viirs_equiv = 1: this is the true VIIRS record (2012 on).
  MODIS  each detection as FIRMS returned it. viirs_equiv = that month's back-cast factor,
         so a month's MODIS viirs_equiv sums to the model's predicted VIIRS count
         (consolidate.add_backcast; this is what fills 2000-2011).

Outputs app/data/hotspots.csv.gz and app/data/meta.json (sampled days, schema, model).
Run after build.py:  python3 hotspots.py
"""
import csv, gzip, io, json, re, statistics
from datetime import date, datetime
from pathlib import Path

import consolidate
import fetch

ROOT = Path(__file__).parent
OUT = ROOT / "app" / "data"
COLUMNS = ["date", "time", "lat", "lon", "sensor", "satellite", "daynight", "confidence", "frp", "pixels", "viirs_equiv"]
CONF_RANK = {"low": 0, "nominal": 1, "high": 2}
SCHEMA = {
    "date": "acquisition date, UTC (YYYY-MM-DD)",
    "time": "acquisition time, UTC (HHMM); earliest pixel for merged VIIRS hotspots",
    "lat": "latitude, degrees", "lon": "longitude, degrees",
    "sensor": "MODIS or VIIRS",
    "satellite": "Terra / Aqua (MODIS) or N = Suomi NPP (VIIRS)",
    "daynight": "D or N",
    "confidence": "low / nominal / high (MODIS 0-100 mapped <30 / 30-79 / 80+, as FIRMS recommends)",
    "frp": "fire radiative power, MW (summed over merged VIIRS pixels)",
    "pixels": "FIRMS detections in this hotspot (1 for MODIS)",
    "viirs_equiv": "weight in consolidated-VIIRS fires: 1 for VIIRS, the model's back-cast factor for MODIS; sum it for harmonized counts",
}


def conf(r):
    c = r["confidence"].strip()
    if c in ("l", "n", "h"):
        return {"l": "low", "n": "nominal", "h": "high"}[c]
    c = float(c)
    return "low" if c < 30 else "nominal" if c < 80 else "high"


def main():
    html = (ROOT / "fires.html").read_text(encoding="utf8")
    D = json.loads(re.search(r"const D = (\{.*?\});\n", html, re.S).group(1))
    west, south, east, north = D["bbox"]
    T = D["to_viirs"]
    mod = {r[0]: r for r in D["sensors"]["MODIS"]["monthly"]}
    factor = {}
    for mo, pred, lo, hi in T["monthly"]:
        if pred is not None and mod[mo][1]:
            factor[mo] = pred * mod[mo][2] / mod[mo][1]   # predicted VIIRS fires per MODIS detection
    avail = {k: (date.fromisoformat(a), date.fromisoformat(b)) for k, (a, b) in json.loads((ROOT / "data" / "availability.json").read_text()).items()}

    out, sampled, counts = [], {}, {}
    for sensor, feeds in fetch.SENSORS.items():
        sdays = set()
        for feed in feeds:
            for path in sorted((fetch.RAW / feed).glob("*.csv.gz")):
                start = date.fromisoformat(path.name[:10])
                sdays.update(d.isoformat() for d in fetch.snippet_days(start, avail[feed][1]))
                with gzip.open(path, "rt") as f:
                    rows = list(csv.DictReader(f))
                if sensor == "MODIS":
                    for r in rows:
                        out.append([r["acq_date"], r["acq_time"].zfill(4), round(float(r["latitude"]), 4), round(float(r["longitude"]), 4),
                                    "MODIS", r["satellite"], r["daynight"], conf(r), round(float(r["frp"]), 1), 1, round(factor.get(r["acq_date"][:7], 0), 4)])
                else:
                    for g in consolidate.merged_groups(rows, west, south, north).values():
                        r0 = min(g, key=lambda r: int(r["acq_time"]))
                        cf = max((conf(r) for r in g), key=CONF_RANK.get)
                        out.append([r0["acq_date"], r0["acq_time"].zfill(4),
                                    round(statistics.fmean(float(r["latitude"]) for r in g), 4),
                                    round(statistics.fmean(float(r["longitude"]) for r in g), 4),
                                    "VIIRS", r0["satellite"], r0["daynight"], cf, round(sum(float(r["frp"]) for r in g), 1),
                                    len(g), 1])
        sampled[sensor] = sorted(sdays)
        counts[sensor] = sum(1 for r in out if r[4] == sensor)
    out.sort(key=lambda r: (r[0], r[1], r[4]))

    OUT.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO(); w = csv.writer(buf, lineterminator="\n"); w.writerow(COLUMNS); w.writerows(out)
    with gzip.GzipFile(OUT / "hotspots.csv.gz", "wb", mtime=0) as f:
        f.write(buf.getvalue().encode())
    meta = {
        "name": D["name"], "bbox": D["bbox"], "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "source": "NASA FIRMS area API: MODIS C6.1 (SP + NRT) and VIIRS S-NPP 375 m (SP + NRT)",
        "sampling": fetch.SAMPLING,
        "sampled_days": sampled, "hotspots": counts,
        "harmonization": {"reference": "VIIRS (consolidated)", "merge_km": T["km"], "model": T["fit"]["chosen"],
                          "trained": T["fit"]["years"], "r2": round(T["fit"]["r2"], 3),
                          "typical_error": round(2.718281828 ** T["fit"]["scores"][T["fit"]["chosen"]]["mae_log"], 3),
                          "interval90_log": [round(v, 4) for v in T["fit"]["interval90_log"]],
                          "coverage90": round(T["fit"]["coverage90"], 3),
                          "note": "The back-cast is fitted on region-wide monthly totals; applying it to a sub-region assumes the same relationship holds there."},
        "columns": SCHEMA,
    }
    (OUT / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
    size = (OUT / "hotspots.csv.gz").stat().st_size
    print(f"wrote app/data/hotspots.csv.gz ({size / 1e6:.1f} MB): " + ", ".join(f"{k} {v:,}" for k, v in counts.items()))


if __name__ == "__main__":
    main()
