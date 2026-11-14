#!/usr/bin/env python3
"""Aggregate the cached FIRMS downloads into fires.html (self-contained page).

The data is a sample (5 days per month), so everything is reported as a rate:
detections per day, averaged over the days actually downloaded.
"""
import csv, gzip, json, calendar, math, statistics
from collections import defaultdict
from datetime import date, timedelta, datetime
from pathlib import Path

import consolidate
import fetch
import harmonize
import weather

ROOT = Path(__file__).parent
DATA = ROOT / "data"
CELL = 0.05  # map grid size in degrees (~5 km)
HEAD = 10    # raw rows shown per snippet in the DataFrame view
NUMERIC = {"latitude", "longitude", "brightness", "bright_ti4", "scan", "track", "bright_t31", "bright_ti5", "frp"}


def describe(rows, cols):
    """pandas-style describe() for the numeric columns: count, mean, std, min, 25%, 50%, 75%, max."""
    out = {}
    for c in cols:
        if c not in NUMERIC:
            continue
        v = sorted(float(r[c]) for r in rows if r[c] not in ("", None))
        if not v:
            continue
        q = statistics.quantiles(v, n=4, method="inclusive") if len(v) > 1 else [v[0]] * 3
        std = statistics.stdev(v) if len(v) > 1 else None
        out[c] = [len(v), statistics.fmean(v), std, v[0], q[0], q[1], q[2], v[-1]]
    return out


def availability():
    cache = DATA / "availability.json"
    try:
        av = {k: [str(a), str(b)] for k, (a, b) in fetch.availability(fetch.map_key()).items()}
        cache.write_text(json.dumps(av))
    except Exception:
        av = json.loads(cache.read_text())
    return {k: (date.fromisoformat(a), date.fromisoformat(b)) for k, (a, b) in av.items()}


def main():
    meta = json.loads((DATA / "meta.json").read_text())
    west, south, east, north = map(float, meta["bbox"].split(","))
    avail = availability()
    out = {"name": meta["name"], "bbox": [west, south, east, north], "cell": CELL,
           "generated": datetime.now().strftime("%Y-%m-%d %H:%M"), "sensors": {}}
    table = []
    merged = {km: defaultdict(int) for km in consolidate.SENSITIVITY}   # VIIRS merged to ~MODIS pixels

    for sensor, feeds in fetch.SENSORS.items():
        first = min(avail[f][0] for f in feeds)
        last = max(avail[f][1] for f in feeds)
        m_count, m_days = defaultdict(int), defaultdict(int)       # keyed (year, month)
        y_cells = defaultdict(lambda: defaultdict(int))             # year -> (ix, iy) -> count
        samples = {}                                                # month -> raw rows + describe()
        for feed in feeds:
            for path in sorted((fetch.RAW / feed).glob("*.csv.gz")):
                start = date.fromisoformat(path.name[:10])
                for i in range(fetch.WINDOW):                       # days this snippet covers
                    d = start + timedelta(days=i)
                    if d <= avail[feed][1] and d.day < fetch.SAMPLE_DAY + fetch.WINDOW:
                        m_days[(d.year, d.month)] += 1
                with gzip.open(path, "rt") as f:
                    rows = list(csv.DictReader(f))
                if rows:
                    cols = list(rows[0].keys())
                    samples[path.name[:7]] = {"feed": feed, "start": str(start), "n": len(rows), "cols": cols,
                                              "head": [[r[c] for c in cols] for r in rows[:HEAD]],
                                              "desc": describe(rows, cols)}
                if True:
                    for r in rows:
                        y, m = int(r["acq_date"][:4]), int(r["acq_date"][5:7])
                        m_count[(y, m)] += 1
                        ix = int((float(r["longitude"]) - west) / CELL)
                        iy = int((float(r["latitude"]) - south) / CELL)
                        y_cells[y][(ix, iy)] += 1
                    if sensor == "VIIRS":
                        for km in merged:
                            for k in consolidate.merged_keys(rows, west, south, north, km):
                                merged[km][k[0][:7]] += 1

        monthly, years = [], {}
        y, m = first.year, first.month
        while (y, m) <= (last.year, last.month):
            lo = max(first, date(y, m, 1))
            hi = min(last, date(y, m, calendar.monthrange(y, m)[1]))
            possible = (hi - lo).days + 1
            monthly.append([f"{y}-{m:02d}", m_count[(y, m)], m_days[(y, m)], possible])
            table.append([sensor, f"{y}-{m:02d}", m_count[(y, m)], m_days[(y, m)], possible,
                          merged[consolidate.KM][f"{y}-{m:02d}"] if sensor == "VIIRS" else ""])
            yr = years.setdefault(y, {"count": 0, "days": 0, "possible": 0})
            yr["count"] += m_count[(y, m)]; yr["days"] += m_days[(y, m)]; yr["possible"] += possible
            y, m = (y + 1, 1) if m == 12 else (y, m + 1)
        for y, yr in years.items():
            yr["cells"] = [v for (ix, iy), c in sorted(y_cells[y].items()) for v in (ix, iy, c)]
        out["sensors"][sensor] = {"first": str(first), "last": str(last), "monthly": monthly, "years": years, "samples": samples}
        got, poss = sum(r[2] for r in monthly), sum(r[3] for r in monthly)
        print(f"{sensor}: {sum(r[1] for r in monthly):,} detections in {got} sampled days ({sum(1 for r in monthly if r[2])}/{len(monthly)} months)")

    # raw-value comparison only makes sense where both sensors returned rows
    both = set(out["sensors"]["MODIS"]["samples"]) & set(out["sensors"]["VIIRS"]["samples"])
    for s in out["sensors"].values():
        s["samples"] = {k: v for k, v in sorted(s["samples"].items()) if k in both}

    with open(ROOT / "monthly.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sensor", "month", "detections", "days_downloaded", "days_in_record", f"detections_merged_{consolidate.KM:g}km"])
        w.writerows(table)
    harmonize.add_harmonized(out)
    consolidate.add_to_modis(out, {km: dict(v) for km, v in merged.items()}, weather.monthly())
    f = out["to_modis"]["fit"]
    print(f"VIIRS per MODIS: raw {f['ratio_raw']:.2f}x, merged {f['ratio_merged']:.2f}x; "
          f"VIIRS->MODIS model: {f['chosen']}, 90% band coverage {f['coverage90']:.0%}")
    for k in f["order"]:
        s = f["scores"][k]
        print(f"  {k:38s} held-out error x{math.exp(s['mae_log']):.3f}  quiet months x{math.exp(s['mae_log_quiet']):.3f}")
    html = (ROOT / "template.html").read_text().replace("/*DATA*/null", json.dumps(out, separators=(",", ":")))
    (ROOT / "fires.html").write_text(html)
    print(f"wrote fires.html ({len(html) / 1e6:.1f} MB) and monthly.csv")


if __name__ == "__main__":
    main()
