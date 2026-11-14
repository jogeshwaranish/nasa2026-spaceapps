#!/usr/bin/env python3
"""Download a sample of MODIS + VIIRS active-fire detections from NASA FIRMS.

The area API allows 5 days per request and 5000 transactions / 10 min (a 5-day
request costs 10), so instead of the full record this takes three 5-day snippets
per month (days 1-5, 11-15 and 21-25, about half of every month) across each
sensor's whole lifetime: ~1,500 requests, about 35 minutes. Snippets are cached in data/raw/<SOURCE>/<start>.csv.gz and
skipped on re-run.

    python3 fetch.py                       # California, MODIS + VIIRS S-NPP
    python3 fetch.py --bbox -10,35,5,44    # west,south,east,north
"""
import argparse, csv, gzip, io, itertools, json, os, sys, threading, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import date, timedelta
from pathlib import Path

API = "https://firms.modaps.eosdis.nasa.gov/api"
ROOT = Path(__file__).parent
RAW = ROOT / "data" / "raw"
WINDOW = 5                      # API maximum day range
SAMPLE_DAYS = (1, 11, 21)       # each month is sampled on days 1-5, 11-15 and 21-25 (~50%)
SAMPLING = "15 days per month (the 1st-5th, 11th-15th and 21st-25th, about half of every month)"


def in_sample(d):
    return any(s <= d.day < s + WINDOW for s in SAMPLE_DAYS)


def snippet_days(start, last):
    """Days a snippet starting on `start` covers (its window, clipped to the feed's last day)."""
    days = itertools.takewhile(lambda d: d <= last and in_sample(d), (start + timedelta(days=i) for i in range(WINDOW)))
    return list(days)
SENSORS = {                     # sensor -> its standard (SP) and near-real-time (NRT) feeds
    "MODIS": ["MODIS_SP", "MODIS_NRT"],
    "VIIRS": ["VIIRS_SNPP_SP", "VIIRS_SNPP_NRT"],
}


def map_key():
    key = os.environ.get("FIRMS_MAP_KEY")
    env = ROOT / ".env"
    if not key and env.exists():
        for line in env.read_text().splitlines():
            if line.startswith("FIRMS_MAP_KEY="):
                key = line.split("=", 1)[1].strip()
    if not key:
        sys.exit("Set FIRMS_MAP_KEY (env var or .env file)")
    return key


def get(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def availability(key):
    rows = csv.DictReader(io.StringIO(get(f"{API}/data_availability/csv/{key}/all")))
    return {r["data_id"]: (date.fromisoformat(r["min_date"]), date.fromisoformat(r["max_date"])) for r in rows}


def windows(avail):
    """Three 5-day snippets per month (SAMPLE_DAYS) for every feed: (source, start, ndays)."""
    out = []
    for feeds in SENSORS.values():
        for src in feeds:
            lo, hi = avail[src]
            y, m = lo.year, lo.month
            while (y, m) <= (hi.year, hi.month):
                for sd in SAMPLE_DAYS:
                    start = max(lo, date(y, m, sd))
                    end = min(hi, date(y, m, sd + WINDOW - 1))
                    if start <= end:
                        out.append((src, start, (end - start).days + 1))
                y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


class Pacer:
    """Space request starts so we stay under the transaction limit."""
    def __init__(self, interval):
        self.interval, self.next, self.lock = interval, 0.0, threading.Lock()

    def wait(self):
        with self.lock:
            now = time.monotonic()
            start = max(now, self.next)
            self.next = start + self.interval
        time.sleep(max(0, start - now))


def fetch(key, bbox, pacer, w):
    src, start, n = w
    path = RAW / src / f"{start}.csv.gz"
    url = f"{API}/area/csv/{key}/{src}/{bbox}/{n}/{start}"
    for attempt in range(8):
        pacer.wait()
        try:
            body = get(url)
        except Exception as e:
            body = f"ERR {e}"
        if body.startswith("latitude"):
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".tmp")
            with gzip.open(tmp, "wt") as f:
                f.write(body)
            tmp.rename(path)
            return body.count("\n") - 1
        # rate limited or transient error: back off and retry
        time.sleep(min(60, 5 * 2 ** attempt))
    print(f"\nFAILED {src} {start}: {body[:120]!r}", file=sys.stderr)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", default="-125,32,-114,42.5", help="west,south,east,north (default: California)")
    ap.add_argument("--name", default="California")
    ap.add_argument("--interval", type=float, default=1.35, help="seconds between requests (1.35 = ~4400 transactions/10min)")
    a = ap.parse_args()

    key = map_key()
    meta_path = ROOT / "data" / "meta.json"
    if meta_path.exists() and json.loads(meta_path.read_text())["bbox"] != a.bbox:
        sys.exit(f"data/ holds a different bbox ({json.loads(meta_path.read_text())['bbox']}); move or delete data/ first")
    RAW.mkdir(parents=True, exist_ok=True)
    meta_path.write_text(json.dumps({"bbox": a.bbox, "name": a.name}))

    todo = [w for w in windows(availability(key)) if not (RAW / w[0] / f"{w[1]}.csv.gz").exists()]
    print(f"{len(todo)} windows to fetch (~{len(todo) * a.interval / 60:.0f} min)", flush=True)
    pacer, done, rows, t0 = Pacer(a.interval), 0, 0, time.time()
    with ThreadPoolExecutor(6) as pool:
        for n in pool.map(lambda w: fetch(key, a.bbox, pacer, w), todo):
            done += 1
            rows += n or 0
            if done % 25 == 0 or done == len(todo):
                print(f"{done}/{len(todo)} windows, {rows:,} detections, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
