#!/usr/bin/env python3
"""Named fires for the web app: the recorded wildfires that Vaangmaya's
viirs_from_modis/fire_spike_matcher.py matched to spikes in the *observed* detections,
plus a Wikipedia article for the large ones.

Input : viirs_from_modis/out/validation/matches_{MODIS,VIIRS}_observed_*.csv
        (best match per spike event; FPA FOD, US Forest Service, 1992-2020)
Output: app/data/fires.json
Wikipedia: uses wiki_lookup() from fire_spike_matcher.py (its own query, title check and
rate limit), loaded on its own so pandas/numpy aren't needed. Results are cached in
data/wiki_cache.json, so re-runs make no new requests.

Run: python3 known_fires.py            (about 1 s per uncached fire >= WIKI_MIN_ACRES)
     python3 known_fires.py --no-wiki
"""
import ast, csv, json, sys, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).parent
VAL = ROOT / "viirs_from_modis" / "out" / "validation"
SERIES = ["MODIS_observed_2001_2011", "MODIS_observed_2012", "VIIRS_observed_2012"]
MIN_ACRES = 100
WIKI_MIN_ACRES = 5000
CACHE = ROOT / "data" / "wiki_cache.json"


def load_wiki_lookup():
    """Pull wiki_lookup out of fire_spike_matcher.py without importing pandas/numpy."""
    src = (ROOT / "viirs_from_modis" / "fire_spike_matcher.py").read_text()
    fn = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "wiki_lookup")
    ns = {"json": json, "time": time, "urllib": urllib}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), "fire_spike_matcher.py", "exec"), ns)
    return ns["wiki_lookup"]


def names_this_fire(title, name):
    """Stricter than wiki_lookup's first-word check, which accepts e.g. 'Kelly Rowland' for the
    Rowland fire: the title must be '<name> Fire...' or '<name> Complex...' (capitalised, as
    Wikipedia names fire articles), so generic pages like 'Oil well fire' are rejected too."""
    t = title.split(" (")[0]
    rest = t[len(name):]
    return t.lower().startswith(name.lower()) and (rest.startswith((" Fire", " Complex")) or
                                                   (name.lower().endswith("complex") and (not rest or rest.lower().startswith(" fire"))))


def main():
    fires = {}
    for s in SERIES:
        for r in csv.DictReader(open(VAL / f"matches_{s}.csv")):
            if r["status"] != "matched" or r["rank"] not in ("1", "1.0"):
                continue
            acres = float(r["fire_size_acres"] or 0)
            if acres < MIN_ACRES:
                continue
            k = (r["fire_name"].strip().upper(), r["fire_start"])
            f = fires.setdefault(k, {"name": r["fire_name"].strip().title(), "start": r["fire_start"], "end": r["fire_end"] or r["fire_start"],
                                     "acres": round(acres), "cause": r["fire_cause"], "state": r["fire_state"],
                                     "lat": round(float(r["event_lat"]), 3), "lon": round(float(r["event_lon"]), 3),
                                     "km": float(r["distance_km"] or 0), "seen_by": []})
            sensor = s.split("_")[0]
            if sensor not in f["seen_by"]:
                f["seen_by"].append(sensor)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    if "--no-wiki" not in sys.argv:
        lookup = load_wiki_lookup()
        todo = [f for f in fires.values() if f["acres"] >= WIKI_MIN_ACRES and f"{f['name']}|{f['start'][:4]}" not in cache]
        for i, f in enumerate(sorted(todo, key=lambda f: -f["acres"])):
            title, url = lookup(f["name"], f["start"][:4])
            cache[f"{f['name']}|{f['start'][:4]}"] = [title, url]
            if i % 20 == 19:
                CACHE.parent.mkdir(exist_ok=True); CACHE.write_text(json.dumps(cache)); print(f"  {i + 1}/{len(todo)} looked up")
        CACHE.parent.mkdir(exist_ok=True); CACHE.write_text(json.dumps(cache))
    out = []
    for f in sorted(fires.values(), key=lambda f: (f["start"], -f["acres"])):
        title, url = cache.get(f"{f['name']}|{f['start'][:4]}", [None, None])
        if title and not names_this_fire(title, f["name"]):
            title = url = None
        out.append({**f, "wiki_title": title, "wiki_url": url})
    meta = {"source": "FPA FOD (US Forest Service RDS-2013-0009.6, 1992-2020), matched to observed MODIS/VIIRS spikes by fire_spike_matcher.py",
            "location": "lat/lon is the centre of the 0.25° cell where the detection spike was; the fire itself is within buffer km of it",
            "wikipedia": f"article found by fire_spike_matcher.wiki_lookup for fires of {WIKI_MIN_ACRES:,}+ acres", "fires": out}
    (ROOT / "app" / "data" / "fires.json").write_text(json.dumps(meta, separators=(",", ":")))
    print(f"wrote app/data/fires.json: {len(out)} fires, {sum(1 for f in out if f['wiki_url'])} with a Wikipedia article")


if __name__ == "__main__":
    main()
