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
import ast, csv, json, re, sys, time, urllib.error, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).parent
VAL = ROOT / "viirs_from_modis" / "out" / "validation"
SERIES = ["MODIS_observed_2001_2011", "MODIS_observed_2012", "VIIRS_observed_2012"]
MIN_ACRES = 100
WIKI_MIN_ACRES = 5000
CACHE = ROOT / "data" / "wiki_cache.json"
UA = {"User-Agent": "nasa-spaceapps-fire-explorer/1.0 (student project)"}
STATES = {"CA": "California", "OR": "Oregon", "NV": "Nevada", "AZ": "Arizona", "ID": "Idaho", "UT": "Utah"}
ACRONYMS = {"SCU", "LNU", "CZU", "SQF", "TP", "MDF", "DH", "WNA", "XSC", "MM", "BTU", "SHU", "TCU", "HWY"}
# FPA FOD files some complexes under one component fire. (record name, year) -> (title shown, Wikipedia article)
ALIASES = {
    ("DOE", "2020"): ("August Complex (Doe Fire)", "August Complex fire"),
    ("CLAREMONT", "2020"): ("North Complex (Claremont Fire)", "North Complex fire"),
    ("RANCH", "2018"): ("Mendocino Complex (Ranch Fire)", "Mendocino Complex Fire"),
    ("OBSERVATIO", "2001"): ("Observation Fire", None),   # truncated in the source
}


def clean_name(raw):
    """FPA FOD names are upper case and sometimes codes: 'SCU LIGHTNING COMPLEX' -> 'SCU Lightning Complex',
    'MCNALLY' -> 'McNally'; codes, blanks and 'UNNAMED ...' -> None."""
    s = " ".join(raw.split())
    if not s or s.upper().startswith("UNNAMED") or re.search(r"\d{5,}", s):
        return None
    def word(w):
        if w.upper() in ACRONYMS:
            return w.upper()
        w = "-".join(x.capitalize() for x in w.split("-"))
        return "Mc" + w[2:].capitalize() if len(w) > 3 and w.startswith("Mc") else w
    return " ".join(word(w) for w in s.split(" "))


def fire_title(name, raw, year, state):
    if (raw.upper(), year) in ALIASES:
        t = ALIASES[(raw.upper(), year)][0]
    elif name is None:
        t = "Unnamed fire"
    else:
        t = name if name.lower().endswith("complex") else name + " Fire"
    return t if state == "CA" else f"{t} ({STATES.get(state, state)})"


def exact_article(f, raw):
    """Look for an article titled exactly after the fire and accept it only if its opening paragraph
    mentions the fire's year and state (so 'Creek Fire' for a 2008 fire doesn't get the 2020 article)."""
    alias = ALIASES.get((raw.upper(), f["start"][:4]))
    titles = [alias[1]] if alias and alias[1] else [] if alias else \
        [f"{f['name']} Fire ({f['start'][:4]})", f"{f['name']} Fire", f"{f['name']} fire"] if f["name"] else []
    if not titles:
        return None, None
    q = urllib.parse.urlencode({"action": "query", "format": "json", "redirects": 1, "prop": "extracts|pageprops|info", "inprop": "url",
                                "exintro": 1, "explaintext": 1, "exlimit": "max", "titles": "|".join(titles)})
    req = urllib.request.Request("https://en.wikipedia.org/w/api.php?" + q, headers=UA)
    for attempt in range(5):                      # back off when Wikipedia rate-limits (HTTP 429)
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                pages = json.load(r)["query"].get("pages", {})
            break
        except urllib.error.HTTPError as e:
            if e.code != 429 or attempt == 4:
                raise
            time.sleep(15 * (attempt + 1))
    for pg in pages.values():
        intro = (pg.get("extract") or "")[:600]
        if "missing" in pg or "disambiguation" in pg.get("pageprops", {}):
            continue
        if f["start"][:4] in intro and STATES.get(f["state"], "") in intro and "fire" in intro.lower():
            return pg["title"], pg["fullurl"]
    return None, None


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
            f = fires.setdefault(k, {"raw": r["fire_name"].strip(), "name": clean_name(r["fire_name"]), "start": r["fire_start"], "end": r["fire_end"] or r["fire_start"],
                                     "acres": round(acres), "cause": r["fire_cause"], "state": r["fire_state"],
                                     "lat": round(float(r["event_lat"]), 3), "lon": round(float(r["event_lon"]), 3),
                                     "km": float(r["distance_km"] or 0), "seen_by": []})
            sensor = s.split("_")[0]
            if sensor not in f["seen_by"]:
                f["seen_by"].append(sensor)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    if "--no-wiki" not in sys.argv:
        lookup = load_wiki_lookup()
        todo = [f for f in fires.values() if f["acres"] >= WIKI_MIN_ACRES and f["name"] and f"{f['name']}|{f['start'][:4]}" not in cache]
        for i, f in enumerate(sorted(todo, key=lambda f: -f["acres"])):
            title, url = lookup(f["name"], f["start"][:4])
            cache[f"{f['name']}|{f['start'][:4]}"] = [title, url]
            if i % 20 == 19:
                CACHE.parent.mkdir(exist_ok=True); CACHE.write_text(json.dumps(cache)); print(f"  {i + 1}/{len(todo)} looked up")
        CACHE.parent.mkdir(exist_ok=True); CACHE.write_text(json.dumps(cache))
        # second pass: exact article titles (and the complex aliases) for large fires the search didn't settle
        todo = [f for f in fires.values() if f["acres"] >= WIKI_MIN_ACRES and f"exact|{f['raw'].upper()}|{f['start'][:4]}" not in cache]
        for i, f in enumerate(sorted(todo, key=lambda f: -f["acres"])):
            try:
                cache[f"exact|{f['raw'].upper()}|{f['start'][:4]}"] = list(exact_article(f, f["raw"]))
            except Exception as e:
                print(f"  {f['raw']}: {e}"); continue
            time.sleep(2)
            if i % 20 == 19:
                CACHE.write_text(json.dumps(cache)); print(f"  exact titles: {i + 1}/{len(todo)}")
        CACHE.write_text(json.dumps(cache))
    out = []
    for f in sorted(fires.values(), key=lambda f: (f["start"], -f["acres"])):
        year = f["start"][:4]
        title, url = cache.get(f"exact|{f['raw'].upper()}|{year}", [None, None])
        if not title and f["name"]:
            title, url = cache.get(f"{f['name']}|{year}", [None, None]) or [None, None]
            if title and not names_this_fire(title, f["name"]):
                title = url = None
        raw = f.pop("raw")
        out.append({**f, "name": f["name"] or "Unnamed", "title": fire_title(f["name"], raw, year, f["state"]), "wiki_title": title, "wiki_url": url})
    meta = {"source": "FPA FOD (US Forest Service RDS-2013-0009.6, 1992-2020), matched to observed MODIS/VIIRS spikes by fire_spike_matcher.py",
            "location": "lat/lon is the centre of the 0.25° cell where the detection spike was; the fire itself is within buffer km of it",
            "wikipedia": f"article found by fire_spike_matcher.wiki_lookup for fires of {WIKI_MIN_ACRES:,}+ acres", "fires": out}
    (ROOT / "app" / "data" / "fires.json").write_text(json.dumps(meta, separators=(",", ":")))
    print(f"wrote app/data/fires.json: {len(out)} fires, {sum(1 for f in out if f['wiki_url'])} with a Wikipedia article")


if __name__ == "__main__":
    main()
