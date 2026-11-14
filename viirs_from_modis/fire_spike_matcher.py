#!/usr/bin/env python3
"""
fire_spike_matcher.py
Find spikes in MODIS / VIIRS active-fire detections and match them to known fires.

PIPELINE
  1. Load FIRMS archive CSVs (MODIS and/or VIIRS), filter low-confidence rows.
  2. Count detections per grid cell, per day, per sensor.
  3. Flag "spike" cell-days: counts far above that cell's normal for that calendar month.
  4. Merge neighbouring spike cell-days (space + time) into spike EVENTS.
  5. Match each event to fire records (FPA FOD SQLite, or any CSV of fires) by place + date.
  6. Optionally look up a Wikipedia article for each matched fire (--wiki).

SECOND MODE: --cells-from fires.html  (validate the MODIS->VIIRS back-cast against incidents)
  The repo only holds the 5-days-a-month sample as per-year 0.05-degree cell counts (not daily
  rows), so spikes are found per CELL-YEAR instead of per cell-day, for three series:
     MODIS_observed     what MODIS actually detected
     VIIRS_observed     what VIIRS actually detected (2012+)
     VIIRS_predicted    the model's VIIRS estimate from MODIS: out-of-year predictions for 2012+
                        (out/cells_predicted_loyo.json), the back-cast for 2001-2011
                        (out/cells_predicted.json)
  Each series is compared with FPA FOD incidents, counting only the part of a fire that burned
  on the sampled days (the 13th-17th of each month), because that is all the detections cover.
  Reported: Spearman correlation per cell-year (zeros included), how often spikes land on a real
  fire vs. chance, and how many large fires were caught. See validation.json / report printed.

INPUTS
  FIRMS CSVs : https://firms.modaps.eosdis.nasa.gov/download/  (country = United States)
               Expected columns: latitude, longitude, acq_date, confidence, frp,
               and either `instrument` or bright_ti4 (used to tell MODIS from VIIRS).
  Fire records, one of:
    --fod  FPA_FOD.sqlite  (Forest Service RDS-2013-0009, table "Fires")
    --fires-csv file.csv   with columns: name, start_date, lat, lon
                           and optionally end_date, size_acres, cause, state

OUTPUTS (in --out directory)
  spikes.csv   one row per spike cell-day per sensor
  events.csv   one row per merged spike event (bbox, dates, counts per sensor)
  matches.csv  candidate fires per event, or status=unmatched
  (cell-year mode: spikes_<series>.csv, matches_<series>.csv, validation.json)

EXAMPLES
  python fire_spike_matcher.py --demo
  python fire_spike_matcher.py --firms modis_usa.csv viirs_usa.csv --fod FPA_FOD.sqlite
  python fire_spike_matcher.py --firms firms/*.csv --fires-csv my_fires.csv --wiki
  python fire_spike_matcher.py --cells-from ../fires.html --fod FPA_FOD.sqlite --out out/validation

Requires: pandas, numpy. Everything else is standard library.
"""
import argparse
import json
import re
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

KM_PER_DEG = 111.0


# ----------------------------------------------------------------------------
# 1. Load detections
# ----------------------------------------------------------------------------
def load_firms(paths, modis_min_conf=30, keep_low_viirs=False):
    """Read FIRMS CSVs into one frame: lat, lon, date, sensor, frp."""
    frames = []
    for p in paths:
        df = pd.read_csv(p, dtype={"confidence": str}, low_memory=False)
        df.columns = [c.lower() for c in df.columns]
        if "instrument" in df.columns:
            sensor = df["instrument"].astype(str).str.upper().str.strip()
        elif "bright_ti4" in df.columns:
            sensor = pd.Series("VIIRS", index=df.index)
        else:
            sensor = pd.Series("MODIS", index=df.index)
        df["sensor"] = np.where(sensor.str.contains("VIIRS"), "VIIRS", "MODIS")

        # FIRMS archive `type`: 0 = presumed vegetation fire; 1-3 = volcano,
        # static land source (e.g. gas flare), offshore. Keep only 0 when present.
        if "type" in df.columns:
            df = df[pd.to_numeric(df["type"], errors="coerce").fillna(0) == 0]

        # Confidence: MODIS is 0-100, VIIRS is l/n/h (or low/nominal/high).
        conf = df["confidence"].astype(str).str.lower().str.strip()
        num = pd.to_numeric(conf, errors="coerce")
        is_modis = df["sensor"] == "MODIS"
        keep = np.where(
            is_modis,
            num.fillna(100) >= modis_min_conf,
            True if keep_low_viirs else ~conf.str.startswith("l"),
        )
        df = df[keep]

        frames.append(pd.DataFrame({
            "lat": pd.to_numeric(df["latitude"], errors="coerce"),
            "lon": pd.to_numeric(df["longitude"], errors="coerce"),
            "date": pd.to_datetime(df["acq_date"], errors="coerce").dt.normalize(),
            "sensor": df["sensor"].values,
            "frp": pd.to_numeric(df.get("frp", np.nan), errors="coerce"),
        }))
    out = pd.concat(frames, ignore_index=True).dropna(subset=["lat", "lon", "date"])
    if out.empty:
        sys.exit("No detections left after filtering. Check the input files.")
    return out


# ----------------------------------------------------------------------------
# 2-3. Grid, daily counts, spike detection
# ----------------------------------------------------------------------------
def find_spikes(det, cell_deg=0.25, z_thresh=4.0, min_count=10):
    """
    A cell-day is a spike when its count is >= min_count AND its z-score is
    >= z_thresh, measured against that cell's own history for the same calendar
    month (days with no detections count as zero). Done separately per sensor,
    because VIIRS sees several times more fire pixels than MODIS.
    """
    det = det.copy()
    det["ix"] = np.floor(det["lon"] / cell_deg).astype(int)
    det["iy"] = np.floor(det["lat"] / cell_deg).astype(int)
    daily = (det.groupby(["sensor", "ix", "iy", "date"])
                .agg(n=("lat", "size"), frp=("frp", "sum"))
                .reset_index())
    daily["month"] = daily["date"].dt.month
    daily["n2"] = daily["n"] ** 2

    # Number of calendar days per month in each sensor's record (the zeros).
    span = det.groupby("sensor")["date"].agg(["min", "max"])
    ndays = []
    for sensor, (dmin, dmax) in span.iterrows():
        months = pd.date_range(dmin, dmax, freq="D").month
        vc = pd.Series(months).value_counts()
        ndays.append(pd.DataFrame({"sensor": sensor, "month": vc.index, "N": vc.values}))
    ndays = pd.concat(ndays)

    stats = (daily.groupby(["sensor", "ix", "iy", "month"])
                  .agg(s=("n", "sum"), ss=("n2", "sum")).reset_index()
                  .merge(ndays, on=["sensor", "month"]))
    stats["mean"] = stats["s"] / stats["N"]
    var = stats["ss"] / stats["N"] - stats["mean"] ** 2
    stats["std"] = np.sqrt(np.clip(var, 1e-9, None))

    daily = daily.merge(stats[["sensor", "ix", "iy", "month", "mean", "std"]],
                        on=["sensor", "ix", "iy", "month"])
    daily["z"] = (daily["n"] - daily["mean"]) / daily["std"]
    spikes = daily[(daily["n"] >= min_count) & (daily["z"] >= z_thresh)].copy()
    spikes["cell_lat"] = (spikes["iy"] + 0.5) * cell_deg
    spikes["cell_lon"] = (spikes["ix"] + 0.5) * cell_deg
    return spikes.drop(columns=["n2", "month"]).rename(
        columns={"n": "detections", "mean": "baseline_mean"})


# ----------------------------------------------------------------------------
# 4. Merge spike cell-days into events (connected components, 3D neighbours)
# ----------------------------------------------------------------------------
def build_events(spikes, cell_deg=0.25):
    if spikes.empty:
        return pd.DataFrame(), spikes.assign(event_id=pd.Series(dtype=int))
    day0 = spikes["date"].min()
    spikes = spikes.copy()
    spikes["t"] = (spikes["date"] - day0).dt.days
    keys = list({(r.ix, r.iy, r.t) for r in spikes.itertuples()})
    parent = {k: k for k in keys}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    offsets = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    for k in keys:
        for a, b, c in offsets:
            nb = (k[0] + a, k[1] + b, k[2] + c)
            if nb in parent:
                ra, rb = find(k), find(nb)
                if ra != rb:
                    parent[ra] = rb
    roots = {}
    spikes["event_id"] = [
        roots.setdefault(find((r.ix, r.iy, r.t)), len(roots) + 1)
        for r in spikes.itertuples()
    ]

    rows = []
    for eid, g in spikes.groupby("event_id"):
        by_day = g.groupby("date")["detections"].sum()
        per_sensor = g.groupby("sensor")["detections"].sum()
        rows.append({
            "event_id": eid,
            "start_date": g["date"].min().date(),
            "end_date": g["date"].max().date(),
            "peak_date": by_day.idxmax().date(),
            "lat_min": g["iy"].min() * cell_deg,
            "lat_max": (g["iy"].max() + 1) * cell_deg,
            "lon_min": g["ix"].min() * cell_deg,
            "lon_max": (g["ix"].max() + 1) * cell_deg,
            "modis_detections": int(per_sensor.get("MODIS", 0)),
            "viirs_detections": int(per_sensor.get("VIIRS", 0)),
            "sensors": "+".join(sorted(per_sensor.index)),
            "max_z": round(float(g["z"].max()), 1),
            "frp_sum_mw": round(float(g["frp"].sum()), 1),
        })
    events = pd.DataFrame(rows)
    events["lat"] = (events["lat_min"] + events["lat_max"]) / 2
    events["lon"] = (events["lon_min"] + events["lon_max"]) / 2
    return events, spikes.drop(columns=["t"])


# ----------------------------------------------------------------------------
# 5. Fire records and matching
# ----------------------------------------------------------------------------
def _fod_date(series):
    """FPA FOD dates are Julian day numbers in older editions, text in newer."""
    num = pd.to_numeric(series, errors="coerce")
    if num.notna().mean() > 0.9 and num.median() > 2_000_000:
        return pd.to_datetime(num - 2440587.5, unit="D", errors="coerce").dt.normalize()
    return pd.to_datetime(series, errors="coerce").dt.normalize()


def load_fod(path, min_acres=10.0, bbox=None):
    """bbox = (west, south, east, north) limits the query to a region (much faster)."""
    con = sqlite3.connect(path)
    cols = {r[1].upper() for r in con.execute("PRAGMA table_info(Fires)")}
    if not cols:
        sys.exit(f'No table named "Fires" in {path}. Is this the FPA FOD SQLite file?')
    cause = next((c for c in ("NWCG_GENERAL_CAUSE", "STAT_CAUSE_DESCR") if c in cols), None)
    want = ["FIRE_NAME", "DISCOVERY_DATE", "CONT_DATE", "FIRE_SIZE",
            "LATITUDE", "LONGITUDE", "STATE"] + ([cause] if cause else [])
    have = [c for c in want if c in cols]
    sql, params = f"SELECT {', '.join(have)} FROM Fires WHERE FIRE_SIZE >= ?", [min_acres]
    if bbox:
        w, s, e, n = bbox
        sql += " AND LONGITUDE BETWEEN ? AND ? AND LATITUDE BETWEEN ? AND ?"
        params += [w, e, s, n]
    df = pd.read_sql_query(sql, con, params=params)
    con.close()
    return pd.DataFrame({
        "name": df.get("FIRE_NAME"),
        "start_date": _fod_date(df["DISCOVERY_DATE"]),
        "end_date": _fod_date(df["CONT_DATE"]) if "CONT_DATE" in df else pd.NaT,
        "lat": df["LATITUDE"], "lon": df["LONGITUDE"],
        "size_acres": df["FIRE_SIZE"],
        "cause": df[cause] if cause else None,
        "state": df.get("STATE"),
    })


def load_fires_csv(path, min_acres=0.0):
    df = pd.read_csv(path)
    df.columns = [c.lower() for c in df.columns]
    for need in ("name", "start_date", "lat", "lon"):
        if need not in df.columns:
            sys.exit(f"--fires-csv is missing required column: {need}")
    df["start_date"] = pd.to_datetime(df["start_date"], errors="coerce").dt.normalize()
    df["end_date"] = pd.to_datetime(df.get("end_date"), errors="coerce").dt.normalize()
    for opt in ("size_acres", "cause", "state"):
        if opt not in df.columns:
            df[opt] = np.nan
    if min_acres:
        df = df[df["size_acres"].fillna(0) >= min_acres]
    return df[["name", "start_date", "end_date", "lat", "lon", "size_acres", "cause", "state"]]


def match_events(events, fires, buffer_km=15.0, pad_days=2, default_days=7, top_k=3):
    """A fire is a candidate if it lies in the event's (buffered) box and was
    burning during the event's (padded) date window. Largest fires rank first."""
    fires = fires.dropna(subset=["start_date", "lat", "lon"]).reset_index(drop=True)
    f_start = fires["start_date"].values
    f_end = fires["end_date"].fillna(
        fires["start_date"] + pd.Timedelta(days=default_days)).values
    f_lat, f_lon = fires["lat"].values, fires["lon"].values
    pad = np.timedelta64(pad_days, "D")
    rows = []
    for ev in events.itertuples():
        e0, e1 = np.datetime64(ev.start_date), np.datetime64(ev.end_date)
        dlat = buffer_km / KM_PER_DEG
        dlon = dlat / max(np.cos(np.radians(ev.lat)), 0.2)
        m = ((f_start <= e1 + pad) & (f_end >= e0 - pad) &
             (f_lat >= ev.lat_min - dlat) & (f_lat <= ev.lat_max + dlat) &
             (f_lon >= ev.lon_min - dlon) & (f_lon <= ev.lon_max + dlon))
        cand = fires[m].sort_values("size_acres", ascending=False).head(top_k)
        base = {"event_id": ev.event_id, "event_start": ev.start_date,
                "event_end": ev.end_date, "event_lat": round(ev.lat, 3),
                "event_lon": round(ev.lon, 3), "sensors": ev.sensors}
        if cand.empty:
            rows.append({**base, "status": "unmatched", "rank": None})
            continue
        for rank, f in enumerate(cand.itertuples(), 1):
            dist = np.hypot((f.lat - ev.lat) * KM_PER_DEG,
                            (f.lon - ev.lon) * KM_PER_DEG * np.cos(np.radians(ev.lat)))
            rows.append({**base, "status": "matched", "rank": rank,
                         "fire_name": f.name, "fire_start": f.start_date.date(),
                         "fire_end": None if pd.isna(f.end_date) else f.end_date.date(),
                         "fire_size_acres": f.size_acres, "fire_cause": f.cause,
                         "fire_state": f.state, "distance_km": round(float(dist), 1)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# 6. Optional Wikipedia lookup (official API, one request per fire, rate limited)
# ----------------------------------------------------------------------------
def wiki_lookup(name, year, pause=1.0):
    if not isinstance(name, str) or not name.strip():
        return None, None
    q = f"{name.title()} Fire {year}"
    url = ("https://en.wikipedia.org/w/api.php?action=query&list=search&format=json"
           "&srlimit=1&srsearch=" + urllib.parse.quote(q))
    req = urllib.request.Request(
        url, headers={"User-Agent": "fire-spike-matcher/1.0 (hackathon research script)"})
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            hits = json.load(r)["query"]["search"]
    except Exception:
        return None, None
    finally:
        time.sleep(pause)
    if not hits:
        return None, None
    title = hits[0]["title"]
    # Only accept the hit if it plausibly names this fire.
    if name.split()[0].lower() not in title.lower():
        return None, None
    return title, "https://en.wikipedia.org/wiki/" + urllib.parse.quote(title.replace(" ", "_"))


def add_wiki(matches, limit=50):
    top = matches[(matches["status"] == "matched") & (matches["rank"] == 1)].head(limit)
    titles, urls = {}, {}
    for r in top.itertuples():
        titles[r.Index], urls[r.Index] = wiki_lookup(
            r.fire_name, pd.Timestamp(r.fire_start).year)
    matches["wiki_title"] = pd.Series(titles)
    matches["wiki_url"] = pd.Series(urls)
    return matches


# ----------------------------------------------------------------------------
# CELL-YEAR MODE: does predicted VIIRS line up with real incidents?
# ----------------------------------------------------------------------------
SAMPLE_FIRST_DAY, SAMPLE_LAST_DAY = 13, 17      # fetch.py samples days 13-17 of each month


def load_cell_series(html_path, pred_json, loyo_json, cell_deg=0.25):
    """Per-year detection counts on a coarse grid for MODIS_observed, VIIRS_observed and
    VIIRS_predicted. Columns: series, year, cx, cy (cell index from the bbox corner), n."""
    h = Path(html_path).read_text(encoding="utf8")
    D = json.loads(re.search(r"const D = (\{.*?\});\n", h, re.S).group(1))
    k = max(1, round(cell_deg / D["cell"]))
    modis_days = {int(y): v["days"] for y, v in D["sensors"]["MODIS"]["years"].items()}
    rows = []

    def add(series, year, flat, scale):
        a = np.asarray(flat, float).reshape(-1, 3)
        df = pd.DataFrame({"cx": (a[:, 0] // k).astype(int), "cy": (a[:, 1] // k).astype(int), "n": a[:, 2] * scale})
        df = df.groupby(["cx", "cy"], as_index=False)["n"].sum()
        df.insert(0, "year", year); df.insert(0, "series", series)
        rows.append(df)

    for sensor in ("MODIS", "VIIRS"):
        for y, v in D["sensors"][sensor]["years"].items():
            if v["days"] and v["cells"]:
                add(f"{sensor}_observed", int(y), v["cells"], 1.0)
    first_viirs = min(int(y) for y in D["sensors"]["VIIRS"]["years"])
    for path, keep in ((pred_json, lambda y: y < first_viirs), (loyo_json, lambda y: y >= first_viirs)):
        P = json.loads(Path(path).read_text())
        for y, flat in P["years"].items():
            y = int(y)
            if keep(y) and flat and modis_days.get(y):
                add("VIIRS_predicted", y, flat, modis_days[y])     # rate/day -> detections in the sample
    out = pd.concat(rows, ignore_index=True)
    nx, ny = int(np.ceil((D["bbox"][2] - D["bbox"][0]) / cell_deg)), int(np.ceil((D["bbox"][3] - D["bbox"][1]) / cell_deg))
    return out, D["bbox"], (nx, ny), first_viirs


def sample_days_burning(start, end):
    """How many of a fire's days fall on the sampled days (13th-17th), and its length in days."""
    end = start if pd.isna(end) or end < start else end
    days = pd.date_range(start, end, freq="D")
    return int(((days.day >= SAMPLE_FIRST_DAY) & (days.day <= SAMPLE_LAST_DAY)).sum()), len(days)


def fires_in_sample(fires, bbox, cell_deg, years):
    """Keep fires that burned on at least one sampled day; add the cell, the year, and
    acres_in_sample = acres x (share of the fire's days that were sampled days)."""
    f = fires.dropna(subset=["start_date", "lat", "lon"]).copy()
    f = f[f["start_date"].dt.year.isin(years)]
    ov = [sample_days_burning(s, e) for s, e in zip(f["start_date"], f["end_date"])]
    f["sample_days"] = [a for a, _ in ov]; f["fire_days"] = [b for _, b in ov]
    f = f[f["sample_days"] > 0].copy()
    f["acres_in_sample"] = f["size_acres"] * f["sample_days"] / f["fire_days"]
    f["year"] = f["start_date"].dt.year
    f["cx"] = np.floor((f["lon"] - bbox[0]) / cell_deg).astype(int)
    f["cy"] = np.floor((f["lat"] - bbox[1]) / cell_deg).astype(int)
    return f


def find_spikes_cellyear(df, years, z_thresh=4.0, min_count=10.0):
    """Cell-year version of find_spikes. z compares a cell-year with that cell's OTHER years
    (leave-one-out mean and std, absent years = 0, std floored at 1 detection), so one huge
    year can't hide itself by inflating its own baseline."""
    d = df[df["year"].isin(years)]
    N = len(years)
    if d.empty or N < 3:
        return d.iloc[0:0].assign(z=[], baseline_mean=[])
    piv = d.pivot_table(index=["cx", "cy"], columns="year", values="n", fill_value=0.0).reindex(columns=years, fill_value=0.0)
    X = piv.values
    S, SS = X.sum(1, keepdims=True), (X ** 2).sum(1, keepdims=True)
    mean = (S - X) / (N - 1)
    std = np.sqrt(np.clip((SS - X ** 2) / (N - 1) - mean ** 2, 1.0, None))
    Z = (X - mean) / std
    ci, yi = np.nonzero((X >= min_count) & (Z >= z_thresh))
    idx = piv.index.to_frame(index=False)
    out = pd.DataFrame({"cx": idx["cx"].values[ci], "cy": idx["cy"].values[ci], "year": np.asarray(years)[yi],
                        "detections": X[ci, yi], "z": Z[ci, yi], "baseline_mean": mean[ci, yi]})
    return out


def spikes_to_events(sp, bbox, cell_deg):
    """Reuse build_events(): put each year on Jan 1 so only same-year neighbours merge, then
    give each event the whole year as its date window."""
    if sp.empty:
        return pd.DataFrame(), sp
    ox, oy = round(bbox[0] / cell_deg), round(bbox[1] / cell_deg)       # absolute cell index of the bbox corner
    s = pd.DataFrame({"ix": sp["cx"] + ox, "iy": sp["cy"] + oy, "date": pd.to_datetime(sp["year"].astype(str) + "-01-01"),
                      "sensor": "X", "detections": sp["detections"], "z": sp["z"], "frp": np.nan})
    ev, _ = build_events(s, cell_deg)
    ev["start_date"] = [pd.Timestamp(d).replace(month=1, day=1).date() for d in ev["start_date"]]
    ev["end_date"] = [pd.Timestamp(d).replace(month=12, day=31).date() for d in ev["end_date"]]
    return ev, sp


def spearman(x, y):
    return float(pd.Series(x).rank().corr(pd.Series(y).rank()))


def validate_series(name, df, fires_s, years, grid, a):
    """All numbers for one series over the given years."""
    nx, ny = grid
    full = pd.MultiIndex.from_product([years, range(nx), range(ny)], names=["year", "cx", "cy"])
    n = df[df["year"].isin(years)].groupby(["year", "cx", "cy"])["n"].sum().reindex(full, fill_value=0.0)
    acres = fires_s[fires_s["year"].isin(years)].groupby(["year", "cx", "cy"])["acres_in_sample"].sum().reindex(full, fill_value=0.0)
    big = (acres >= a.big_acres)
    res = {"years": f"{years[0]}-{years[-1]}", "n_years": len(years), "cell_years": int(len(n)),
           "cell_years_with_detections": int((n > 0).sum()), "cell_years_with_big_fire": int(big.sum())}
    res["spearman_all_cells"] = spearman(n.values, acres.values)
    both = (n.values > 0) | (acres.values > 0)
    res["spearman_active_cells"] = spearman(n.values[both], acres.values[both])
    # year-block bootstrap: resample whole years so spatial/temporal clustering is respected
    rng = np.random.default_rng(0)
    by_year = {y: (n.xs(y, level="year").values, acres.xs(y, level="year").values) for y in years}
    boots = []
    for _ in range(a.boot):
        pick = rng.choice(years, len(years))
        boots.append(spearman(np.concatenate([by_year[y][0] for y in pick]), np.concatenate([by_year[y][1] for y in pick])))
    res["spearman_ci95"] = [float(np.nanpercentile(boots, 2.5)), float(np.nanpercentile(boots, 97.5))]
    # yearly totals, for the scatter on the page
    tot = pd.DataFrame({"n": n.groupby("year").sum(), "acres": acres.groupby("year").sum()})
    res["yearly_totals"] = {int(y): [float(r.n), float(r.acres)] for y, r in tot.iterrows()}
    res["spearman_yearly_totals"] = spearman(tot["n"], tot["acres"]) if len(tot) > 3 else None
    # spikes
    scale = a._scale.get(name, 1.0)
    sp = find_spikes_cellyear(df, list(years), a.z, max(3.0, a.min_count * scale))
    sp_idx = pd.MultiIndex.from_frame(sp[["year", "cx", "cy"]]) if len(sp) else pd.MultiIndex.from_arrays([[], [], []], names=["year", "cx", "cy"])
    hit = big.reindex(sp_idx, fill_value=False).values if len(sp) else np.array([], bool)
    base = float(big.mean())
    prec = float(hit.mean()) if len(hit) else None
    # chance level: same number of spikes per year, placed uniformly at random among that year's cells
    per_year = sp.groupby("year").size().to_dict()
    bigm = {y: big.xs(y, level="year").values for y in years}
    null = [sum(rng.choice(bigm[y], k, replace=False).sum() for y, k in per_year.items()) for _ in range(a.boot)] if len(sp) else [0]
    obs = int(hit.sum())
    res.update({"min_count_used": round(max(3.0, a.min_count * scale), 1), "spike_cell_years": int(len(sp)),
                "spikes_on_big_fire": obs, "precision": prec, "base_rate": base,
                "lift": (prec / base) if prec is not None and base else None,
                "chance_expected_hits": float(np.mean(null)), "perm_p": float((np.sum(np.asarray(null) >= obs) + 1) / (len(null) + 1)),
                "big_fire_cell_years": int(big.sum()), "recall_of_big_fire_cells": (obs / int(big.sum())) if big.sum() else None})
    huge = (acres >= a.huge_acres)
    if huge.sum():
        res["huge_fire_cell_years"] = int(huge.sum())
        res["recall_of_huge_fire_cells"] = float(huge.reindex(sp_idx, fill_value=False).values.sum() / huge.sum()) if len(sp) else 0.0
    return res, sp


def run_cell_year(a, out):
    series, bbox, grid, first_viirs = load_cell_series(a.cells_from, a.pred_json, a.loyo_json, a.cell_deg)
    lo, hi = a.years
    fires = load_fod(a.fod, a.min_acres, bbox) if a.fod else load_fires_csv(a.fires_csv, a.min_acres)
    fs = fires_in_sample(fires, bbox, a.cell_deg, list(range(lo, hi + 1)))
    print(f"Incidents: {len(fires):,} fires >= {a.min_acres:g} acres in the region; "
          f"{len(fs):,} burned on a sampled day ({lo}-{hi}). 'Big' = >= {a.big_acres:g} acres-in-sample in a cell-year.")
    # spike thresholds are in detections-in-sample; scale MODIS/VIIRS min_count to each series' volume
    ref = series[(series["series"] == "VIIRS_predicted") & series["year"].between(lo, hi)]["n"].sum()
    a._scale = {s: min(1.0, series[(series["series"] == s) & series["year"].between(lo, hi)]["n"].sum() / ref) for s in series["series"].unique()} if ref else {}
    plans = [("VIIRS_predicted", "VIIRS_predicted (back-cast, 2001-%d)" % (first_viirs - 1), range(max(lo, 2001), first_viirs)),
             ("VIIRS_predicted", "VIIRS_predicted (held out, %d+)" % first_viirs, range(first_viirs, hi + 1)),
             ("VIIRS_observed", "VIIRS_observed (%d+)" % first_viirs, range(first_viirs, hi + 1)),
             ("MODIS_observed", "MODIS_observed (2001-%d)" % (first_viirs - 1), range(max(lo, 2001), first_viirs)),
             ("MODIS_observed", "MODIS_observed (%d+)" % first_viirs, range(first_viirs, hi + 1))]
    report = {"settings": {"cell_deg": a.cell_deg, "z": a.z, "min_count": a.min_count, "big_acres": a.big_acres,
                           "huge_acres": a.huge_acres, "min_acres": a.min_acres, "buffer_km": a.buffer_km, "years": [lo, hi],
                           "sampled_days": f"{SAMPLE_FIRST_DAY}th-{SAMPLE_LAST_DAY}th of each month"}, "results": {}}
    for key, label, yrs in plans:
        yrs = [y for y in yrs if y in set(series[series["series"] == key]["year"])]
        if len(yrs) < 3:
            continue
        res, sp = validate_series(key, series[series["series"] == key], fs, yrs, grid, a)
        report["results"][label] = res
        ev, _ = spikes_to_events(sp, bbox, a.cell_deg)
        tag = re.sub(r"\W+", "_", label).strip("_")
        sp.to_csv(out / f"spikes_{tag}.csv", index=False)
        if len(ev):
            m = match_events(ev, fs.assign(end_date=fs["end_date"].fillna(fs["start_date"])), a.buffer_km, 0, top_k=a.top_k)
            m.to_csv(out / f"matches_{tag}.csv", index=False)
            res["events"] = int(len(ev)); res["events_matched_to_a_fire"] = int(m.loc[m["status"] == "matched", "event_id"].nunique())
    (out / "validation.json").write_text(json.dumps(report, indent=1, default=str))
    print(f"\n{'series':42s} {'rho(all)':>8s} {'95% CI':>15s} {'spikes':>6s} {'on fire':>7s} {'prec':>6s} {'base':>6s} {'lift':>5s} {'p':>6s} {'recall':>6s}")
    for label, r in report["results"].items():
        f = lambda v, p="{:.2f}": "-" if v is None else p.format(v)
        print(f"{label:42s} {r['spearman_all_cells']:8.2f} {'[%.2f, %.2f]' % tuple(r['spearman_ci95']):>15s} {r['spike_cell_years']:6d} "
              f"{r['spikes_on_big_fire']:7d} {f(r['precision']):>6s} {f(r['base_rate']):>6s} {f(r['lift'], '{:.1f}x'):>5s} {f(r['perm_p'], '{:.3f}'):>6s} {f(r['recall_of_big_fire_cells']):>6s}")
    print(f"\nWrote validation.json, spikes_*.csv, matches_*.csv to {out}/")


# ----------------------------------------------------------------------------
# Demo data (synthetic, so the script can be tried with no downloads)
# ----------------------------------------------------------------------------
def make_demo(outdir):
    rng = np.random.default_rng(7)
    days = pd.date_range("2019-01-01", "2021-12-31", freq="D")
    rows = []
    for sensor, rate in (("MODIS", 60), ("VIIRS", 220)):          # background scatter
        n = rate * len(days)
        rows.append(pd.DataFrame({
            "latitude": rng.uniform(32, 42, n), "longitude": rng.uniform(-124, -114, n),
            "acq_date": rng.choice(days, n), "instrument": sensor,
            "confidence": "80" if sensor == "MODIS" else "n",
            "frp": rng.gamma(2, 8, n)}))
    fires = [("CEDAR RIDGE", "2020-08-17", 12, 39.80, -121.40, 95000, "Natural"),
             ("DRY CREEK", "2021-07-04", 6, 34.60, -118.90, 14000, "Equipment and vehicle use")]
    for name, start, dur, lat, lon, acres, cause in fires:
        for d in range(dur):
            day = pd.Timestamp(start) + pd.Timedelta(days=d)
            for sensor, k in (("MODIS", 45), ("VIIRS", 260)):
                n = int(k * np.exp(-((d - dur / 3) ** 2) / dur))
                rows.append(pd.DataFrame({
                    "latitude": lat + rng.normal(0, 0.08, n),
                    "longitude": lon + rng.normal(0, 0.08, n),
                    "acq_date": day, "instrument": sensor,
                    "confidence": "90" if sensor == "MODIS" else "h",
                    "frp": rng.gamma(3, 30, n)}))
    firms = outdir / "demo_firms.csv"
    pd.concat(rows).to_csv(firms, index=False)
    fcsv = outdir / "demo_fires.csv"
    pd.DataFrame([{"name": n, "start_date": s,
                   "end_date": (pd.Timestamp(s) + pd.Timedelta(days=d)).date(),
                   "lat": la, "lon": lo, "size_acres": a, "cause": c, "state": "CA"}
                  for n, s, d, la, lo, a, c in fires]).to_csv(fcsv, index=False)
    return [firms], fcsv


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("PIPELINE")[0].strip())
    ap.add_argument("--firms", nargs="+", help="FIRMS archive CSV file(s)")
    ap.add_argument("--fod", help="FPA FOD SQLite file")
    ap.add_argument("--fires-csv", help="CSV of fire records (see module docstring)")
    ap.add_argument("--out", default="out", help="output directory")
    ap.add_argument("--cell-deg", type=float, default=0.25, help="grid cell size in degrees")
    ap.add_argument("--z", type=float, default=4.0, help="z-score threshold for a spike")
    ap.add_argument("--min-count", type=int, default=10, help="min detections per cell-day")
    ap.add_argument("--modis-min-conf", type=int, default=30)
    ap.add_argument("--keep-low-viirs", action="store_true")
    ap.add_argument("--min-acres", type=float, default=10.0, help="ignore smaller fires")
    ap.add_argument("--buffer-km", type=float, default=15.0)
    ap.add_argument("--pad-days", type=int, default=2)
    ap.add_argument("--top-k", type=int, default=3, help="candidate fires kept per event")
    ap.add_argument("--wiki", action="store_true", help="look up Wikipedia for top matches")
    ap.add_argument("--wiki-limit", type=int, default=50)
    ap.add_argument("--demo", action="store_true", help="run on synthetic data")
    # cell-year validation mode
    ap.add_argument("--cells-from", help="fires.html with the embedded cell-year data -> validate predicted VIIRS against incidents")
    ap.add_argument("--pred-json", default="out/cells_predicted.json", help="back-cast cells from train.py")
    ap.add_argument("--loyo-json", default="out/cells_predicted_loyo.json", help="out-of-year cells from train.py")
    ap.add_argument("--years", nargs=2, type=int, default=[2001, 2020], metavar=("FROM", "TO"), help="years to validate (FPA FOD ends 2020)")
    ap.add_argument("--big-acres", type=float, default=100.0, help="a cell-year has a 'big fire' at this many acres burned on sampled days")
    ap.add_argument("--huge-acres", type=float, default=5000.0, help="recall is also reported for cell-years this big")
    ap.add_argument("--boot", type=int, default=1000, help="bootstrap / permutation repetitions")
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.cells_from:
        if not (a.fod or a.fires_csv):
            ap.error("--cells-from needs --fod or --fires-csv")
        return run_cell_year(a, out)
    if a.demo:
        a.firms, a.fires_csv = make_demo(out)
        a.min_acres = 0
    if not a.firms:
        ap.error("give --firms files, or use --demo")

    det = load_firms(a.firms, a.modis_min_conf, a.keep_low_viirs)
    print(f"Detections kept: {len(det):,} "
          f"({det['sensor'].value_counts().to_dict()}), "
          f"{det['date'].min().date()} to {det['date'].max().date()}")

    spikes = find_spikes(det, a.cell_deg, a.z, a.min_count)
    events, spikes = build_events(spikes, a.cell_deg)
    spikes.to_csv(out / "spikes.csv", index=False)
    events.to_csv(out / "events.csv", index=False)
    print(f"Spike cell-days: {len(spikes):,}   merged events: {len(events):,}")
    if events.empty:
        print("No spikes found. Try lowering --z or --min-count.")
        return

    if a.fod:
        fires = load_fod(a.fod, a.min_acres)
    elif a.fires_csv:
        fires = load_fires_csv(a.fires_csv, a.min_acres)
    else:
        print("No fire records given (--fod or --fires-csv); skipping matching.")
        return
    print(f"Fire records loaded: {len(fires):,}")

    matches = match_events(events, fires, a.buffer_km, a.pad_days, top_k=a.top_k)
    if a.wiki:
        matches = add_wiki(matches, a.wiki_limit)
    matches.to_csv(out / "matches.csv", index=False)
    n_matched = matches.loc[matches["status"] == "matched", "event_id"].nunique()
    print(f"Events matched to a recorded fire: {n_matched} of {len(events)}")
    print(f"Unmatched events (possibly agricultural or prescribed burning, or a gap "
          f"in the fire records): {len(events) - n_matched}")
    print(f"Wrote spikes.csv, events.csv, matches.csv to {out}/")


if __name__ == "__main__":
    main()
