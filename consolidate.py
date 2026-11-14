#!/usr/bin/env python3
"""Express VIIRS in MODIS units (stdlib only).

MODIS is the reference: it has the long record (Nov 2000 onward), so 2000-2011 stays
as measured and only VIIRS is converted. Two steps:

1. Physical merge. A VIIRS pixel is ~375 m and a MODIS pixel ~1 km, so one fire lights
   up several VIIRS pixels where MODIS lights up one. Within each overpass, VIIRS
   detections that fall in the same KM x KM cell are merged into one MODIS-pixel-sized
   detection. Then  merged VIIRS ≈ the count VIIRS would give with 1 km pixels.

2. Statistical correction for what merging can't remove (small fires only VIIRS sees,
   MODIS flying two satellites vs one VIIRS satellite, off-nadir pixel growth):
       log(1 + MODIS/day) = a + b * log(1 + mergedVIIRS/day)
   fitted on the overlap months (MODIS as the target, so it's fitted in this
   direction, not inverted from a VIIRS-from-MODIS fit).

Validation is leave-one-year-out: every score is on a year the fit never saw, and the
90% band comes from those out-of-year errors, so its coverage is honest.

build.py calls merged_keys() while it reads the raw files and add_to_modis() before
writing fires.html.
"""
import math, statistics
from collections import defaultdict

KM = 1.0                    # merge cell size: MODIS nadir pixel
SENSITIVITY = (0.75, 1.0, 1.5)
PASS_GAP = 30               # minutes; granules closer than this belong to one overpass


def _minutes(hhmm):
    t = int(hhmm)
    return t // 100 * 60 + t % 100


def merged_keys(rows, west, south, north, km=KM):
    """One key per merged detection: (date, satellite, overpass, cell_x, cell_y).

    acq_time is per granule, so one overpass across California can carry several
    times a few minutes apart; those are grouped into one overpass first.
    """
    kx, ky = 111.32 * math.cos(math.radians((south + north) / 2)), 110.57
    times = defaultdict(set)
    for r in rows:
        times[(r["acq_date"], r["satellite"])].add(_minutes(r["acq_time"]))
    pass_of = {}
    for k, ts in times.items():
        n, prev = 0, None
        for t in sorted(ts):
            if prev is not None and t - prev > PASS_GAP:
                n += 1
            pass_of[(k, t)] = n; prev = t
    keys = set()
    for r in rows:
        k = (r["acq_date"], r["satellite"])
        x = (float(r["longitude"]) - west) * kx / km
        y = (float(r["latitude"]) - south) * ky / km
        keys.add(k + (pass_of[(k, _minutes(r["acq_time"]))], int(x), int(y)))
    return keys


def _ols(x, y):
    mx, my = statistics.fmean(x), statistics.fmean(y)
    b = sum((a - mx) * (c - my) for a, c in zip(x, y)) / sum((a - mx) ** 2 for a in x)
    return my - b * mx, b


def _quantile(v, q):
    v = sorted(v); i = q * (len(v) - 1); lo = int(i)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (i - lo)


def _loyo(P, xi, kind):
    """Leave-one-year-out predictions of log1p(MODIS/day). P rows: (month, modis, viirs_raw, viirs_merged)."""
    pred = []
    for i, p in enumerate(P):
        tr = [q for q in P if q[0][:4] != p[0][:4]]
        if kind == "fit":
            a, b = _ols([math.log1p(q[xi]) for q in tr], [math.log1p(q[1]) for q in tr])
            pred.append(a + b * math.log1p(p[xi]))
        elif kind == "ratio":
            k = sum(q[1] for q in tr) / sum(q[xi] for q in tr)
            pred.append(math.log1p(k * p[xi]))
        else:                                   # "none": take the VIIRS number as-is
            pred.append(math.log1p(p[xi]))
    return pred


def add_to_modis(D, merged_by_km):
    """merged_by_km: {km: {"YYYY-MM": merged VIIRS count}}; KM is the one used, the rest are a sensitivity check."""
    merged_monthly = merged_by_km[KM]
    mod = {r[0]: r for r in D["sensors"]["MODIS"]["monthly"]}
    P = []
    for r in D["sensors"]["VIIRS"]["monthly"]:
        m = mod.get(r[0])
        if r[2] and m and m[2]:
            P.append((r[0], m[1] / m[2], r[1] / r[2], merged_monthly.get(r[0], 0) / r[2]))
    y = [math.log1p(p[1]) for p in P]
    err = lambda pr: [t - q for t, q in zip(y, pr)]
    cands = {
        "raw VIIRS as-is": ("none", 2), "raw VIIRS ÷ constant ratio": ("ratio", 2), "raw VIIRS, log-log fit": ("fit", 2),
        "merged VIIRS as-is": ("none", 3), "merged VIIRS ÷ constant ratio": ("ratio", 3), "merged VIIRS, log-log fit": ("fit", 3),
    }
    scores, chosen_res = {}, None
    for name, (kind, xi) in cands.items():
        e = err(_loyo(P, xi, kind))
        scores[name] = {"mae_log": statistics.fmean(abs(v) for v in e), "bias_log": statistics.fmean(e)}
        if name == "merged VIIRS, log-log fit":
            chosen_res = e
    q_lo, q_hi = _quantile(chosen_res, 0.05), _quantile(chosen_res, 0.95)
    a, b = _ols([math.log1p(p[3]) for p in P], y)
    res_full = [t - (a + b * math.log1p(p[3])) for t, p in zip(y, P)]
    r2 = 1 - sum(v * v for v in res_full) / sum((t - statistics.fmean(y)) ** 2 for t in y)

    conv = lambda v: [round(math.expm1(a + b * math.log1p(v) + d), 3) for d in (0, q_lo, q_hi)]
    monthly = []
    for r in D["sensors"]["VIIRS"]["monthly"]:
        monthly.append([r[0]] + (conv(merged_monthly.get(r[0], 0) / r[2]) if r[2] else [None] * 3))
    years = defaultdict(list)
    for mo, *v in monthly:
        if v[0] is not None:
            years[mo[:4]].append(v)
    yearly = {k: [round(statistics.fmean(c[i] for c in v), 3) for i in range(3)] for k, v in years.items()}
    raw_tot = sum(p[2] for p in P); mer_tot = sum(p[3] for p in P); mod_tot = sum(p[1] for p in P)
    days = {r[0]: r[2] for r in D["sensors"]["VIIRS"]["monthly"]}
    sens = {str(k): sum(mm.get(p[0], 0) / days[p[0]] for p in P) / mod_tot for k, mm in sorted(merged_by_km.items())}
    D["to_modis"] = {
        "km": KM, "monthly": monthly, "yearly": yearly,
        "merged_counts": {k: merged_monthly[k] for k in sorted(merged_monthly)},
        "fit": {"a": a, "b": b, "r2": r2, "n": len(P), "years": f"{P[0][0][:4]}–{P[-1][0][:4]}",
                "interval90_log": [q_lo, q_hi], "coverage90": sum(q_lo <= v <= q_hi for v in chosen_res) / len(chosen_res),
                "ratio_raw": raw_tot / mod_tot, "ratio_merged": mer_tot / mod_tot, "scores": scores,
                "sensitivity": sens},
    }
    return D
