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

3. Optional extra terms, kept only if they cut the held-out error: "fire clustering"
   (VIIRS pixels per merged cell: big clumped fires vs scattered small burns),
   season, and weather on the sampled days from NASA POWER (weather.py).

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
    """One key per merged detection: (date, satellite, overpass, cell_x, cell_y)."""
    return set(merged_groups(rows, west, south, north, km))


def merged_groups(rows, west, south, north, km=KM):
    """{merge key: [rows merged into it]}; key = (date, satellite, overpass, cell_x, cell_y).

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
    groups = defaultdict(list)
    for r in rows:
        k = (r["acq_date"], r["satellite"])
        x = (float(r["longitude"]) - west) * kx / km
        y = (float(r["latitude"]) - south) * ky / km
        groups[k + (pass_of[(k, _minutes(r["acq_time"]))], int(x), int(y))].append(r)
    return groups


def _solve(A, b):
    n = len(A); M = [row[:] + [v] for row, v in zip(A, b)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(M[r][i])); M[i], M[piv] = M[piv], M[i]
        for r in range(n):
            if r != i:
                f = M[r][i] / M[i][i]; M[r] = [a - f * c for a, c in zip(M[r], M[i])]
    return [M[i][n] / M[i][i] for i in range(n)]


def _fit(X, y):
    """OLS on standardized columns; returns (intercept, coefs) in original units."""
    k = len(X[0])
    mu = [statistics.fmean(x[j] for x in X) for j in range(k)]
    sd = [statistics.pstdev([x[j] for x in X]) or 1 for j in range(k)]
    Z = [[(x[j] - mu[j]) / sd[j] for j in range(k)] for x in X]
    my = statistics.fmean(y)
    A = [[sum(z[i] * z[j] for z in Z) + (1e-6 if i == j else 0) for j in range(k)] for i in range(k)]
    w = _solve(A, [sum(z[i] * (t - my) for z, t in zip(Z, y)) for i in range(k)])
    coef = [w[j] / sd[j] for j in range(k)]
    return my - sum(c * m for c, m in zip(coef, mu)), coef


def _quantile(v, q):
    v = sorted(v); i = q * (len(v) - 1); lo = int(i)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (i - lo)


# Candidate models for log1p(MODIS/day). Each "fit" lists its terms; a term maps a month row to a number.
TERMS = {
    "log(1 + raw VIIRS/day)": lambda p: math.log1p(p["raw"]),
    "log(1 + merged VIIRS/day)": lambda p: math.log1p(p["merged"]),
    "fire clustering": lambda p: math.log((p["raw"] + 1) / (p["merged"] + 1)),   # VIIRS pixels per merged cell
    "season sin": lambda p: math.sin(2 * math.pi * int(p["month"][5:7]) / 12),
    "season cos": lambda p: math.cos(2 * math.pi * int(p["month"][5:7]) / 12),
    "cloud %": lambda p: p["w"]["cloud"],
    "log(1 + rain mm/day)": lambda p: math.log1p(p["w"]["precip"]),
    "humidity %": lambda p: p["w"]["humidity"],
    "max temp °C": lambda p: p["w"]["tmax"],
    "wind m/s": lambda p: p["w"]["wind"],
}
M_, CL = "log(1 + merged VIIRS/day)", "fire clustering"
CANDIDATES = [  # (name, kind, terms or column); merged-VIIRS fits are listed simplest first
    ("raw VIIRS ÷ constant ratio", "ratio", "raw"),
    ("merged VIIRS ÷ constant ratio", "ratio", "merged"),
    ("raw VIIRS, log-log fit", "fit", ["log(1 + raw VIIRS/day)"]),
    ("merged VIIRS, log-log fit", "fit", [M_]),
    ("merged + fire clustering", "fit", [M_, CL]),
    ("merged + clustering + season", "fit", [M_, CL, "season sin", "season cos"]),
    ("merged + clustering + cloud & rain", "fit", [M_, CL, "cloud %", "log(1 + rain mm/day)"]),
    ("merged + clustering + all weather", "fit", [M_, CL, "cloud %", "log(1 + rain mm/day)", "humidity %", "max temp °C", "wind m/s"]),
]
MIN_GAIN = 0.005   # a more complex merged model must cut the held-out log error by this much to be chosen


def _loyo(P, kind, spec):
    """Leave-one-year-out predictions of log1p(MODIS/day)."""
    pred = []
    for p in P:
        tr = [q for q in P if q["month"][:4] != p["month"][:4]]
        if kind == "ratio":
            k = sum(q["modis"] for q in tr) / sum(q[spec] for q in tr)
            pred.append(math.log1p(k * p[spec]))
        else:
            a, c = _fit([[TERMS[t](q) for t in spec] for q in tr], [math.log1p(q["modis"]) for q in tr])
            pred.append(a + sum(ci * TERMS[t](p) for ci, t in zip(c, spec)))
    return pred


def _fill_weather(W):
    """Months missing a variable (POWER cloud lags a few weeks) get that calendar month's average."""
    clim = {}
    for m, v in W.items():
        for k, x in v.items():
            clim.setdefault((m[5:7], k), []).append(x)
    clim = {k: statistics.fmean(v) for k, v in clim.items()}
    keys = {k for v in W.values() for k in v}
    return lambda m: {k: W.get(m, {}).get(k, clim.get((m[5:7], k))) for k in keys}


def add_to_modis(D, merged_by_km, W=None):
    """merged_by_km: {km: {"YYYY-MM": merged VIIRS count}}; KM is the one used, the rest are a sensitivity check.
    W: weather.monthly() or None (weather candidates are then skipped)."""
    merged_monthly = merged_by_km[KM]
    wx = _fill_weather(W) if W else (lambda m: None)
    mod = {r[0]: r for r in D["sensors"]["MODIS"]["monthly"]}
    row = lambda r: {"month": r[0], "raw": r[1] / r[2], "merged": merged_monthly.get(r[0], 0) / r[2], "w": wx(r[0])}
    P = []
    for r in D["sensors"]["VIIRS"]["monthly"]:
        m = mod.get(r[0])
        if r[2] and m and m[2]:
            P.append(dict(row(r), modis=m[1] / m[2]))
    y = [math.log1p(p["modis"]) for p in P]
    cands = [c for c in CANDIDATES if W or not any(t in ("cloud %",) for t in (c[2] if c[1] == "fit" else []))]
    scores, res = {}, {}
    for name, kind, spec in cands:
        e = [t - q for t, q in zip(y, _loyo(P, kind, spec))]
        scores[name] = {"mae_log": statistics.fmean(abs(v) for v in e), "bias_log": statistics.fmean(e),
                        "mae_log_quiet": statistics.fmean(abs(v) for v, p in zip(e, P) if p["modis"] < 10)}
        res[name] = e
    chosen = "merged VIIRS, log-log fit"
    for name, kind, spec in cands:
        if kind == "fit" and spec[0] == M_ and scores[name]["mae_log"] <= scores[chosen]["mae_log"] - MIN_GAIN:
            chosen = name
    spec = next(c[2] for c in cands if c[0] == chosen)
    q_lo, q_hi = _quantile(res[chosen], 0.05), _quantile(res[chosen], 0.95)
    a, coef = _fit([[TERMS[t](p) for t in spec] for p in P], y)
    f_ = lambda p: a + sum(c * TERMS[t](p) for c, t in zip(coef, spec))
    r2 = 1 - sum((t - f_(p)) ** 2 for t, p in zip(y, P)) / sum((t - statistics.fmean(y)) ** 2 for t in y)

    monthly = []
    for r in D["sensors"]["VIIRS"]["monthly"]:
        if not r[2]:
            monthly.append([r[0], None, None, None]); continue
        v = f_(row(r)); monthly.append([r[0]] + [round(math.expm1(v + d), 3) for d in (0, q_lo, q_hi)])
    years = defaultdict(list)
    for mo, *v in monthly:
        if v[0] is not None:
            years[mo[:4]].append(v)
    yearly = {k: [round(statistics.fmean(c[i] for c in v), 3) for i in range(3)] for k, v in years.items()}
    raw_tot = sum(p["raw"] for p in P); mer_tot = sum(p["merged"] for p in P); mod_tot = sum(p["modis"] for p in P)
    days = {r[0]: r[2] for r in D["sensors"]["VIIRS"]["monthly"]}
    sens = {str(k): sum(mm.get(p["month"], 0) / days[p["month"]] for p in P) / mod_tot for k, mm in sorted(merged_by_km.items())}
    D["to_modis"] = {
        "km": KM, "monthly": monthly, "yearly": yearly,
        "merged_counts": {k: merged_monthly[k] for k in sorted(merged_monthly)},
        "fit": {"chosen": chosen, "intercept": a, "terms": [[t, c] for t, c in zip(spec, coef)], "r2": r2, "n": len(P),
                "years": f"{P[0]['month'][:4]}–{P[-1]['month'][:4]}", "order": [c[0] for c in cands], "min_gain": MIN_GAIN,
                "weather": {"source": "NASA POWER daily, mean of sampled days at 8 fire-prone points"} if W else None,
                "interval90_log": [q_lo, q_hi], "coverage90": sum(q_lo <= v <= q_hi for v in res[chosen]) / len(res[chosen]),
                "ratio_raw": raw_tot / mod_tot, "ratio_merged": mer_tot / mod_tot, "scores": scores,
                "sensitivity": sens},
    }
    return D
