#!/usr/bin/env python3
"""Consolidate VIIRS fires, then back-cast VIIRS to 2000-2011 from MODIS (stdlib only).

VIIRS is the reference, after one physical step:

1. Consolidate. A VIIRS pixel is ~375 m and a MODIS pixel ~1 km, so one fire lights up
   several VIIRS pixels. Within each overpass, VIIRS detections in the same KM x KM cell
   are merged into one fire. This consolidated VIIRS count is the "true VIIRS" record.

2. Back-cast. VIIRS only exists from 2012, MODIS from Nov 2000. On the months both
   flew, a model learns consolidated VIIRS/day from MODIS/day:
       log(1 + VIIRS/day) = a + b * log(1 + MODIS/day) [+ season, neighbouring months]
   and is then applied to every MODIS month, giving the VIIRS that would probably have
   been seen in 2000-2011 (and a prediction to compare with the true VIIRS from 2012).

Validation is leave-one-year-out: every score is on a year the fit never saw, and the
90% band comes from those out-of-year errors. A more complex candidate is used only if
it cuts the held-out error by MIN_GAIN.

build.py calls merged_keys() while it reads the raw files and add_backcast() before
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


# Candidate back-cast models for log1p(consolidated VIIRS/day); each term maps a month row to a number.
TERMS = {
    "log(1 + MODIS/day)": lambda p: math.log1p(p["modis"]),
    "MODIS previous month": lambda p: math.log1p(p["prev"]),
    "MODIS next month": lambda p: math.log1p(p["next"]),
    "season sin": lambda p: math.sin(2 * math.pi * int(p["month"][5:7]) / 12),
    "season cos": lambda p: math.cos(2 * math.pi * int(p["month"][5:7]) / 12),
}
L_ = "log(1 + MODIS/day)"
CANDIDATES = [  # (name, kind, terms), simplest first
    ("constant ratio", "ratio", None),
    ("log-log fit", "fit", [L_]),
    ("+ season", "fit", [L_, "season sin", "season cos"]),
    ("+ neighbouring months", "fit", [L_, "MODIS previous month", "MODIS next month"]),
    ("+ neighbouring months + season", "fit", [L_, "MODIS previous month", "MODIS next month", "season sin", "season cos"]),
]
MIN_GAIN = 0.005   # a more complex model must cut the held-out log error by this much to be chosen


def _loyo(P, kind, spec):
    """Leave-one-year-out predictions of log1p(VIIRS/day)."""
    pred = []
    for p in P:
        tr = [q for q in P if q["month"][:4] != p["month"][:4]]
        if kind == "ratio":
            k = sum(q["viirs"] for q in tr) / sum(q["modis"] for q in tr)
            pred.append(math.log1p(k * p["modis"]))
        else:
            a, c = _fit([[TERMS[t](q) for t in spec] for q in tr], [math.log1p(q["viirs"]) for q in tr])
            pred.append(a + sum(ci * TERMS[t](p) for ci, t in zip(c, spec)))
    return pred


def add_backcast(D, merged_by_km):
    """merged_by_km: {km: {"YYYY-MM": consolidated VIIRS count}}; KM is the one used, the rest are a sensitivity check.
    Adds D["to_viirs"]: the back-cast for every MODIS month, the fit, and its held-out scores."""
    merged = merged_by_km[KM]
    mrows = D["sensors"]["MODIS"]["monthly"]
    rate = {r[0]: r[1] / r[2] for r in mrows if r[2]}
    order = [r[0] for r in mrows]
    def row(mo):
        i = order.index(mo)
        nb = lambda k: rate.get(order[i + k]) if 0 <= i + k < len(order) and order[i + k] in rate else rate[mo]
        return {"month": mo, "modis": rate[mo], "prev": nb(-1), "next": nb(1)}
    vdays = {r[0]: r[2] for r in D["sensors"]["VIIRS"]["monthly"]}
    vraw = {r[0]: r[1] for r in D["sensors"]["VIIRS"]["monthly"]}
    P = [dict(row(mo), viirs=merged.get(mo, 0) / vdays[mo], raw=vraw[mo] / vdays[mo]) for mo in order if mo in rate and vdays.get(mo)]
    y = [math.log1p(p["viirs"]) for p in P]
    scores, res = {}, {}
    for name, kind, spec in CANDIDATES:
        e = [t - q for t, q in zip(y, _loyo(P, kind, spec))]
        scores[name] = {"mae_log": statistics.fmean(abs(v) for v in e), "bias_log": statistics.fmean(e),
                        "mae_log_quiet": statistics.fmean(abs(v) for v, p in zip(e, P) if p["modis"] < 10)}
        res[name] = e
    chosen = "log-log fit"
    for name, kind, spec in CANDIDATES:
        if kind == "fit" and scores[name]["mae_log"] <= scores[chosen]["mae_log"] - MIN_GAIN:
            chosen = name
    spec = next(c[2] for c in CANDIDATES if c[0] == chosen)
    q_lo, q_hi = _quantile(res[chosen], 0.05), _quantile(res[chosen], 0.95)
    a, coef = _fit([[TERMS[t](p) for t in spec] for p in P], y)
    f_ = lambda p: a + sum(c * TERMS[t](p) for c, t in zip(coef, spec))
    r2 = 1 - sum((t - f_(p)) ** 2 for t, p in zip(y, P)) / sum((t - statistics.fmean(y)) ** 2 for t in y)
    monthly = [[mo] + ([round(math.expm1(f_(row(mo)) + d), 3) for d in (0, q_lo, q_hi)] if mo in rate else [None] * 3) for mo in order]
    years = defaultdict(list)
    for mo, *v in monthly:
        if v[0] is not None:
            years[mo[:4]].append(v)
    yearly = {k: [round(statistics.fmean(c[i] for c in v), 3) for i in range(3)] for k, v in years.items()}
    mod_tot = sum(p["modis"] for p in P)
    sens = {str(k): sum(mm.get(p["month"], 0) / vdays[p["month"]] for p in P) / mod_tot for k, mm in sorted(merged_by_km.items())}
    D["to_viirs"] = {
        "km": KM, "monthly": monthly, "yearly": yearly,
        "merged_counts": {k: merged[k] for k in sorted(merged)},
        "fit": {"chosen": chosen, "intercept": a, "terms": [[t, c] for t, c in zip(spec, coef)], "r2": r2, "n": len(P),
                "years": f"{P[0]['month'][:4]}–{P[-1]['month'][:4]}", "order": [c[0] for c in CANDIDATES], "min_gain": MIN_GAIN,
                "interval90_log": [q_lo, q_hi], "coverage90": sum(q_lo <= v <= q_hi for v in res[chosen]) / len(res[chosen]),
                "ratio_raw": sum(p["raw"] for p in P) / mod_tot, "ratio_merged": sum(p["viirs"] for p in P) / mod_tot,
                "scores": scores, "sensitivity": sens},
    }
    return D
