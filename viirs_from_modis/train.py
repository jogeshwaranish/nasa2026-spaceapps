#!/usr/bin/env python3
"""Train models that predict VIIRS detections from MODIS, then back-cast VIIRS for 2000-2011.

Inputs : ../fires.html (data embedded by build.py: monthly rates + per-year 0.05-degree cell counts)
Outputs: out/monthly_predictions.csv, out/cells_predicted.json, out/cells_predicted_loyo.json, out/metrics.json
Run    : python train.py          (needs numpy, scikit-learn)

Two models, both trained only on 2012+ (when VIIRS exists) and applied to MODIS-only years:
  A. monthly  : log1p(VIIRS/day) ~ MODIS this/prev/next month + season
  B. spatial  : log1p(VIIRS/day in a cell-year) ~ MODIS in the cell and its neighbourhood + location
Validation is leave-one-year-out, so every score is on a year the model never saw.
"""
import json, re, math
from pathlib import Path
import numpy as np
from scipy.ndimage import uniform_filter
from sklearn.ensemble import HistGradientBoostingRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge

HERE = Path(__file__).parent
OUT = HERE / "out"; OUT.mkdir(exist_ok=True)
html = (HERE.parent / "fires.html").read_text(encoding="utf8")
D = json.loads(re.search(r"const D = (\{.*?\});\n", html, re.S).group(1))
SM, SV = D["sensors"]["MODIS"], D["sensors"]["VIIRS"]
west, south, east, north = D["bbox"]; CELL = D["cell"]
NX, NY = math.ceil((east - west) / CELL), math.ceil((north - south) / CELL)
VIIRS_YEARS = sorted(int(y) for y in SV["years"])
MODIS_ONLY = [y for y in sorted(int(y) for y in SM["years"]) if y < VIIRS_YEARS[0]]
rng = np.random.default_rng(0)


def mae(a, b): return float(np.mean(np.abs(np.asarray(a) - np.asarray(b))))
def r2(y, p): y, p = np.asarray(y), np.asarray(p); return float(1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum())


# ============================ A. monthly model ============================
mon = {r[0]: r for r in SM["monthly"]}; vmon = {r[0]: r for r in SV["monthly"]}
months = [r[0] for r in SM["monthly"]]
mrate = {m: (mon[m][1] / mon[m][2] if mon[m][2] else None) for m in months}
def nb(m, k):                       # MODIS rate k months away (fallback to this month)
    i = months.index(m) + k
    v = mrate[months[i]] if 0 <= i < len(months) else None
    return v if v is not None else mrate[m]
def feat(m):
    t = mrate[m]; mo = int(m[5:7])
    return [math.log1p(t), math.log1p(nb(m, -1)), math.log1p(nb(m, 1)), math.sin(2 * math.pi * mo / 12), math.cos(2 * math.pi * mo / 12)]
tr_m = [m for m in months if m in vmon and vmon[m][2] and mrate[m] is not None]
X = np.array([feat(m) for m in tr_m]); y = np.array([math.log1p(vmon[m][1] / vmon[m][2]) for m in tr_m])
yr = np.array([int(m[:4]) for m in tr_m])

MODELS = {
    "log-log line (previous prototype)": (lambda: Ridge(alpha=1e-6), [0]),
    "ridge: MODIS + neighbours + season": (lambda: Ridge(alpha=1.0), [0, 1, 2, 3, 4]),
    "gradient boosting": (lambda: GradientBoostingRegressor(n_estimators=150, max_depth=2, learning_rate=0.05, subsample=0.8, random_state=0), [0, 1, 2, 3, 4]),
}
ratio_cv = np.zeros(len(y))
cv = {}
for name, (mk, cols) in MODELS.items():
    p = np.zeros(len(y))
    for yy in sorted(set(yr)):
        te = yr == yy
        p[te] = mk().fit(X[~te][:, cols], y[~te]).predict(X[te][:, cols])
    cv[name] = p
for yy in sorted(set(yr)):           # constant-ratio baseline
    te = yr == yy
    ratio = sum(vmon[m][1] / vmon[m][2] for m, t in zip(tr_m, te) if not t) / sum(mrate[m] for m, t in zip(tr_m, te) if not t)
    ratio_cv[te] = np.log1p(ratio * np.expm1(X[te][:, 0]))
cv["constant ratio"] = ratio_cv
mscore = {k: {"mae_log": mae(v, y), "r2_log": r2(y, v), "typical_factor": math.exp(mae(v, y))} for k, v in cv.items()}
best = min((k for k in MODELS), key=lambda k: mscore[k]["mae_log"])
res = y - cv[best]
q_lo, q_hi = np.quantile(res, [0.05, 0.95])                 # empirical 90% interval from out-of-year errors
cover = float(np.mean((res >= q_lo) & (res <= q_hi)))
mk, cols = MODELS[best]; final = mk().fit(X[:, cols], y)

rows = []
for m in months:
    if mrate[m] is None: rows.append([m, None, None, None, None, None]); continue
    f = np.array([feat(m)])[:, cols]; p = float(final.predict(f)[0])
    rows.append([m, round(mrate[m], 4), round(vmon[m][1] / vmon[m][2], 4) if m in vmon and vmon[m][2] else None,
                 round(math.expm1(p), 4), round(math.expm1(p + q_lo), 4), round(math.expm1(p + q_hi), 4)])
with open(OUT / "monthly_predictions.csv", "w") as f:
    f.write("month,modis_per_day,viirs_actual_per_day,viirs_predicted_per_day,pred_lo90,pred_hi90\n")
    for r in rows: f.write(",".join("" if v is None else str(v) for v in r) + "\n")

# ============================ B. spatial model ============================
def grid(sensor, year):
    s = D["sensors"][sensor]["years"].get(str(year)); g = np.zeros((NY, NX))
    if not s or not s["days"]: return None
    c = s["cells"]
    for i in range(0, len(c), 3):
        if c[i] < NX and c[i + 1] < NY: g[c[i + 1], c[i]] += c[i + 2]
    return g / s["days"]
LAT = (south + (np.arange(NY) + .5) * CELL)[:, None] * np.ones((1, NX)); LON = (west + (np.arange(NX) + .5) * CELL)[None, :] * np.ones((NY, 1))
def cell_features(gm):
    L = np.log1p(gm)
    return np.stack([L, uniform_filter(L, 3), uniform_filter(L, 7), uniform_filter(L, 21), np.full_like(L, math.log1p(gm.sum())), LAT, LON], -1)
def universe(gm, gv=None):
    u = uniform_filter((gm > 0).astype(float), 5) > 0
    return u | (gv > 0) if gv is not None else u
sets = {}
for yy in VIIRS_YEARS:
    gm, gv = grid("MODIS", yy), grid("VIIRS", yy)
    if gm is None or gv is None: continue
    u = universe(gm, gv); sets[yy] = (cell_features(gm)[u], np.log1p(gv)[u], u, gm, gv)
HGB = lambda: HistGradientBoostingRegressor(max_iter=200, learning_rate=0.06, max_leaf_nodes=15, min_samples_leaf=30, random_state=0)
sp = {"per_year": {}}; allp, ally, allb = [], [], []
loyo = {}
for yy in sets:
    Xtr = np.vstack([sets[k][0] for k in sets if k != yy]); ytr = np.concatenate([sets[k][1] for k in sets if k != yy])
    p = HGB().fit(Xtr, ytr).predict(sets[yy][0]); t = sets[yy][1]; Xte = sets[yy][0]
    ratio = np.expm1(ytr).sum() / np.expm1(np.vstack([sets[k][0] for k in sets if k != yy])[:, 0]).sum()
    base = np.log1p(ratio * np.expm1(Xte[:, 0]))            # constant-ratio baseline: VIIRS = k * MODIS in the cell
    act = (np.expm1(t) > 0) | (np.expm1(Xte[:, 0]) > 0)      # only cells where either sensor saw something
    sp["per_year"][yy] = {"active_cells": int(act.sum()), "mae_log": mae(p[act], t[act]), "mae_log_ratio_baseline": mae(base[act], t[act]),
                          "total_pred": float(np.expm1(p).sum()), "total_actual": float(np.expm1(t).sum())}
    allp.append(p[act]); ally.append(t[act]); allb.append(base[act])
    iy_, ix_ = np.nonzero(sets[yy][2]); r_ = np.expm1(p); r_[r_ < 0.005] = 0; kk = r_ > 0   # out-of-year prediction for this year
    loyo[yy] = [round(float(v), 3) if j == 2 else int(v) for a_, b_, c_ in zip(ix_[kk], iy_[kk], r_[kk]) for j, v in enumerate((a_, b_, c_))]
allp, ally, allb = map(np.concatenate, (allp, ally, allb))
sp.update({"scored_on": "cells where MODIS or VIIRS detected fire that year, leave-one-year-out",
           "mae_log": mae(allp, ally), "r2_log": r2(ally, allp), "mae_log_ratio_baseline": mae(allb, ally), "r2_log_ratio_baseline": r2(ally, allb),
           "corr_cell_rate": float(np.corrcoef(np.expm1(allp), np.expm1(ally))[0, 1])})
final_sp = HGB().fit(np.vstack([sets[k][0] for k in sets]), np.concatenate([sets[k][1] for k in sets]))
cells_out = {}
for yy in MODIS_ONLY + VIIRS_YEARS:                     # predicted VIIRS cell rates for every MODIS year
    gm = grid("MODIS", yy)
    if gm is None: continue
    u = universe(gm); p = np.expm1(final_sp.predict(cell_features(gm)[u])); p[p < 0.005] = 0
    iy, ix = np.nonzero(u); keep = p > 0
    cells_out[yy] = [round(float(v), 3) if j == 2 else int(v) for a, b, c in zip(ix[keep], iy[keep], p[keep]) for j, v in enumerate((a, b, c))]
(OUT / "cells_predicted.json").write_text(json.dumps({"cell": CELL, "bbox": D["bbox"], "years": cells_out}, separators=(",", ":")))

(OUT / "cells_predicted_loyo.json").write_text(json.dumps({"cell": CELL, "bbox": D["bbox"], "years": loyo}, separators=(",", ":")))   # honest 2012+ predictions
metrics = {"monthly": {"n_months": len(y), "validation": "leave-one-year-out, 2012-2026", "scores": mscore, "chosen": best,
                       "interval90_log": [float(q_lo), float(q_hi)], "interval_coverage_oos": cover},
           "spatial": sp, "viirs_years": VIIRS_YEARS, "backcast_years": MODIS_ONLY}
(OUT / "metrics.json").write_text(json.dumps(metrics, indent=1))
print(json.dumps({"monthly": metrics["monthly"], "spatial": {k: v for k, v in sp.items() if k != "per_year"}}, indent=1))
