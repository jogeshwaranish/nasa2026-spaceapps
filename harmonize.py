#!/usr/bin/env python3
"""Prototype MODIS -> VIIRS-equivalent harmonization (stdlib only).

Model: in the overlap (months where both sensors have sampled days),
    log(1 + VIIRS/day) = a + b * log(1 + MODIS/day) + e,   e ~ N(0, s^2)
fit by OLS, then applied to every MODIS month (incl. 2000-2011) to give a
"VIIRS-equivalent" rate with a ~90% prediction band. Validated by fitting on
the early overlap years and testing on the later ones.

Run `python harmonize.py` to regenerate fires.html from template.html + the data already embedded in it (no raw data needed);
build.py also calls add_harmonized() when it regenerates the page.
"""
import json, math, re, statistics
from pathlib import Path

ROOT = Path(__file__).parent
Z90 = 1.645


def ols(x, y):
    mx, my = statistics.fmean(x), statistics.fmean(y)
    sxx = sum((a - mx) ** 2 for a in x)
    b = sum((a - mx) * (c - my) for a, c in zip(x, y)) / sxx
    a = my - b * mx
    res = [c - (a + b * p) for p, c in zip(x, y)]
    s = math.sqrt(sum(r * r for r in res) / (len(x) - 2))
    ss_tot = sum((c - my) ** 2 for c in y)
    return a, b, s, 1 - sum(r * r for r in res) / ss_tot


def rate(r):                       # monthly row [month, count, days_sampled, days_possible] -> /day
    return r[1] / r[2] if r[2] else None


def pairs(D):
    v = {r[0]: rate(r) for r in D["sensors"]["VIIRS"]["monthly"]}
    return [(m[0], rate(m), v[m[0]]) for m in D["sensors"]["MODIS"]["monthly"]
            if m[0] in v and rate(m) is not None and v[m[0]] is not None]


def add_harmonized(D):
    P = pairs(D)
    x = [math.log1p(m) for _, m, _ in P]; y = [math.log1p(v) for _, _, v in P]
    a, b, s, r2 = ols(x, y)
    # holdout: fit through 2018, test 2019+
    tr = [i for i, p in enumerate(P) if p[0] < "2019"]; te = [i for i, p in enumerate(P) if p[0] >= "2019"]
    ha, hb, hs, _ = ols([x[i] for i in tr], [y[i] for i in tr])
    err = [y[i] - (ha + hb * x[i]) for i in te]
    # naive baseline: constant VIIRS/MODIS ratio from training years
    ratio = sum(P[i][2] for i in tr) / sum(P[i][1] for i in tr)
    base = [y[i] - math.log1p(ratio * P[i][1]) for i in te]
    cover = sum(abs(e) <= Z90 * hs for e in err) / len(err)
    pred = lambda m: (math.expm1(a + b * math.log1p(m)), math.expm1(a + b * math.log1p(m) - Z90 * s),
                      math.expm1(a + b * math.log1p(m) + Z90 * s))
    monthly = []
    for r in D["sensors"]["MODIS"]["monthly"]:
        m = rate(r)
        monthly.append([r[0]] + ([round(t, 3) for t in pred(m)] if m is not None else [None] * 3))
    years = {}
    for mo, *p in monthly:
        if p[0] is not None:
            years.setdefault(mo[:4], []).append(p)
    yearly = {y_: [round(statistics.fmean(c[i] for c in v), 3) for i in range(3)] for y_, v in years.items()}
    D["harmonized"] = {
        "monthly": monthly, "yearly": yearly,
        "fit": {"a": a, "b": b, "sigma": s, "r2": r2, "n": len(P), "ratio_all": sum(p[2] for p in P) / sum(p[1] for p in P),
                "holdout": {"train": f"{P[tr[0]][0]}–{P[tr[-1]][0]}", "test": f"{P[te[0]][0]}–{P[te[-1]][0]}", "n_test": len(te),
                            "mae_log": statistics.fmean(abs(e) for e in err), "mae_log_ratio": statistics.fmean(abs(e) for e in base),
                            "bias_log": statistics.fmean(err), "coverage90": cover}},
    }
    return D


def main():
    """Rebuild fires.html from template.html using the data already embedded in fires.html (no raw downloads needed)."""
    html = (ROOT / "fires.html").read_text(encoding="utf8")
    D = json.loads(re.search(r"const D = (\{.*?\});\n", html, re.S).group(1))
    D = add_harmonized(D)
    out = (ROOT / "template.html").read_text(encoding="utf8").replace("/*DATA*/null", json.dumps(D, separators=(",", ":")))
    (ROOT / "fires.html").write_text(out, encoding="utf8")
    print(json.dumps(D["harmonized"]["fit"], indent=1))


if __name__ == "__main__":
    main()
