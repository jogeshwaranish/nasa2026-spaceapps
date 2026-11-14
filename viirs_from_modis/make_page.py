#!/usr/bin/env python3
"""Build index.html from out/ (run after train.py)."""
import json, re, csv
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "out"
D = json.loads(re.search(r"const D = (\{.*?\});\n", (HERE.parent / "fires.html").read_text(encoding="utf8"), re.S).group(1))
P = json.loads((OUT / "cells_predicted.json").read_text())


def rates(sensor):  # year -> flat [ix, iy, detections/day]
    return {y: [v for i in range(0, len(s["cells"]), 3)
                for v in (s["cells"][i], s["cells"][i + 1], round(s["cells"][i + 2] / s["days"], 3))]
            for y, s in D["sensors"][sensor]["years"].items() if s["days"]}


data = {"name": D["name"], "bbox": D["bbox"], "cell": D["cell"],
        "metrics": json.loads((OUT / "metrics.json").read_text()),
        "to_modis": D.get("to_modis"),   # VIIRS-in-MODIS-units series from ../consolidate.py (via build.py)
        "monthly": list(csv.DictReader(open(OUT / "monthly_predictions.csv"))),
        "modis": rates("MODIS"), "viirs": rates("VIIRS"), "pred": {str(k): v for k, v in P["years"].items()}}
html = (HERE / "page_template.html").read_text(encoding="utf8").replace("/*DATA*/null", json.dumps(data, separators=(",", ":")))
(HERE / "index.html").write_text(html, encoding="utf8")
print(f"wrote index.html ({len(html) / 1e6:.1f} MB)")
