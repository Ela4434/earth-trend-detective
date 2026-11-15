"""Known-answer and false-positive tests for the statistics engine.   python tests/test_stats.py

Each check prints PASS/FAIL with the numbers, so the stats person can show them to judges.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.stats import trend as TR, relations as REL  # noqa: E402
from backend.sources import giovanni_csv, power  # noqa: E402

warnings.filterwarnings("ignore")
rng = np.random.default_rng(42)
IDX = pd.date_range("1981-01-01", "2025-12-01", freq="MS")
N = len(IDX)
SEASON = 8 * np.cos(2 * np.pi * (IDX.month.values - 7) / 12)
T = np.arange(N) / 12
fails = 0


def ar1(phi, sd=1.0, n=N):
    e = rng.standard_normal(n) * sd
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


def check(name, cond, detail):
    global fails
    print(f"{'PASS' if cond else 'FAIL'}  {name}: {detail}")
    fails += 0 if cond else 1


# 1. a known trend is recovered
s = pd.Series(10 + SEASON + 0.03 * T + ar1(0.5), index=IDX)
r = TR.analyse(s)
lo, hi = r["annual"]["sen_ci_dec"]
check("trend recovered", r["verdict"] == "PROVEN" and lo <= 0.3 <= hi,
      f"true +0.300/decade, Sen {r['sen_dec']:+.3f} (95% CI {lo:+.3f}…{hi:+.3f}), OLS {r['ols']['slope_dec']:+.3f}, verdict {r['verdict']}")

# 2. false positives under no trend (autocorrelated noise, seasonal cycle)
M = 300
for phi in (0.3, 0.6):
    proven = plain_ols = plain_sk = 0
    for _ in range(M):
        x = pd.Series(10 + SEASON + ar1(phi), index=IDX)
        rr = TR.analyse(x)
        proven += rr["verdict"] == "PROVEN"
        plain_ols += rr["ols"]["p"] < 0.05
        plain_sk += rr["seasonal_kendall_plain"]["p"] < 0.05
    check(f"no-trend false-positive rate (AR1 φ={phi})", proven / M <= 0.08,
          f"PROVEN in {100 * proven / M:.1f}% of {M} trendless series (target ≈5%); "
          f"naive OLS t-test: {100 * plain_ols / M:.1f}%, plain Seasonal Kendall: {100 * plain_sk / M:.1f}%")

# 2b. power: a modest real trend is still found
hits = 0
for _ in range(100):
    x = pd.Series(10 + SEASON + 0.03 * T + ar1(0.6), index=IDX)
    hits += TR.analyse(x)["verdict"] == "PROVEN"
check("power for +0.3/decade (noise SD≈1.25)", hits >= 80, f"PROVEN in {hits}% of 100 series with a real trend")

# 3. correlation of independent autocorrelated series: N_eff keeps false positives near 5%
naive = adj = 0
for _ in range(M):
    a = pd.Series(ar1(0.8), index=IDX)
    b = pd.Series(ar1(0.8), index=IDX)
    p = REL.pair(a, b)
    from scipy import stats as st
    naive += st.pearsonr(a, b)[1] < 0.05
    adj += p["pearson"]["p"] < 0.05
check("independent series: effective-N correction", adj / M <= 0.08,
      f"significant in {100 * adj / M:.1f}% with N_eff (target ≈5%), naive Pearson p says {100 * naive / M:.1f}%")

# 4. two unrelated trending series: 'trend-only', not 'linked'
df = pd.DataFrame({"A": 5 + SEASON + 0.05 * T + ar1(0.3, 0.5), "B": 2 - 0.04 * T + ar1(0.3, 0.5),
                   "C": None}, index=IDX)
df["C"] = 0.8 * (df["A"] - SEASON - 0.05 * T) + ar1(0.3, 0.4)  # truly linked to A year-to-year, no trend
d = REL.discover(df, ["A", "B", "C"], family=lambda v: v)
kinds = {(f["x"], f["y"]): f["kind"] for f in d["found"]}
check("shared trend vs real link", kinds.get(("A", "B")) == "trend-only" and kinds.get(("A", "C")) == "linked",
      f"A×B (both trend, unrelated) → {kinds.get(('A', 'B'))}; A×C (co-vary, no shared trend) → {kinds.get(('A', 'C'))}")

# 5. Pettitt + model comparison: step vs smooth trend
step = pd.Series(10 + SEASON + np.where(IDX.year >= 2003, 1.5, 0) + ar1(0.3, 0.6), index=IDX)
smooth = pd.Series(10 + SEASON + 0.035 * T + ar1(0.3, 0.6), index=IDX)
rs, rl = TR.analyse(step), TR.analyse(smooth)
check("step detected only for a step", rs["pettitt"]["step_better"] and not rl["pettitt"]["step_better"],
      f"step series: break {rs['pettitt']['year']} step_better={rs['pettitt']['step_better']}; "
      f"smooth trend: break {rl['pettitt']['year']} step_better={rl['pettitt']['step_better']}")

# 6. Weatherhead years-to-detect shrinks as the trend grows
res = np.array(ar1(0.5))
y1, y2 = TR.years_to_detect(res, 0.01), TR.years_to_detect(res, 0.04)
check("years-to-detect monotone", y1 > y2 > 0, f"0.01/yr → {y1:.1f} years, 0.04/yr → {y2:.1f} years")

# 7. multiple regression recovers known standardised effects
x1, x2 = ar1(0.3), ar1(0.3)
yv = 0.7 * x1 - 0.4 * x2 + 0.5 * rng.standard_normal(N)
dfr = pd.DataFrame({"y": yv, "x1": x1, "x2": x2}, index=IDX)
rg = REL.regression(dfr, "y", ["x1", "x2"])
b = {c["name"]: c["beta"] for c in rg["coefs"]}
check("regression signs and order", b["x1"] > 0 > b["x2"] and abs(b["x1"]) > abs(b["x2"]),
      f"β x1 {b['x1']:+.2f}, β x2 {b['x2']:+.2f}, R² {rg['r2']:.2f}")

# 8. BH-FDR matches a hand-computed example
q = TR.bh_fdr([0.01, 0.04, 0.03, 0.20])
# sorted p 0.01,0.03,0.04,0.20 -> p*m/rank 0.04,0.06,0.0533,0.20 -> running min from the top 0.04,0.0533,0.0533,0.20
check("BH-FDR", np.allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.20], atol=1e-9),
      f"q = {np.round(q, 4).tolist()} (expected [0.04, 0.0533, 0.0533, 0.2])")

# 9. Giovanni CSV parser on a file in Giovanni's layout
fx = Path(__file__).parent / "fixtures" / "g4.areaAvgTimeSeries.M2TMNXSLV_5_12_4_T2M.19800101-19891231.16E_45N_23E_48N.csv"
g = giovanni_csv.parse(fx)
check("Giovanni CSV parser", g["bbox"] == (16.0, 45.0, 23.0, 48.5) and g["series"].notna().sum() == 119 and "temperature" in g["label"].lower(),
      f"label '{g['label']}', unit '{g['unit']}', bbox {g['bbox']}, {g['series'].notna().sum()} valid months (1 fill value dropped)")

# 10. NASA POWER JSON parser + cos-lat weighting
import json
pj = json.loads((Path(__file__).parent / "fixtures" / "power_point_sample.json").read_text())
df1, units = power.point_frame(pj)
check("POWER JSON parser", list(df1.columns) == ["T2M", "PRECTOTCORR"] and df1["T2M"].isna().sum() == 1 and len(df1) == 24 and units["T2M"] == "C",
      f"{len(df1)} months, annual '13' rows dropped, fill -999 → NaN ({int(df1['T2M'].isna().sum())}), units {units}")

print("\nALL PASS" if fails == 0 else f"\n{fails} FAILED")
sys.exit(1 if fails else 0)
