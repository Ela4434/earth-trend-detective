"""Step 1 on a new laptop: is NASA POWER reachable, and does it return what the app expects?

    python tools/check_power.py

Downloads ONE point (Madrid) and prints, per variable, the months found, the first/last month
and the unit. Paste the output to the team chat if anything says MISSING or ERROR.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.catalog import VARIABLES  # noqa: E402
from backend.sources import power  # noqa: E402

lat, lon = 40.5, -3.5
print(f"Requesting NASA POWER monthly point {lat}, {lon} ...")
try:
    j = power._fetch_point_remote(lat, lon)
except Exception as e:
    print("ERROR:", e)
    sys.exit(1)
df, units = power.point_frame(j)
for v in VARIABLES:
    if v in df and df[v].notna().any():
        s = df[v].dropna()
        print(f"  OK       {v:<18} {len(s):>4} months  {s.index[0]:%Y-%m} … {s.index[-1]:%Y-%m}  unit '{units.get(v, '')}'"
              f"  mean {s.mean():.3g}")
    else:
        print(f"  MISSING  {v}")
print("requests made:", power.status()["requests"])
