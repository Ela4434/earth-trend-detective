"""SYNTHETIC data in NASA POWER's JSON format, for offline development and UI testing only.

Enabled with TREND_GLOBE_DEMO=1. Every response carries "synthetic": true and the app shows
a red banner. Never present these numbers as findings.
Built-in structure (so the statistics have something to find): warming trend, soil wetness
falling when it is hot and dry, clear-sky sunlight falling when aerosols rise, aerosols
trending down north of 35°N and up south of it (a "mirror trend").
"""
from __future__ import annotations
import zlib

import numpy as np
import pandas as pd

UNITS = {"T2M": "C", "T2M_MAX": "C", "T2M_MIN": "C", "TS": "C", "PRECTOTCORR": "mm/day", "RH2M": "%",
         "QV2M": "g/kg", "GWETROOT": "1", "GWETTOP": "1", "EVPTRNS": "mm/day", "ALLSKY_SFC_SW_DWN": "kW-hr/m^2/day",
         "CLRSKY_SFC_SW_DWN": "kW-hr/m^2/day", "CLOUD_AMT": "%", "AOD_55": "1", "WS10M": "m/s", "PS": "kPa"}


def _ar1(rng, n, phi, sd):
    e = rng.standard_normal(n) * sd
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


def point_json(lat, lon, start=1981, end=2025):
    seed = zlib.crc32(f"{lat:.1f},{lon:.1f}".encode())
    rng = np.random.default_rng(seed)
    idx = pd.date_range(f"{start}-01-01", f"{end}-12-01", freq="MS")
    n = len(idx)
    t = (np.arange(n) / 12.0)
    m = idx.month.values
    hemi = 1 if lat >= 0 else -1
    season = np.cos(2 * np.pi * (m - 7) / 12) * hemi
    regional = 0.5 * np.sin(np.deg2rad(lon)) + 0.3 * np.cos(np.deg2rad(lat * 2))
    heat = _ar1(rng, n, 0.5, 0.8)
    wet = _ar1(rng, n, 0.4, 1.0)
    aer = _ar1(rng, n, 0.6, 1.0)
    tbase = 27 - 0.45 * abs(lat)
    t2m = tbase + (6 + 0.15 * abs(lat)) * season + (0.03 + 0.01 * regional) * t + heat
    pr = np.clip(2.0 + 1.0 * regional + 0.8 * season * (-hemi * np.sign(lat - 30)) + 0.6 * wet - 0.004 * t * regional, 0.05, None)
    aod_trend = -0.0025 if lat > 35 else 0.002
    aod = np.clip(0.18 + 0.05 * season + aod_trend * t + 0.03 * aer, 0.02, None)
    clr = 6.0 + 2.2 * season - 4.0 * (aod - 0.18) + 0.05 * rng.standard_normal(n)
    cloud = np.clip(55 - 10 * season + 6 * wet + 2 * rng.standard_normal(n), 2, 98)
    allsky = clr * (1 - 0.0065 * cloud) + 0.05 * rng.standard_normal(n)
    gw = np.clip(0.55 + 0.06 * wet - 0.05 * heat - 0.0015 * t - 0.05 * season, 0.05, 0.95)
    gwt = np.clip(gw + 0.05 * rng.standard_normal(n) + 0.03 * wet, 0.03, 0.99)
    rh = np.clip(70 - 8 * season + 3 * wet - 2.5 * heat + rng.standard_normal(n), 10, 100)
    qv = np.clip(6 + 3 * season + 0.4 * heat + 0.3 * wet + 0.01 * t, 0.5, None)
    evp = np.clip(1.5 + 1.2 * season + 0.3 * wet + 0.1 * heat, 0.0, None)
    ws = np.clip(4 + 0.5 * season * -1 - 0.004 * t + 0.3 * rng.standard_normal(n), 0.3, None)
    ps = 101.3 - 0.0115 * max(0, 300 * (1 + regional)) / 10 + 0.3 * rng.standard_normal(n)
    vals = {"T2M": t2m, "T2M_MAX": t2m + 5 + 0.3 * rng.standard_normal(n), "T2M_MIN": t2m - 5 + 0.3 * rng.standard_normal(n),
            "TS": t2m + 0.5 + 0.4 * rng.standard_normal(n), "PRECTOTCORR": pr, "RH2M": rh, "QV2M": qv,
            "GWETROOT": gw, "GWETTOP": gwt, "EVPTRNS": evp, "ALLSKY_SFC_SW_DWN": allsky, "CLRSKY_SFC_SW_DWN": clr,
            "CLOUD_AMT": cloud, "AOD_55": aod, "WS10M": ws, "PS": ps}
    param = {}
    for k, v in vals.items():
        d = {f"{ts.year}{ts.month:02d}": round(float(x), 4) for ts, x in zip(idx, v)}
        if k in ("ALLSKY_SFC_SW_DWN", "CLRSKY_SFC_SW_DWN", "CLOUD_AMT"):
            for y in range(start, 1984):
                for mm in range(1, 13):
                    d[f"{y}{mm:02d}"] = -999.0
        param[k] = d
    return {"synthetic": True, "header": {"fill_value": -999.0},
            "properties": {"parameter": param},
            "parameters": {k: {"units": u} for k, u in UNITS.items()}}
