"""NASA POWER monthly point API -> region-mean monthly series.

API: https://power.larc.nasa.gov/api/temporal/monthly/point
     ?parameters=T2M,...&community=RE&longitude=..&latitude=..&start=YYYY&end=YYYY&format=JSON
Response: properties.parameter.<PARAM>.<YYYYMM> (MM=13 is the annual value, dropped),
          parameters.<PARAM>.units, header.fill_value (-999).
No login needed. Every point response is cached on disk forever (monthly history doesn't change;
delete data/cache/ to refresh the latest months).
"""
from __future__ import annotations
import datetime as dt
import json
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from ..catalog import VARIABLES

CACHE = Path(__file__).resolve().parents[2] / "data" / "cache"
URL = "https://power.larc.nasa.gov/api/temporal/monthly/point"
START_YEAR = 1981
PARAMS = list(VARIABLES)
WORKERS = int(os.environ.get("POWER_WORKERS", "4"))
_lock = threading.Lock()
_status = {"requests": 0, "cache_hits": 0, "errors": []}


class PowerError(RuntimeError):
    pass


def _cache_path(lat, lon):
    return CACHE / f"power_{lat:+06.1f}_{lon:+07.1f}.json"


def _get(lat, lon, params, end, tries=4):
    """GET with retry on rate limits / server hiccups (429, 5xx, timeouts)."""
    q = dict(parameters=",".join(params), community="RE", longitude=lon, latitude=lat,
             start=START_YEAR, end=end, format="JSON")
    for k in range(tries):
        try:
            r = requests.get(URL, params=q, timeout=90)
        except requests.RequestException as e:
            if k == tries - 1:
                raise PowerError(f"cannot reach NASA POWER ({type(e).__name__}: {e})")
            time.sleep(2 * 2 ** k)
            continue
        with _lock:
            _status["requests"] += 1
        if r.status_code in (429, 500, 502, 503, 504) and k < tries - 1:
            time.sleep(2 * 2 ** k)
            continue
        return r
    return r


def _fetch_point_remote(lat, lon) -> dict:
    """Try all parameters at once; on a parameter error fall back to one request per parameter."""
    this_year = dt.date.today().year
    last_err = None
    for end in (this_year, this_year - 1, this_year - 2):
        r = _get(lat, lon, PARAMS, end)
        if r.status_code == 200:
            return r.json()
        last_err = f"HTTP {r.status_code}: {r.text[:300]}"
        if r.status_code != 422:
            break
        # 422 can mean "end year not available" or "unknown parameter": test with one safe parameter
        probe = _get(lat, lon, ["T2M"], end)
        if probe.status_code != 200:
            continue  # the year is the problem -> try an earlier end year
        merged = probe.json()
        for p in PARAMS:
            if p == "T2M":
                continue
            rp = _get(lat, lon, [p], end)
            if rp.status_code == 200:
                j = rp.json()
                merged["properties"]["parameter"].update(j["properties"]["parameter"])
                merged.setdefault("parameters", {}).update(j.get("parameters", {}))
        return merged
    raise PowerError(f"NASA POWER request failed for {lat},{lon}: {last_err}")


def fetch_point(lat, lon, offline_fallback=None) -> dict:
    p = _cache_path(lat, lon)
    if p.exists():
        with _lock:
            _status["cache_hits"] += 1
        return json.loads(p.read_text())
    if offline_fallback is not None:
        return offline_fallback(lat, lon)
    j = _fetch_point_remote(lat, lon)
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(j))
    tmp.replace(p)
    return j


def point_frame(j: dict) -> tuple[pd.DataFrame, dict]:
    fill = j.get("header", {}).get("fill_value", -999)
    cols, units = {}, {}
    for name, series in j["properties"]["parameter"].items():
        items = [(k, v) for k, v in series.items() if not k.endswith("13")]
        idx = pd.to_datetime([k[:4] + "-" + k[4:6] + "-01" for k, _ in items])
        vals = np.array([np.nan if (v is None or v == fill or v <= -990) else v for _, v in items], float)
        cols[name] = pd.Series(vals, index=idx)
        units[name] = j.get("parameters", {}).get(name, {}).get("units", "")
    df = pd.DataFrame(cols).sort_index()
    return df, units


def region_frame(points, offline_fallback=None):
    """cos(lat)-weighted mean over points. Returns (DataFrame, units, info)."""
    with ThreadPoolExecutor(WORKERS) as ex:
        js = list(ex.map(lambda pt: fetch_point(pt[0], pt[1], offline_fallback), points))
    frames, units = [], {}
    for j in js:
        df, u = point_frame(j)
        frames.append(df)
        units.update({k: v for k, v in u.items() if v})
    w = np.cos(np.deg2rad([pt[0] for pt in points]))
    allcols = sorted(set().union(*[f.columns for f in frames]))
    idx = frames[0].index
    for f in frames[1:]:
        idx = idx.union(f.index)
    out = {}
    for c in allcols:
        M = np.vstack([f.reindex(idx)[c].values if c in f else np.full(len(idx), np.nan) for f in frames])
        W = np.where(np.isnan(M), 0, w[:, None])
        with np.errstate(invalid="ignore"):
            out[c] = np.where(W.sum(0) > 0, np.nansum(M * W, 0) / np.where(W.sum(0) == 0, 1, W.sum(0)), np.nan)
    df = pd.DataFrame(out, index=idx)
    # drop trailing months where everything is missing (current year not yet published)
    df = df.loc[: df.dropna(how="all").index.max()] if df.notna().any().any() else df
    synthetic = any(j.get("synthetic") for j in js)
    return df, units, dict(points=[list(p) for p in points], synthetic=synthetic)


def status():
    return dict(_status)
