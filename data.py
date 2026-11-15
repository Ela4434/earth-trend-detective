"""Assemble all variables for a region: NASA POWER (+ Giovanni CSVs that match the region)."""
from __future__ import annotations
import os
import threading

import pandas as pd

from . import catalog, regions as R
from .sources import power, giovanni_csv, demo

DEMO = os.environ.get("TREND_GLOBE_DEMO", "0") == "1"
POINT_CAP = int(os.environ.get("POINT_CAP", "8"))   # countries (also used by the globe scan)
BOX_CAP = int(os.environ.get("BOX_CAP", "16"))      # study boxes / Giovanni boxes / points: denser, closer to Giovanni
_mem: dict = {}
_lock = threading.Lock()


def region_data(rid: str):
    """Return a data frame with explicit provenance.

    Giovanni boxes use their exact, area-averaged Giovanni series as the source.
    Other regions use NASA POWER grid-point averages, with matching Giovanni
    series added as separate columns. When POWER is unavailable but a Giovanni
    file matches the region, we still allow analysis of the imported file.
    """
    region = R.get(rid)
    with _lock:
        if rid in _mem:
            return _mem[rid]
    gfiles = [g for g in giovanni_csv.for_region(region) if g.get("series") is not None]
    only_giovanni = rid.startswith("g-") and len(gfiles) > 0
    points = []
    power_warning = None
    if only_giovanni:
        df, units = pd.DataFrame(), {}
        info = {"points": [], "synthetic": False, "power_used": False}
    else:
        points = R.sample_points(region, POINT_CAP if region["kind"] == "country" else BOX_CAP)
        try:
            df, units, info = power.region_frame(points, offline_fallback=demo.point_json if DEMO else None)
            info["power_used"] = True
        except power.PowerError as e:
            if not gfiles:
                raise
            df, units = pd.DataFrame(), {}
            power_warning = "NASA POWER is temporarily unavailable. Showing your imported Giovanni data only."
            info = {"points": [], "synthetic": False, "power_used": False}

    labels = {k: v["label"] for k, v in catalog.VARIABLES.items()}
    pct_ok = {k: v["pct"] for k, v in catalog.VARIABLES.items()}
    cats = {k: v["category"] for k, v in catalog.VARIABLES.items()}
    keep = [c for c in catalog.VARIABLES if c in df.columns and df[c].notna().sum() > 0]
    df = df[keep]
    attached = []
    for g in gfiles:
        s, u = g["series"], g["unit"]
        if u == "K":
            s, u = s - 273.15, "°C"
        vid = g["id"]
        # Different products and exports may use the same original column name.
        if vid in df.columns:
            vid = vid + ":" + str(len(attached) + 1)
        df = df.join(s.rename(vid), how="outer")
        labels[vid] = g["label"]
        units[vid] = u
        pct_ok[vid] = not (u in ("K", "C", "°C") or "temperature" in g["label"].lower())
        cats[vid] = "Giovanni"
        attached.append(g["file"])
    if df.empty:
        raise power.PowerError("No observations were found for this region")
    info.update(region=R.public(region), giovanni_files=attached, categories=cats,
                demo=DEMO and bool(info.get("synthetic")), power_warning=power_warning,
                source_type="Giovanni CSV" if not info["power_used"] else
                            ("NASA POWER + Giovanni CSV" if attached else "NASA POWER"))
    out = (df.sort_index(), labels, units, pct_ok, info)
    with _lock:
        _mem[rid] = out
    return out


def clear():
    with _lock:
        _mem.clear()
