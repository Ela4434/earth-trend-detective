"""Plain-language text for the left sidebar ("Case notes") and result cards."""
from __future__ import annotations
import math

UNIT_MAP = {"C": "°C", "1": "", "kW-hr/m^2/day": "kWh/m²/day", "MJ/m^2/day": "MJ/m²/day", "W/m^2": "W/m²",
            "mm/day": "mm/day", "%": "%", "g/kg": "g/kg", "m/s": "m/s", "kPa": "kPa", "dimensionless": ""}


def unit(u: str | None) -> str:
    if not u:
        return ""
    return UNIT_MAP.get(u, u)


def fnum(x, sig=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "–"
    if x == 0:
        return "0"
    mag = math.floor(math.log10(abs(x)))
    dec = max(0, sig - 1 - mag)
    return f"{x:,.{min(dec, 6)}f}"


def fsigned(x, sig=3):
    s = fnum(x, sig)
    return s if s.startswith("-") or s in ("–", "0") else "+" + s


def fp(p):
    if p is None or not math.isfinite(p):
        return "–"
    if p < 0.001:
        return "< 0.001"
    return f"= {p:.3f}" if p < 0.1 else f"= {p:.2f}"


def amount(res, u, pct_ok):
    if pct_ok and res.get("pct") is not None and abs(res["pct"]) >= 0.5:
        return f"{abs(res['pct']):.0f} %"
    return f"{fnum(abs(res['total']))} {u}".strip()


def trend_sentence(label, res, u, pct_ok):
    if not res.get("ok"):
        return f"{label}: not enough data ({res.get('reason', '')})."
    y0 = res["start"][:4]
    p = max(res["seasonal_kendall"]["p"], res["hamed_rao"]["p"])
    if res["verdict"] == "NO TREND":
        return f"{label}: no clear trend since {y0} (p {fp(p)}). There is not enough evidence of a steady long-term increase or decrease."
    verb = "rose" if res["sen_dec"] > 0 else "fell"
    core = (f"{label} {verb} by {amount(res, u, pct_ok)} since {y0} "
            f"({fsigned(res['sen_dec'])} {u} per decade, p {fp(p)})").replace("  ", " ")
    if res["verdict"] == "PROVEN":
        return core + "."
    extra = ""
    if res.get("more_years_needed") is not None and res["more_years_needed"] > 0:
        extra = f" At this pace about {res['more_years_needed']:.0f} more years of data would settle it."
    return core + ", but the evidence is not strong enough yet." + extra


def relation_sentence(xl, yl, f):
    tau = f["tau_detrended"] if f["kind"] == "linked" else f["tau"]
    q = f["q_detrended"] if f["kind"] == "linked" else f["q"]
    direction = "higher" if tau > 0 else "lower"
    s = (f"When {xl.lower()} is higher than usual, {yl.lower()} tends to be {direction} "
         f"(Kendall τ = {tau:+.2f}, FDR q {fp(q)}).")
    if f["kind"] == "linked":
        s += " This holds after removing both long-term trends, so they really move together."
    else:
        s += " Only true with the long-term trends kept: they share a trend, but don't co-vary year to year."
    return s


def strength(tau):
    a = abs(tau)
    return "strong" if a >= 0.5 else "moderate" if a >= 0.3 else "weak" if a >= 0.1 else "negligible"


def case_notes(region, period, results, labels, units, pct_ok, info, discover=None, label_of=None):
    ok = {k: r for k, r in results.items() if r.get("ok")}
    proven = sorted([k for k, r in ok.items() if r["verdict"] == "PROVEN"], key=lambda k: -abs(ok[k].get("snr") or 0))
    susp = [k for k, r in ok.items() if r["verdict"] == "SUSPICIOUS"]
    none = [k for k, r in ok.items() if r["verdict"] == "NO TREND"]
    head = (f"{region['name']}, {period}: {len(proven)} of {len(ok)} measurements show clear trends, "
            f"{len(susp)} show possible trends, and {len(none)} have no clear trend.")
    notes = dict(headline=head,
                 what=[trend_sentence(labels[k], ok[k], unit(units.get(k)), pct_ok.get(k, True)) for k in proven],
                 ongoing=[trend_sentence(labels[k], ok[k], unit(units.get(k)), pct_ok.get(k, True)) for k in susp],
                 stable=[labels[k] for k in none],
                 where=(f"{region['name']}: the series was exported from NASA Giovanni for this selected area."
                        if info.get("source_type") == "Giovanni CSV" else
                        f"{region['name']}: NASA POWER was sampled at {len(info.get('points', []))} locations within this area."
                        + (" Uploaded Giovanni series are additional measurements for their own selected boxes." if info.get("giovanni_files") else "")),
                 significance=("A clear trend means several statistical tests agree and the evidence is strong. "
                               "A possible trend means the pattern is interesting but the evidence is weaker. "
                               "No clear trend means the tests cannot confirm a steady direction. "
                               "We also account for seasons, values that are similar from month to month, and the number of tests performed."),
                 relations=[], caveats=[])
    if discover:
        lab = label_of or (lambda v: labels.get(v, v))
        notes["relations"] = [relation_sentence(lab(f["x"]), lab(f["y"]), f) for f in discover["found"][:5]]
    if info.get("synthetic"):
        notes["caveats"].append("SYNTHETIC DEMO DATA — the numbers are invented for testing the app. Turn demo mode off before presenting.")
    ur = [labels[k] for k, r in ok.items() if r.get("unit_root_warning")]
    if ur:
        notes["caveats"].append("Unit-root warning (ADF test) for " + ", ".join(ur[:4]) +
                                ": these series wander like a random walk, so trend p-values may be too optimistic.")
    cps = [f"{labels[k]} ({r['pettitt']['year']})" for k, r in ok.items()
           if r["pettitt"]["year"] and r["pettitt"].get("step_better") and r["verdict"] != "NO TREND"]
    if cps:
        notes["caveats"].append("Looks more like a sudden step than a smooth trend (Pettitt test + model comparison): " + ", ".join(cps[:4]) +
                                ". Check for a data or instrument change in that year.")
    if len(info.get("points", [])) == 1:
        notes["caveats"].append("Small region: represented by a single grid point.")
    notes["caveats"].append("A trend says what changed, not why. Correlations between variables are not proof of cause.")
    return notes
