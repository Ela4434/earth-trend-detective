"""Relationships between variables in one region.

Everything works on monthly ANOMALIES (each calendar month's long-term mean removed), so
the shared seasonal cycle can't fake a relationship (summer is hot AND sunny every year).
Two modes:
  anomalies  - "do they move together, including their long-term trends?"
  detrended  - each anomaly series also has its own linear trend removed:
               "do they move together year to year?" Two unrelated variables that both
               trend over 40 years correlate strongly in the first mode but not here.
Tests: Pearson r, Spearman rho, Kendall tau-b. Monthly anomalies are autocorrelated, so
p-values use an effective sample size N_eff = N (1 - r1 r2) / (1 + r1 r2)
(Bretherton et al. 1999, J. Climate 12), r1/r2 = lag-1 autocorrelations.
Multiple testing: Benjamini-Hochberg q-values and Bonferroni across all pairs tested
(Statsig guide: "multiple comparisons ... Bonferroni correction").
Multiple regression: OLS with Newey-West (HAC) errors, standardised coefficients, VIF.
Correlation is not causation; the sentences say "move together", never "causes".
"""
from __future__ import annotations
import itertools
import math
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as st
from statsmodels.stats.outliers_influence import variance_inflation_factor

from .trend import ALPHA, bh_fdr, deseasonalise, decimal_years, hac_lags, lag1, monthly

METHODS = ("pearson", "spearman", "kendall")
MIN_N = 36


def prepare(df: pd.DataFrame, mode="anomalies") -> pd.DataFrame:
    out = {}
    for c in df.columns:
        s = monthly(df[c])
        if s.notna().sum() < MIN_N:
            continue
        a = deseasonalise(s)
        if mode == "detrended":
            t = decimal_years(a.index)
            ok = a.notna().values
            b = np.polyfit(t[ok], a.values[ok], 1)
            a = a - np.polyval(b, t)
        out[c] = a
    return pd.DataFrame(out)


def n_eff(x, y):
    n = len(x)
    r1, r2 = lag1(x), lag1(y)
    ne = n * (1 - r1 * r2) / (1 + r1 * r2) if (1 + r1 * r2) > 0 else n
    return float(min(n, max(3.0, ne)))


def _adj(p, n, ne):
    """Shrink a test's z-score by sqrt(N_eff/N) and recompute the two-sided p."""
    if not np.isfinite(p) or p <= 0:
        p = max(p, 1e-300) if np.isfinite(p) else 1.0
    z = st.norm.isf(p / 2) * math.sqrt(ne / n)
    return float(2 * st.norm.sf(abs(z)))


def pair(x: pd.Series, y: pd.Series) -> dict | None:
    d = pd.concat([x, y], axis=1).dropna()
    n = len(d)
    if n < MIN_N:
        return None
    a, b = d.iloc[:, 0].values, d.iloc[:, 1].values
    if np.std(a) == 0 or np.std(b) == 0:
        return None
    ne = n_eff(a, b)
    r, _ = st.pearsonr(a, b)
    tt = r * math.sqrt(max(ne - 2, 1) / max(1e-12, 1 - r * r))
    p_pear = float(2 * st.t.sf(abs(tt), df=max(ne - 2, 1)))
    rho, p_rho = st.spearmanr(a, b)
    tau, p_tau = st.kendalltau(a, b)
    return dict(n=n, n_eff=round(ne, 1),
                pearson=dict(r=float(r), p=p_pear),
                spearman=dict(r=float(rho), p=_adj(float(p_rho), n, ne)),
                kendall=dict(r=float(tau), p=_adj(float(p_tau), n, ne)))


def lagged(x: pd.Series, y: pd.Series, max_lag=6):
    """Pearson r of x(t) with y(t+lag); positive lag = x leads y. Bonferroni over all lags."""
    best = None
    rows = []
    for L in range(-max_lag, max_lag + 1):
        d = pd.concat([x, y.shift(-L)], axis=1).dropna()
        if len(d) < MIN_N:
            continue
        a, b = d.iloc[:, 0].values, d.iloc[:, 1].values
        r = float(np.corrcoef(a, b)[0, 1])
        ne = n_eff(a, b)
        tt = r * math.sqrt(max(ne - 2, 1) / max(1e-12, 1 - r * r))
        p = float(2 * st.t.sf(abs(tt), df=max(ne - 2, 1)))
        rows.append(dict(lag=L, r=r, p=p))
        if best is None or abs(r) > abs(best["r"]):
            best = dict(lag=L, r=r, p=p)
    if best:
        best["p_bonf"] = min(1.0, best["p"] * len(rows))
    return dict(best=best, all=rows)


def matrix(df: pd.DataFrame, cols, method="kendall", mode="anomalies", skip_same_family=None, with_lags=False):
    A = prepare(df[cols], mode)
    cols = [c for c in cols if c in A.columns]
    rows = []
    for x, y in itertools.combinations(cols, 2):
        if skip_same_family and skip_same_family(x) == skip_same_family(y):
            continue
        r = pair(A[x], A[y])
        if r is None:
            continue
        r.update(x=x, y=y)
        if with_lags:
            r["lag"] = lagged(A[x], A[y])["best"]
        rows.append(r)
    m = len(rows)
    for meth in METHODS:
        q = bh_fdr([r[meth]["p"] for r in rows])
        for r, qq in zip(rows, q):
            r[meth]["q"] = float(qq)
            r[meth]["p_bonf"] = min(1.0, r[meth]["p"] * m)
    for r in rows:
        h = r[method]
        r["significant"] = bool(h["q"] < ALPHA)
    return dict(cols=cols, method=method, mode=mode, tests=m, pairs=rows)


def regression(df: pd.DataFrame, y: str, xs: list[str], mode="anomalies", add_trend=False) -> dict:
    A = prepare(df[[y] + xs], mode).dropna()
    xs = [x for x in xs if x in A.columns]
    if y not in A.columns or not xs or len(A) < MIN_N:
        return dict(ok=False, reason="not enough overlapping months for these variables")
    Z = (A - A.mean()) / A.std(ddof=0)
    X = Z[xs].copy()
    if add_trend:
        X["time trend"] = (decimal_years(Z.index) - decimal_years(Z.index).mean()) / 10
    Xc = sm.add_constant(X)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = sm.OLS(Z[y], Xc).fit(cov_type="HAC", cov_kwds={"maxlags": hac_lags(len(Z))})
        plain = sm.OLS(Z[y], Xc).fit()
    vif = {}
    if X.shape[1] > 1:
        for i, c in enumerate(X.columns):
            vif[c] = float(variance_inflation_factor(X.values, i))
    coefs = []
    for c in X.columns:
        lo, hi = fit.conf_int(ALPHA).loc[c]
        coefs.append(dict(name=c, beta=float(fit.params[c]), p=float(fit.pvalues[c]), ci=[float(lo), float(hi)],
                          vif=vif.get(c)))
    return dict(ok=True, y=y, xs=xs, mode=mode, n=int(len(Z)), r2=float(plain.rsquared), r2_adj=float(plain.rsquared_adj),
                F=float(plain.fvalue), p_F=float(plain.f_pvalue), p_F_hac=float(fit.f_pvalue) if fit.f_pvalue is not None else None,
                coefs=coefs, add_trend=add_trend,
                resid_lag1=lag1(plain.resid.values))


def scatter(df: pd.DataFrame, x: str, y: str, mode="anomalies"):
    A = prepare(df[[x, y]], mode).dropna()
    return dict(x=[round(float(v), 6) for v in A[x].values], y=[round(float(v), 6) for v in A[y].values],
                t=[d.strftime("%Y-%m") for d in A.index])


def discover(df: pd.DataFrame, cols, family, top=8):
    """Auto-detective: every pair, both modes, Kendall tau as the headline test.
    Classifies each significant pair:
      'linked'      significant after removing both trends -> they co-vary month to month / year to year
      'trend-only'  significant only with trends kept -> they share a long-term trend, little more
    """
    a = matrix(df, cols, "kendall", "anomalies", family)
    d = matrix(df, cols, "kendall", "detrended", family, with_lags=False)
    dmap = {(r["x"], r["y"]): r for r in d["pairs"]}
    found = []
    for r in a["pairs"]:
        rd = dmap.get((r["x"], r["y"]))
        sig_a = r["kendall"]["q"] < ALPHA
        sig_d = rd is not None and rd["kendall"]["q"] < ALPHA
        if not (sig_a or sig_d):
            continue
        kind = "linked" if sig_d else "trend-only"
        found.append(dict(x=r["x"], y=r["y"], kind=kind, tau=r["kendall"]["r"], q=r["kendall"]["q"],
                          tau_detrended=rd["kendall"]["r"] if rd else None, q_detrended=rd["kendall"]["q"] if rd else None,
                          pearson=r["pearson"]["r"], n=r["n"], n_eff=r["n_eff"]))
    found.sort(key=lambda f: (f["kind"] != "linked", -abs(f["tau_detrended"] if f["kind"] == "linked" else f["tau"])))
    return dict(found=found[:top], total_found=len(found), tests=a["tests"] + d["tests"],
                matrix_anomalies=a, matrix_detrended=d)
