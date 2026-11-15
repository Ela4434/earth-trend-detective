"""Trend statistics for one monthly series.

Methods, and where they come from in the challenge resources:
  * OLS linear trend with slope t-test (Medium article "Statistical test methods to detect trends")
    and the ANOVA F-test + R² (SeaWiFS chlorophyll paper, NASA NTRS 20090042747). We add
    Newey-West (HAC) standard errors because monthly anomalies are autocorrelated.
  * Yearly anomalies regressed on time: subtract each calendar month's long-term mean,
    average within each year, regress on year (SeaWiFS paper, following Gregg et al. 2005).
  * Robust linear regression to limit outlier influence (SeaWiFS paper). Here: Huber RLM.
  * Mann-Kendall test and Seasonal Kendall test, period 12 (Medium article). The verdict uses the
    Hirsch-Slack version of the Seasonal Kendall test (corrects for correlation between months)
    and the Hamed-Rao modified Mann-Kendall on yearly anomalies (autocorrelation-corrected).
  * Sen's slope (with the seasonal Kendall) and its 95 % CI (Gilbert 1987) on annual means.
  * Augmented Dickey-Fuller unit-root test (Medium article): warns when a series may wander
    like a random walk, where trend tests can be misleading.
  * Seasonal split, warm vs cold half-year (SeaWiFS paper splits winter/summer).
  * Event sensitivity: exclude years and refit (SeaWiFS paper excludes the 2000 Po flood).
  * Pettitt change-point test, and Weatherhead et al. (1998) years-to-detect.
Significance level α = 0.05, confidence intervals 95 % (Statsig guide).
"""
from __future__ import annotations
import math
import warnings

import numpy as np
import pandas as pd
try:
    import pymannkendall as mk
    MK_FALLBACK = False
except ImportError:
    from . import mk_fallback as mk
    MK_FALLBACK = True
import statsmodels.api as sm
from scipy import stats as st
from statsmodels.tsa.stattools import adfuller

ALPHA = 0.05
MIN_MONTHS = 60


def monthly(s: pd.Series) -> pd.Series:
    """Regular month-start index (gaps become NaN)."""
    s = s.dropna()
    if s.empty:
        return s
    idx = pd.date_range(s.index.min(), s.index.max(), freq="MS")
    return s.reindex(idx)


def deseasonalise(s: pd.Series) -> pd.Series:
    return s - s.groupby(s.index.month).transform("mean")


def decimal_years(idx) -> np.ndarray:
    idx = pd.DatetimeIndex(idx)
    return idx.year.values + (idx.month.values - 0.5) / 12.0


def lag1(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 4:
        return 0.0
    x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 else 0.0


def hac_lags(n: int) -> int:
    return max(1, int(math.floor(4 * (n / 100.0) ** (2 / 9))))


def ols(y: np.ndarray, t: np.ndarray) -> dict:
    ok = ~np.isnan(y)
    y, t = y[ok], t[ok]
    X = sm.add_constant(t)
    plain = sm.OLS(y, X).fit()
    hac = sm.OLS(y, X).fit(cov_type="HAC", cov_kwds={"maxlags": hac_lags(len(y))})
    ci = plain.conf_int(ALPHA)[1]
    hci = hac.conf_int(ALPHA)[1]
    return dict(slope=float(plain.params[1]), intercept=float(plain.params[0]), se=float(plain.bse[1]),
                t=float(plain.tvalues[1]), p=float(plain.pvalues[1]), r2=float(plain.rsquared),
                F=float(plain.fvalue), p_F=float(plain.f_pvalue), ci=[float(ci[0]), float(ci[1])],
                p_hac=float(hac.pvalues[1]), ci_hac=[float(hci[0]), float(hci[1])], n=int(len(y)),
                resid=plain.resid)


def sen_ci(y: np.ndarray, alpha=ALPHA):
    """Sen's slope with confidence interval (Gilbert 1987), per step."""
    y = np.asarray(y, float)
    y = y[~np.isnan(y)]
    n = len(y)
    i, j = np.triu_indices(n, 1)
    d = np.sort((y[j] - y[i]) / (j - i))
    _, cnt = np.unique(y, return_counts=True)
    var_s = (n * (n - 1) * (2 * n + 5) - np.sum(cnt * (cnt - 1) * (2 * cnt + 5))) / 18
    c = st.norm.ppf(1 - alpha / 2) * math.sqrt(var_s)
    N = len(d)
    lo, hi = int(np.floor((N - c) / 2)), int(np.ceil((N + c) / 2))
    return float(np.median(d)), float(d[max(lo, 0)]), float(d[min(hi, N - 1)])


def _hamed_rao(x):
    x = x[~np.isnan(x)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = mk.hamed_rao_modification_test(x)
    p = float(r.p)
    if not np.isfinite(p):  # rare negative corrected variance -> fall back to plain MK
        p = float(mk.original_test(x).p)
    return p, float(r.z) if np.isfinite(r.z) else float("nan")


def pettitt(x: np.ndarray):
    """Pettitt (1979) change point. Returns (index of last point before change, approx p)."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 8:
        return None, 1.0
    sgn = np.sign(x[None, :] - x[:, None])  # sgn[i,j] = sign(x_j - x_i)
    U = np.array([sgn[: k + 1, k + 1:].sum() for k in range(n - 1)])
    k = int(np.argmax(np.abs(U)))
    K = abs(U[k])
    p = min(1.0, 2 * math.exp(-6 * K ** 2 / (n ** 3 + n ** 2)))
    return k, p


def years_to_detect(resid: np.ndarray, slope_per_year: float) -> float | None:
    """Weatherhead et al. (1998), J. Geophys. Res. 103(D14): years of monthly data needed to detect
    a trend of this size (5 % level, 90 % probability), given the noise size and its autocorrelation."""
    r = np.asarray(resid, float)
    r = r[~np.isnan(r)]
    if len(r) < 24 or slope_per_year == 0:
        return None
    sd = float(np.std(r, ddof=1))
    phi = min(max(lag1(r), -0.95), 0.95)
    return float((3.3 * sd / abs(slope_per_year) * math.sqrt((1 + phi) / (1 - phi))) ** (2 / 3))


def annual_means(s: pd.Series, min_months=9) -> pd.Series:
    g = s.groupby(s.index.year)
    a = g.mean()
    return a[g.count() >= min_months]


def verdict(p_sk, p_hr, p_ols_hac, agree, q=None):
    p_main = max(p_sk, p_hr)
    if q is not None:
        p_main = max(p_main, q)
    if p_main < ALPHA and agree:
        return "PROVEN"
    if min(p_sk, p_hr, p_ols_hac) < 0.10:
        return "SUSPICIOUS"
    return "NO TREND"


def analyse(s: pd.Series, start_year=None, end_year=None, exclude=None, pct_ok=True, detail=False) -> dict:
    """All trend statistics for one monthly series. exclude: list of [y0, y1] year ranges to drop."""
    s = monthly(s)
    if start_year:
        s = s[s.index.year >= int(start_year)]
    if end_year:
        s = s[s.index.year <= int(end_year)]
    if exclude:
        for y0, y1 in exclude:
            s = s.where(~((s.index.year >= int(y0)) & (s.index.year <= int(y1))))
    s = monthly(s) if s.notna().any() else s
    nvalid = int(s.notna().sum())
    if nvalid < MIN_MONTHS:
        return dict(ok=False, reason=f"only {nvalid} months of data in the chosen period (need {MIN_MONTHS})", n=nvalid)
    anom = deseasonalise(s)
    t = decimal_years(s.index)
    years = (s.index[-1] - s.index[0]).days / 365.25 + 1 / 12

    o = ols(anom.values, t)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sk_plain = mk.seasonal_test(s.values, period=12)
        sk = mk.correlated_seasonal_test(s.values, period=12)  # Hirsch & Slack (1984)
        mko = mk.original_test(anom.dropna().values)
    p_hr_m, _ = _hamed_rao(anom.values)
    sen = float(sk_plain.slope)  # seasonal Sen slope, per year

    ann_anom = annual_means(anom)
    # second main test: Hamed-Rao Mann-Kendall on yearly anomalies. In simulations (tests/test_stats.py)
    # the monthly Hamed-Rao and plain Seasonal Kendall gave 9-23 % false positives on autocorrelated
    # trendless data; this pair stays near 5 %.
    p_hr, z_hr = _hamed_rao(ann_anom.values) if len(ann_anom) >= 8 else (1.0, float("nan"))
    ann_res = None
    if len(ann_anom) >= 5:
        oy = ols(ann_anom.values, ann_anom.index.values.astype(float))
        _, lo, hi = sen_ci(ann_anom.values)
        ann_res = dict(slope=oy["slope"], p=oy["p"], r2=oy["r2"], p_F=oy["p_F"], n=oy["n"], sen_lo=lo, sen_hi=hi)

    ok = ~np.isnan(anom.values)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rlm = sm.RLM(anom.values[ok], sm.add_constant(t[ok]), M=sm.robust.norms.HuberT()).fit()
        try:
            adf_p = float(adfuller(anom.dropna().values, regression="ct", autolag="AIC")[1])
        except Exception:
            adf_p = float("nan")
    k, p_pet = pettitt(ann_anom.values) if len(ann_anom) >= 8 else (None, 1.0)
    cp_year = int(ann_anom.index[k + 1]) if k is not None and p_pet < ALPHA else None
    # Pettitt also fires on any smooth trend. Only call it a STEP if a two-level step fits the
    # annual anomalies clearly better than a straight line (BIC lower by > 2).
    step_better = False
    if cp_year is not None:
        y = ann_anom.values
        n_a = len(y)
        xx = ann_anom.index.values.astype(float)
        sse_lin = float(np.sum((y - np.polyval(np.polyfit(xx, y, 1), xx)) ** 2))
        a1, a2 = y[: k + 1], y[k + 1:]
        sse_step = float(np.sum((a1 - a1.mean()) ** 2) + np.sum((a2 - a2.mean()) ** 2))
        bic_lin = n_a * math.log(max(sse_lin, 1e-300) / n_a) + 2 * math.log(n_a)
        bic_step = n_a * math.log(max(sse_step, 1e-300) / n_a) + 3 * math.log(n_a)
        step_better = bic_step < bic_lin - 2

    agree = np.sign(sen) == np.sign(o["slope"]) and sen != 0
    total = sen * years
    mean = float(s.mean())
    base = mean - sen * years / 2
    pct = 100 * total / base if (pct_ok and base > 0) else None
    noise = float(np.nanstd(anom.values, ddof=1))
    ytd = years_to_detect(o["resid"], sen if sen != 0 else o["slope"])
    last_year = annual_means(s, 12)
    clim = s.groupby(s.index.month).mean()
    res = dict(
        ok=True, n=nvalid, start=str(s.index[0].date()), end=str(s.index[-1].date()), years=round(years, 2),
        mean=mean, latest_month=dict(date=str(s.dropna().index[-1].date())[:7], value=float(s.dropna().iloc[-1])),
        latest_year=(dict(year=int(last_year.index[-1]), value=float(last_year.iloc[-1])) if len(last_year) else None),
        first_year=(dict(year=int(last_year.index[0]), value=float(last_year.iloc[0])) if len(last_year) else None),
        sen_dec=sen * 10, total=total, pct=pct, noise_sd=noise, snr=(total / noise if noise > 0 else None),
        ols=dict(slope_dec=o["slope"] * 10, p=o["p"], t=o["t"], r2=o["r2"], F=o["F"], p_F=o["p_F"],
                 ci_dec=[c * 10 for c in o["ci"]], p_hac=o["p_hac"], ci_hac_dec=[c * 10 for c in o["ci_hac"]]),
        mk=dict(p=float(mko.p), tau=float(mko.Tau), z=float(mko.z)),
        hamed_rao=dict(p=p_hr, z=z_hr, on="yearly anomalies"),
        hamed_rao_monthly=dict(p=p_hr_m),
        seasonal_kendall=dict(p=float(sk.p), z=float(sk.z), tau=float(sk.Tau), variant="Hirsch-Slack"),
        seasonal_kendall_plain=dict(p=float(sk_plain.p), z=float(sk_plain.z)),
        annual=(dict(slope_dec=ann_res["slope"] * 10, p=ann_res["p"], r2=ann_res["r2"], p_F=ann_res["p_F"], n=ann_res["n"],
                     sen_ci_dec=[ann_res["sen_lo"] * 10, ann_res["sen_hi"] * 10]) if ann_res else None),
        robust=dict(slope_dec=float(rlm.params[1]) * 10, p=float(rlm.pvalues[1])),
        adf_p=adf_p, unit_root_warning=bool(np.isfinite(adf_p) and adf_p > ALPHA),
        pettitt=dict(p=p_pet, year=cp_year, step_better=bool(step_better)),
        lag1=lag1(o["resid"]), years_to_detect=ytd,
        agree=bool(agree),
        spark=dict(year=[int(y) for y in ann_anom.index], anomaly=_clean(ann_anom.values)),
    )
    res["verdict"] = verdict(res["seasonal_kendall"]["p"], p_hr, o["p_hac"], agree)
    if res["verdict"] != "PROVEN" and ytd is not None:
        res["more_years_needed"] = max(0.0, ytd - years)
    if detail:
        warm = s[s.index.month.isin([5, 6, 7, 8, 9, 10])]
        cold = s[~s.index.month.isin([5, 6, 7, 8, 9, 10])]
        split = {}
        for name, part in [("May–Oct", warm), ("Nov–Apr", cold)]:
            a = part.groupby(part.index.year).mean()[part.groupby(part.index.year).count() >= 5]
            if len(a) >= 8:
                oy = ols(a.values, a.index.values.astype(float))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    mm = mk.original_test(a.values)
                split[name] = dict(slope_dec=oy["slope"] * 10, p_ols=oy["p"], p_mk=float(mm.p), sen_dec=float(mm.slope) * 10)
        res["seasonal_split"] = split
        ts = pd.DatetimeIndex(s.index)
        b = sen
        a0 = float(np.nanmedian(anom.values - b * t))
        res["chart"] = dict(
            t=[d.strftime("%Y-%m") for d in ts], value=_clean(s.values), anomaly=_clean(anom.values),
            anomaly_12m=_clean(anom.rolling(12, center=True, min_periods=9).mean().values),
            sen_line=_clean(a0 + b * t), ols_line=_clean(o["intercept"] + o["slope"] * t),
            annual=dict(year=[int(y) for y in ann_anom.index], anomaly=_clean(ann_anom.values)),
            climatology=_clean(clim.reindex(range(1, 13)).values),
        )
    return res


def _clean(a):
    return [None if (x is None or not np.isfinite(x)) else round(float(x), 6) for x in np.asarray(a, float)]


def bh_fdr(p) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q
