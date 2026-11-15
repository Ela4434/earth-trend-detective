"""Minimal, dependency-free implementations of Mann–Kendall routines.

Used only if optional pymannkendall cannot be imported. The correlated seasonal
approximation uses an annual-mean rank test; it is intentionally conservative
for persistent seasonal signals but is NOT a Hirsch–Slack exact estimator.
Install pymannkendall for the full published methods and challenge reporting.
"""
from collections import namedtuple
import numpy as np
from scipy.stats import norm, kendalltau

Result = namedtuple('MKResult', 'trend h p z Tau s var_s slope intercept')

def _calc(a):
    a=np.asarray(a, dtype=float)
    a=a[np.isfinite(a)]
    n=len(a)
    if n<4:
        return Result('no trend',False,1.0,0.0,0.0,0.0,0.0,0.0,float(np.nanmedian(a)) if n else float('nan'))
    i,j=np.triu_indices(n,1)
    diff=np.sign(a[j]-a[i]); S=float(diff.sum())
    unique, counts=np.unique(a,return_counts=True)
    var=(n*(n-1)*(2*n+5)-np.sum(counts*(counts-1)*(2*counts+5)))/18
    z=(S-np.sign(S))/np.sqrt(var) if var>0 else 0.0
    p=float(2*norm.sf(abs(z)))
    tau=float(kendalltau(np.arange(n),a).statistic)
    slopes=(a[j]-a[i])/(j-i)
    slope=float(np.median(slopes))
    intercept=float(np.median(a)-np.median(np.arange(n))*slope)
    return Result(('increasing' if slope>0 else 'decreasing') if p<.05 else 'no trend',p<.05,p,float(z),tau,S,float(var),slope,intercept)

def original_test(x):
    return _calc(x)

def hamed_rao_modification_test(x):
    # Conservative trend test on ranks with an AR(1) effective sample-size correction.
    a=np.asarray(x,float)
    a=a[np.isfinite(a)]
    r=_calc(a)
    if len(a)<6:return r
    t=np.arange(len(a),dtype=float)
    fit=np.polyval(np.polyfit(t,a,1),t)
    res=a-fit
    phi=np.corrcoef(res[:-1],res[1:])[0,1] if np.std(res)>0 else 0.0
    phi=float(np.clip(np.nan_to_num(phi),-.8,.95))
    inflation=max(1.0,(1+phi)/(1-phi))
    z=r.z/np.sqrt(inflation)
    p=float(2*norm.sf(abs(z)))
    return r._replace(z=float(z),p=p,h=p<.05)

def seasonal_test(x,period=12):
    a=np.asarray(x,float)
    s=0.;var=0.;slopes=[];j=0
    for month in range(period):
        vals=a[month::period]
        r=_calc(vals)
        s+=r.s;var+=r.var_s
        v=vals[np.isfinite(vals)]
        if len(v)>1:
            i0,i1=np.triu_indices(len(v),1)
            slopes.extend(((v[i1]-v[i0])/(i1-i0)).tolist())
    z=(s-np.sign(s))/np.sqrt(var) if var>0 else 0
    p=float(2*norm.sf(abs(z)))
    slope=float(np.nanmedian(slopes)) if slopes else 0.
    return Result(('increasing' if slope>0 else 'decreasing') if p<.05 else 'no trend',p<.05,p,float(z),float('nan'),s,var,slope,0.)

def correlated_seasonal_test(x,period=12):
    # Fallback only: compare annual means. Explicitly not the exact Hirsch–Slack statistic.
    a=np.asarray(x,float)
    arr=np.array([np.nanmean(a[k:k+period]) for k in range(0,len(a),period)])
    r=hamed_rao_modification_test(arr)
    return r._replace(slope=r.slope)
