# Trend Globe — Be an Earth System Trend Detective

NASA Space Apps Challenge 2026 · *“create a tool that can select variables measured by Earth science missions or produced by Earth science models, visualize how they change over time, and determine what is changing, where it is changing, how much it is changing, and if these changes are what scientists call ‘significant.’”*

A 3-D globe web app. Search any country, a study box or a point → the app pulls 16 NASA variables for that region (plus any Giovanni CSVs you add), tests every variable for a trend, writes plain-language case notes in the left sidebar, and finds which variables move together, all with autocorrelation and multiple-testing corrections running in the Python backend.

## Run it (5 minutes)

```bash
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python tools/check_power.py          # 1) can this laptop reach NASA POWER? every line should say OK
uvicorn backend.app:app --port 8000  # 2) start the app  (or: ./run.sh  /  run.bat)
```
Open **http://localhost:8000**. Type “Spain” and press Enter.

* The first visit to a region downloads its NASA data (about 10–30 s); after that everything is cached in `data/cache/` and instant, also offline.
* **Offline / UI testing:** `TREND_GLOBE_DEMO=1 uvicorn backend.app:app --port 8000` uses **synthetic** data and shows a red banner. Never present those numbers.
* `http://localhost:8000/#Hungary` opens a region directly (handy for the demo video).

## What it does

| Area | Feature |
|---|---|
| **Globe** | 3-D globe with all countries. Click a country, or search: a country (“Spain”, “usa”), a study box (“Sahel”, “Carpathian Basin”), a point (`40.4, -3.7`) or a box in Giovanni order (`-9,36,3,43`). Yellow dots = the NASA grid points averaged for the region. |
| **Scan the globe** | Pick a variable → trend for every country, coloured red (rising) / blue (falling), faded when not significant (FDR across countries). Shows mirror trends at a glance. First scan downloads ~1,000 points (15–40 min, run it once before the demo); after that it's ~2 min. |
| **Case notes (left)** | The four challenge questions in plain words: what is changing and how much, investigation ongoing (suspicious + years still needed), no trend, where, is it significant, what moves together, caveats. |
| **Variables tab** | One card per variable: latest year, long-term mean, trend per decade, total change, verdict stamp (PROVEN / SUSPICIOUS / NO TREND), p- and q-values, sparkline. Click → full investigation: raw series, anomalies with Sen and OLS lines, all tests, seasonal split, years-to-detect meter. |
| **Relationships tab** | **Auto-detective:** every pair of variables, Kendall τ, with and without long-term trends → “linked year to year” vs “shared trend only”. **Own selection:** pick variables, choose Kendall / Spearman / Pearson, anomalies or detrended, optional lags ±6 months → sentences, heatmap with significance stars, table with Bonferroni and FDR, scatter plot. **Multiple OLS** with standardised β, HAC errors and VIF. |
| **Compare regions** | Same variable in two regions + test of the difference → flags **MIRROR TREND** (rising here, falling there). |
| **Header** | Period (from/to), **exclude years** (e.g. `1991-1993` for Pinatubo: does the trend survive without the event?), **＋ Giovanni CSV** upload. |
| **Method tab** | The challenge text, every test with its source in the challenge resources. |

## Data

* **NASA POWER** monthly API (`power.larc.nasa.gov`, no login): MERRA-2 meteorology (1981–) and CERES/GEWEX satellite radiation (1984–). 16 variables: temperature (mean, max, min, skin), precipitation, humidity (relative, specific), soil wetness (root zone, surface), evapotranspiration, sunlight (all-sky, clear-sky), cloud amount, aerosol optical depth, wind, surface pressure. Region value = cos(latitude)-weighted mean of up to 8 grid points (countries) or 16 points (boxes) inside the region.
* **Giovanni** (the challenge's main tool): download *Time Series, Area-Averaged* → *CSV* and drop it into `data/giovanni/` or use **＋ Giovanni CSV**. The app reads the bounding box from the file, creates a region “Giovanni box W,S,E,N”, and shows the file as an extra variable there and in any matching country/box. Kelvin is converted to °C. Every POWER card also names the matching Giovanni variable, so you can cross-check numbers for the same box.
* **Hydrology Time Series / AppEEARS / NOAA CDO:** any CSV with a `time,value` table and a `Bounding Box:,"W,S,E,N"` line in the header is read the same way (add the line by hand if the tool doesn't write it).

## Statistics (backend/stats)

| Question | Method | Source in the challenge resources |
|---|---|---|
| Is there a trend? (verdict) | **Seasonal Kendall** (Hirsch-Slack, period 12) **and** **Hamed-Rao Mann-Kendall** on yearly anomalies, both p < 0.05 after **Benjamini-Hochberg FDR** across the region's variables, and Sen and OLS agree on the direction | Medium trend-test article (Mann-Kendall, Seasonal Kendall); Statsig guide (multiple comparisons) |
| How much? | **Sen's slope** per decade + 95 % CI, total and % change, change ÷ noise | |
| Linear trend | **OLS slope t-test**, **ANOVA F-test + R²**, Newey-West (HAC) errors, **yearly-anomaly regression**, **robust (Huber) regression** | Medium article; SeaWiFS chlorophyll paper (NTRS 20090042747) |
| Robustness | **Seasonal split** (May–Oct vs Nov–Apr), **exclude years**, **ADF** unit-root warning, **Pettitt** break + step-vs-trend model comparison | SeaWiFS paper; Medium article |
| Not significant yet? | **Weatherhead et al. (1998)** years of data needed → “investigation ongoing” meter | |
| Relationships | **Pearson, Spearman, Kendall** on anomalies or detrended anomalies; p-values with **effective sample size** (Bretherton 1999); **FDR** + **Bonferroni**; lagged correlation; **multiple OLS** with HAC + VIF | Statsig guide |

**Calibration (run `python tests/test_stats.py`):** on 300 simulated trendless, autocorrelated series, the verdict says PROVEN in ~5 % of cases, as it should, while a naive OLS t-test says “significant” in 17–30 %. Independent autocorrelated series: ~5 % significant correlations with the effective-N correction vs ~38 % naively. Use these numbers on the validity slide.

## Tests

```bash
python tests/test_stats.py   # known answers + false-positive rates (≈1.5 min)
python tests/test_api.py     # every endpoint incl. Giovanni upload, on synthetic data
```

## Project layout

```
backend/app.py              FastAPI: API + serves the frontend
backend/catalog.py          the 16 variables: label, plain-language definition, Giovanni equivalent
backend/regions.py          countries (Natural Earth), study boxes, search, point-in-polygon, grid sampling
backend/sources/power.py    NASA POWER client with disk cache and retries
backend/sources/giovanni_csv.py   Giovanni CSV parser + region matching
backend/sources/demo.py     synthetic data (TREND_GLOBE_DEMO=1 only)
backend/stats/trend.py      trend tests        backend/stats/relations.py  correlation / regression / auto-detective
backend/narrative.py        plain-language sentences for the case notes
frontend/                   index.html, app.js, style.css; vendor/ = globe.gl, Plotly, topojson (local, no CDN)
data/cache/                 downloaded NASA POWER points (safe to delete)   data/giovanni/  your Giovanni CSVs
tools/check_power.py        connectivity check
```

## Known limits (say them before the judges do)

* Region values are averages of grid points, not of every pixel; for a big country 8 points is coarse. Boxes use up to 16.
* POWER's meteorology comes from MERRA-2, a reanalysis (model + observations). Trends in reanalyses can contain artefacts from changes in the assimilated observations (e.g. MERRA-2 aerosol inputs around 2000) — use **exclude years** / a later start year to check.
* Correlation is not causation. The app never says “causes”.
* The globe scan's FDR is across countries; regions of very different size are treated equally.
