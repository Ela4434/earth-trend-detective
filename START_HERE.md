# EarthTrend Explorer — getting started

This is a real, locally runnable website. It has a **3D globe**, simple explanations, time-series charts and statistical trend tests.

## Start on a Mac (recommended)

1. Unzip the project.
2. Right-click `START_ON_MAC.command` and select **Open**. If macOS warns that it is from an unidentified developer, open it using Terminal instead (see below). You need Python 3 installed and an internet connection for the first installation.
3. Open **http://localhost:8000** in Chrome. You can drag the globe to rotate it or search for a country.

If the Mac script doesn't run, open the Terminal application, move to the extracted `earthtrend-explorer` folder, and run:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m uvicorn backend.app:app --port 8000
```

Then open **http://localhost:8000**. Stop the server with **Control+C** in Terminal.

Open **Findings** to browse notebook chapters or search for a topic. **Overview** shows the selected region. Open **Country profile** for geographic facts and choose a variable to reveal its chart. Both sidebars can be closed to give the globe more space.

Population and surface area are from the World Bank, with observation years shown. Capital data are from [mledoze/countries](https://github.com/mledoze/countries) under the accompanying ODbL notice.

## Use real NASA data

**Automatic:** Select a country or search for a location to load NASA POWER climate data. The first request requires internet and can take a while. NASA POWER combines satellite observations, reanalysis and models, depending on the variable. Its data are **not** Giovanni exports.

**NASA Giovanni:** Open <https://giovanni.gsfc.nasa.gov/giovanni/>. Choose a variable, a date range and a bounding box. Set **Select Plot → Time Series, Area-Averaged**. Plot, download the resulting **CSV**, and use the site's **+ Giovanni data** button. The CSV is analyzed as a separate Giovanni measurement and its bounding box appears in search.

You can also place downloaded Giovanni CSVs into `data/giovanni/` and restart the app. CSVs need `time,<variable>` rows plus a Giovanni bounding box (either in metadata or the file name).

**Important:** The website does not directly authenticate to Giovanni's cloud API. It supports Giovanni CSV imports and a separate automatic NASA POWER source. Eight original Giovanni exports are bundled for the Carpathian Basin and Iberian study areas. **Try a study area** opens the Carpathian Basin. The `tests/fixtures` folder is for tests.

## If you are offline

For UI testing only (data are **invented**, not scientific results):

```bash
TREND_GLOBE_DEMO=1 python -m uvicorn backend.app:app --port 8000
```

A red warning banner labels the test data. **Never present these numbers as a real finding.**

## What the results mean

- **Clear trend:** Multiple statistical checks agree the data show a consistent change.
- **Possible trend:** A pattern appears, but the evidence is weaker.
- **No clear trend:** The data don't establish a steady rise or fall.
- **Change per decade:** The estimated amount the variable increases or decreases in ten years.
- **The graph:** The solid line shows the smoothed pattern; the dashed line shows the longer-term trend.

These tests identify patterns, not causes. Data quality, measurement methods, missing records, and region size matter when interpreting results.

## For developers

The app uses FastAPI/Python, Plotly, globe.gl and a local world atlas. The backend calculates trends using regression, Seasonal Kendall, a Hamed–Rao test, Sen's slope and multiple-testing corrections. Country averages use a limited set of NASA POWER grid points and are **approximate**, not detailed pixel-by-pixel national averages. Read `README.md` for the existing technical documentation.

Run integration checks:

```bash
python tests/test_api.py
python tests/test_stats.py
```

The second test performs slow simulations.

## Publish with Render

Push the contents of this project folder to the root of the GitHub repository,
including `render.yaml`. In the [Render Dashboard](https://dashboard.render.com/),
choose **New → Blueprint**, connect that repository, and deploy the free web
service. Render provides the public `onrender.com` URL. Check the homepage and
**Try a study area** in a private browser window. The free service can take
longer to respond after inactivity.

On the hosted site, visitors can use the bundled Giovanni data and NASA POWER
region analyses. CSV upload and globe-wide batch scans remain available locally
because those actions can modify shared data or start a server-wide job.
