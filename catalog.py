"""Variable catalog: everything the app knows about a variable lives here.

Each variable is fetched from NASA POWER (monthly, global, no login, derived from the
NASA GMAO MERRA-2 reanalysis and NASA CERES/GEWEX satellite radiation products).
`giovanni` names the matching Giovanni variable so the team can cross-check any number
in Giovanni (the challenge's main resource) with the same region box.

family: variables in the same family are near-identical by construction (e.g. mean,
max and min temperature). Their correlations are trivial, so the auto-detective skips
those pairs instead of "discovering" them.
"""

VARIABLES = {
    # ---------------- heat
    "T2M": dict(
        label="Air temperature (2 m)", category="Heat", family="temp",
        definition="Air temperature at 2 m above the ground, the height of a weather station. Monthly mean of day and night.",
        giovanni="MERRA-2 M2TMNXSLV · 2-meter air temperature (T2M)", pct=False),
    "T2M_MAX": dict(
        label="Daily maximum temperature", category="Heat", family="temp",
        definition="Average of each day's highest 2 m air temperature in the month. Tracks hot afternoons.",
        giovanni="MERRA-2 M2SDNXSLV · T2MMAX (daily statistics)", pct=False),
    "T2M_MIN": dict(
        label="Daily minimum temperature", category="Heat", family="temp",
        definition="Average of each day's lowest 2 m air temperature in the month. Tracks warm nights.",
        giovanni="MERRA-2 M2SDNXSLV · T2MMIN (daily statistics)", pct=False),
    "TS": dict(
        label="Earth skin temperature", category="Heat", family="temp",
        definition="Temperature of the very surface of the land or sea, as seen by a satellite's infrared eye.",
        giovanni="MERRA-2 M2TMNXSLV · surface skin temperature (TS)", pct=False),
    # ---------------- water
    "PRECTOTCORR": dict(
        label="Precipitation", category="Water", family="precip",
        definition="Rain plus snow (as water), bias-corrected with gauge data. Average per day in the month.",
        giovanni="MERRA-2 M2TMNXFLX · PRECTOTCORR  (or GPM IMERG monthly for 2000+)", pct=True),
    "RH2M": dict(
        label="Relative humidity (2 m)", category="Water", family="humid",
        definition="How close the air is to saturation, in percent. 100 % means fog or dew.",
        giovanni="MERRA-2 (derived from QV2M, T2M, PS)", pct=False),
    "QV2M": dict(
        label="Specific humidity (2 m)", category="Water", family="humid2",
        definition="How much water vapour the air actually holds, grams of water per kilogram of air.",
        giovanni="MERRA-2 M2TMNXSLV · QV2M", pct=True),
    "GWETROOT": dict(
        label="Root-zone soil wetness", category="Water", family="soil",
        definition="How wet the soil is where plant roots live (top ~1 m). 0 = bone dry, 1 = saturated.",
        giovanni="MERRA-2 M2TMNXLND · GWETROOT  (or GLDAS soil moisture via Hydrology Time Series)", pct=True),
    "GWETTOP": dict(
        label="Surface soil wetness", category="Water", family="soil",
        definition="Wetness of the top ~5 cm of soil. Reacts within days to rain and heat.",
        giovanni="MERRA-2 M2TMNXLND · GWETTOP", pct=True),
    "EVPTRNS": dict(
        label="Evapotranspiration", category="Water", family="evap",
        definition="Water returned to the air by evaporation from soil and water plus transpiration by plants.",
        giovanni="MERRA-2 M2TMNXLND · EVLAND  (or GLDAS Evap_tavg)", pct=True),
    # ---------------- sun, clouds, air
    "ALLSKY_SFC_SW_DWN": dict(
        label="Sunlight at the surface (all sky)", category="Sun & air", family="sw",
        definition="Solar energy reaching the ground including the effect of clouds and haze.",
        giovanni="MERRA-2 M2TMNXRAD · SWGDN  (CERES for 2000+)", pct=True),
    "CLRSKY_SFC_SW_DWN": dict(
        label="Sunlight at the surface (clear sky)", category="Sun & air", family="sw",
        definition="Solar energy that would reach the ground if there were no clouds. Changes here come from the air itself, mostly aerosols.",
        giovanni="MERRA-2 M2TMNXRAD · SWGDNCLR", pct=True),
    "CLOUD_AMT": dict(
        label="Cloud amount", category="Sun & air", family="cloud",
        definition="Share of the sky covered by clouds, in percent.",
        giovanni="MERRA-2 M2TMNXRAD · CLDTOT  (or MODIS cloud fraction)", pct=True),
    "AOD_55": dict(
        label="Aerosol optical depth (550 nm)", category="Sun & air", family="aod",
        definition="How hazy the air is: how much sunlight is blocked by particles (soot, sulfate, dust, sea salt). 0.1 is clean, 0.5 is heavy haze. No unit.",
        giovanni="MERRA-2 M2TMNXAER · TOTEXTTAU  (or MODIS AOD)", pct=True),
    # ---------------- wind & pressure
    "WS10M": dict(
        label="Wind speed (10 m)", category="Wind", family="wind",
        definition="Average wind speed 10 m above the ground. Matters for wind energy and for dust and pollution transport.",
        giovanni="MERRA-2 M2TMNXSLV · U10M/V10M", pct=True),
    "PS": dict(
        label="Surface pressure", category="Wind", family="ps",
        definition="Weight of the air column above the ground. Mostly follows altitude; its changes reflect shifting weather patterns.",
        giovanni="MERRA-2 M2TMNXSLV · PS", pct=False),
}

CATEGORY_ORDER = ["Heat", "Water", "Sun & air", "Wind", "Giovanni"]

SOURCE_NOTE = (
    "Data: NASA POWER monthly (power.larc.nasa.gov), built from the NASA GMAO MERRA-2 reanalysis "
    "(meteorology, 0.5° × 0.625°, 1981–present) and NASA CERES / GEWEX SRB satellite radiation (1984–present). "
    "Region values are cos(latitude)-weighted means of POWER grid points inside the region. "
    "Giovanni area-averaged CSVs you add are shown as their own variables."
)


def family_of(var_id: str) -> str:
    v = VARIABLES.get(var_id)
    if v:
        return v["family"]
    return "giovanni:" + var_id  # every Giovanni file is its own family


def meta(var_id: str, units: dict | None = None) -> dict:
    v = dict(VARIABLES.get(var_id, {}))
    v["id"] = var_id
    if units and var_id in units:
        v["unit"] = units[var_id]
    return v
