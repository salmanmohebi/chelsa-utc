"""Daily ERA5 from the local hourly mirror (``era5_hourly``), for a chosen day window.

Day windows, in hours after 00 UTC of day d:

    window        forecast fields          analysis fields
    "utc"         valid hours 01..24       00..23        the 00-24 UTC day (CHELSA-UTC)
    "chelsa0606"  valid hours 07..30       00..23        CHELSA V2.1 inputs (06 UTC d -> 06 UTC d+1)

A forecast value at valid hour t covers (t - 1 h, t]. Forecast fields: mn2t, mx2t, mtpr, ssrd, and
cbh: cbh is instantaneous, but CHELSA V2.1 averages it over the forecast window (its cbase.nc =
NaN-mean over valid hours 07..30, checked on two days). Analysis fields: t2m, sp, tcc, blh, u10, v10.

Daily variables:

    name       hourly   daily value                   unit
    t2m_min    mn2t     min                           K
    t2m_max    mx2t     max                           K
    tp_mean    mtpr     sum x 3600                    mm/day
    ssrd       ssrd     sum / 86400                   W m-2
    cbh        cbh      NaN-mean (NaN = no cloud)     m
    t2m_mean   t2m      mean                          K
    sp         sp       mean / 100                    hPa
    tcc        tcc      mean                          1
    blh        blh      mean                          m
    u10        u10      mean                          m s-1
    v10        v10      mean                          m s-1

Grid: latitude 90 -> -90 (721), longitude -180 -> 179.75 (1440); the hourly mirror keeps ARCO's
0 -> 359.75, so every field is rolled by 720 columns.

Cache: the daily store of a window (config ``paths.era5_daily.<window>``) holds ``daily(d, window)``
for every day already built. ``daily`` reads a day from the store when it is there and computes it
from the hourly mirror otherwise. The store values are bit-identical to the computed ones
(``docs/validation.md``).
"""
from __future__ import annotations

import concurrent.futures as _cf
import multiprocessing as _mp
import warnings
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import zarr
from zarr.codecs import BloscCodec

from chelsa_utc import era5_hourly
from chelsa_utc.config import paths

WINDOWS = {                      # (first, last + 1) hour after 00 UTC of day d
    "utc":        {"forecast": (1, 25), "analysis": (0, 24)},
    "chelsa0606": {"forecast": (7, 31), "analysis": (0, 24)},
}
FORECAST = {"mn2t", "mx2t", "mtpr", "ssrd", "cbh"}
DAILY = {                        # name: (hourly variable, aggregation, factor, unit)
    "t2m_min":  ("mn2t", "min", 1.0, "K"),
    "t2m_max":  ("mx2t", "max", 1.0, "K"),
    "tp_mean":  ("mtpr", "sum", 3600.0, "mm/day"),
    "ssrd":     ("ssrd", "sum", 1 / 86400, "W m-2"),
    "cbh":      ("cbh", "nanmean", 1.0, "m"),
    "t2m_mean": ("t2m", "mean", 1.0, "K"),
    "sp":       ("sp", "mean", 0.01, "hPa"),
    "tcc":      ("tcc", "mean", 1.0, "1"),
    "blh":      ("blh", "mean", 1.0, "m"),
    "u10":      ("u10", "mean", 1.0, "m s-1"),
    "v10":      ("v10", "mean", 1.0, "m s-1"),
}
AGG = {"min": np.min, "max": np.max, "sum": np.sum, "mean": np.mean, "nanmean": np.nanmean}


def _hourly() -> zarr.Group:
    return zarr.open_group(str(era5_hourly.store_path()), mode="r")


def from_hourly(day, window: str = "utc", variables: Iterable[str] | None = None) -> dict[str, np.ndarray]:
    """The daily fields of one day computed from the hourly mirror:
    {name: float32 (721, 1440), lat 90 -> -90, lon -180 -> 179.75}."""
    g = _hourly()
    h0 = int((pd.Timestamp(day) - era5_hourly.T0) / pd.Timedelta(hours=1))
    fields = {}
    for name in variables or DAILY:
        var, how, factor, _ = DAILY[name]
        a, b = WINDOWS[window]["forecast" if var in FORECAST else "analysis"]
        x = g[var][h0 + a:h0 + b].astype(np.float64)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)      # cbh: all-NaN pixel = never cloudy
            v = AGG[how](x, axis=0) * factor
        fields[name] = np.roll(v.astype(np.float32), -720, axis=1)
    return fields


# ── daily store (cache) ───────────────────────────────────────────────────────

def store_path(window: str) -> Path:
    return paths().era5_daily[window]


def _store_index(g: zarr.Group, day: pd.Timestamp) -> int | None:
    """Index of ``day`` on the time axis of a daily store, or None."""
    epoch = pd.Timestamp(g["time"].attrs["units"].split("since")[1].strip())
    k = int((day - epoch).days)
    t = g["time"]
    if 0 <= k < t.shape[0] and int(t[k]) == k:                  # contiguous axis (the usual case)
        return k
    hit = np.nonzero(t[:] == k)[0]
    return int(hit[0]) if hit.size else None


def from_store(day, window: str, variables: Iterable[str] | None = None) -> dict[str, np.ndarray] | None:
    """The daily fields of one day from the window's daily store, or None if the store, the day or
    one of the fields is not there (a field that was never written is all NaN)."""
    p = store_path(window)
    if not p.exists():
        return None
    g = zarr.open_group(str(p), mode="r")
    t = _store_index(g, pd.Timestamp(day))
    if t is None:
        return None
    fields = {}
    for name in variables or DAILY:
        if name not in g:
            return None
        x = np.asarray(g[name][t], dtype=np.float32)
        if name != "cbh" and not np.isfinite(x).any():
            return None
        fields[name] = x
    if "cbh" in fields and not np.isfinite(fields["cbh"]).any():
        return None
    return fields


# The daily ERA5 input files of CHELSA V2.1 (one netCDF per variable and day, 06->06 forecasts and
# 00-23 analyses), if available: <paths.v21_era5_inputs>/<year>/<yyyymmdd>/<stem>.nc.
# Window "v21files" reads them instead of the hourly mirror (the "mirror test" of docs/validation.md).
V21_FILES = {               # name: (file stem, factor to the daily unit)
    "t2m_min": ("tmin", 1.0), "t2m_max": ("tmax", 1.0), "t2m_mean": ("tmean", 1.0),
    "tp_mean": ("ps", 1000.0),          # ps.nc = precipitation in m
    "ssrd": ("ssrd", 1 / 86400), "sp": ("ps0", 0.01), "blh": ("blh", 1.0),
    "cbh": ("cbase", 1.0), "tcc": ("tcc", 1.0), "u10": ("u", 1.0), "v10": ("v", 1.0),
}


def from_v21_files(day, variables: Iterable[str] | None = None) -> dict[str, np.ndarray]:
    """The daily fields of one day from the V2.1 daily input files (same grid and units as ``daily``).
    A file with several time steps (e.g. 24 hourly fields on 1958-05-31) gives their mean, as V2.1 used."""
    import xarray as xr
    day = pd.Timestamp(day)
    base = paths().v21_era5_inputs
    if base is None:
        raise FileNotFoundError("paths.v21_era5_inputs is not set in the config")
    fields = {}
    for name in variables or DAILY:
        stem, factor = V21_FILES[name]
        with xr.open_dataset(base / str(day.year) / f"{day:%Y%m%d}" / f"{stem}.nc") as ds:
            da = ds[list(ds.data_vars)[0]]
            a = da.values.astype(np.float64)
            if a.ndim == 3:
                a = a.mean(axis=0)
            lat = ds.latitude.values; lon = ds.longitude.values
        lon = np.where(lon >= 180, lon - 360, lon)
        fields[name] = a[np.argsort(-lat, kind="stable")][:, np.argsort(lon, kind="stable")] * factor   # float64
    return fields


def daily(day, window: str = "utc", variables: Iterable[str] | None = None,
          source: str = "auto") -> dict[str, np.ndarray]:
    """The daily fields of one day: {name: float32 (721, 1440), lat 90 -> -90, lon -180 -> 179.75}.
    source "auto": the daily store if it has the day, else the hourly mirror; "store" / "hourly": only that.
    window "v21files": the V2.1 daily input files (``from_v21_files``)."""
    if window == "v21files":
        return from_v21_files(day, variables)
    if source in ("auto", "store"):
        f = from_store(day, window, variables)
        if f is not None:
            return f
        if source == "store":
            raise KeyError(f"{day} not in {store_path(window)}")
    return from_hourly(day, window, variables)


def _create(path: Path, days: pd.DatetimeIndex, window: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    g = zarr.open_group(str(path), mode="w", zarr_format=3)
    t = g.create_array("time", shape=(len(days),), dtype="int64", chunks=(len(days),), dimension_names=["time"],
                       attributes={"units": f"days since {days[0]:%Y-%m-%d}", "calendar": "proleptic_gregorian"})
    t[:] = np.arange(len(days))
    lat = g.create_array("latitude", shape=(721,), dtype="float64", chunks=(721,), dimension_names=["latitude"],
                         attributes={"units": "degrees_north"})
    lat[:] = np.linspace(90.0, -90.0, 721)
    lon = g.create_array("longitude", shape=(1440,), dtype="float64", chunks=(1440,), dimension_names=["longitude"],
                         attributes={"units": "degrees_east"})
    lon[:] = -180.0 + 0.25 * np.arange(1440)
    comp = BloscCodec(cname="zstd", clevel=3, shuffle="bitshuffle")
    for name, (var, how, _, unit) in DAILY.items():
        g.create_array(name, shape=(len(days), 721, 1440), dtype="float32", chunks=(1, 721, 1440),
                       fill_value=np.float32("nan"), compressors=comp,
                       dimension_names=["time", "latitude", "longitude"],
                       attributes={"units": unit, "cell_methods": f"time: {how} of hourly {var}"})
    g.attrs.update({"window": window, "hours": WINDOWS[window], "source": "hourly ARCO-ERA5 mirror (era5_hourly)",
                    "built_by": "chelsa_utc.era5_daily"})
    zarr.consolidate_metadata(str(path))


def _write_day(path: str, window: str, t: int, day: str) -> None:
    g = zarr.open_group(path, mode="r+", use_consolidated=False)
    for name, x in from_hourly(day, window).items():
        g[name][t] = x


def build(window: str = "utc", first: str = "1941-01-01", last: str = "1978-12-31", workers: int = 16,
          redo: bool = False) -> Path:
    """Fill the daily store of one window for the days first..last (one process per day).
    The store is created if it does not exist (time axis first..last); days already in it are skipped."""
    days = pd.date_range(first, last, freq="D")
    path = store_path(window)
    if not path.exists():
        _create(path, days, window)
    g = zarr.open_group(str(path), mode="r")
    todo = []
    for d in days:
        t = _store_index(g, d)
        if t is None:
            print(f"[era5-daily] {d:%Y-%m-%d} is not on the time axis of {path}, skipped")
        elif redo or from_store(d, window, ["t2m_mean"]) is None:
            todo.append((t, f"{d:%Y-%m-%d}"))
    print(f"[era5-daily] {window}: {len(days) - len(todo)} of {len(days)} days already in the store, {len(todo)} to compute")
    if todo:
        with _cf.ProcessPoolExecutor(workers, mp_context=_mp.get_context("spawn")) as ex:
            list(ex.map(_write_day, [str(path)] * len(todo), [window] * len(todo), *zip(*todo)))
    return path
