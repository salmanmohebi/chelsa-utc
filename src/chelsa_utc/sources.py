"""Published CHELSA V2.1 files: local paths, readers and download from the public buckets.

CHELSA-UTC reuses these published files (they are already on the 00-23 UTC day or do not depend on
the day):

    daily helper maps   tz  (0.25 deg lapse rate, K m-1)          -> tasmin, tasmax, tas
                        we  (0.025 deg wind effect)                -> prec
                        clt (0.025 deg total cloud cover)          -> rsds
    unchanged outputs   hurs, sfcWind (30")                        -> copied into the store
    static grids        dem_latlong, expocor, patch, orography (SAGA .sdat + .sgrd)
    clear sky           stot: CHELSA_stot_pj_<doy>_V.2.1.tif, doy 1..366 (rsds)

Local layout (config ``paths.chelsa_v21`` and ``paths.stot``) = the layout of the public bucket:

    <chelsa_v21>/daily/<var>/<year>/CHELSA_<var>_<dd>_<mm>_<yyyy>_V.2.1.tif
    <chelsa_v21>/input/static/<name>.sdat (+ .sgrd)
    <stot>/CHELSA_stot_pj_<doy>_V.2.1.tif

``download_*`` fetch only the files that are missing locally (or have another size).
"""
from __future__ import annotations

import concurrent.futures as _cf
import time as _time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

from chelsa_utc.config import paths
from chelsa_utc.saga_ops import CS, NX, NY, X30, XMIN, NodeGrid, rows_y

PUBLIC = "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global"
STOT_URL = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/daily_normals/rsdscs/CHELSA_stot_pj_{doy}_V.2.1.tif"
STATIC = {"dem": "dem_latlong", "expocor": "expocor", "patch": "patch", "orog": "orography"}
HELPER = ("tz", "we", "clt")
UNCHANGED = ("hurs", "sfcWind")


# ── paths ─────────────────────────────────────────────────────────────────────

def daily_name(var: str, day) -> str:
    return f"CHELSA_{var}_{pd.Timestamp(day):%d_%m_%Y}_V.2.1.tif"


def daily_tif(var: str, day) -> Path:
    day = pd.Timestamp(day)
    return paths().chelsa_v21 / "daily" / var / f"{day.year}" / daily_name(var, day)


def static_path(name: str, ext: str = ".sdat") -> Path:
    return paths().chelsa_v21 / "input" / "static" / f"{STATIC.get(name, name)}{ext}"


def stot_path(day) -> Path:
    return paths().stot / f"CHELSA_stot_pj_{pd.Timestamp(day).dayofyear}_V.2.1.tif"


# ── readers (cached per process) ─────────────────────────────────────────────

_CACHE: dict = {}


def static(name: str):
    """Static SAGA grids, memory-mapped once per process."""
    if name not in _CACHE:
        if name == "dem":
            _CACHE[name] = np.memmap(static_path("dem"), dtype="<i2", mode="r", shape=(NY, NX))
        elif name == "expocor":
            _CACHE[name] = np.memmap(static_path("expocor"), dtype="<f4", mode="r", shape=(NY, NX))
        elif name == "patch":
            _CACHE[name] = np.memmap(static_path("patch"), dtype="<f4", mode="r", shape=(7201, 14401))
        elif name == "orog":   # (721, 1440) lat ascending, lon 0..359.75 -> -180..179.75
            a = np.fromfile(static_path("orog"), dtype="<f4").reshape(721, 1440).astype(np.float64)
            _CACHE[name] = np.roll(a, -720, axis=1)
        else:
            raise KeyError(name)
    return _CACHE[name]


def sdat_rows(name: str, nodata: float, r0: int, r1: int) -> np.ndarray:
    """30" SAGA grid rows r0..r1-1 (tif order, north -> south) as float64, NoData -> NaN."""
    a = np.array(static(name)[NY - r1:NY - r0], dtype=np.float64)[::-1]
    a[a == nodata] = np.nan
    return a


def patch_rows(r0: int, r1: int) -> np.ndarray:
    """patch.sdat (0.025 deg) at the 30" pixels of rows r0..r1-1, nearest node."""
    pxmin, pymin, pcs = -180.0123559950, -90.0124999950, 0.0249999900
    pix = np.floor((X30 - pxmin) / pcs + 0.5).astype(np.int64)
    piy = np.floor((rows_y(r0, r1) - pymin) / pcs + 0.5).astype(np.int64)
    a = np.array(static("patch")[piy.min():piy.max() + 1], dtype=np.float64)
    return a[piy - piy.min()][:, pix]


def tif_grid(path: Path, ylo: float, yhi: float) -> NodeGrid:
    """Full-width window of a north-up tif (the 0.025 deg we / clt) covering latitudes ylo..yhi
    with a 3-node margin, as a NodeGrid (rows ascending)."""
    with rasterio.open(path) as s:
        t = s.transform; W, H = s.width, s.height
        cs = t.a; x0 = t.c + cs / 2; ytop = t.f + t.e / 2; y0 = ytop + (H - 1) * t.e
        rt0 = max(0, int(np.floor((ytop - yhi) / (-t.e))) - 3)
        rt1 = min(H - 1, int(np.ceil((ytop - ylo) / (-t.e))) + 3)
        raw = s.read(1, window=Window(0, rt0, W, rt1 - rt0 + 1))
        a = raw.astype(np.float64) * s.scales[0] + s.offsets[0]
        if s.nodata is not None:
            a[raw == s.nodata] = np.nan
    return NodeGrid(a[::-1], 0, (H - 1) - rt1, W, H, x0, y0, cs)


def tz(day) -> np.ndarray:
    """Published daily lapse rate (0.25 deg), lat ascending, lon -180..179.75."""
    with rasterio.open(daily_tif("tz", day)) as s:
        a = s.read(1).astype(np.float64)
        if s.nodata is not None:
            a[a == s.nodata] = np.nan
    return a[::-1]


def stot_rows(day, r0: int, r1: int) -> np.ndarray:
    """Clear-sky radiation of the day of year, tif rows r0..r1-1, NoData -> NaN."""
    with rasterio.open(stot_path(day)) as s:
        st = s.read(1, window=Window(0, r0, NX, r1 - r0)).astype(np.float64)
        if s.nodata is not None:
            st[st == s.nodata] = np.nan
    return st


# ── download ─────────────────────────────────────────────────────────────────

def _fetch(url: str, dest: Path, tries: int = 6) -> str:
    """Download url to dest unless dest exists with the same size. Returns a status line."""
    import requests
    err: Exception | None = None
    for k in range(tries):
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                if r.status_code == 404:
                    return f"missing {url}"
                r.raise_for_status()
                size = int(r.headers.get("Content-Length", 0))
                if dest.exists() and size and dest.stat().st_size == size:
                    return f"cached  {dest.name}"
                dest.parent.mkdir(parents=True, exist_ok=True)
                tmp = dest.with_name(dest.name + ".part")
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 24):
                        fh.write(chunk)
                if size and tmp.stat().st_size != size:
                    raise IOError(f"size {tmp.stat().st_size} != {size}")
                tmp.rename(dest)
            return f"ok      {dest.name} ({dest.stat().st_size / 1e6:.0f} MB)"
        except Exception as e:                                # retry with back-off
            err = e
            _time.sleep(min(120, 5 * 2 ** k))
    return f"FAIL    {url}: {type(err).__name__}: {err}"


def _run(jobs: list[tuple[str, Path]], threads: int, skip_existing: bool) -> int:
    if skip_existing:
        jobs = [(u, d) for u, d in jobs if not d.exists()]
    with _cf.ThreadPoolExecutor(threads) as ex:
        msgs = list(ex.map(lambda j: _fetch(*j), jobs))
    for m in msgs:
        print(m, flush=True)
    fails = [m for m in msgs if m.startswith("FAIL")]
    print(f"[download] {len(jobs) - len(fails)} of {len(jobs)} ok")
    return 1 if fails else 0


def download_static(threads: int = 4, skip_existing: bool = True) -> int:
    """The four static grids (.sdat + .sgrd, ~5.8 GB)."""
    jobs = [(f"{PUBLIC}/input/static/{n}{ext}", static_path(n, ext)) for n in STATIC.values() for ext in (".sdat", ".sgrd")]
    return _run(jobs, threads, skip_existing)


def download_stot(threads: int = 8, skip_existing: bool = True) -> int:
    """The 366 clear-sky files (~1.1 GB each, ~393 GB)."""
    jobs = [(STOT_URL.format(doy=d), paths().stot / f"CHELSA_stot_pj_{d}_V.2.1.tif") for d in range(1, 367)]
    return _run(jobs, threads, skip_existing)


def download_daily(variables: Iterable[str], days: Iterable, threads: int = 8, skip_existing: bool = True) -> int:
    """Published daily files of the given variables and days (tz ~4 MB, we / clt ~100 MB,
    hurs / sfcWind ~300-450 MB per file)."""
    jobs = []
    for d in pd.DatetimeIndex(list(days)):
        for v in variables:
            jobs.append((f"{PUBLIC}/daily/{v}/{d.year}/{daily_name(v, d)}", daily_tif(v, d)))
    return _run(jobs, threads, skip_existing)


def inputs_present(day, variables=HELPER + UNCHANGED) -> bool:
    """True if the published daily files that CHELSA-UTC needs exist locally for this day."""
    return all(daily_tif(v, day).exists() for v in variables)


def check_grid_30s(path: Path) -> None:
    """Assert that a 30" tif is on the CHELSA grid."""
    with rasterio.open(path) as s:
        assert (s.width, s.height) == (NX, NY), (path, s.width, s.height)
        assert abs(s.transform.c - (XMIN - CS / 2)) < 1e-6 and abs(s.transform.a - CS) < 1e-9, (path, s.transform)
