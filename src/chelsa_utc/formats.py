"""The two file formats and the lossless conversion between them.

Native CHELSA format (as published, and as on the public bucket):

    <root>/daily/<var>/<year>/CHELSA_<var>_<dd>_<mm>_<yyyy>_V.2.1.tif
    Cloud Optimized GeoTIFF, uint16, nodata 65535, DEFLATE (predictor 2), 512 x 512 tiles,
    overviews 2..128, EPSG:4326, 43200 x 20880 pixels, value = raw * scale + offset
    (scale 0.1; hurs 0.001; sfcWind 0.01; offset 0, except V2.1 tas of 1941-1950: -273.15),
    CHELSA metadata tags (variable, units, citation, datetime, ...)

Zarr format: ``store`` (one shard per (variable, day), the same uint16 integers and scale).

Functions:
    read_tif / write_tif   one native file
    read_raw               one (variable, day) from either format
    open_chelsa            decoded values (xarray) from either format, for a region and days
    tif_to_zarr            native files -> zarr store
    zarr_to_tif            zarr store -> native files (e.g. CHELSA-UTC in the original file format)
"""
from __future__ import annotations

import concurrent.futures as _cf
import multiprocessing as _mp
from pathlib import Path
from typing import Iterable, NamedTuple

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.transform import Affine
from rasterio.windows import Window

from chelsa_utc import store
from chelsa_utc.saga_ops import NX, NY

CRS_WKT = ('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],'
           'AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,'
           'AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]')
NAME = "CHELSA_{var}_{day:%d_%m_%Y}_V.2.1.tif"
# CHELSA V2.1 metadata per variable: (cf_standard_name, variable_long_name, variable_unit), as in the published tags
TAGS = {
    "tasmin":  ("air_temperature", "Daily Minimum Near-Surface Air Temperature", "K"),
    "tasmax":  ("air_temperature", "Daily Maximum Near-Surface Air Temperature", "K"),
    "tas":     ("air_temperature", "Daily Mean Near-Surface Air Temperature", "K"),
    "prec":    ("precipitation_flux", "Precipitation", "kg m-2 day-1"),
    "rsds":    ("surface_downwelling_shortwave_flux_in_air", "Surface Downwelling Shortwave Flux in Air", "W m-2"),
    "ps":      ("surface_air_pressure", "Surface Air Pressure", "hPa"),
    "sfcWind": ("wind_speed", "Near-Surface Wind Speed", "m s-1"),
    "hurs":    ("relative_humidity", "Near-Surface Relative Humidity", "percent"),
}
COMMON_TAGS = {
    "activity_id": "observation", "experiment_id": "observation", "forcing": "ERA5", "frequency": "day",
    "institution_id": "WSL", "product": "model-output", "project": "CHELSA", "source_id": "CHELSA_V2.1",
    "version": "2.1",
    "citation1": "Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., "
                 "Linder, P., Kessler, M. (2017). Climatologies at high resolution for the Earth land surface areas. "
                 "Scientific Data. 4 170122. https://doi.org/10.1038/sdata.2017.122",
    "data_cite": "Karger, D. N. (2025). CHELSA-daily climate data at high resolution. EnviDat. "
                 "https://www.doi.org/10.16904/envidat.687",
}


class Raw(NamedTuple):
    raw: np.ndarray          # uint16, 65535 = missing
    scale: float
    offset: float
    tags: dict


# ── one native file ──────────────────────────────────────────────────────────

def read_tif(path: str | Path, rows: slice | None = None) -> Raw:
    """Integers, scale, offset and tags of a native CHELSA GeoTIFF (``rows``: a band of tif rows)."""
    with rasterio.open(path) as s:
        win = None if rows is None else Window(0, rows.start, s.width, rows.stop - rows.start)
        return Raw(s.read(1, window=win), float(s.scales[0]), float(s.offsets[0]) + 0.0, s.tags())


def write_tif(path: str | Path, raw: np.ndarray, scale: float, offset: float = 0.0, tags: dict | None = None) -> Path:
    """Write one global day as a native CHELSA GeoTIFF (COG, the layout of the published files).
    Overviews use nearest-neighbour resampling, so they hold valid integers."""
    path = Path(path)
    assert raw.shape == (NY, NX) and raw.dtype == np.uint16, (raw.shape, raw.dtype)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    profile = dict(driver="COG", dtype="uint16", nodata=65535, width=NX, height=NY, count=1, crs=CRS_WKT,
                   transform=Affine(*store.TRANSFORM), compress="DEFLATE", predictor=2, blocksize=512,
                   overview_resampling="NEAREST", BIGTIFF="IF_SAFER")
    with rasterio.open(tmp, "w", **profile) as d:
        d.write(raw, 1)
        d.scales = (scale,)
        d.offsets = (offset,)
        d.update_tags(**{k: str(v) for k, v in (tags or {}).items()}, AREA_OR_POINT="Area")
    tmp.rename(path)
    return path


def tags_for(var: str, day, extra: dict | None = None) -> dict:
    """CHELSA V2.1 tags of one (variable, day), plus ``extra`` (e.g. the day window of CHELSA-UTC)."""
    std, long_name, unit = TAGS[var]
    return {**COMMON_TAGS, "cf_standard_name": std, "variable": var, "variable_long_name": long_name,
            "variable_unit": unit, "datetime": f"{pd.Timestamp(day):%Y-%m-%d}", **(extra or {})}


def tif_path(root: str | Path, var: str, day, name: str = NAME) -> Path:
    day = pd.Timestamp(day)
    return Path(root) / "daily" / var / f"{day.year}" / name.format(var=var, day=day)


# ── either format ─────────────────────────────────────────────────────────────

def read_raw(source: str | Path, var: str, day) -> tuple[np.ndarray, float, float]:
    """(uint16, scale, offset) of one global (variable, day) from a zarr store or a native tree."""
    if store.is_store(source):
        g = store.open_store(source)
        t = store.index(g).get(pd.Timestamp(day))
        if t is None:
            raise KeyError(f"{day} not on the time axis of {source}")
        return store.raw_day(g, var, t)
    r = read_tif(tif_path(source, var, day))
    return r.raw, r.scale, r.offset


def _bbox_rows_cols(bbox) -> tuple[slice, slice]:
    """Tif rows / columns of a (west, south, east, north) box on the 30" grid."""
    a, _, c, _, e, f = store.TRANSFORM
    w, s, ea, n = bbox
    c0 = max(0, int(np.floor((w - c) / a))); c1 = min(NX, int(np.ceil((ea - c) / a)))
    r0 = max(0, int(np.floor((n - f) / e))); r1 = min(NY, int(np.ceil((s - f) / e)))
    return slice(r0, r1), slice(c0, c1)


def open_chelsa(source: str | Path, variables: Iterable[str], days=None, bbox=None) -> xr.Dataset:
    """Decoded float32 values (NaN = missing) of a zarr store or a native tree, for the given days
    (default: all days of the store; required for a native tree) and a (west, south, east, north) box.
    The result is in memory: one decoded global day is ~3.6 GB per variable."""
    rows, cols = _bbox_rows_cols(bbox) if bbox else (slice(0, NY), slice(0, NX))
    a, _, c, _, e, f = store.TRANSFORM
    y = f + (np.arange(rows.start, rows.stop) + 0.5) * e
    x = c + (np.arange(cols.start, cols.stop) + 0.5) * a
    if store.is_store(source):
        g = store.open_store(source)
        idx = store.index(g)
        days = pd.DatetimeIndex(days) if days is not None else store.dates(g)
        ts = [idx[d] for d in days]
    elif days is None:
        raise ValueError("days are required for a native tree")
    else:
        days = pd.DatetimeIndex(days)
    out = {}
    for v in variables:
        a3 = np.full((len(days), len(y), len(x)), np.nan, dtype=np.float32)
        for i, d in enumerate(days):
            if store.is_store(source):
                raw, sc, off = store.raw_day(g, v, ts[i], rows)
                raw = raw[:, cols]
            else:
                p = tif_path(source, v, d)
                if not p.exists():
                    continue
                with rasterio.open(p) as s:
                    raw = s.read(1, window=Window(cols.start, rows.start, cols.stop - cols.start, rows.stop - rows.start))
                    sc, off = float(s.scales[0]), float(s.offsets[0])
            a3[i] = store.decode(raw, sc, off)
        units = store.VARIABLES[v][1] if v in store.VARIABLES else ""
        out[v] = xr.DataArray(a3, dims=("time", "y", "x"), coords={"time": days, "y": y, "x": x},
                              attrs={"units": units})
    return xr.Dataset(out, attrs={"source": str(source)})


# ── conversions ───────────────────────────────────────────────────────────────

def _import_item(src: str, dst: str, var: str, t: int, day: str) -> str | None:
    p = tif_path(src, var, day)
    if not p.exists():
        return f"missing {p}"
    r = read_tif(p)
    g = store.open_store(dst, "r+")
    a = g[var]
    if abs(r.scale - float(a.attrs["scale_factor"])) > 1e-12:
        return f"{p.name}: scale {r.scale} != store {a.attrs['scale_factor']}"
    if r.offset != 0.0:
        name = a.attrs.get("add_offset_per_time")
        if not name:
            return f"{p.name}: offset {r.offset} but the store has no per-day offset for {var}"
        g[name][t] = r.offset
    a[t] = r.raw
    return None


def tif_to_zarr(src: str | Path, dst: str | Path, variables: Iterable[str], first: str, last: str,
                workers: int = 8, redo: bool = False) -> int:
    """Copy native files (the published layout under ``src``) into a zarr store, integers unchanged.
    The store is created if needed (time axis first..last; temperatures get a per-day offset array,
    as some published tas files carry offset -273.15). Pairs already written are skipped."""
    dst = Path(dst); variables = list(variables)
    days = pd.date_range(first, last, freq="D")
    if not dst.exists():
        store.create(dst, days, variables, {"source": f"native CHELSA GeoTIFFs ({Path(src).name})",
                                            "built_by": "chelsa_utc.formats.tif_to_zarr"},
                     per_day_offset=[v for v in variables if v in ("tas", "tasmin", "tasmax")])
    g = store.open_store(dst)
    idx = store.index(g)
    items = [(v, idx[d], f"{d:%Y-%m-%d}") for d in days if d in idx for v in variables
             if v in g and (redo or not store.written(dst, v, idx[d]))]
    print(f"[tif-to-zarr] {len(items)} (variable, day) pairs to copy -> {dst}")
    with _cf.ProcessPoolExecutor(workers, mp_context=_mp.get_context("spawn")) as ex:
        msgs = [m for m in ex.map(_import_item, [str(src)] * len(items), [str(dst)] * len(items), *zip(*items)) if m] \
            if items else []
    for m in msgs:
        print(f"[tif-to-zarr] {m}")
    print(f"[tif-to-zarr] {len(items) - len(msgs)} copied, {len(msgs)} not copied")
    return 1 if msgs else 0


def _export_item(src: str, dst: str, var: str, t: int, day: str, name: str, extra: dict, overwrite: bool) -> str | None:
    p = tif_path(dst, var, day, name)
    if p.exists() and not overwrite:
        return None
    raw, scale, offset = store.raw_day(store.open_store(src), var, t)
    write_tif(p, raw, scale, offset, tags_for(var, day, extra))
    return str(p)


def zarr_to_tif(src: str | Path, dst: str | Path, variables: Iterable[str] | None = None, first: str | None = None,
                last: str | None = None, name: str = NAME, workers: int = 4, overwrite: bool = False) -> int:
    """Write a zarr store as native files under ``dst`` (``<dst>/daily/<var>/<year>/<name>``), integers,
    scale and offset unchanged. ``name`` is a format string with {var} and {day}; the default is the
    V2.1 file name, so tools that read the published files read the output too. The day window of the
    store is written into the tags (``day_window``), and ``source_id`` says which product it is."""
    g = store.open_store(src)
    variables = [v for v in (variables or g.attrs.get("variables", list(store.VARIABLES))) if v in g]
    dates = store.dates(g)
    sel = (dates >= pd.Timestamp(first or dates[0])) & (dates <= pd.Timestamp(last or dates[-1]))
    window = g.attrs.get("day_window", "")
    extra = {"day_window": window} if window else {}
    if "CHELSA-UTC" in str(g.attrs.get("title", "")):
        extra.update({"source_id": "CHELSA-UTC", "version": str(g.attrs.get("version", "")),
                      "comment": "CHELSA V2.1 algorithm on the 00-24 UTC day (chelsa-utc)"})
    items = [(v, int(t), f"{d:%Y-%m-%d}") for t, d in zip(np.nonzero(sel)[0], dates[sel]) for v in variables
             if store.written(src, v, int(t))]
    print(f"[zarr-to-tif] {len(items)} files -> {dst}")
    n = 0
    with _cf.ProcessPoolExecutor(workers, mp_context=_mp.get_context("spawn")) as ex:
        futs = [ex.submit(_export_item, str(src), str(dst), v, t, d, name, extra, overwrite) for v, t, d in items]
        for f in _cf.as_completed(futs):
            if f.result():
                n += 1
    print(f"[zarr-to-tif] {n} written, {len(items) - n} existed already")
    return 0
