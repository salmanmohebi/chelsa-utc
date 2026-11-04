"""The zarr layout of the 30" CHELSA stores (CHELSA-UTC, rebuilt V2.1, and imported published files).

    <store>.zarr/                zarr v3
        time  (time,)            int64, "days since <first day>"; may have gaps (missing days)
        y     (20880,)           pixel-centre latitudes, 84 N -> 90 S (from the tif transform)
        x     (43200,)           pixel-centre longitudes, 180 W -> 180 E
        <var> (time, y, x)       uint16, value = raw * scale_factor (+ offset), 65535 = missing;
                                 inner chunk 720 x 720, one shard (= one file) per (day, variable):
                                 <store>.zarr/<var>/c/<t>/0/0
        <var>_add_offset (time,) only in stores imported from published tifs with a non-zero offset
                                 on some days (V2.1 tas of 1941-1950 stores degrees C with offset
                                 -273.15); the array name is in the attribute `add_offset_per_time`

This is the layout of the files that the published V2.1 GeoTIFFs give one to one (same integers,
same scale), so the conversion between the two formats is lossless (``formats``).
A (variable, day) counts as written when its shard file exists.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from zarr.codecs import BloscCodec

from chelsa_utc.saga_ops import NX, NY

FILL = 65535
CHUNK = 720             # divides 20880 (29x) and 43200 (60x)
TRANSFORM = (0.0083333333, 0.0, -180.00013888885002, 0.0, -0.0083333333, 83.99986041515001)   # of the published tifs
VARIABLES = {           # name: (scale_factor, units, long_name)
    "tasmin":  (0.1, "K", "daily minimum 2 m air temperature"),
    "tasmax":  (0.1, "K", "daily maximum 2 m air temperature"),
    "tas":     (0.1, "K", "daily mean 2 m air temperature"),
    "prec":    (0.1, "mm day-1", "precipitation"),
    "rsds":    (0.1, "W m-2", "surface downwelling shortwave flux"),
    "ps":      (0.1, "hPa", "surface air pressure"),
    "hurs":    (0.001, "1", "near-surface relative humidity (fraction)"),
    "sfcWind": (0.01, "m s-1", "near-surface wind speed"),
}


def is_store(path: str | Path) -> bool:
    return (Path(path) / "zarr.json").exists()


def create(path: Path, days: pd.DatetimeIndex, variables, attrs: dict, per_day_offset=()) -> zarr.Group:
    """A new, empty store with the time axis ``days`` (may have gaps)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    days = pd.DatetimeIndex(days)
    g = zarr.open_group(str(path), mode="w-", zarr_format=3)
    t = g.create_array("time", shape=(len(days),), dtype="int64", chunks=(len(days),), dimension_names=["time"],
                       attributes={"units": f"days since {days[0]:%Y-%m-%d}", "calendar": "proleptic_gregorian",
                                   "standard_name": "time"})
    t[:] = (days - days[0]).days.to_numpy()
    a, _, c, _, e, f = TRANSFORM
    y = g.create_array("y", shape=(NY,), dtype="float64", chunks=(NY,), dimension_names=["y"],
                       attributes={"units": "degrees_north", "standard_name": "latitude"})
    y[:] = f + (np.arange(NY, dtype="float64") + 0.5) * e
    x = g.create_array("x", shape=(NX,), dtype="float64", chunks=(NX,), dimension_names=["x"],
                       attributes={"units": "degrees_east", "standard_name": "longitude"})
    x[:] = c + (np.arange(NX, dtype="float64") + 0.5) * a
    comp = BloscCodec(cname="zstd", clevel=3, shuffle="shuffle")
    for v in variables:
        scale, units, long_name = VARIABLES[v]
        va = {"units": units, "long_name": long_name, "scale_factor": scale, "add_offset": 0.0, "_FillValue": FILL}
        if v in per_day_offset:
            off = g.create_array(f"{v}_add_offset", shape=(len(days),), dtype="float64", chunks=(len(days),),
                                 fill_value=0.0, dimension_names=["time"],
                                 attributes={"long_name": f"per-day add_offset of the source {v} tifs", "units": units})
            off[:] = 0.0
            va["add_offset_per_time"] = f"{v}_add_offset"
        g.create_array(v, shape=(len(days), NY, NX), dtype="uint16", chunks=(1, CHUNK, CHUNK), shards=(1, NY, NX),
                       fill_value=FILL, compressors=comp, dimension_names=["time", "y", "x"], attributes=va)
    g.attrs.update({"grid": "30arcsec", "transform": list(TRANSFORM), "variables": list(variables), **attrs})
    zarr.consolidate_metadata(str(path))
    return zarr.open_group(str(path), mode="r+", use_consolidated=False)


def open_store(path: str | Path, mode: str = "r") -> zarr.Group:
    return zarr.open_group(str(path), mode=mode, use_consolidated=False)


def dates(g: zarr.Group) -> pd.DatetimeIndex:
    epoch = pd.Timestamp(g["time"].attrs["units"].split("since")[1].strip())
    return epoch + pd.to_timedelta(g["time"][:], unit="D")


def index(g: zarr.Group) -> dict[pd.Timestamp, int]:
    return {d: i for i, d in enumerate(dates(g))}


def written(path: str | Path, var: str, t: int) -> bool:
    """True if the shard file of (var, day index t) exists."""
    return (Path(path) / var / "c" / str(t) / "0" / "0").exists()


def raw_day(g: zarr.Group, var: str, t: int, rows: slice = slice(None)) -> tuple[np.ndarray, float, float]:
    """(uint16 raw, scale, offset) of one day; value = raw * scale + offset, 65535 = missing."""
    a = g[var]
    off = float(a.attrs.get("add_offset", 0.0))
    name = a.attrs.get("add_offset_per_time")
    if name:
        off = float(g[name][t])
    return np.asarray(a[t, rows]), float(a.attrs["scale_factor"]), off


def decode(raw: np.ndarray, scale: float, offset: float = 0.0) -> np.ndarray:
    """float32(raw) * float32(scale) + float32(offset), 65535 -> NaN (the float32 of the source files)."""
    v = raw.astype(np.float32) * np.float32(scale) + np.float32(offset)
    v[raw == FILL] = np.nan
    return v


def relink(src: str | Path, dst: str | Path, days, note: str = "") -> zarr.Group:
    """A new store ``dst`` with the time axis ``days`` and the arrays, chunking and attributes of ``src``.
    Every (variable, day) already written in ``src`` is hard-linked into ``dst`` (no copy, no extra disk
    space; both stores must be on the same file system); the other days stay empty for ``build``.
    zarr writes a chunk by writing a new file and renaming it, so a later write into ``dst`` replaces
    the link and never changes ``src``."""
    import os
    src, dst = Path(src), Path(dst)
    if dst.exists():
        raise FileExistsError(dst)
    g = open_store(src)
    old = index(g)
    days = pd.DatetimeIndex(sorted(set(pd.DatetimeIndex(days))))
    cube = [k for k in g.array_keys() if g[k].ndim == 3]
    per_time = [k for k in g.array_keys() if g[k].ndim == 1 and k not in ("time", "y", "x")]
    h = zarr.open_group(str(dst), mode="w-", zarr_format=3)
    t = h.create_array("time", shape=(len(days),), dtype="int64", chunks=(len(days),), dimension_names=["time"],
                       attributes=dict(g["time"].attrs) | {"units": f"days since {days[0]:%Y-%m-%d}"})
    t[:] = (days - days[0]).days.to_numpy()
    for c in ("y", "x"):
        a = h.create_array(c, shape=g[c].shape, dtype=g[c].dtype, chunks=g[c].shape, dimension_names=[c],
                           attributes=dict(g[c].attrs))
        a[:] = g[c][:]
    for k in per_time:                                    # e.g. <var>_add_offset(time)
        a = g[k]; vals = a[:]
        b = h.create_array(k, shape=(len(days),), dtype=a.dtype, chunks=(len(days),), fill_value=a.fill_value,
                           dimension_names=["time"], attributes=dict(a.attrs))
        b[:] = np.array([vals[old[d]] if d in old else a.fill_value for d in days], dtype=a.dtype)
    for k in cube:
        a = g[k]
        h.create_array(k, shape=(len(days),) + a.shape[1:], dtype=a.dtype, chunks=a.chunks, shards=a.shards,
                       fill_value=a.fill_value, compressors=a.compressors, dimension_names=a.metadata.dimension_names,
                       attributes=dict(a.attrs))
    h.attrs.update(dict(g.attrs) | {"relinked_from": src.name, "relink_note": note})
    zarr.consolidate_metadata(str(dst))
    n = 0
    for k in cube:
        for i, d in enumerate(days):
            if d in old and written(src, k, old[d]):
                p = dst / k / "c" / str(i) / "0" / "0"
                p.parent.mkdir(parents=True, exist_ok=True)
                os.link(src / k / "c" / str(old[d]) / "0" / "0", p)
                n += 1
    print(f"[relink] {dst}: {len(days)} days ({len(days) - len(set(days) & set(old))} new), {n} shards linked")
    return open_store(dst, "r+")
