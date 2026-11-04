"""Hourly ERA5 single-level fields, 1940-1979, mirrored from ARCO-ERA5 (lossless, bit-identical).

    <era5_hourly>                    zarr v3, one array per variable (config ``paths.era5_hourly``)
    <era5_hourly stem>.manifest/     <shard>.parquet: per (variable, hour) QC row
    <era5_hourly stem>.check.parquet result of `check`

Why: the local source of every daily ERA5 window (00-24 UTC, CHELSA 06->06, station windows),
so that no window needs a new download. Variables = the 11 single-level ERA5 inputs of CHELSA
V2.1 (its daily input files, one netCDF per variable and day), which are also the 11 variables
of the daily stores of ``era5_daily``:

    store   ARCO name                                               V2.1 file     units (ARCO, unchanged)
    mn2t    minimum_2m_temperature_since_previous_post_processing   tmin.nc       K
    mx2t    maximum_2m_temperature_since_previous_post_processing   tmax.nc       K
    t2m     2m_temperature                                          tmean.nc      K
    mtpr    mean_total_precipitation_rate                           ps.nc         kg m-2 s-1
    sp      surface_pressure                                        ps0.nc        Pa
    ssrd    surface_solar_radiation_downwards                       ssrd.nc       J m-2 (per hour)
    blh     boundary_layer_height                                   blh.nc        m
    cbh     cloud_base_height                                       cbase.nc      m (NaN = no cloud)
    tcc     total_cloud_cover                                       tcc.nc        1
    u10     10m_u_component_of_wind                                 u.nc          m s-1
    v10     10m_v_component_of_wind                                 v.nc          m s-1

The pressure-level inputs (taz, geoz, hursz, uz, vz, ccz; 37 levels) are not mirrored: ~520 MB per
hour in ARCO (~180 TB for 1940-1979), and ARCO has no relative humidity on levels.

Source: https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3
(zarr v2, one chunk = one global hourly field, float32 blosc-lz4), read over plain HTTPS. ARCO's
root attributes (``last_updated``, ...) and the download date are kept in the store attributes.

Layout: time = hours since 1940-01-01 00 UTC (350 640 hours, to 1979-12-31 23 UTC), latitude
90 -> -90 (721), longitude 0 -> 359.75 (1440): the ARCO grid, not rolled. Hour convention (ARCO):
analyses (t2m, sp, blh, cbh, tcc, u10, v10) are instantaneous at hour t; forecast fields (mn2t,
mx2t, mtpr, ssrd) at valid hour t cover (t - 1 h, t]. float32, NaN = missing (also a missing ARCO
chunk), blosc zstd (lossless). Inner chunk = 1 day (24 h); shard (= one file) = SHARD_H = 744 h
(31 days, not calendar months: zarr shards are regular), so 472 shards x 11 variables = 5 192 files.
Shard s holds hours s*744 .. s*744+743, i.e. days s*31 .. s*31+30 after 1940-01-01.

Steps (restartable; a shard is committed by its manifest file and skipped when run again):
    init       create the store (coordinates, empty arrays, ARCO metadata)
    download   fetch and write the given shards, one process per (shard, variable)
    check      completeness, missing chunks / NaN summary from the manifests, and a bit-for-bit
               comparison with ARCO on random hours
"""
from __future__ import annotations

import concurrent.futures as _cf
import json
import multiprocessing as _mp
import threading
import time as _time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import zarr
from zarr.codecs import BloscCodec

from chelsa_utc.config import paths

BASE = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
STORE = "era5-arco-hourly"
VARS = {
    "mn2t": "minimum_2m_temperature_since_previous_post_processing",
    "mx2t": "maximum_2m_temperature_since_previous_post_processing",
    "t2m": "2m_temperature",
    "mtpr": "mean_total_precipitation_rate",
    "sp": "surface_pressure",
    "ssrd": "surface_solar_radiation_downwards",
    "blh": "boundary_layer_height",
    "cbh": "cloud_base_height",
    "tcc": "total_cloud_cover",
    "u10": "10m_u_component_of_wind",
    "v10": "10m_v_component_of_wind",
}
FORECAST = ("mn2t", "mx2t", "mtpr", "ssrd")
T0 = pd.Timestamp("1940-01-01")
T1 = pd.Timestamp("1980-01-01")                       # exclusive
NT = int((T1 - T0) / pd.Timedelta(hours=1))           # 350 640
NLAT, NLON = 721, 1440
CHUNK_H = 24
SHARD_H = 744
NSHARD = -(-NT // SHARD_H)                             # 472
MAX_ROUNDS = 3          # fresh process pools for failed items before a shard counts as failed
EPOCH = pd.Timestamp("1900-01-01")                     # ARCO time axis: hours since EPOCH
ARCO_OFFSET = int((T0 - EPOCH) / pd.Timedelta(hours=1))   # ARCO index of 1940-01-01 00 UTC

_local = threading.local()


# ── ARCO access ──────────────────────────────────────────────────────────────

def _get(url: str, tries: int = 8) -> bytes | None:
    """Bytes of one ARCO object; None if it does not exist (404)."""
    import requests
    if not hasattr(_local, "s"):
        _local.s = requests.Session()
    for k in range(tries):
        try:
            r = _local.s.get(url, timeout=120)
            if r.status_code == 200:
                return r.content
            if r.status_code == 404:
                return None
        except Exception:                                   # network hiccup: back off and retry
            pass
        _time.sleep(min(60, 2 ** k))
    raise RuntimeError(f"failed after {tries} tries: {url}")


def _decode(raw: bytes, n: int) -> np.ndarray:
    from numcodecs import Blosc
    return np.frombuffer(Blosc().decode(raw), dtype="<f4", count=n)


def arco_field(arco_name: str, hour: int) -> np.ndarray | None:
    """One hourly global field (721, 1440) float32 as stored in ARCO; None if the chunk is missing.
    ``hour`` = index in our time axis (hours since 1940-01-01 00 UTC)."""
    raw = _get(f"{BASE}/{arco_name}/{hour + ARCO_OFFSET}.0.0")
    return None if raw is None else _decode(raw, NLAT * NLON).reshape(NLAT, NLON)



def _arco_meta() -> dict:
    m = json.loads(_get(f"{BASE}/.zmetadata"))["metadata"]
    assert m["time/.zattrs"]["units"].startswith("hours since 1900-01-01"), m["time/.zattrs"]
    t = _decode(_get(f"{BASE}/time/0"), m["time/.zarray"]["chunks"][0] * 2).view("<i8")
    assert t[0] == 0 and np.all(np.diff(t) == 1), "ARCO time axis is not hourly from 1900-01-01"
    for v in VARS.values():
        a = m[f"{v}/.zarray"]
        assert a["chunks"] == [1, NLAT, NLON] and a["dtype"] == "<f4", (v, a)
    return m


# ── store ────────────────────────────────────────────────────────────────────

def store_path(out=None) -> Path:
    """The mirror (``out``: another path instead of the configured one)."""
    return Path(out) if out else paths().era5_hourly


def _sibling(out, suffix: str) -> Path:
    p = store_path(out)
    return p.with_name(p.name.removesuffix(".zarr") + suffix)


def _manifest_dir(out=None) -> Path:
    return _sibling(out, ".manifest")


def shard_hours(s: int) -> tuple[int, int]:
    return s * SHARD_H, min((s + 1) * SHARD_H, NT)


def hours_to_time(h) -> pd.DatetimeIndex:
    return T0 + pd.to_timedelta(np.asarray(h), unit="h")


def init(out=None, overwrite: bool = False) -> Path:
    store = store_path(out)
    if store.exists() and not overwrite:
        print(f"[init] {store} exists, reused")
        return store
    if store.exists():
        import shutil
        shutil.rmtree(store); shutil.rmtree(_manifest_dir(out), ignore_errors=True)
    m = _arco_meta()
    store.parent.mkdir(parents=True, exist_ok=True)
    g = zarr.open_group(str(store), mode="w", zarr_format=3)
    t = g.create_array("time", shape=(NT,), dtype="int64", chunks=(NT,), dimension_names=["time"],
                       attributes={"units": f"hours since {T0:%Y-%m-%d %H:%M:%S}",
                                   "calendar": "proleptic_gregorian", "standard_name": "time"})
    t[:] = np.arange(NT, dtype="int64")
    for c, n in (("latitude", NLAT), ("longitude", NLON)):
        a = g.create_array(c, shape=(n,), dtype="float32", chunks=(n,), dimension_names=[c],
                           attributes={k: v for k, v in m[f"{c}/.zattrs"].items() if not k.startswith("_")})
        a[:] = _decode(_get(f"{BASE}/{c}/0"), n)
    comp = BloscCodec(cname="zstd", clevel=3, shuffle="bitshuffle")
    for short, arco in VARS.items():
        attrs = {k: v for k, v in m[f"{arco}/.zattrs"].items() if not k.startswith("_")}
        attrs.update({"arco_name": arco,
                      "hour_convention": "value at valid hour t covers (t-1h, t]" if short in FORECAST
                      else "instantaneous at hour t"})
        g.create_array(short, shape=(NT, NLAT, NLON), dtype="float32", chunks=(CHUNK_H, NLAT, NLON),
                       shards=(SHARD_H, NLAT, NLON), fill_value=np.float32("nan"), compressors=comp,
                       dimension_names=["time", "latitude", "longitude"], attributes=attrs)
    g.attrs.update({
        "title": "ERA5 hourly single-level fields, mirror of ARCO-ERA5",
        "source": BASE,
        "arco_root_attrs": m.get(".zattrs", {}),
        "created": pd.Timestamp.now(tz="UTC").isoformat(),
        "variables": VARS,
        "layout": f"inner chunk {CHUNK_H} h, shard {SHARD_H} h (shard s = hours s*{SHARD_H} ..), ARCO grid not rolled",
        "values": "bit-identical to ARCO float32; NaN = missing value or missing ARCO chunk",
        "built_by": "chelsa_utc.era5_hourly",
    })
    zarr.consolidate_metadata(str(store))
    print(f"[init] {store}: {NT} hours ({T0:%Y-%m-%d} .. {T1 - pd.Timedelta(hours=1):%Y-%m-%d %H}), "
          f"{len(VARS)} variables, {NSHARD} shards of {SHARD_H} h")
    return store


# ── download ─────────────────────────────────────────────────────────────────

def _fetch_hours(arco: str, hours: list[int], threads: int) -> tuple[np.ndarray, list[dict]]:
    """Fetch the given hours into one (n, 721, 1440) array; QC row per hour."""
    a = np.empty((len(hours), NLAT, NLON), dtype=np.float32)
    rows: list[dict | None] = [None] * len(hours)

    def one(i: int) -> None:
        f = arco_field(arco, hours[i])
        if f is None:
            a[i] = np.nan
            rows[i] = {"hour": hours[i], "missing": True, "n_nan": NLAT * NLON,
                       "min": np.nan, "max": np.nan, "mean": np.nan}
            return
        a[i] = f
        ok = np.isfinite(f)
        nv = int(ok.sum())
        rows[i] = {"hour": hours[i], "missing": False, "n_nan": NLAT * NLON - nv,
                   "min": float(f.min(where=ok, initial=np.inf)) if nv else np.nan,
                   "max": float(f.max(where=ok, initial=-np.inf)) if nv else np.nan,
                   "mean": float(f.mean(where=ok, dtype="float64")) if nv else np.nan}

    with _cf.ThreadPoolExecutor(threads) as ex:
        list(ex.map(one, range(len(hours))))
    return a, rows


def _shard_item(out: str | None, var: str, shard: int, threads: int, hours: list[int] | None = None) -> list[dict]:
    """Fetch and write one (shard, variable); ``hours`` (smoke test): only these hours of the shard."""
    h0, h1 = shard_hours(shard)
    hrs = hours if hours is not None else list(range(h0, h1))
    a, rows = _fetch_hours(VARS[var], hrs, threads)
    arr = zarr.open_group(str(store_path(out)), mode="r+", use_consolidated=False)[var]
    if hours is None:
        arr[h0:h1] = a                                     # the whole shard in one write
    else:
        for i, h in enumerate(hrs):
            arr[h] = a[i]
    return [{"variable": var, **r} for r in rows]


def _worker_init() -> None:
    import numcodecs.blosc
    numcodecs.blosc.set_nthreads(1)


def _commit(out, shard: int, rows: list[dict]) -> None:
    d = _manifest_dir(out); d.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows).sort_values(["variable", "hour"])
    df["time"] = hours_to_time(df["hour"].to_numpy())
    f = d / f"{shard:04d}.parquet"
    tmp = f.with_suffix(".tmp"); df.to_parquet(tmp, index=False); tmp.rename(f)


def download(shards: Iterable[int] | None = None, variables=None, out=None, procs: int = 8,
             threads: int = 24, redo: bool = False, days=None) -> int:
    """Write the given shards (default: all). ``days`` (smoke test): only the hours of these
    dates, written but NOT committed."""
    variables = list(variables or VARS)
    outs = str(store_path(out)) if out else None
    ctx = _mp.get_context("spawn")
    if days:
        hrs = sorted({int((pd.Timestamp(d) - T0) / pd.Timedelta(hours=1)) + k for d in days for k in range(24)})
        by_shard: dict[int, list[int]] = {}
        for h in hrs:
            by_shard.setdefault(h // SHARD_H, []).append(h)
        with _cf.ProcessPoolExecutor(procs, mp_context=ctx, initializer=_worker_init) as ex:
            futs = [ex.submit(_shard_item, outs, v, s, threads, hh) for s, hh in by_shard.items() for v in variables]
            rows = [r for f in futs for r in f.result()]
        df = pd.DataFrame(rows)
        print(f"[download] smoke: {len(hrs)} hours x {len(variables)} variables written (not committed)")
        print(df.groupby("variable")[["missing", "n_nan"]].sum().to_string())
        return 0
    shards = sorted(set(shards if shards is not None else range(NSHARD)))
    if set(variables) != set(VARS):
        raise ValueError("a committed shard needs all variables; use --days for a subset test")
    mdir = _manifest_dir(out)
    todo_shards = [s for s in shards if redo or not (mdir / f"{s:04d}.parquet").exists()]
    for s in sorted(set(shards) - set(todo_shards)):
        print(f"[download] shard {s}: committed before, skipped")
    if not todo_shards:
        return 0
    t0 = _time.time()
    rows: dict[int, list[dict]] = {s: [] for s in todo_shards}
    left: dict[int, set] = {s: set(variables) for s in todo_shards}
    todo = [(s, v) for s in todo_shards for v in variables]
    # A killed worker breaks the whole pool (BrokenProcessPool for every pending item), so failed
    # items get up to MAX_ROUNDS fresh pools before their shard counts as failed.
    for rnd in range(1, MAX_ROUNDS + 1):
        failed = []
        with _cf.ProcessPoolExecutor(procs, mp_context=ctx, initializer=_worker_init) as ex:
            futs = {ex.submit(_shard_item, outs, v, s, threads): (s, v) for s, v in todo}
            for f in _cf.as_completed(futs):
                s, v = futs[f]
                try:
                    rows[s].extend(f.result())
                except Exception as e:                       # retried in the next round
                    failed.append((s, v))
                    if len(failed) <= 3:
                        print(f"[download] round {rnd} FAIL shard {s} {v}: {type(e).__name__}: {e}", flush=True)
                    continue
                left[s].discard(v)
                if not left[s]:
                    _commit(out, s, rows.pop(s))
                    h0, h1 = shard_hours(s)
                    print(f"[download] shard {s} ({hours_to_time([h0])[0]:%Y-%m-%d} .. "
                          f"{hours_to_time([h1 - 1])[0]:%Y-%m-%d}) committed, "
                          f"{(_time.time() - t0) / 60:.1f} min", flush=True)
        if not failed:
            break
        print(f"[download] round {rnd}: {len(failed)} items failed"
              + (", retrying in a fresh process pool" if rnd < MAX_ROUNDS else ""), flush=True)
        todo = failed
    bad = sorted(s for s, v in left.items() if v)
    if bad:
        print(f"[download] NOT committed: shards {bad}; submit again")
        return 1
    print(f"[download] {len(todo_shards)} shards in {(_time.time() - t0) / 60:.1f} min")
    return 0


# ── check ────────────────────────────────────────────────────────────────────

def _check_item(out: str | None, var: str, hour: int) -> dict:
    ref = arco_field(VARS[var], hour)
    got = zarr.open_group(str(store_path(out)), mode="r", use_consolidated=False)[var][hour]
    if ref is None:
        same = bool(np.isnan(got).all())
    else:
        same = bool(np.array_equal(ref.view("<u4"), np.ascontiguousarray(got).view("<u4")))
    return {"variable": var, "hour": hour, "time": hours_to_time([hour])[0], "arco_missing": ref is None,
            "bit_identical": same}


def check(out=None, n_sample: int = 20, seed: int = 0, workers: int = 16, days=None, variables=None) -> int:
    """Completeness and NaN summary from the manifests; bit-for-bit comparison with ARCO on
    ``n_sample`` random hours per variable (from committed shards, or from ``days``)."""
    mdir = _manifest_dir(out)
    done = sorted(int(f.stem) for f in mdir.glob("[0-9]*.parquet")) if mdir.exists() else []
    print(f"[check] committed shards: {len(done)} / {NSHARD}"
          + (f"; missing: {sorted(set(range(NSHARD)) - set(done))[:20]}" if len(done) < NSHARD else ""))
    if done:
        man = pd.concat([pd.read_parquet(mdir / f"{s:04d}.parquet") for s in done], ignore_index=True)
        n_exp = sum(shard_hours(s)[1] - shard_hours(s)[0] for s in done) * len(VARS)
        print(f"[check] manifest rows {len(man)} (expected {n_exp})")
        s = man.groupby("variable").agg(hours=("hour", "size"), missing_chunks=("missing", "sum"),
                                        hours_with_nan=("n_nan", lambda x: int((x > 0).sum())),
                                        min=("min", "min"), max=("max", "max"))
        print(s.to_string())
    rng = np.random.default_rng(seed)
    if days:
        pool = np.array([int((pd.Timestamp(d) - T0) / pd.Timedelta(hours=1)) + k for d in days for k in range(24)])
    else:
        pool = np.concatenate([np.arange(*shard_hours(s)) for s in done]) if done else np.array([], int)
    if not len(pool):
        print("[check] nothing to compare")
        return 1
    items = [(v, int(h)) for v in (variables or VARS) for h in rng.choice(pool, size=min(n_sample, len(pool)), replace=False)]
    ctx = _mp.get_context("spawn")
    outs = str(store_path(out)) if out else None
    with _cf.ProcessPoolExecutor(workers, mp_context=ctx, initializer=_worker_init) as ex:
        res = pd.DataFrame(list(ex.map(_check_item, [outs] * len(items), *zip(*items))))
    if not days:
        res.to_parquet(_sibling(out, ".check.parquet"), index=False)
    bad = res[~res.bit_identical]
    print(f"[check] bit-for-bit vs ARCO: {len(res) - len(bad)} / {len(res)} identical")
    if len(bad):
        print(bad.head(20).to_string())
    return 1 if len(bad) else 0
