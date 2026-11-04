"""CHELSA V2.1 daily on the 30" grid from daily ERA5 of a chosen day window.

The last steps of the CHELSA V2.1 algorithm (``chelsa_main.py`` of CHELSA V2.1): daily 0.25 deg ERA5
fields are brought to the 30" grid with the published daily helper maps and the static grids. The
SAGA GIS operators are re-implemented exactly in ``saga_ops`` (K = "Bicubic Spline", B = "B-Spline").

Inputs (``sources``, ``era5_daily``):
    ERA5          era5_daily.daily(day, window)  (the hourly mirror, aggregated over the window)
    helper maps   published daily tz (0.25 deg lapse rate), we (0.025 deg wind effect),
                  clt (0.025 deg cloud cover)
    static        dem_latlong, expocor, patch, orography (.sdat)
    clear sky     stot (one file per day of year)

Variables (uint16 x 0.1, 65535 = missing):
    tasmin, tasmax, tas   K(T) + K(tz) (DEM - K(orog))          T = t2m_min, t2m_max, t2m_mean
    prec                  wind_cor x B(max(P, 0)) / B(mean wind_cor)    P = tp_mean;
                          wind_cor = expocor x B(we / (1 - |DEM - BLH| / 9000)), patched with `patch`
    rsds                  srad / B((S + 0.001) / (ssrd + 0.001));
                          srad = stot(doy) x 0.01157407407 x (1 - 0.75 B(clt)^3.4), S = 0.25 deg mean of srad
    ps                    sp (1 - 0.0065 dz / T)^5.2559, dz = DEM - B(orog), T = t2m_mean (K)
hurs and sfcWind are not recomputed (already on the 00-23 UTC day): the stores hold the published files.

``v21=True`` reproduces the published V2.1 (use it with window "chelsa0606"): the V2.1 ps formula
sp (1 - L dz / (273.15 + T))^5.255 with L = -tz (T in K in SAGA's Celsius slot), and tz without the
+180 column (V2.1 temperature is NoData east of 179.875 E). ``v21=False`` (CHELSA-UTC) fixes both.
"""
from __future__ import annotations

import concurrent.futures as _cf
import json
import multiprocessing as _mp
import time as _time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from chelsa_utc import era5_daily, flags, formats, sources, store
from chelsa_utc import saga_ops as so
from chelsa_utc.config import paths
from chelsa_utc.saga_ops import NX, NY, X30, NodeGrid, coarse, rows_y

COMPUTED = ["prec", "rsds", "tasmin", "tasmax", "tas", "ps"]       # heavy first
UNCHANGED = list(sources.UNCHANGED)                                # copied from the published files
SCALE = 0.1
BAND = 720                  # tif rows per band (= chunk rows)
VERSION = "1.0"
MAX_ROUNDS = 3              # fresh process pools for failed items before giving up


# ── inputs ────────────────────────────────────────────────────────────────────

def _era5(day: pd.Timestamp, window: str, names: list[str]) -> dict[str, np.ndarray]:
    """Daily ERA5 fields on the 0.25 deg node grid (lat ascending, lon -180..179.75), float64."""
    return {k: v.astype(np.float64)[::-1] for k, v in era5_daily.daily(day, window, names).items()}


# ── the CHELSA operators, one band of tif rows at a time ─────────────────────

def _temperature(T: NodeGrid, TZ: NodeGrid, OR: NodeGrid, r0: int, r1: int) -> np.ndarray:
    y = rows_y(r0, r1)
    dem = sources.sdat_rows("dem", -32768, r0, r1)
    return T.sample(X30, y, "bicubic") + TZ.sample(X30, y, "bicubic") * (dem - OR.sample(X30, y, "bicubic"))


def _prec(P: NodeGrid, PBLH: NodeGrid, day: pd.Timestamp, r0: int, r1: int) -> np.ndarray:
    """chelsa_main.py:425-461. d2b domain +-2.25 deg, wind_cor domain +-1.25 deg."""
    D0, D1 = max(0, r0 - 270), min(NY, r1 + 270)
    W0, W1 = max(0, r0 - 150), min(NY, r1 + 150)
    d2b = np.abs(sources.sdat_rows("dem", -32768, D0, D1) - PBLH.sample(X30, rows_y(D0, D1), "bspline"))
    maxd = so.block_max(d2b, D0, D1)
    yW = rows_y(W0, W1)
    maxd2 = np.abs(maxd.sample(X30, yW, "bspline") - d2b[W0 - D0:W1 - D0])
    del d2b
    we = sources.tif_grid(sources.daily_tif("we", day), yW.min(), yW.max())
    wc = so.same_grid_bspline(we.sample(X30, yW, "bspline"), W0, W1) / (1 - maxd2 / 9000)
    del maxd2
    wc = np.where(np.isfinite(wc), wc, sources.patch_rows(W0, W1))
    wind_cor = sources.sdat_rows("expocor", -9999, W0, W1) * so.same_grid_bspline(wc, W0, W1)
    del wc
    wcoarse = so.area_mean(wind_cor, W0, W1)
    y = rows_y(r0, r1)
    return wind_cor[r0 - W0:r1 - W0] * P.sample(X30, y, "bspline") / wcoarse.sample(X30, y, "bspline")


def _rsds(ssrd: np.ndarray, day: pd.Timestamp, r0: int, r1: int) -> np.ndarray:
    """chelsa_main.py:570-627. srad domain +-1 deg."""
    W0, W1 = max(0, r0 - 120), min(NY, r1 + 120)
    yW = rows_y(W0, W1)
    st = sources.stot_rows(day, W0, W1)
    clt = sources.tif_grid(sources.daily_tif("clt", day), yW.min(), yW.max()).sample(X30, yW, "bspline")
    srad = st * 0.01157407407 * (1 - 0.75 * clt ** 3.4)
    del st, clt
    S = so.area_mean(srad, W0, W1)
    r = NodeGrid((S.G + 0.001) / (so.patch(ssrd) + 0.001), 0, 0, 1441, 721, -180.0, -90.0, 0.25)
    return srad[r0 - W0:r1 - W0] / r.sample(X30, rows_y(r0, r1), "bspline")


def _ps(SP: NodeGrid, T: NodeGrid, L: NodeGrid | None, OR: NodeGrid, r0: int, r1: int) -> np.ndarray:
    """L None: sp (1 - 0.0065 dz / T)^5.2559. L given: the V2.1 formula."""
    y = rows_y(r0, r1)
    dz = sources.sdat_rows("dem", -32768, r0, r1) - OR.sample(X30, y, "bspline")
    p = SP.sample(X30, y, "bspline"); t = T.sample(X30, y, "bspline")
    if L is None:
        return p * (1 - 0.0065 * dz / t) ** 5.2559
    return p * (1 - L.sample(X30, y, "bspline") * dz / (273.15 + t)) ** 5.255


def compute_day(var: str, day, window: str = "utc", v21: bool = False,
                rows: tuple[int, int] = (0, NY)) -> np.ndarray:
    """uint16 field of tif rows rows[0]..rows[1]-1 (default: the whole day, NY x NX)."""
    day = pd.Timestamp(day)
    orog = sources.static("orog")
    OR = coarse(orog)
    if var in ("tasmin", "tasmax", "tas"):
        name = {"tasmin": "t2m_min", "tasmax": "t2m_max", "tas": "t2m_mean"}[var]
        T, TZ = coarse(_era5(day, window, [name])[name]), coarse(sources.tz(day), patched=not v21)
        f = lambda a, b: _temperature(T, TZ, OR, a, b)
    elif var == "prec":
        e = _era5(day, window, ["tp_mean", "blh"])
        P, PBLH = coarse(np.maximum(e["tp_mean"], 0)), coarse(e["blh"] + orog)
        f = lambda a, b: _prec(P, PBLH, day, a, b)
    elif var == "rsds":
        ssrd = _era5(day, window, ["ssrd"])["ssrd"]
        f = lambda a, b: _rsds(ssrd, day, a, b)
    elif var == "ps":
        e = _era5(day, window, ["sp", "t2m_mean"])
        SP, T = coarse(e["sp"]), coarse(e["t2m_mean"])
        L = coarse(-sources.tz(day), patched=False) if v21 else None
        f = lambda a, b: _ps(SP, T, L, OR, a, b)
    else:
        raise ValueError(var)
    r_lo, r_hi = rows
    u = np.full((r_hi - r_lo, NX), 65535, dtype=np.uint16)
    for b0 in range(r_lo - r_lo % BAND, r_hi, BAND):
        a, b = max(b0, r_lo), min(b0 + BAND, r_hi)
        u[a - r_lo:b - r_lo] = so.to_u16(f(a, b), SCALE)
    return u


def day_raw(var: str, day, window: str = "utc", v21: bool = False) -> tuple[np.ndarray, float]:
    """(uint16, scale) of one global day: computed, or the published file for hurs / sfcWind."""
    if var in UNCHANGED:
        raw, scale, offset, _ = formats.read_tif(sources.daily_tif(var, day))
        if offset != 0.0 or scale != store.VARIABLES[var][0]:
            raise ValueError(f"{var} {day}: unexpected scale {scale} / offset {offset} in the published file")
        return raw, scale
    return compute_day(var, day, window, v21), SCALE


# ── the output store ──────────────────────────────────────────────────────────

def store_path(window: str) -> Path:
    return paths().chelsa[window]


def _attrs(window: str, v21: bool) -> dict:
    if window == "utc" and not v21:
        return {
            "title": "CHELSA-UTC: CHELSA V2.1 daily on the 00-24 UTC day", "version": VERSION,
            "day_window": "00-24 UTC for every variable (forecast fields over valid hours (00,24]; "
                          "analysis fields as means of 00..23 UTC)",
            "source": "CHELSA V2.1 algorithm (last steps re-implemented, SAGA-exact) on hourly ARCO-ERA5; "
                      "published daily tz, we, clt, hurs, sfcWind and static grids of CHELSA V2.1",
            "computed": COMPUTED, "copied_from_v21": UNCHANGED,
            "storage": "uint16 * scale_factor, 65535 = missing; no per-day offsets",
            "fixes_vs_v21": ["window of tasmin/tasmax/prec/rsds", "ps formula (T in K in SAGA's Celsius slot)",
                             "tz +180 column (temperature east of 179.875 E)", "broken ERA5 input days repaired"],
            "known_limits": ["rsds: clt keeps the 06->06 cloud base and V2.1's one-column cloud shift",
                             "prec/rsds: flagged on days where the published we/clt were built from broken inputs",
                             "north of 83.75 N (ocean only): SAGA NoData fill approximated"],
            "notes": flags.NOTES, "built_by": "chelsa_utc.chelsa",
        }
    return {"title": f"CHELSA V2.1 daily algorithm on the {window} day window" + (" (V2.1 formulas)" if v21 else ""),
            "day_window": window, "v21_formulas": v21, "computed": COMPUTED, "copied_from_v21": UNCHANGED,
            "storage": "uint16 * scale_factor, 65535 = missing; no per-day offsets",
            "built_by": "chelsa_utc.chelsa"}


def _qc(u: np.ndarray, var: str, scale: float) -> dict:
    lo, hi = flags.PIXEL_RANGE[var]
    n_valid = n_out = 0; s = 0.0; vmin, vmax = np.inf, -np.inf
    for r0 in range(0, u.shape[0], BAND):
        b = u[r0:r0 + BAND]
        ok = b != 65535
        nv = int(ok.sum())
        if nv:
            xv = b[ok].astype(np.float64) * scale
            n_valid += nv; s += float(xv.sum()); vmin = min(vmin, float(xv.min())); vmax = max(vmax, float(xv.max()))
            n_out += int(((xv < lo) | (xv > hi)).sum())
    return {"n_valid": n_valid, "n_total": int(u.size), "n_out": n_out, "min": vmin if n_valid else None,
            "max": vmax if n_valid else None, "mean": s / n_valid if n_valid else None}


def _worker_init() -> None:
    import numcodecs.blosc
    import zarr
    numcodecs.blosc.set_nthreads(1)
    zarr.config.set({"threading.max_workers": 2})


def _item(path: str, window: str, v21: bool, var: str, t: int, day: str) -> dict:
    t0 = _time.time()
    u, scale = day_raw(var, day, window, v21)
    store.open_store(path, "r+")[var][t] = u
    return {"date": day, "variable": var, **_qc(u, var, scale), "seconds": round(_time.time() - t0, 1)}


def plan(window: str, first: str, last: str, variables: Iterable[str] | None = None, redo: bool = False):
    """(items to compute, days not on the store's time axis); items = (var, t, "YYYY-MM-DD")."""
    path = store_path(window)
    days = pd.date_range(first, last, freq="D")
    variables = list(variables or COMPUTED + UNCHANGED)
    if not path.exists():
        return [(v, t, f"{d:%Y-%m-%d}") for t, d in enumerate(days) for v in variables], []
    idx = store.index(store.open_store(path))
    items, off_axis = [], []
    for d in days:
        if d not in idx:
            off_axis.append(d)
            continue
        t = idx[d]
        items += [(v, t, f"{d:%Y-%m-%d}") for v in variables if redo or not store.written(path, v, t)]
    return items, off_axis


def build(window: str = "utc", first: str = "1941-01-01", last: str = "1978-12-31",
          variables: Iterable[str] | None = None, v21: bool | None = None, workers: int = 8,
          redo: bool = False, dry_run: bool = False) -> int:
    """Fill the store of one window for the days first..last, one process per (variable, day).

    The store (config ``paths.chelsa.<window>``) is created if it does not exist, with the time axis
    first..last. An existing store is a cache: (variable, day) pairs already written are skipped, days
    off its time axis are reported. Memory: prec needs ~10 GB per process, the others 4-6 GB.
    ``v21`` default: True for window "chelsa0606" (rebuild V2.1), False for "utc" (CHELSA-UTC)."""
    v21 = (window == "chelsa0606") if v21 is None else v21
    path = store_path(window)
    variables = list(variables or COMPUTED + UNCHANGED)
    items, off_axis = plan(window, first, last, variables, redo)
    n_days = len(pd.date_range(first, last, freq="D"))
    print(f"[build] {window} (v21 formulas: {v21}) -> {path}")
    print(f"[build] {n_days} days x {len(variables)} variables: {len(items)} to compute, "
          f"{n_days * len(variables) - len(items) - len(off_axis) * len(variables)} already written, "
          f"{len(off_axis)} days not on the time axis")
    for d in off_axis[:20]:
        print(f"[build]   not on the time axis: {d:%Y-%m-%d}")
    if dry_run or not items:
        return 0
    if not path.exists():
        store.create(path, pd.date_range(first, last, freq="D"), list(store.VARIABLES), _attrs(window, v21))
    qc_dir = path.with_name(path.name.removesuffix(".zarr") + ".qc"); qc_dir.mkdir(exist_ok=True)
    ctx = _mp.get_context("spawn")
    t0 = _time.time(); rows = []; todo = list(items)
    # A killed worker breaks the whole pool (BrokenProcessPool for every pending item), so
    # failed items get up to MAX_ROUNDS fresh pools.
    for rnd in range(1, MAX_ROUNDS + 1):
        failed = []
        with _cf.ProcessPoolExecutor(workers, mp_context=ctx, initializer=_worker_init) as ex:
            futs = {ex.submit(_item, str(path), window, v21, v, t, d): (v, t, d) for v, t, d in todo}
            for i, f in enumerate(_cf.as_completed(futs), 1):
                v, t, d = futs[f]
                try:
                    rows.append(f.result())
                except Exception as e:                       # retried in the next round
                    failed.append((v, t, d))
                    if len(failed) <= 5:
                        print(f"[build] round {rnd} FAIL {v} {d}: {type(e).__name__}: {e}", flush=True)
                if i % 50 == 0:
                    print(f"[build] round {rnd}: {i}/{len(todo)} items, {(_time.time() - t0) / 60:.1f} min", flush=True)
        if not failed:
            break
        print(f"[build] round {rnd}: {len(failed)} items failed"
              + (", retrying in a fresh process pool" if rnd < MAX_ROUNDS else ""), flush=True)
        todo = failed
    if rows:
        f = qc_dir / f"{pd.Timestamp.now():%Y%m%dT%H%M%S}.parquet"
        pd.DataFrame(rows).sort_values(["date", "variable"]).to_parquet(f, index=False)
        print(f"[build] QC of {len(rows)} items -> {f}")
    print(f"[build] {len(rows)} items in {(_time.time() - t0) / 60:.1f} min, {len(failed)} failed"
          + (" (missing inputs? see the messages above; run again to retry)" if failed else ""))
    return 1 if failed else 0


def finalize(window: str = "utc") -> Path:
    """Flag file ``<store>.flags.json`` from the QC of all builds and the known source issues."""
    path = store_path(window)
    qc_dir = path.with_name(path.name.removesuffix(".zarr") + ".qc")
    files = sorted(qc_dir.glob("*.parquet")) if qc_dir.exists() else []
    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True) if files else pd.DataFrame()
    if len(df):
        df = df.drop_duplicates(["date", "variable"], keep="last")
    on_axis = {f"{d:%Y-%m-%d}" for d in store.dates(store.open_store(path))}
    records = {}
    for r in df.to_dict("records"):
        reason = flags.flag(r)
        if reason:
            records[(r["date"], r["variable"])] = {"variable": r["variable"], "date": r["date"], "reason": reason,
                                                   "min": r.get("min"), "max": r.get("max"), "mean": r.get("mean")}
    if window == "utc":
        for (d, v), why in flags.KNOWN_ISSUES.items():
            if d in on_axis and (d, v) not in records:
                records[(d, v)] = {"variable": v, "date": d, "reason": f"known source issue: {why}"}
    recs = sorted(records.values(), key=lambda r: (r["date"], r["variable"]))
    out = path.with_name(path.name.removesuffix(".zarr") + ".flags.json")
    out.write_text(json.dumps({"store": path.name, "qc_items": len(df), "dates": sorted({r["date"] for r in recs}),
                               "records": recs}, indent=1, default=lambda x: None))
    print(f"[finalize] {len(recs)} flagged (variable, day) pairs from {len(df)} QC items -> {out}")
    return out


# ── compare with a reference (published V2.1 or another store) ────────────────

def compare(var: str, day, window: str, v21: bool, reference: str | Path) -> dict:
    """Compute one (variable, day) and compare it with ``reference``: a zarr store, or the root of a
    published V2.1 tree (``<root>/daily/<var>/<year>/...``). Integers are compared where the
    reference has no offset; with an offset (V2.1 tas of 1941-1950) physical values within half a step."""
    t0 = _time.time()
    u, scale = day_raw(var, day, window, v21)
    ref, rscale, roff = formats.read_raw(reference, var, day)
    if rscale == scale and roff == 0.0:
        both = (u != 65535) & (ref != 65535)
        dd = np.abs(u[both].astype(np.int64) - ref[both].astype(np.int64))
        eq = dd == 0
        dmax = float(dd.max()) if dd.size else 0.0
    else:
        x = store.decode(u, scale).astype(np.float64); o = store.decode(ref, rscale, roff).astype(np.float64)
        both = np.isfinite(x) & np.isfinite(o)
        d = np.abs(x[both] - o[both])
        eq = d <= 0.5 * scale + 1e-6
        dmax = float(d.max() / scale) if d.size else 0.0
    nd = int(((u == 65535) != (ref == 65535)).sum())
    n = int(both.sum())
    return {"date": f"{pd.Timestamp(day):%Y-%m-%d}", "variable": var, "window": window, "v21": v21,
            "reference": str(reference), "n": n, "identical": float(eq.sum()) / n if n else None,
            "max_abs_steps": dmax, "nodata_mismatch": nd, "seconds": round(_time.time() - t0, 1)}
