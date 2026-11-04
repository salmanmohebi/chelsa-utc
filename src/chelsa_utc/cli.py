"""Command line: ``chelsa-utc <command> ...`` (``python -m chelsa_utc`` works too). See README."""
from __future__ import annotations

import argparse
import json
import sys

import pandas as pd

from chelsa_utc import config


def _days(a) -> pd.DatetimeIndex:
    if getattr(a, "dates", None):
        return pd.DatetimeIndex(a.dates.split(","))
    return pd.date_range(a.first, a.last, freq="D")


def _csv(s: str | None) -> list[str] | None:
    return s.split(",") if s else None


def cmd_info(a) -> int:
    from chelsa_utc import era5_hourly, store
    p = config.paths()
    print(f"config: {p.source}")
    rows = [("era5_hourly", p.era5_hourly), ("chelsa_v21", p.chelsa_v21), ("stot", p.stot)]
    rows += [(f"era5_daily.{w}", q) for w, q in p.era5_daily.items()]
    rows += [(f"chelsa.{w}", q) for w, q in p.chelsa.items()]
    for k, q in rows:
        state = "missing"
        if q.exists():
            state = "exists"
            if store.is_store(q) and (q / "time").exists():
                g = store.open_store(q)
                if g["time"].attrs["units"].startswith("days"):
                    d = store.dates(g)
                    state = f"zarr, {len(d)} days {d[0]:%Y-%m-%d} .. {d[-1]:%Y-%m-%d}"
                else:
                    state = f"zarr, {g['time'].shape[0]} steps ({g['time'].attrs['units']})"
        print(f"  {k:22s} {state:44s} {q}")
    if p.era5_hourly.exists():
        m = era5_hourly._manifest_dir()
        n = len(list(m.glob("[0-9]*.parquet"))) if m.exists() else 0
        print(f"  hourly mirror: {n} / {era5_hourly.NSHARD} shards committed")
    if p.stot.exists():
        print(f"  stot: {len(list(p.stot.glob('CHELSA_stot_pj_*_V.2.1.tif')))} / 366 files")
    return 0


def cmd_download(a) -> int:
    from chelsa_utc import sources
    if a.what == "static":
        return sources.download_static(a.threads)
    if a.what == "stot":
        return sources.download_stot(a.threads)
    return sources.download_daily(_csv(a.vars) or list(sources.HELPER), _days(a), a.threads)


def cmd_era5_hourly(a) -> int:
    from chelsa_utc import era5_hourly
    if a.step == "init":
        era5_hourly.init()
        return 0
    if a.step == "download":
        shards = None
        if a.shards:
            shards = [int(x) for x in a.shards.split(",")]
        elif a.first:
            h0 = int((pd.Timestamp(a.first) - era5_hourly.T0) / pd.Timedelta(hours=1))
            h1 = int((pd.Timestamp(a.last) + pd.Timedelta(days=2) - era5_hourly.T0) / pd.Timedelta(hours=1))
            shards = list(range(max(0, h0 - 24) // era5_hourly.SHARD_H, min(era5_hourly.NT - 1, h1) // era5_hourly.SHARD_H + 1))
        return era5_hourly.download(shards, procs=a.workers, threads=a.threads)
    return era5_hourly.check(workers=a.workers)


def cmd_era5_daily(a) -> int:
    from chelsa_utc import era5_daily
    era5_daily.build(a.window, a.first, a.last, workers=a.workers, redo=a.redo)
    return 0


def cmd_build(a) -> int:
    from chelsa_utc import chelsa
    return chelsa.build(a.window, a.first, a.last, _csv(a.vars), v21=a.v21, workers=a.workers,
                        redo=a.redo, dry_run=a.dry_run)


def cmd_finalize(a) -> int:
    from chelsa_utc import chelsa
    chelsa.finalize(a.window)
    return 0


def cmd_day(a) -> int:
    from chelsa_utc import chelsa, formats
    v21 = (a.window != "utc") if a.v21 is None else a.v21
    u, scale = chelsa.day_raw(a.var, a.date, a.window, v21)
    extra = {"day_window": a.window, "source_id": "CHELSA-UTC" if a.window == "utc" else "CHELSA_V2.1 (rebuilt)"}
    print(formats.write_tif(a.out, u, scale, 0.0, formats.tags_for(a.var, a.date, extra)))
    return 0


def cmd_compare(a) -> int:
    from chelsa_utc import chelsa
    v21 = (a.window != "utc") if a.v21 is None else a.v21
    ref = config.paths().chelsa_v21 if a.reference == "v21" else a.reference
    rows = []
    for d in _days(a):
        for v in _csv(a.vars) or chelsa.COMPUTED:
            r = chelsa.compare(v, d, a.window, v21, ref)
            rows.append(r)
            print(f"[compare] {v:7s} {r['date']} identical {r['identical']:.6f}  max {r['max_abs_steps']:.0f} steps  "
                  f"nodata mismatch {r['nodata_mismatch']}  ({r['seconds']:.0f} s)", flush=True)
    if a.out:
        pd.DataFrame(rows).to_parquet(a.out, index=False)
        print(f"[compare] -> {a.out}")
    return 0


def cmd_tif_to_zarr(a) -> int:
    from chelsa_utc import formats
    src = a.src or config.paths().chelsa_v21
    return formats.tif_to_zarr(src, a.dst, a.vars.split(","), a.first, a.last, a.workers, a.redo)


def cmd_zarr_to_tif(a) -> int:
    from chelsa_utc import formats
    src = a.src or config.paths().chelsa[a.window]
    return formats.zarr_to_tif(src, a.dst, _csv(a.vars), a.first, a.last, a.name, a.workers, a.overwrite)


def cmd_relink(a) -> int:
    from chelsa_utc import sources, store
    src = a.src or config.paths().chelsa[a.window]
    g = store.open_store(src)
    days = set(store.dates(g))
    span = pd.date_range(a.first or min(days), a.last or max(days), freq="D")
    add = [d for d in span if d not in days and (a.all_days or sources.inputs_present(d))]
    print(f"[relink] adding {len(add)} days: " + " ".join(f"{d:%Y-%m-%d}" for d in add))
    if a.dry_run:
        return 0
    store.relink(src, a.dst, sorted(days | set(add)), note=f"added {len(add)} days with inputs")
    return 0


def cmd_flags(a) -> int:
    from chelsa_utc import flags
    print(json.dumps(flags.released_flags(), indent=1))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="chelsa-utc", description=__doc__)
    ap.add_argument("--config", help="settings file (default: $CHELSA_UTC_CONFIG or ./config.toml)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def days_args(p, required=True):
        p.add_argument("--first", default="1941-01-01" if not required else None, required=required)
        p.add_argument("--last", default="1978-12-31" if not required else None, required=required)

    p = sub.add_parser("info", help="resolved paths and what exists already (cache state)")
    p.set_defaults(f=cmd_info)

    p = sub.add_parser("download", help="the CHELSA V2.1 inputs (not CHELSA-UTC) from the public CHELSA buckets (missing files only)")
    p.add_argument("what", choices=["static", "stot", "daily"])
    p.add_argument("--vars", help="daily: comma list (default tz,we,clt; add hurs,sfcWind for the unchanged outputs)")
    p.add_argument("--first"); p.add_argument("--last"); p.add_argument("--dates", help="comma list instead of first/last")
    p.add_argument("--threads", type=int, default=8)
    p.set_defaults(f=cmd_download)

    p = sub.add_parser("era5-hourly", help="hourly ERA5 (input) from the public ARCO-ERA5 store: init | download | check")
    p.add_argument("step", choices=["init", "download", "check"])
    p.add_argument("--shards", help="comma list of shard numbers (744 h each)")
    p.add_argument("--first", help="download: the shards that cover first..last (with the window margins)")
    p.add_argument("--last")
    p.add_argument("--workers", type=int, default=8); p.add_argument("--threads", type=int, default=24)
    p.set_defaults(f=cmd_era5_hourly)

    p = sub.add_parser("era5-daily", help="daily ERA5 store of one window (cache of the daily aggregation)")
    p.add_argument("--window", choices=config.WINDOWS, default="utc"); days_args(p, False)
    p.add_argument("--workers", type=int, default=16); p.add_argument("--redo", action="store_true")
    p.set_defaults(f=cmd_era5_daily)

    p = sub.add_parser("build", help="compute the 30\" store of one window (skips what is written already)")
    p.add_argument("--window", choices=config.WINDOWS, default="utc"); days_args(p, False)
    p.add_argument("--vars", help="comma list (default: all 8)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--v21", dest="v21", action="store_true", default=None, help="V2.1 formulas (default for chelsa0606)")
    g.add_argument("--no-v21", dest="v21", action="store_false", help="CHELSA-UTC formulas (default for utc)")
    p.add_argument("--workers", type=int, default=8, help="processes; prec needs ~10 GB each")
    p.add_argument("--redo", action="store_true"); p.add_argument("--dry-run", action="store_true")
    p.set_defaults(f=cmd_build)

    p = sub.add_parser("finalize", help="flag file of a store from the build QC and the known issues")
    p.add_argument("--window", choices=config.WINDOWS, default="utc")
    p.set_defaults(f=cmd_finalize)

    p = sub.add_parser("day", help="compute one (variable, day) and write it as a native GeoTIFF")
    p.add_argument("var"); p.add_argument("date"); p.add_argument("out")
    p.add_argument("--window", choices=config.WINDOWS + ("v21files",), default="utc",
                   help="v21files: ERA5 from the V2.1 daily input files (paths.v21_era5_inputs)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--v21", dest="v21", action="store_true", default=None)
    g.add_argument("--no-v21", dest="v21", action="store_false")
    p.set_defaults(f=cmd_day)

    p = sub.add_parser("compare", help="compute days and compare them with published V2.1 or a store")
    p.add_argument("--window", choices=config.WINDOWS + ("v21files",), default="utc",
                   help="v21files: ERA5 from the V2.1 daily input files (paths.v21_era5_inputs)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--v21", dest="v21", action="store_true", default=None)
    g.add_argument("--no-v21", dest="v21", action="store_false")
    p.add_argument("--reference", required=True, help="'v21' (the published files) or a zarr store / native tree")
    p.add_argument("--dates", required=True, help="comma list")
    p.add_argument("--vars", help="comma list (default: the 6 computed variables)")
    p.add_argument("--out", help="parquet file for the results")
    p.set_defaults(f=cmd_compare)

    p = sub.add_parser("tif-to-zarr", help="native CHELSA GeoTIFFs -> zarr store (lossless)")
    p.add_argument("dst"); p.add_argument("--src", help="root of the native tree (default: paths.chelsa_v21)")
    p.add_argument("--vars", required=True); days_args(p)
    p.add_argument("--workers", type=int, default=8); p.add_argument("--redo", action="store_true")
    p.set_defaults(f=cmd_tif_to_zarr)

    p = sub.add_parser("zarr-to-tif", help="zarr store -> native CHELSA GeoTIFFs (lossless)")
    p.add_argument("dst"); p.add_argument("--src", help="zarr store (default: the store of --window)")
    p.add_argument("--window", choices=config.WINDOWS, default="utc")
    p.add_argument("--vars"); p.add_argument("--first"); p.add_argument("--last")
    p.add_argument("--name", default="CHELSA_{var}_{day:%d_%m_%Y}_V.2.1.tif", help="file name format")
    p.add_argument("--workers", type=int, default=4); p.add_argument("--overwrite", action="store_true")
    p.set_defaults(f=cmd_zarr_to_tif)

    p = sub.add_parser("relink", help="new store with more days on the time axis; written days are hard-linked")
    p.add_argument("dst"); p.add_argument("--src", help="zarr store (default: the store of --window)")
    p.add_argument("--window", choices=config.WINDOWS, default="utc")
    p.add_argument("--first"); p.add_argument("--last")
    p.add_argument("--all-days", action="store_true", help="add every missing day, not only days whose inputs exist")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(f=cmd_relink)

    p = sub.add_parser("flags", help="print the flag file of the released CHELSA-UTC v1.0")
    p.set_defaults(f=cmd_flags)

    a = ap.parse_args(argv)
    config.configure(a.config)
    return a.f(a)


if __name__ == "__main__":
    sys.exit(main())
