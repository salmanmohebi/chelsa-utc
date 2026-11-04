# chelsa-utc

CHELSA V2.1 daily climate data recomputed on the **00–24 UTC day** (CHELSA-UTC), with tools to
convert between CHELSA's native GeoTIFF files and zarr in both directions.

Published CHELSA V2.1 daily mixes two day windows inside one "day":

| Variables | Day window in CHELSA V2.1 | In CHELSA-UTC |
|---|---|---|
| tasmin, tasmax, prec, rsds | 06 UTC *d* → 06 UTC *d+1* (ERA5 forecast fields) | 00–24 UTC *d* |
| tas, ps, hurs, sfcWind | 00–23 UTC *d* (ERA5 analyses) | 00–23 UTC *d* (unchanged) |

This repository re-implements the last step of the CHELSA V2.1 algorithm in numpy, with operators
that copy SAGA GIS exactly, and runs it on ERA5 aggregated over the 00–24 UTC day from hourly
ERA5. With the original 06→06 inputs, the same code reproduces the published V2.1 files: more than
99.99 % of pixels have the identical stored integer, and the rest differ by one storage step
([docs/validation.md](docs/validation.md)). CHELSA-UTC is therefore a deterministic recomputation,
not a statistical correction.

Period: 1941-01-01 to 1978-12-31 (13 864 of 13 879 days; 15 days have no published `tz` / `clt`).
Inputs: the CHELSA V2.1 files on the public bucket on 2026-11-04 (checked file by file, see
[docs/known-issues.md](docs/known-issues.md)).
Grid: global 30″ (43200 × 20880 pixels, 84° N to 90° S).

## What changes in CHELSA-UTC

| Variable | CHELSA-UTC |
|---|---|
| tasmin, tasmax | ERA5 mn2t / mx2t over 00–24 UTC; the lapse rate `tz` gets the +180° column (V2.1 is NoData east of 179.875° E) |
| prec | ERA5 precipitation over 00–24 UTC; the published daily wind effect `we` |
| rsds | ERA5 ssrd over 00–24 UTC; the published daily cloud cover `clt` |
| tas | same window, recomputed (all days use one encoding; broken V2.1 input days repaired) |
| ps | standard barometric formula; V2.1 passes the temperature in K where SAGA expects °C (+56 hPa error per 1000 m of height difference) |
| hurs, sfcWind | the published files, unchanged (already on the 00–23 UTC day) |

All ERA5 values come from hourly ERA5, so the broken daily ERA5 input files of V2.1 are repaired.
Days where a published helper map (`tz`, `we`, `clt`) is broken cannot be repaired and are
flagged ([docs/known-issues.md](docs/known-issues.md)). One 06→06 piece is left: the published
cloud cover `clt` used for rsds was built with the cloud base height of the 06→06 window (issue W1,
not fixed; it changes rsds by about 0.3 W m⁻² RMS). Mean changes are near zero; single days
change: tasmin 0.79 K RMS, tasmax 0.73 K, prec 2.36 mm/day, rsds 19 W m⁻², ps 5.35 hPa
([docs/validation.md](docs/validation.md)).

## Getting CHELSA-UTC

The CHELSA-UTC data (1941–1978, 8 variables, ~23 TB as zarr) will be published soon; the link will be added here.
Until then, compute it with this repository: download the **inputs** (published CHELSA V2.1 files
and hourly ERA5, all public) and run `build` (see the quick start). No command of this repository
downloads CHELSA-UTC itself yet.

## Documentation

| File | Content |
|---|---|
| [docs/method.md](docs/method.md) | the day windows, the SAGA operators and the formulas, what CHELSA-UTC changes |
| [docs/data-and-formats.md](docs/data-and-formats.md) | inputs and where to get them, the zarr and GeoTIFF formats, conversion in both directions, caching |
| [docs/validation.md](docs/validation.md) | how the code was checked and the results |
| [docs/known-issues.md](docs/known-issues.md) | problems found in CHELSA V2.1 daily, and which are repaired or flagged in CHELSA-UTC |

## Install

Python ≥ 3.11, CPU only (numpy, scipy, zarr, rasterio).

```bash
git clone https://github.com/salmanmohebi/chelsa-utc.git && cd chelsa-utc
uv sync                      # or: pip install -e .
cp config.example.toml config.toml   # then edit the paths
uv run chelsa-utc info       # shows the resolved paths and what already exists
```

## Quick start

**1. Download the inputs.** These are published CHELSA V2.1 files and ERA5, not CHELSA-UTC. All
are public, no account needed:

- CHELSA V2.1 static grids, clear-sky radiation and daily helper maps, from the public CHELSA
  buckets of WSL (`https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/`,
  `https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/daily_normals/rsdscs/`);
- hourly ERA5 from ARCO-ERA5 on Google Cloud
  (`https://storage.googleapis.com/gcp-public-data-arco-era5/`).

```bash
chelsa-utc download static                                         # V2.1: 4 static grids, ~5.8 GB
chelsa-utc download stot                                           # V2.1: clear-sky radiation, 366 files, ~393 GB
chelsa-utc download daily --vars tz,we,clt,hurs,sfcWind --first 1976-07-01 --last 1976-07-31   # V2.1 daily files
chelsa-utc era5-hourly init
chelsa-utc era5-hourly download --first 1976-07-01 --last 1976-07-31   # ERA5: ~19 GB per 31 days
```

`tz`, `we` and `clt` are the V2.1 helper maps the computation needs; `hurs` and `sfcWind` are V2.1
outputs that CHELSA-UTC keeps unchanged.

**2. Compute CHELSA-UTC.**

One day as a GeoTIFF in the original CHELSA format:

```bash
chelsa-utc day tasmin 1976-07-15 out/CHELSA-UTC_tasmin_15_07_1976.tif            # CHELSA-UTC
chelsa-utc day tasmin 1976-07-15 out/V21_tasmin_15_07_1976.tif --window chelsa0606  # V2.1 rebuilt
```

The same in Python:

```python
from chelsa_utc import chelsa, config

config.configure("config.toml")
u = chelsa.compute_day("tasmin", "1976-07-15", window="utc")                  # CHELSA-UTC
v = chelsa.compute_day("tasmin", "1976-07-15", window="chelsa0606", v21=True)  # V2.1 rebuilt
# uint16 (20880, 43200); value = 0.1 * u, 65535 = NoData
```

**3. A whole period** into a zarr store (one process per variable and day; prec needs ~10 GB per process):

```bash
chelsa-utc build --window utc --first 1976-01-01 --last 1976-12-31 --workers 16 --dry-run
chelsa-utc build --window utc --first 1976-01-01 --last 1976-12-31 --workers 16
chelsa-utc finalize --window utc        # flag file from the build QC and the known issues
```

## Back to the original CHELSA format

Two directions are supported:

1. **File format.** Any store converts losslessly to native CHELSA GeoTIFFs (Cloud Optimized
   GeoTIFF, uint16, the V2.1 file names and folder layout), and native files convert to zarr:

   ```bash
   chelsa-utc zarr-to-tif out/chelsa-utc-tif --window utc --vars tasmin,prec --first 1976-07-01 --last 1976-07-31
   chelsa-utc tif-to-zarr out/v21.zarr --vars tasmin,tas --first 1976-07-01 --last 1976-07-31   # published files -> zarr
   ```

   The output tree is `<dst>/daily/<var>/<year>/CHELSA_<var>_<dd>_<mm>_<yyyy>_V.2.1.tif`, the same
   as the published tree, so tools that read the published files also read these. The GeoTIFF tags
   record the day window (`day_window`) and the product (`source_id`). Use `--name` for another
   file name pattern.

2. **Day window.** `--window chelsa0606` runs the same code on the 06→06 window with the V2.1
   formulas, which rebuilds the published V2.1 from hourly ERA5 (99.5–99.9 % identical integers,
   the rest ±1 step).

To read either format in Python (decoded values, NaN = missing):

```python
from chelsa_utc import formats
ds = formats.open_chelsa("out/chelsa-utc.zarr", ["tasmin"], ["1976-07-15"], bbox=(5, 45, 11, 48))
ds = formats.open_chelsa("data/chelsa-v2.1", ["tasmin"], ["1976-07-15"], bbox=(5, 45, 11, 48))  # native tree
```

## Reusing existing files

Every path in `config.toml` can point to files that already exist. They are used as a cache:

- `era5_hourly`: an existing hourly mirror; `era5-hourly download` skips committed shards.
- `era5_daily.<window>`: an existing daily ERA5 store; a day in it is read, not recomputed (the
  values are bit-identical to a recomputation, [docs/validation.md](docs/validation.md)).
- `chelsa_v21`, `stot`: existing copies of the published files; `download` fetches only missing files.
- `chelsa.<window>`: an existing output store; `build` computes only the (variable, day) pairs that
  are not written yet (a pair is written when its shard file exists). Days that are not on the
  store's time axis are reported, not added. Run `build --dry-run` first to see what would be computed.

To add days to an existing store (its time axis is fixed), make a new store with `relink`: it
hard-links every written (variable, day) of the old store (no copy, no extra disk space, the old
store stays unchanged) and adds the missing days whose inputs exist; `build` then computes only
the new days:

```bash
chelsa-utc relink out/chelsa-utc-v1.1.zarr --src out/chelsa-utc.zarr --dry-run   # lists the days to add
chelsa-utc relink out/chelsa-utc-v1.1.zarr --src out/chelsa-utc.zarr
chelsa-utc build --window utc      # with paths.chelsa.utc = out/chelsa-utc-v1.1.zarr
```

A (variable, day) whose result is all NoData (e.g. rsds on the 10 days with an empty published
`clt`) writes no shard, so `build` computes it again on each run; this costs a few minutes.

## Commands

| Command | What it does |
|---|---|
| `info` | resolved paths, what exists, how complete the caches are |
| `download static\|stot\|daily` | the CHELSA V2.1 **inputs** (not CHELSA-UTC) from the public CHELSA buckets (missing files only) |
| `era5-hourly init\|download\|check` | hourly ERA5 mirror from ARCO-ERA5 (bit-identical, restartable) |
| `era5-daily` | daily ERA5 store of a window (optional cache) |
| `build` | the 30″ store of a window |
| `finalize` | flag file of a store |
| `day` | one (variable, day) as a native GeoTIFF |
| `compare` | compute days and compare them with the published files or a store |
| `tif-to-zarr`, `zarr-to-tif` | lossless format conversion |
| `relink` | new store with more days; written days are hard-linked from the old store |
| `flags` | the flag file of the released CHELSA-UTC v1.0 |

Run `chelsa-utc <command> --help` for the options.

## Cost

On one CPU core, one (variable, day) takes 1–2.5 min (temperatures, ps, rsds) or 3–7 min (prec),
depending on the CPU and the file system. The full 1941–1978 period (6 variables × 13 879 days) is about 175 node-hours on
128-core, 256 GB nodes. The store is about 23 TB: 15 TB for the 6 computed variables and 8 TB for the
copied hurs and sfcWind. The hourly ERA5 mirror for 1940–1979 is 8.4 TB.

## Tests

```bash
uv run pytest
```

The tests need no data. The checks against real data are the `compare` command
([docs/validation.md](docs/validation.md)).

## Credits and license

CHELSA: Karger, D. N. et al. (2017), *Climatologies at high resolution for the Earth's land surface
areas*, Scientific Data 4, 170122, https://doi.org/10.1038/sdata.2017.122; CHELSA V2.1 daily:
Karger, D. N. (2025), *CHELSA-daily climate data at high resolution*, EnviDat,
https://doi.org/10.16904/envidat.687. The CHELSA V2.1 code is published under GPL v3. ERA5:
Hersbach, H. et al. (2020), https://doi.org/10.1002/qj.3803, read from the public ARCO-ERA5 store
(Google Research). SAGA GIS: Conrad, O. et al. (2015), https://doi.org/10.5194/gmd-8-1991-2015.

This repository re-implements parts of the CHELSA V2.1 algorithm and is released under the
GNU General Public License v3.0 or later ([LICENSE](LICENSE)).
