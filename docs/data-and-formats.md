# Data and formats

## 1. Inputs

These are the inputs from which CHELSA-UTC is computed (published CHELSA V2.1 files and ERA5), not
CHELSA-UTC itself; the CHELSA-UTC data will be published separately. All inputs are public. The paths are set in `config.toml` (format: `config.example.toml`).

| Input | Config key | Source | Size | Used for |
|---|---|---|---|---|
| Hourly ERA5, 11 single-level fields, 1940–1979 | `era5_hourly` | ARCO-ERA5, `https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3` (`chelsa-utc era5-hourly`) | 8.4 TB (19 GB per 31 days) | all ERA5 fields |
| Static grids `dem_latlong`, `expocor`, `patch`, `orography` (SAGA `.sdat` + `.sgrd`) | `chelsa_v21` → `input/static/` | `https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/input/static/` (`download static`) | 5.8 GB | all variables |
| Daily lapse rate `tz` (0.25°) | `chelsa_v21` → `daily/tz/` | `.../chelsa/global/daily/tz/<year>/` (`download daily --vars tz`) | ~4 MB per day | tasmin, tasmax, tas |
| Daily wind effect `we` (0.025°) | `chelsa_v21` → `daily/we/` | `.../daily/we/<year>/` | ~100 MB per day | prec |
| Daily cloud cover `clt` (0.025°) | `chelsa_v21` → `daily/clt/` | `.../daily/clt/<year>/` | ~100 MB per day | rsds |
| Daily hurs, sfcWind (30″) | `chelsa_v21` → `daily/{hurs,sfcWind}/` | `.../daily/<var>/<year>/` | 300–450 MB per day | copied unchanged |
| Clear-sky radiation `stot`, one file per day of year | `stot` | `https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/daily_normals/rsdscs/CHELSA_stot_pj_<doy>_V.2.1.tif` (`download stot`) | 393 GB (366 files) | rsds |
| Published V2.1 outputs (optional) | `chelsa_v21` → `daily/<var>/` | `.../daily/<var>/<year>/` | 150–450 MB per file | `compare`, `tif-to-zarr` |

For a test of the method, one month of inputs is enough: the static grids, `stot` for the days of
year needed, `tz`/`we`/`clt` for the days, and the hourly ERA5 shards that cover the days plus one
day after (the 06→06 window and the 01–24 UTC forecast hours reach into the next day).

**Public files change.** WSL maintains the public bucket and has replaced some files after their
first release (mostly on broken days, see [known-issues.md](known-issues.md)); some files are
missing for single days (e.g. `clt` on 15 days, the `tz` of some days). A missing input makes the
(variable, day) fail in `build`; it is reported and stays empty.

## 2. Native CHELSA format (GeoTIFF)

The published layout, used for the inputs and for `zarr-to-tif`:

```
<root>/daily/<var>/<year>/CHELSA_<var>_<dd>_<mm>_<yyyy>_V.2.1.tif
<root>/input/static/<name>.sdat
```

| Property | Value |
|---|---|
| Format | Cloud Optimized GeoTIFF, one band |
| Data type, NoData | uint16, 65535 |
| Value | raw × scale + offset; scale 0.1 (temperatures K, prec mm day⁻¹ = kg m⁻² day⁻¹, rsds W m⁻², ps hPa), 0.01 (sfcWind m s⁻¹), 0.001 (hurs, fraction although the tag says "percent"); offset 0, except V2.1 tas files of 1941–1950 (711 files: °C with offset −273.15) |
| Grid | 43200 × 20880, EPSG:4326, transform (0.0083333333, 0, −180.00013888885, 0, −0.0083333333, 83.99986041515) |
| Compression, tiles | DEFLATE with predictor 2, 512 × 512 |
| Overviews | 2, 4, …, 128 (written with nearest-neighbour resampling by `zarr-to-tif`) |
| Tags | CHELSA metadata: `variable`, `variable_long_name`, `variable_unit`, `cf_standard_name`, `datetime`, `source_id`, `version`, citations, … |

`zarr-to-tif` adds `day_window` to the tags and sets `source_id = CHELSA-UTC` for CHELSA-UTC stores.

## 3. Zarr format

One store per product (`store.py`):

```
<store>.zarr/                zarr v3
    time  (time,)            int64 "days since <first day>"; may have gaps
    y     (20880,)           pixel-centre latitudes, 84 N -> 90 S
    x     (43200,)           pixel-centre longitudes, 180 W -> 180 E
    <var> (time, y, x)       uint16, value = raw * scale_factor (+ offset), 65535 = missing
                             inner chunk 720 x 720, one shard (= one file) per (variable, day):
                             <store>.zarr/<var>/c/<t>/0/0
    <var>_add_offset (time,) only in stores made from native files with a non-zero offset on some
                             days; named in the variable attribute `add_offset_per_time`
<store>.qc/<timestamp>.parquet   per (variable, day) of each build: valid pixels, min, max, mean, seconds
<store>.flags.json               flagged (variable, day) pairs (`finalize`)
```

The integers and scales are those of the GeoTIFFs, so the conversion is lossless both ways.
Decode as `float32(raw) * float32(scale_factor) + float32(offset)` with 65535 as NaN
(`store.decode`; `formats.open_chelsa` does this for a region and days). Plain `xarray.open_zarr`
applies only the constant `add_offset` attribute, so use the decoder for stores with
`<var>_add_offset`.

Daily ERA5 stores (`era5_daily`): `time` (days), `latitude` 90 → −90 (721), `longitude` −180 →
179.75 (1440), float32 arrays `t2m_min`, `t2m_max`, `tp_mean`, `ssrd`, `cbh`, `t2m_mean`, `sp`,
`tcc`, `blh`, `u10`, `v10` (units in [method.md](method.md)); one chunk per (variable, day).

## 4. Conversion

| From → to | Command | Function |
|---|---|---|
| native tree → zarr | `chelsa-utc tif-to-zarr <dst.zarr> --src <root> --vars ... --first ... --last ...` | `formats.tif_to_zarr` |
| zarr → native tree | `chelsa-utc zarr-to-tif <dst root> --src <store.zarr> [--name ...]` | `formats.zarr_to_tif` |
| one file | | `formats.read_tif`, `formats.write_tif` |
| either → decoded values | | `formats.open_chelsa(source, variables, days, bbox)` |
| either → integers | | `formats.read_raw(source, var, day)` |

Both conversions skip what exists already (`--redo` / `--overwrite` to rewrite) and run one process
per file. Writing one global GeoTIFF takes about one minute.

## 5. Caching and restarts

- Hourly mirror: a shard is committed by its manifest file (`<stem>.manifest/<shard>.parquet`);
  `era5-hourly download` skips committed shards.
- Daily ERA5: `era5_daily.daily` reads a day from the daily store of the window if it is there,
  otherwise it aggregates the hourly mirror. `era5-daily` fills only missing days.
- 30″ stores: `build` skips (variable, day) pairs whose shard file exists. A failed item (missing
  input, killed process) gets up to three fresh process pools, then it is reported; run `build` again.
- Downloads: files that exist with the same size are not fetched again.
- More days: a store's time axis is fixed. `relink` makes a new store with more days and hard-links
  the shard files of all written (variable, day) pairs of the old one (same file system needed).
  zarr writes a chunk to a new file and renames it, so writing into the new store never changes
  the old one.
