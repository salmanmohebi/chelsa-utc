# Method

CHELSA V2.1 daily downscales ERA5 (0.25°) to 30″. Its last step combines three kinds of input:

- **daily ERA5 fields:** minimum, maximum and mean 2 m temperature, precipitation, solar radiation,
  surface pressure, boundary-layer height (BLH);
- **daily helper maps that CHELSA computes itself:** lapse rate `tz` (0.25°), wind effect `we`
  (0.025°), total cloud cover `clt` (0.025°);
- **static grids:** 30″ DEM (`dem_latlong`), ERA5 orography, `expocor`, `patch`, and the clear-sky
  radiation `stot` (one file per day of year).

CHELSA-UTC keeps the helper maps and the static grids, takes every ERA5 field from hourly ERA5 over
the 00–24 UTC day, and re-runs the last step with a numpy copy of the SAGA operators.

| CHELSA step | In CHELSA-UTC |
|---|---|
| 1. hourly → daily ERA5 | **replaced**: hourly ARCO-ERA5 aggregated over the chosen window (`era5_hourly`, `era5_daily`) |
| 2. daily helper maps `tz`, `we`, `clt` | **reused** from the published files (already on the 00–23 UTC day) |
| 3. final 30″ formulas | **re-implemented** (`saga_ops`, `chelsa`) |

## 1. Hourly ERA5 and the day windows

**Source.** Hourly ERA5 single-level fields from the public ARCO-ERA5 zarr store
(`https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3`),
the same 11 variables as the daily input files of CHELSA V2.1. `era5_hourly` keeps a local copy
for 1940-01-01 to 1979-12-31 (8.4 TB, zarr v3, one shard = 31 days); every value is bit-identical
to the source.

| V2.1 input file | ERA5 field | Type | Daily value | Unit |
|---|---|---|---|---|
| `tmin.nc` | `mn2t` | forecast | minimum | K |
| `tmax.nc` | `mx2t` | forecast | maximum | K |
| `ps.nc` | `mtpr` (precipitation rate) | forecast | sum × 3600 | mm day⁻¹ |
| `ssrd.nc` | `ssrd` | forecast | sum / 86400 | W m⁻² |
| `cbase.nc` | `cbh` | instantaneous | mean, NaN = no cloud | m |
| `tmean.nc` | `t2m` | analysis | mean | K |
| `ps0.nc` | `sp` | analysis | mean / 100 | hPa |
| `tcc.nc`, `blh.nc`, `u.nc`, `v.nc` | `tcc`, `blh`, `u10`, `v10` | analysis | mean | 1, m, m s⁻¹ |

**Hour convention.** A forecast value at valid hour *t* covers (*t* − 1 h, *t*]. Valid hours 01–06
of day *d* come from the 18 UTC run of *d* − 1 (steps 7–12), hours 07–18 from the 06 UTC run of *d*
(steps 1–12), hours 19–24 from the 18 UTC run of *d* (steps 1–6). Analyses are instantaneous at
hour *t*.

**Day windows** (`era5_daily.WINDOWS`, hours after 00 UTC of day *d*):

| Window | Forecast fields | Analyses | Use |
|---|---|---|---|
| `chelsa0606` | valid hours 07–30 (06 UTC *d* → 06 UTC *d+1*) | 00–23 | CHELSA V2.1 |
| `utc` | valid hours 01–24 (00–24 UTC *d*) | 00–23 (unchanged) | CHELSA-UTC |

**Cloud base height.** `cbh` is instantaneous, but V2.1 uses it like a forecast field. Of four
windows tested against the published `cbase.nc`, only the NaN-aware mean over valid hours 07–30
matches (≤ 0.13 m, identical NaN pattern). On the UTC day we use valid hours 01–24. (`cbh` is an
input of `clt`, which CHELSA-UTC reuses; it is in the daily stores for completeness.)

## 2. SAGA operators (`saga_ops.py`)

The last CHELSA step runs SAGA GIS tools. To get the same integers, the code copies what SAGA does,
not a textbook interpolation.

- **Grids.** The 30″ grid has 43200 × 20880 pixel centres *x*ᵢ = −179.9959722222 + *i*·Δ with SAGA's
  cell size Δ = 0.0083333333 (not exactly 1/120). The 0.25° ERA5 grid is used as 721 × 1441 nodes
  from −90° and −180°, where column 1441 is a copy of column 1 (SAGA `change_latlong` with `PATCH`).
- **Sampling** (SAGA `CSG_Grid::Get_Value`). The target must lie inside the cell-edge extent and its
  nearest node must be valid, otherwise NoData. The value uses a 4 × 4 stencil; nodes outside the
  grid or NoData are filled by SAGA's iterative 3 × 3 neighbour mean. There is no longitude wrap.
- **Kernels.** "Bicubic Spline" (Catmull–Rom, *a* = −0.5) for temperature, written K(·) below;
  "B-Spline" (cubic B-spline without prefilter) for all other fields, written B(·).
- **Same-grid sampling.** SAGA samples with the B-spline even when the source and the target grid
  are the same, which is a (1/6, 2/3, 1/6) smoothing; the code applies it where V2.1 does.
- **Block statistics.** Block maximum and area mean per 0.25° cell use SAGA's nearest-node assignment
  and the true pixel geometry: the CHELSA grid is 0.5″ off the ERA5 cell edges, so a plain 30 × 30
  mean gives only 87–99.9 % identical rsds.
- **Storage.** round(*v* / 0.1) as uint16, 65535 = NoData.

The fast path is a separable gather; pixels whose stencil leaves the grid or touches NoData are
recomputed with the exact per-stencil fill. The code works on bands of 720 tif rows, with margins
large enough for each operator.

## 3. The formulas (`chelsa.py`)

Line numbers refer to `chelsa_main.py` of the public CHELSA V2.1 code.

**Temperature** (lines 742–881). For *T* = daily minimum, maximum or mean 2 m temperature:

    T_30" = K(T) + K(tz) · (DEM − K(orog))

**Precipitation** (lines 419–483). With *P* = max(*P*_ERA5, 0), the daily sum on the 0.25° grid,
and S(·) the same-grid B-spline:

1. boundary-layer top above sea level *h* = B(BLH + orog); *d* = |DEM − *h*| on the 30″ grid;
   *m* = block maximum of *d* per 0.25° cell;
2. *w* = S(B(we)) / (1 − |B(*m*) − *d*| / 9000), NoData filled from `patch`;
3. *c* = expocor · S(*w*), and *c̄* = area mean of *c* per 0.25° cell;
4. prec = *c* · B(*P*) / B(*c̄*).

**Solar radiation** (lines 570–627):

    s    = stot_doy · 0.01157407407 · (1 − 0.75 · B(clt)^3.4)
    rsds = s / B((s̄ + 0.001) / (ssrd + 0.001))

with *s̄* the area mean of *s* per 0.25° cell and ssrd the daily ERA5 radiation in W m⁻².

**Surface pressure** (lines 696–712). V2.1 (`v21=True`):

    ps_V2.1 = B(sp) · (1 − B(L) · Δz / (273.15 + B(T)))^5.255,   L = −tz,   Δz = DEM − B(orog)

with *T* in kelvin. SAGA's `Air_Pressure` expects *T* in °C and adds 273.15 itself, so the
denominator is about 560 K instead of 288 K and the height correction is roughly halved
(issue B1 in [known-issues.md](known-issues.md)). CHELSA-UTC (`v21=False`) uses the standard
barometric formula with the standard lapse rate:

    ps_UTC = B(sp) · (1 − 0.0065 · Δz / B(T))^5.2559,   T in K

This is SAGA's `Air_Pressure` with *T* in °C and *L* = 0.0065 K m⁻¹; the exponent differs from
5.255 by at most 0.05 hPa for |Δz| ≤ 3 km, below the storage step.

**The +180° column of `tz`** (`v21=False`). Every 0.25° ERA5 grid gets a copy of its first column at
+180° before interpolation. The published `tz` has 1440 columns without it, so in V2.1 tasmin,
tasmax and tas are NoData in the last 15 columns of the 30″ grid (0.0096 % of land: far eastern
Chukotka, Wrangel Island, some Fiji islands). CHELSA-UTC adds the column to `tz` like to every other
grid.

## 4. CHELSA-UTC in one table

| Variable | V2.1 | CHELSA-UTC |
|---|---|---|
| tasmin, tasmax | mn2t / mx2t over 06 UTC *d* → 06 UTC *d+1* | over 00–24 UTC *d*; `tz` with the +180° column |
| prec | precipitation over 06 UTC *d* → 06 UTC *d+1* | over 00–24 UTC *d*; published `we` |
| rsds | ssrd over 06 UTC *d* → 06 UTC *d+1* | over 00–24 UTC *d*; published `clt` |
| tas | 00–23 UTC mean | same window, recomputed (one encoding) |
| ps | formula with *T* in K in the °C slot | standard formula |
| hurs, sfcWind | 00–23 UTC | the published files, unchanged |
| all | daily ERA5 input files (with broken days) | hourly ERA5, aggregated (broken input days repaired) |

**Known limits.**
- rsds reuses the published `clt`, which keeps the 06→06 cloud base inside it (issue W1: zero mean,
  RMS 0.33 W m⁻² globally, max 67 W m⁻²; not fixed) and V2.1's one-column shift of the cloud fields
  (issue B2, RMS 2–3 W m⁻²).
- prec and rsds inherit the published `we` / `clt` on days where those were built from broken
  inputs (flagged).
- North of 83.75° N (ocean only), SAGA's NoData fill is approximated.

**Why the 00–24 UTC day.** The ERA5 inputs of a downscaling model, station observations converted
to UTC, and the target should describe the same 24 hours. In V2.1, tas (00–24 UTC) and
tasmin/tasmax (06→06) describe different hours, so tas can lie outside its own extremes; in
CHELSA-UTC such pixels are 3–5 times fewer. The cost of the 00–24 UTC day: where the coldest hour
falls near 00 UTC (e.g. 60–90° E), the same night can enter two consecutive days, and rain gauges
read at 06 UTC (much of Europe) match the 06→06 window better.
