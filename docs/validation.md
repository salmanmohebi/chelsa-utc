# Validation

Four questions, each answered by a test of the code in this repository:

1. Is the daily ERA5 aggregation correct?
2. Does the 30″ code reproduce the published CHELSA V2.1?
3. Is the released CHELSA-UTC store the output of this code?
4. What changes against V2.1, and is the new day right?

Integers are compared as stored (1 step = 0.1 K, 0.1 mm day⁻¹, 0.1 W m⁻², 0.1 hPa).

## 1. Daily ERA5

- **Hourly mirror = ARCO-ERA5.** 472/472 shards complete, no missing hour; 220 random hours
  bit-identical to the source (`era5-hourly check`). NaN only in `cbh` (no cloud, expected) and in
  the forecast fields of 1940-01-01 00–06 UTC (ERA5 forecasts start with the 06 UTC run; not used).
- **06→06 window = the V2.1 daily input files** (6 days, `chelsa0606`). Largest differences:
  t2m_min/t2m_max 0.0010 K, t2m_mean 0.0005 K, tp 0.0077 mm, ssrd 0.0052 W m⁻², sp 0.0025 hPa,
  cbh 0.13 m, blh 0.028 m, u/v 0.0003 m s⁻¹. This is rounding: the V2.1 files reached CHELSA by
  another route than ARCO.
- **Regression over all days.** The 06→06 window rebuilt from hourly ERA5 agrees with the V2.1
  daily inputs on 99.9 % of the 13 879 days to 0.001 K / 0.011 mm / 0.007 W m⁻²; every day over
  the tolerance is a known broken input day (D1, D2, D7 in [known-issues.md](known-issues.md)).
- **UTC window = the released ERA5-UTC daily store**: bit-identical on 43 days (one per year
  1941–1978 and 5 special days) × 11 variables. So the released daily store is used as a cache.

## 2. The 30″ code reproduces V2.1

**Mirror test.** The code on the original V2.1 daily ERA5 input files (`--window v21files`), every
61st day 1941–1978 (227 days), against the published files:

| Variable | Identical, median day | Identical, worst day | Worst day | Max. difference (steps) |
|---|---|---|---|---|
| tasmin | 99.9928 % | 65.32 % | 1974-04-30 | 35 |
| tasmax | 99.9926 % | 99.99 % | 1946-05-26 | 1 |
| tas | 99.9927 % | 99.99 % | 1965-06-19 | 5.5 |
| prec | 99.9999 % | 93.72 % | 1955-08-06 | 422 |
| rsds | 99.9943 % | 99.99 % | 1943-11-22 | 1 |
| ps (V2.1 formula) | 99.9614 % | 99.96 % | 1954-06-05 | 1 |

On a typical day more than 99.99 % of the pixels are identical and the rest differ by one step
(float vs SAGA arithmetic). The worst days are known input problems, not code errors: tasmin
1974-04-30 (D3), prec 1955-08-06 (D6), and the tas files of 1941–1950 that store °C (B4).

**V2.1 rebuilt from hourly ERA5** (`--window chelsa0606`, V2.1 formulas) against the published files,
% identical integers (largest difference in steps):

| Variable | 1941-01-01 | 1958-05-31 | 1959-07-24 | 1966-03-10 | 1976-07-15 | 1978-12-31 |
|---|---|---|---|---|---|---|
| tasmin | 99.72 (1) | 99.64 (1) | 99.63 (1) | 99.67 (H1) | 99.65 (1) | 99.68 (1) |
| tasmax | 99.72 (1) | 99.64 (1) | 99.63 (1) | 99.67 (H1) | 99.65 (1) | 99.68 (1) |
| tas | 99.94 (1) | 99.93 (1) | 99.93 (1) | 99.90 (H1) | 99.93 (1) | 99.94 (1) |
| prec | 99.71 (1) | 99.66 (1) | 99.73 (1) | 99.74 (1) | 99.56 (1) | 99.49 (1) |
| rsds | 99.69 (1) | 99.69 (1) | 99.75 (1) | 99.66 (1) | 99.74 (1) | 99.70 (1) |
| ps | 99.80 (1) | 99.78 (1) | 99.79 (1) | 99.70 (H1) | 99.81 (1) | 99.80 (1) |

The ±1 steps come from the ERA5 rounding of point 1. 1966-03-10 is the day with the broken
published `tz` (H1). prec has 12 NoData-mismatch pixels on every day (grid edge, constant).

Reproduce:

```bash
chelsa-utc compare --window chelsa0606 --reference v21 --dates 1976-07-15 --vars tasmin,prec
chelsa-utc compare --window v21files  --reference v21 --dates 1976-07-15      # needs the V2.1 input files
```

## 3. The released CHELSA-UTC store is the output of this code

One random day per year (1941–1978) and 5 special days (1941-01-01, 1958-05-31, 1966-03-10,
1976-07-15, 1978-12-31), all 6 computed variables, recomputed with this code: **258 of 258 maps
identical** to the stored integers. The released store is therefore used as a cache
(`build` skips what is written). Reproduce with
`chelsa-utc compare --window utc --reference <released store> --dates ...`.

The copied hurs and sfcWind are byte-identical to the published files of the copy they come from
(file counts, sizes and 20 random files per variable checked).

**Re-checked with this repository** (2026-11-04, 1976-07-15, the code exactly as committed):

| Check | Result |
|---|---|
| daily ERA5 from the hourly mirror vs the released ERA5-UTC daily store, 11 variables | bit-identical (cbh: identical with NaN = NaN, 33 206 NaN cells in both) |
| CHELSA-UTC tasmin, ps, rsds, prec vs the released store | 100 % identical integers, 0 NoData mismatch |
| V2.1 rebuilt (`chelsa0606`) vs published: tasmin / ps | 99.65 % / 99.81 % identical, max 1 step, 0 NoData mismatch |
| store → GeoTIFF → store (tasmin, hurs) | identical integers both ways |
| hurs in the store vs the published file | identical |
| exported GeoTIFF vs a published file | same layout: COG, uint16, nodata 65535, DEFLATE + predictor 2, 512 × 512 tiles, overviews 2–128, same CRS, transform and scale (file 179 MB vs 161 MB: other DEFLATE level) |
| `open_chelsa` on the store and on the exported tree (Alps box) | identical values |

Time per (variable, day) on one core of a recent CPU: temperatures and ps 1–1.5 min, rsds 1.5 min,
prec 3–3.5 min.

## 4. What changes, and is the new day right?

**Change against V2.1**, all pixels and days 1941–1978 (*d* = CHELSA-UTC − V2.1, cos-latitude
weighted, flagged days left out):

| Variable | Unit | Day window | Mean *d* | RMS *d* | 99 % of \|*d*\| below | Share \|*d*\| > threshold |
|---|---|---|---|---|---|---|
| tasmin | K | moved | −0.026 | 0.79 | 3.9 | 11.6 % > 1 K |
| tasmax | K | moved | −0.015 | 0.73 | 3.9 | 9.8 % > 1 K |
| tas | K | same | −0.000 | 0.01 | 0.1 | 0.0 % > 0.5 K |
| prec | mm day⁻¹ | moved | +0.000 | 2.36 | 9.0 | 2.8 % > 5 mm |
| rsds | W m⁻² | moved | +0.009 | 19.15 | 82.0 | 11.7 % > 20 W m⁻² |
| ps | hPa | same (formula) | +0.043 | 5.35 | 28.5 | 9.1 % > 5 hPa |

The long-term means do not change (prec: median change of the 1941–1978 mean per 0.5° cell
0.03 %); single days do. tasmin/tasmax change with longitude (whether a 24-hour window holds one
or two nights or afternoons depends on local time); ERA5 alone shows the same pattern, so it comes
from the window, not from the downscaling. ps changes only along steep terrain (B1).

**Physical consistency.** Pixels with tasmin > tas or tas > tasmax are 3–5 times fewer than in V2.1
in every month; the rest comes from ERA5 itself (extremes are forecasts, tas is an analysis).

**Against stations.** Stations with sub-daily data, so that the same station-day can be aggregated
to both windows: ISD hourly or 3-hourly 2 m temperature and station pressure, and HPD hourly
precipitation (US), 1941–1978, flagged days left out (14.98 M station-days at 8 214 stations for
temperature, 11.42 M at 1 499 stations for precipitation). Metric: RMSE after removing each
station's mean error (the point-vs-pixel offset).

| Variable | Product | V2.1 vs UTC truth | UTC vs UTC truth | V2.1 vs 06→06 truth | UTC vs 06→06 truth | Own-window ratio |
|---|---|---|---|---|---|---|
| tasmin | CHELSA | 2.36 | 2.15 | 2.16 | 2.38 | 0.995 |
| tasmin | ERA5 | 2.37 | 2.16 | 2.17 | 2.39 | 0.995 |
| tasmax | CHELSA | 2.14 | 1.91 | 1.90 | 2.23 | 1.002 |
| tasmax | ERA5 | 2.17 | 1.94 | 1.94 | 2.26 | 1.001 |
| prec | CHELSA | 5.78 | 5.41 | 5.40 | 5.92 | 1.002 |
| prec | ERA5 | 5.90 | 5.53 | 5.52 | 6.04 | 1.002 |

(own-window ratio = UTC vs UTC truth / V2.1 vs 06→06 truth; K or mm day⁻¹.)

- Each product is best on its own window, with the same skill (ratio 0.995–1.002): CHELSA-UTC is not
  more accurate than V2.1, it is a consistent day.
- A window mismatch costs 9–11 % (temperature) and 6 % (precipitation) of the error. ERA5 alone gives
  the same numbers, so this measures the day definition, not the downscaling.
- ps: station bias −4.1 → 0.0 hPa (at 2500–6000 m: −13.4 → +0.3 hPa); day-to-day error 2.63 → 1.15 hPa.

The station check is not part of this repository (it needs the ISD and HPD archives).
