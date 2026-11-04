# Known issues of CHELSA V2.1 daily (1941–1978) and how CHELSA-UTC treats them

Found while re-implementing the algorithm and re-running the published days
(B: formula, D: daily ERA5 input files, H: published helper maps). The day lists are also in
`src/chelsa_utc/flags.py`; the flag file of the released CHELSA-UTC v1.0 is
`src/chelsa_utc/data/flags-v1.0.json` (`chelsa-utc flags`).

| Issue | CHELSA-UTC |
|---|---|
| B1 ps formula | fixed |
| B2 one-column shift in `clt` | kept (published `clt` reused) |
| B3 no +180° column in `tz` | fixed |
| B4 two encodings of tas | fixed (one encoding) |
| D1–D7 broken daily ERA5 input days | repaired (hourly ERA5) |
| H1–H3 broken published helper maps | flagged |
| W1 06→06 cloud base inside `clt` | **not fixed** (published `clt` reused); affects rsds |

## B: formula issues

**B1: surface pressure, temperature unit.** `chelsa_main.py:696–712` calls SAGA `Air_Pressure` with
the published lapse rate and the daily mean 2 m temperature. The published ps matches (99th
percentile 0.05 hPa)

    ps_V2.1 = B(sp) · (1 − B(L) · Δz / (273.15 + B(T)))^5.255,   L = −tz,   Δz = DEM − B(orog)

with *T* in kelvin. SAGA expects °C and adds 273.15 itself (`air_pressure.cpp:74, 137`), so the
denominator is about 560 K instead of 288 K and the height correction is roughly halved. Against ISD
station pressure (1976, 976 stations) the error grows by +56.4 hPa per 1000 m of Δz: median
−62.5 hPa in deep valleys, +62 to +67 hPa on summits; |ps_V2.1 − ps_physical| > 10 hPa on
4.4–6.0 % of land and on 42 % of the Alps. CHELSA-UTC uses the standard formula
([method.md](method.md)); over 1941–1978 the station bias goes from −4.1 to 0.0 hPa and no longer
depends on height.

**B2: one-column shift in the cloud field.** In the `clt` chain (`chelsa_main.py:529–539`) the
`change_latlong360` round trips shift the 0.025° cloud field and the 0.25° tcc bias field one column
to the east each. A port of the `clt` chain reproduces the published `clt` only with this shift.
Effect on rsds: RMS 2–3 W m⁻², maximum 113 W m⁻². CHELSA-UTC reuses the published `clt` and keeps it.

**B3: no +180° column in the lapse rate.** The published `tz` has 1440 columns, so tasmin, tasmax and
tas are NoData in the last 15 columns of the 30″ grid (0.0096 % of land). CHELSA-UTC adds the column.

**B4: two encodings for tas.** 711 tas files of 1941–1950 store °C with offset −273.15; all other
days store K. CHELSA-UTC stores all variables as uint16 × scale without offset.

## W1: the 06→06 cloud base inside `clt`, and so inside CHELSA-UTC rsds (not fixed)

V2.1 builds the daily cloud cover `clt` from 00–23 UTC analyses (tcc, pressure-level cloud cover,
geopotential) and from the ERA5 cloud base height `cbh`, which it averages over the forecast window
06 UTC *d* → 06 UTC *d+1* (the V2.1 input `cbase.nc` equals the NaN-mean over valid hours 07–30 to
≤ 0.13 m). `cbh` enters the Cloud Overlap step (the wind correction below the cloud base).
CHELSA-UTC reuses the published `clt`, so its rsds still holds this piece of the 06→06 day. This is
the only 06→06 element left in CHELSA-UTC.

Size, measured with a numpy port of the `clt` chain (1976-07-15, whole globe; the port reproduces
the published `clt` on 99.75 % of pixels, 99.96 % within one step): replacing the 06→06 cloud base
with the 00–24 UTC one changes `clt` by RMS 0.002 (max 0.12) and rsds by mean 0, RMS 0.33 W m⁻²,
p99 0.8 W m⁻², max 67 W m⁻²; 0.37 % of pixels change by ≥ 2 W m⁻² (largest at high northern
latitudes in summer). Earlier box tests (Alps, East Asia, Rockies, 12 days 1950–1977): zero mean,
RMS 0.5–1.5 W m⁻². This is small against the change of rsds from the window of ssrd itself
(19 W m⁻² RMS). A fix needs the `clt` chain re-implemented and run for every day; it is not part
of this repository.

## D: daily ERA5 input days

Found by comparing the V2.1 daily input files with hourly ERA5 and by re-running the published
days. CHELSA-UTC takes every ERA5 value from hourly ERA5, so all are repaired.

| | Variable | Days | Problem |
|---|---|---|---|
| D1 | tasmin | 1943-06-17, 1949-02-27, 1949-08-28, 1950-01-24, 1956-07-31, 1957-03-16, 1957-04-13, 1973-11-04 | `mn2t` is the daily mean of hourly `mn2t` over 00–23 UTC, not the minimum: tasmin 3.0–3.7 K too warm |
| D2 | prec | 1945-08-04, 1945-08-15, 1946-04-05, 1949-03-01, 1957-03-17, 1957-05-02, 1961-01-11, 1965-02-09, 1972-08-15 | precipitation already on the 00–24 UTC day, unlike all other days |
| D3 | tasmin | 1971-10-18, 1972-10-15, 1974-04-30, 1974-08-16 | the published value does not follow from the input file (probably replaced after the run) |
| D4 | rsds | 1972-09-28 | published from another ssrd file (implied ssrd = 0.088 × the real value) |
| D5 | tas, tasmin, tasmax, ps | 1966-03-10 | broken published output (values 0 and 6553.4) |
| D6 | prec | 33 days: 32 in 1950–1955, and 1964-12-25 | the wind factor does not follow from the `we` file (`we` probably regenerated after prec); up to 8.5 mm |
| D7 | several | about 70 (day, variable) pairs | broken input files: `tmax.nc` all zero on 9 days, `ssrd.nc` all zero on 8 days, `tmean.nc` on 8 days, `sp` all zero on several days; published ps a constant field on 13 days, rsds on 9 days |

Note on 1958-05-31: `tmean.nc`, `blh.nc`, `u.nc` and `v.nc` of this day hold 24 hourly fields
instead of one daily mean. The published V2.1 agrees with the mean of the 24 fields, so V2.1 is
correct on this day; only a reader that takes the first field gets a wrong value.

## H: published helper maps

CHELSA-UTC reuses the published `tz`, `we` and `clt`; on these days the output cannot be repaired
and is flagged.

- **H1:** `tz` of 1966-03-10 is about 0.3 K m⁻¹ everywhere (normal ≈ −0.006), so tas, tasmin and
  tasmax reach 0–1459 K.
- **H2:** `clt` is empty (all NaN) on 10 days: 1946-10-06, 1946-10-09, 1946-12-06, 1949-10-25,
  1954-11-16, 1954-11-18, 1954-11-19, 1961-12-06, 1962-03-18, 1972-04-22. CHELSA-UTC rsds of these
  days is missing; the published rsds is a constant field.
- **H3:** `we` built from a broken 10 m wind input on 21 days, and `clt` from a broken tcc input on
  11 days (e.g. `clt` of 1965-07-15 has a global mean of 0.243 against ≈ 0.44 on the days around it;
  `we` of 1963-05-29 is the constant 1.000).

## Flags in the released CHELSA-UTC v1.0

13 864 days. 54 (variable, day) records on 52 dates: prec 19, rsds 26, tas / tasmin / tasmax 1 each
(H1–H3), and the QC flags of the unchanged published hurs (2) and sfcWind (4).

## The public files

CHELSA-UTC v1.0 uses exactly the CHELSA V2.1 files that are on the public bucket on 2026-11-04.
The store was first built from a copy of the published files made on 2026-06-04; on 2026-11-04 the
copy was compared with the bucket file by file (1941–1978; `tz`, `we`, `clt`, hurs, sfcWind, the
outputs, the static grids and stot; 153 246 files):

- Files on both sides with the same size are identical: MD5 against the bucket ETag for all 3 605
  `we` files that WSL re-uploaded after the copy, and for 25 random files of each group.
- 9 sfcWind files were replaced by WSL after the copy (broken files of 7 MB or 1.4 GB, now
  357–440 MB): 1941-01-24, 1954-11-28, 1955-09-22, 1962-09-21, 1965-09-22, 1966-03-27, 1968-03-29,
  1968-04-18, 1969-06-27. The new files were downloaded (ETag checked) and copied into the store;
  these 9 days are no longer flagged.
- WSL also replaced 8 tas and 9 tasmax files (broken days of V2.1) and added files for missing
  days. CHELSA-UTC does not use the published outputs (it recomputes them), so they do not change it.
- 15 days stay missing (13 864 of 13 879 days): 1946-05-20, 1949-07-13, 1949-09-21, 1951-07-05,
  1952-11-09, 1956-01-28, 1956-03-12, 1960-08-15, 1963-04-15, 1963-05-29, 1968-06-12, 1973-06-06,
  1973-06-12, 1973-07-08, 1974-08-30. `clt` is not public for any of them, and `tz` only for 10
  (all except 1946-05-20, 1949-07-13, 1949-09-21, 1952-11-09, 1973-07-08). When WSL publishes `clt`
  for them, `relink` + `build` can add them.
- The 19 days that the first build had left out only because a published output file was missing in
  the copy (18 days in January–March 1941 and 1961-12-07) were added on 2026-11-04; none of them is
  flagged.

WSL still updates the bucket (latest change in 1941–1978: 2026-10-25). A comparison like the one
above shows whether a later change affects CHELSA-UTC: only `tz`, `we`, `clt`, the static grids,
stot, hurs and sfcWind matter.
