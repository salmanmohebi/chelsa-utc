"""Known problems of the published CHELSA V2.1 inputs that CHELSA-UTC cannot repair, and the flag file.

CHELSA-UTC takes every ERA5 value from hourly ERA5, so broken daily ERA5 input files of V2.1 are
repaired. It reuses the published helper maps tz, we and clt, so days where one of those is broken
(or was built from a broken input) stay unreliable and are flagged. Details: ``docs/known-issues.md``.

The flag file of the released store (``data/flags-v1.0.json``, 54 records on 52 dates) also holds
the QC flags of the unchanged published hurs and sfcWind.
"""
from __future__ import annotations

import json
from importlib import resources

# Published `we` (10 m wind effect) built from a broken 10 m wind input (H3) -> prec.
BAD_WIND = ["1946-04-26", "1946-06-04", "1950-03-03", "1951-05-31", "1951-09-27", "1951-09-28",
            "1951-11-23", "1955-07-23", "1956-05-11", "1956-08-10", "1956-08-17", "1957-03-02",
            "1957-09-11", "1957-10-21", "1958-05-31", "1960-05-17", "1962-10-22", "1963-04-15",
            "1963-05-28", "1963-05-29", "1966-11-19"]
# Published `clt` built from a broken total-cloud-cover input (H3) -> rsds.
BAD_TCC = ["1947-01-29", "1948-07-23", "1949-02-21", "1951-07-17", "1955-10-26", "1957-02-07",
           "1962-08-15", "1965-07-15", "1973-12-12", "1973-12-26", "1974-03-31"]
# Published `tz` of this day is ~0.3 K/m everywhere (normal ~ -0.006) (H1) -> tas/tasmin/tasmax 0..1459 K.
BAD_TZ = ["1966-03-10"]
# Published `clt` is empty (all NaN) (H2) -> no rsds (V2.1 rsds of these days is a constant field).
BAD_CLT = ["1946-10-06", "1946-10-09", "1946-12-06", "1949-10-25", "1954-11-16", "1954-11-18",
           "1954-11-19", "1961-12-06", "1962-03-18", "1972-04-22"]

KNOWN_ISSUES: dict[tuple[str, str], str] = {
    **{(d, "prec"): "published `we` (wind effect) was built from a broken 10 m wind input" for d in BAD_WIND},
    **{(d, "rsds"): "published `clt` was built from a broken total-cloud-cover input" for d in BAD_TCC},
    **{(d, v): "published `tz` (lapse rate) of this day is broken (~0.3 K/m)" for d in BAD_TZ
       for v in ("tas", "tasmin", "tasmax")},
    **{(d, "rsds"): "published `clt` of this day is broken (no rsds; V2.1 rsds is a constant field)"
       for d in BAD_CLT},
}

# Physical range per variable; a day with more than MAX_OUT_FRAC of its valid pixels outside is flagged.
PIXEL_RANGE: dict[str, tuple[float, float]] = {
    "tas": (150.0, 350.0), "tasmax": (150.0, 350.0), "tasmin": (150.0, 350.0),
    "prec": (0.0, 2000.0), "rsds": (0.0, 600.0), "ps": (250.0, 1100.0),
    "sfcWind": (0.0, 200.0), "hurs": (0.0, 1.05),
}
MAX_OUT_FRAC = 1e-4

# Days where CHELSA-UTC knowingly differs from published V2.1 beyond the window change.
NOTES = {
    "fixed": {
        "tasmin": "8 days whose V2.1 input mn2t was a daily MEAN (tasmin 3.0-3.7 K too warm): "
                  "1943-06-17 1949-02-27 1949-08-28 1950-01-24 1956-07-31 1957-03-16 1957-04-13 1973-11-04; "
                  "4 days whose published value does not reproduce from its input file: "
                  "1971-10-18 1972-10-15 1974-04-30 1974-08-16",
        "tasmax": "9 days with an all-zero input tmax.nc",
        "tas": "days with broken tmean inputs; one encoding for all days (V2.1 1941-1950 tifs have offset -273.15)",
        "ps": "formula bug (T in K in SAGA's Celsius slot) on every day; days with all-zero sp inputs; "
              "1966-03-10 (does not use tz)",
        "rsds": "8 days with an all-zero input ssrd; 1972-09-28 (published from another ssrd file)",
        "prec": "9 days whose V2.1 input was already on the UTC day: 1945-08-04 1945-08-15 "
                "1946-04-05 1949-03-01 1957-03-17 1957-05-02 1961-01-11 1965-02-09 1972-08-15",
        "temperature": "east of 179.875 E (V2.1: NoData, tz without the +180 column)",
    },
    "not_fixed_flagged": {
        "tas, tasmin, tasmax": "1966-03-10: the published lapse rate tz of that day is broken",
        "rsds": "10 days with a broken published clt (all-NaN): " + " ".join(BAD_CLT),
        "prec / rsds": "days whose published we / clt were built from broken wind / cloud inputs",
    },
    "differs_from_v21": {
        "prec": "33 days where the published `we` file does not reproduce the published prec (probably "
                "regenerated after the V2.1 run): 32 in 1950-1955 and 1964-12-25; CHELSA-UTC is "
                "consistent with the published `we`",
    },
}


def released_flags() -> dict:
    """The flag file of the released CHELSA-UTC v1.0 store."""
    return json.loads(resources.files("chelsa_utc").joinpath("data/flags-v1.0.json").read_text())


def flag(record: dict) -> str | None:
    """The reason to flag one QC record ({variable, date, n_valid, n_out, ...}), or None."""
    k = (record["date"], record["variable"])
    if k in KNOWN_ISSUES:
        return f"known source issue: {KNOWN_ISSUES[k]}"
    if not record.get("n_valid"):
        return "all-NaN"
    if record.get("n_out", 0) > MAX_OUT_FRAC * record["n_valid"]:
        lo, hi = PIXEL_RANGE[record["variable"]]
        return f"{record['n_out']} pixels outside [{lo:g}, {hi:g}]"
    return None
