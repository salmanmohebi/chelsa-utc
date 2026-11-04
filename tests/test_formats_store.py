from pathlib import Path

import numpy as np
import pandas as pd

from chelsa_utc import config, flags, formats, store

ROOT = Path(__file__).resolve().parents[1]


def test_store_layout_and_time_axis_with_gaps(tmp_path):
    days = pd.DatetimeIndex(["1941-01-01", "1941-01-02", "1941-01-05"])
    p = tmp_path / "s.zarr"
    store.create(p, days, ["tas", "prec"], {"title": "test"}, per_day_offset=["tas"])
    g = store.open_store(p)
    assert list(store.dates(g)) == list(days)
    assert store.index(g)[pd.Timestamp("1941-01-05")] == 2
    assert g["tas"].shape == (3, 20880, 43200) and g["tas"].dtype == np.uint16
    assert g["tas"].attrs["add_offset_per_time"] == "tas_add_offset"
    assert not store.written(p, "tas", 0)
    assert np.isclose(g["y"][0], 83.99569444, atol=1e-6) and np.isclose(g["x"][0], -179.99597222, atol=1e-6)


def test_decode():
    raw = np.array([2731, 0, 65535], dtype=np.uint16)
    v = store.decode(raw, 0.1)
    assert np.isclose(v[0], 273.1) and v[1] == 0 and np.isnan(v[2])
    assert np.isclose(store.decode(raw[:1], 0.1, -273.15)[0], -0.05, atol=1e-4)


def test_names_and_tags():
    p = formats.tif_path("/x", "tasmin", "1976-07-15")
    assert p == Path("/x/daily/tasmin/1976/CHELSA_tasmin_15_07_1976_V.2.1.tif")
    t = formats.tags_for("prec", "1976-07-15", {"day_window": "utc"})
    assert t["datetime"] == "1976-07-15" and t["variable_unit"] == "kg m-2 day-1" and t["day_window"] == "utc"


def test_bbox_rows_cols():
    # every pixel that touches the box; the grid is 0.5" off the full degrees, so one extra row / column
    rows, cols = formats._bbox_rows_cols((5.0, 45.0, 11.0, 48.0))
    assert rows.stop - rows.start in (360, 361) and cols.stop - cols.start in (720, 721)


def test_config_example(tmp_path):
    c = config._load(ROOT / "config.example.toml")
    assert c.root == ROOT / "data"
    assert c.chelsa["utc"] == ROOT / "data" / "chelsa-utc.zarr"
    assert c.era5_daily["chelsa0606"].name == "era5-chelsa0606-daily.zarr"


def test_released_flags_match_known_issues():
    rel = flags.released_flags()
    keys = {(r["date"], r["variable"]) for r in rel["records"]}
    assert len(rel["records"]) == 54
    # every known issue on a day of the released store is flagged there (2 BAD_WIND days are not on its time axis)
    missing = {k for k in flags.KNOWN_ISSUES if (k[0], k[1]) not in keys}
    assert missing <= {("1963-04-15", "prec"), ("1963-05-29", "prec")}


def test_flag_rule():
    assert flags.flag({"date": "1966-03-10", "variable": "tas", "n_valid": 10, "n_out": 0}).startswith("known")
    assert flags.flag({"date": "1976-07-15", "variable": "tas", "n_valid": 0, "n_out": 0}) == "all-NaN"
    assert flags.flag({"date": "1976-07-15", "variable": "tas", "n_valid": 10 ** 6, "n_out": 0}) is None


def test_relink_adds_days_and_links_shards(tmp_path):
    import os
    days = pd.DatetimeIndex(["1941-01-01", "1941-01-03"])
    src = tmp_path / "a.zarr"
    store.create(src, days, ["tas"], {"title": "t"}, per_day_offset=["tas"])
    shard = src / "tas" / "c" / "1" / "0" / "0"           # pretend day 1 is written
    shard.parent.mkdir(parents=True); shard.write_bytes(b"x")
    store.open_store(src, "r+")["tas_add_offset"][1] = -273.15
    dst = tmp_path / "b.zarr"
    store.relink(src, dst, list(days) + [pd.Timestamp("1941-01-02")])
    g = store.open_store(dst)
    assert list(store.dates(g)) == list(pd.date_range("1941-01-01", "1941-01-03"))
    linked = dst / "tas" / "c" / "2" / "0" / "0"
    assert linked.exists() and os.path.samefile(linked, shard)
    assert not store.written(dst, "tas", 1)
    assert g["tas_add_offset"][2] == -273.15 and g["tas"].shards == (1, 20880, 43200)
