import numpy as np
import pandas as pd
import pytest

from chelsa_utc import era5_daily, era5_hourly

DAY = pd.Timestamp("1976-07-15")
H0 = int((DAY - era5_hourly.T0) / pd.Timedelta(hours=1))


@pytest.fixture
def fake_hourly(monkeypatch):
    """Every hourly value = its hour index (since 1940-01-01), on a 2 x 1440 grid around DAY."""
    lo = H0 - 24
    hours = np.arange(lo, H0 + 72, dtype=np.float64)

    class Var:
        def __getitem__(self, s):
            h = hours[s.start - lo:s.stop - lo]
            return np.broadcast_to(h[:, None, None], (len(h), 2, 1440)).copy()

    monkeypatch.setattr(era5_daily, "_hourly", lambda: {v: Var() for v in era5_hourly.VARS})


def test_utc_window(fake_hourly):
    f = era5_daily.from_hourly(DAY, "utc", ["t2m_min", "t2m_max", "t2m_mean", "tp_mean"])
    assert f["t2m_min"][0, 0] == H0 + 1           # forecast fields: valid hours 01..24
    assert f["t2m_max"][0, 0] == H0 + 24
    assert f["t2m_mean"][0, 0] == H0 + 11.5       # analyses: 00..23
    assert np.isclose(f["tp_mean"][0, 0], 3600 * sum(range(H0 + 1, H0 + 25)), rtol=1e-6)


def test_chelsa0606_window(fake_hourly):
    f = era5_daily.from_hourly(DAY, "chelsa0606", ["t2m_min", "t2m_max", "t2m_mean"])
    assert f["t2m_min"][0, 0] == H0 + 7           # 06 UTC d -> 06 UTC d+1
    assert f["t2m_max"][0, 0] == H0 + 30
    assert f["t2m_mean"][0, 0] == H0 + 11.5       # analyses unchanged


def test_daily_falls_back_to_hourly(fake_hourly, monkeypatch, tmp_path):
    monkeypatch.setattr(era5_daily, "store_path", lambda w: tmp_path / "none.zarr")
    f = era5_daily.daily(DAY, "utc", ["t2m_min"])
    assert f["t2m_min"][0, 0] == H0 + 1
