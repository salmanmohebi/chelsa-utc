"""Paths of all inputs, caches and outputs.

Settings come from a TOML file (format: ``config.example.toml``). Lookup order:
``configure(path)`` (CLI ``--config``) > ``$CHELSA_UTC_CONFIG`` > ``./config.toml`` > defaults.
``configure`` also sets ``$CHELSA_UTC_CONFIG``, so worker processes (spawn) see the same file.
"""
from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass
from pathlib import Path

WINDOWS = ("utc", "chelsa0606")

DEFAULTS = {
    "root": "data",
    "era5_hourly": "{root}/era5/era5-arco-hourly.zarr",
    "chelsa_v21": "{root}/chelsa-v2.1",
    "stot": "{root}/chelsa-v2.1/rsdscs",
    "era5_daily": {"utc": "{root}/era5/era5-utc-daily.zarr",
                   "chelsa0606": "{root}/era5/era5-chelsa0606-daily.zarr"},
    "chelsa": {"utc": "{root}/chelsa-utc.zarr", "chelsa0606": "{root}/chelsa-chelsa0606.zarr"},
}


@dataclass(frozen=True)
class Paths:
    root: Path
    era5_hourly: Path
    chelsa_v21: Path
    stot: Path
    era5_daily: dict[str, Path]
    chelsa: dict[str, Path]
    v21_era5_inputs: Path | None     # optional: the V2.1 daily ERA5 input files (mirror test)
    source: str                      # the config file, or "defaults"


_PATHS: Paths | None = None


def _resolve(value: str, root: Path | None, base: Path) -> Path:
    s = value.replace("{root}", str(root)) if root is not None else value
    p = Path(os.path.expandvars(s)).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def _load(path: Path | None) -> Paths:
    cfg: dict = {}
    base = Path.cwd()
    if path is not None:
        with open(path, "rb") as f:
            cfg = tomllib.load(f).get("paths", {})
        base = path.parent.resolve()
    merged = {k: (dict(v, **cfg.get(k, {})) if isinstance(v, dict) else cfg.get(k, v)) for k, v in DEFAULTS.items()}
    root = _resolve(merged["root"], None, base)
    r = lambda v: _resolve(v, root, base)
    return Paths(root=root, era5_hourly=r(merged["era5_hourly"]), chelsa_v21=r(merged["chelsa_v21"]),
                 stot=r(merged["stot"]),
                 era5_daily={w: r(merged["era5_daily"][w]) for w in WINDOWS},
                 chelsa={w: r(merged["chelsa"][w]) for w in WINDOWS},
                 v21_era5_inputs=r(cfg["v21_era5_inputs"]) if cfg.get("v21_era5_inputs") else None,
                 source=str(path) if path else "defaults")


def configure(path: str | os.PathLike | None = None) -> Paths:
    """Load the settings (see the module docstring) and make them the current ones."""
    global _PATHS
    if path is None:
        env = os.environ.get("CHELSA_UTC_CONFIG")
        path = env if env else ("config.toml" if Path("config.toml").exists() else None)
    p = Path(path).resolve() if path else None
    if p is not None:
        if not p.exists():
            raise FileNotFoundError(p)
        os.environ["CHELSA_UTC_CONFIG"] = str(p)
    _PATHS = _load(p)
    return _PATHS


def paths() -> Paths:
    """The current settings (loaded on first use)."""
    return _PATHS if _PATHS is not None else configure()
