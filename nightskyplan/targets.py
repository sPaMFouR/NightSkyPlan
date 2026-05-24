from __future__ import annotations

from collections.abc import Iterable

import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from nightskyplan.models import Target


EXAMPLE_TARGETS = pd.DataFrame(
    {
        "name": ["SN2022jli", "SN2020tlf", "SN2018zd"],
        "ra": ["00:34:45.2", "14:40:10.03", "06:18:03.18"],
        "dec": ["-08:23:25.0", "+42:46:39.6", "+78:22:00.9"],
        "exposure_min": [45, 30, 20],
        "priority": [3, 2, 1],
    }
)


def parse_angle(value: object, is_ra: bool) -> u.Quantity:
    """Parse decimal degrees or sexagesimal RA/Dec."""
    text = str(value).strip()
    if not text:
        raise ValueError("Coordinate value cannot be blank")
    if ":" in text or " " in text:
        if is_ra:
            return SkyCoord(text, "0d", unit=(u.hourangle, u.deg)).ra
        return SkyCoord("0h", text, unit=(u.hourangle, u.deg)).dec
    return float(text.removesuffix("deg").strip()) * u.deg


def parse_target_rows(df: pd.DataFrame) -> list[Target]:
    normalized = df.rename(columns={c: str(c).strip().lower() for c in df.columns})
    required = {"name", "ra", "dec"}
    missing = required - set(normalized.columns)
    if missing:
        raise ValueError(f"Missing target columns: {sorted(missing)}")

    if "exposure_min" not in normalized:
        normalized["exposure_min"] = 30.0
    if "priority" not in normalized:
        normalized["priority"] = 1.0

    targets = []
    for _, row in normalized.iterrows():
        coord = SkyCoord(parse_angle(row["ra"], True), parse_angle(row["dec"], False), frame="icrs")
        targets.append(
            Target(
                name=str(row["name"]),
                coord=coord,
                exposure_min=_float_or_default(row["exposure_min"], 30.0),
                priority=_float_or_default(row["priority"], 1.0),
            )
        )
    return targets


def targets_to_dataframe(targets: Iterable[Target]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "name": target.name,
                "coord": target.coord,
                "exposure_min": float(target.exposure_min),
                "priority": float(target.priority),
            }
            for target in targets
        ]
    )


def parse_targets(df: pd.DataFrame) -> pd.DataFrame:
    return targets_to_dataframe(parse_target_rows(df))


def _float_or_default(value: object, default: float) -> float:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(parsed):
        return default
    return float(parsed)
