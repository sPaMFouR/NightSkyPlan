from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord


@dataclass(frozen=True)
class Observatory:
    name: str
    latitude_deg: float
    longitude_deg: float
    elevation_m: float
    timezone: str
    horizon_deg: float = 25.0
    zenith_deg: float = 88.0

    @property
    def location(self) -> EarthLocation:
        return EarthLocation(
            lat=self.latitude_deg * u.deg,
            lon=self.longitude_deg * u.deg,
            height=self.elevation_m * u.m,
        )


@dataclass(frozen=True)
class Target:
    name: str
    coord: SkyCoord
    exposure_min: float = 30.0
    priority: float = 1.0


@dataclass(frozen=True)
class ScheduleBlock:
    target: str
    status: str
    start: datetime | None
    end: datetime | None
    duration_min_with_overhead: float
    mean_alt_deg: float | None = None
    mean_airmass: float | None = None
    mean_moon_sep_deg: float | None = None
