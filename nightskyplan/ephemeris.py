from __future__ import annotations

from datetime import date, datetime, time, timedelta
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, get_body, get_sun
from astropy.time import Time
from astropy.utils import iers

from nightskyplan.models import Observatory

iers.conf.auto_download = False
iers.conf.iers_degraded_accuracy = "warn"


def make_time_grid(obs_date: date, timezone: str, cadence_min: int = 5) -> tuple[Time, pd.DatetimeIndex]:
    """Grid from local noon to next local noon, covering the whole observing night."""
    tz = ZoneInfo(timezone)
    start_local = datetime.combine(obs_date, time(12, 0), tzinfo=tz)
    end_local = start_local + timedelta(days=1)
    local_index = pd.date_range(start_local, end_local, freq=f"{cadence_min}min", inclusive="both")
    return Time(local_index.to_pydatetime()), local_index


def crossing_time(times: Time, y: np.ndarray, threshold: float, direction: str) -> Time | None:
    z = y - threshold
    if direction == "down":
        idx = np.where((z[:-1] > 0) & (z[1:] <= 0))[0]
    else:
        idx = np.where((z[:-1] < 0) & (z[1:] >= 0))[0]
    if len(idx) == 0:
        return None
    i = idx[0] if direction == "down" else idx[-1]
    frac = abs(z[i]) / (abs(z[i]) + abs(z[i + 1]))
    return times[i] + frac * (times[i + 1] - times[i])


def solar_lunar_context(times: Time, obs: Observatory, twilight_alt_deg: float) -> dict[str, object]:
    frame = AltAz(obstime=times, location=obs.location)
    sun_alt = get_sun(times).transform_to(frame).alt.deg
    moon = get_body("moon", times, obs.location)
    moon_alt = moon.transform_to(frame).alt.deg
    dusk = crossing_time(times, sun_alt, twilight_alt_deg, "down")
    dawn = crossing_time(times, sun_alt, twilight_alt_deg, "up")
    return {"sun_alt": sun_alt, "moon": moon, "moon_alt": moon_alt, "dusk": dusk, "dawn": dawn}
