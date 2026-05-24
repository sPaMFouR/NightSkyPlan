from __future__ import annotations

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import AltAz, SkyCoord
from astropy.time import Time

from nightskyplan.models import Observatory


def target_track(target: pd.Series, times: Time, obs: Observatory, context: dict[str, object]) -> pd.DataFrame:
    coord: SkyCoord = target["coord"]
    frame = AltAz(obstime=times, location=obs.location)
    altaz = coord.transform_to(frame)
    alt = altaz.alt.deg
    az = altaz.az.deg
    airmass = np.where(alt > 0, 1.0 / np.cos(np.deg2rad(90.0 - alt)), np.nan)
    lst = times.sidereal_time("apparent", longitude=obs.location.lon)
    ha = (lst - coord.ra).wrap_at(12 * u.hourangle).hour
    moon_sep = coord.separation(context["moon"]).deg
    return pd.DataFrame(
        {
            "name": target["name"],
            "time": times.to_datetime(),
            "alt_deg": alt,
            "az_deg": az,
            "airmass": airmass,
            "ha_hour": ha,
            "moon_sep_deg": moon_sep,
            "priority": float(target["priority"]),
            "exposure_min": float(target["exposure_min"]),
        }
    )


def compute_tracks(targets: pd.DataFrame, times: Time, obs: Observatory, context: dict[str, object]) -> pd.DataFrame:
    if targets.empty:
        return pd.DataFrame()
    return pd.concat([target_track(row, times, obs, context) for _, row in targets.iterrows()], ignore_index=True)


def apply_constraints(
    tracks: pd.DataFrame,
    times: Time,
    context: dict[str, object],
    obs: Observatory,
    twilight_alt_deg: float,
    max_airmass: float,
    min_moon_sep_deg: float,
    ha_limit_hour: float,
) -> pd.DataFrame:
    sun_alt_map = pd.Series(context["sun_alt"], index=times.to_datetime())
    constrained = tracks.copy()
    constrained["sun_alt_deg"] = constrained["time"].map(sun_alt_map)
    constrained["pass_twilight"] = constrained["sun_alt_deg"] < twilight_alt_deg
    constrained["pass_altitude"] = (constrained["alt_deg"] >= obs.horizon_deg) & (constrained["alt_deg"] <= obs.zenith_deg)
    constrained["pass_airmass"] = constrained["airmass"] <= max_airmass
    constrained["pass_moon_sep"] = constrained["moon_sep_deg"] >= min_moon_sep_deg
    constrained["pass_hour_angle"] = constrained["ha_hour"].abs() <= ha_limit_hour
    constrained["valid"] = constrained[
        ["pass_twilight", "pass_altitude", "pass_airmass", "pass_moon_sep", "pass_hour_angle"]
    ].all(axis=1)
    return constrained


def visibility_mask(
    tracks: pd.DataFrame,
    horizon_deg: float,
    zenith_deg: float,
    max_airmass: float,
    min_moon_sep_deg: float,
    ha_limit_hour: float,
) -> pd.Series:
    return (
        (tracks["alt_deg"] >= horizon_deg)
        & (tracks["alt_deg"] <= zenith_deg)
        & (tracks["airmass"] <= max_airmass)
        & (tracks["moon_sep_deg"] >= min_moon_sep_deg)
        & (tracks["ha_hour"].abs() <= ha_limit_hour)
    )
