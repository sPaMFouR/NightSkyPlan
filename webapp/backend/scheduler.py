from __future__ import annotations

import math
from datetime import datetime
from zoneinfo import ZoneInfo

import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from nightskyplan.constraints import apply_constraints, compute_tracks
from nightskyplan.ephemeris import make_time_grid, solar_lunar_context
from nightskyplan.models import Observatory
from nightskyplan.scheduling import build_greedy_plan, contiguous_windows
from webapp.backend.models import (
    DiagnosticOut,
    ObservatoryOut,
    ScheduleBlockOut,
    ScheduleRequest,
    ScheduleResponse,
    ScoreSummaryOut,
    SolarLunarContextOut,
    TargetTrackOut,
    TimeWindowOut,
    TrackSampleOut,
)
from webapp.backend.observatories import get_observatory


def build_schedule_response(request: ScheduleRequest) -> ScheduleResponse:
    observatory = get_observatory(request.observatory_id)
    targets = _targets_to_frame(request)
    times, _ = make_time_grid(request.date, observatory.timezone, request.cadence_min)
    context = solar_lunar_context(times, observatory, request.constraints.twilight_alt_deg)
    tracks = compute_tracks(targets, times, observatory, context)
    constrained = apply_constraints(
        tracks,
        times,
        context,
        observatory,
        request.constraints.twilight_alt_deg,
        request.constraints.max_airmass,
        request.constraints.min_moon_sep_deg,
        request.constraints.ha_limit_hour,
    )
    schedule_df, diagnostics_df = build_greedy_plan(constrained, request.cadence_min, request.overhead_percent)

    return ScheduleResponse(
        observatory=_observatory_out(request.observatory_id, observatory),
        date=request.date.isoformat(),
        cadence_min=request.cadence_min,
        overhead_percent=request.overhead_percent,
        constraints=request.constraints,
        context=_context_out(context, observatory.timezone),
        tracks=_tracks_out(constrained, observatory.timezone, request.cadence_min),
        schedule=_schedule_out(schedule_df, observatory.timezone),
        diagnostics=_diagnostics_out(diagnostics_df, observatory.timezone),
        score=_score_out(schedule_df, len(request.targets)),
    )


def _targets_to_frame(request: ScheduleRequest) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "name": target.name,
                "coord": SkyCoord(target.ra_deg * u.deg, target.dec_deg * u.deg, frame="icrs"),
                "exposure_min": target.exposure_min,
                "priority": target.priority,
            }
            for target in request.targets
        ]
    )


def _observatory_out(observatory_id: str, observatory: Observatory) -> ObservatoryOut:
    return ObservatoryOut(
        id=observatory_id.upper(),
        name=observatory.name,
        latitude_deg=observatory.latitude_deg,
        longitude_deg=observatory.longitude_deg,
        elevation_m=observatory.elevation_m,
        timezone=observatory.timezone,
        horizon_deg=observatory.horizon_deg,
        zenith_deg=observatory.zenith_deg,
    )


def _context_out(context: dict[str, object], timezone: str) -> SolarLunarContextOut:
    dusk = _time_to_datetime(context.get("dusk"))
    dawn = _time_to_datetime(context.get("dawn"))
    dark_window_min = None
    if dusk is not None and dawn is not None:
        dark_window_min = max(0.0, (pd.Timestamp(dawn) - pd.Timestamp(dusk)).total_seconds() / 60.0)
    return SolarLunarContextOut(
        dusk_utc=_utc_iso(dusk),
        dusk_local=_local_iso(dusk, timezone),
        dawn_utc=_utc_iso(dawn),
        dawn_local=_local_iso(dawn, timezone),
        dark_window_min=dark_window_min,
    )


def _tracks_out(tracks: pd.DataFrame, timezone: str, cadence_min: int) -> list[TargetTrackOut]:
    output = []
    for name, group in tracks.sort_values(["name", "time"]).groupby("name", sort=False):
        valid = group["valid"].fillna(False).to_numpy()
        times = group["time"].to_numpy()
        windows = []
        for start, end, slots in contiguous_windows(times, valid):
            windows.append(
                TimeWindowOut(
                    start_utc=_utc_iso(start) or "",
                    end_utc=_utc_iso(end) or "",
                    start_local=_local_iso(start, timezone) or "",
                    end_local=_local_iso(end, timezone) or "",
                    duration_min=float(slots * cadence_min),
                )
            )

        samples = [
            TrackSampleOut(
                time_utc=_utc_iso(row["time"]) or "",
                time_local=_local_iso(row["time"], timezone) or "",
                alt_deg=_clean_float(row.get("alt_deg")),
                az_deg=_clean_float(row.get("az_deg")),
                airmass=_clean_float(row.get("airmass")),
                ha_hour=_clean_float(row.get("ha_hour")),
                moon_sep_deg=_clean_float(row.get("moon_sep_deg")),
                sun_alt_deg=_clean_float(row.get("sun_alt_deg")),
                valid=bool(row.get("valid")),
            )
            for _, row in group.iterrows()
        ]
        output.append(TargetTrackOut(target=str(name), windows=windows, samples=samples))
    return output


def _schedule_out(schedule: pd.DataFrame, timezone: str) -> list[ScheduleBlockOut]:
    rows = []
    for _, row in schedule.iterrows():
        rows.append(
            ScheduleBlockOut(
                target=str(row["target"]),
                status=str(row["status"]),
                start_utc=_utc_iso(row.get("start")),
                end_utc=_utc_iso(row.get("end")),
                start_local=_local_iso(row.get("start"), timezone),
                end_local=_local_iso(row.get("end"), timezone),
                duration_min_with_overhead=float(row["duration_min_with_overhead"]),
                reason=str(row.get("reason") or ""),
                mean_alt_deg=_clean_float(row.get("mean_alt_deg")),
                mean_airmass=_clean_float(row.get("mean_airmass")),
                mean_moon_sep_deg=_clean_float(row.get("mean_moon_sep_deg")),
                mean_abs_ha_hour=_clean_float(row.get("mean_abs_ha_hour")),
            )
        )
    return rows


def _diagnostics_out(diagnostics: pd.DataFrame, timezone: str) -> list[DiagnosticOut]:
    rows = []
    for _, row in diagnostics.iterrows():
        rows.append(
            DiagnosticOut(
                target=str(row["target"]),
                observable_slots=int(row.get("observable_slots") or 0),
                window_count=int(row.get("window_count") or 0),
                first_window_start_utc=_utc_iso(row.get("first_window_start")),
                first_window_start_local=_local_iso(row.get("first_window_start"), timezone),
                last_window_end_utc=_utc_iso(row.get("last_window_end")),
                last_window_end_local=_local_iso(row.get("last_window_end"), timezone),
                longest_window_slots=int(row.get("longest_window_slots") or 0),
                duration_min_with_overhead=float(row.get("duration_min_with_overhead") or 0),
                reason=str(row.get("reason") or ""),
            )
        )
    return rows


def _score_out(schedule: pd.DataFrame, total_targets: int) -> ScoreSummaryOut:
    scheduled = schedule[schedule["status"] == "SCHEDULED"] if not schedule.empty else pd.DataFrame()
    scheduled_count = int(len(scheduled))
    completion = 0.0 if total_targets == 0 else scheduled_count / total_targets
    mean_airmass = _clean_float(scheduled["mean_airmass"].mean()) if not scheduled.empty else None
    mean_moon_sep = _clean_float(scheduled["mean_moon_sep_deg"].mean()) if not scheduled.empty else None
    airmass_component = 0.0 if mean_airmass is None else max(0.0, 1.0 - (mean_airmass - 1.0) / 2.4)
    moon_component = 0.0 if mean_moon_sep is None else min(1.0, mean_moon_sep / 120.0)
    plan_score = round(100 * (0.62 * completion + 0.24 * airmass_component + 0.14 * moon_component))
    return ScoreSummaryOut(
        total_targets=total_targets,
        scheduled_targets=scheduled_count,
        completion_percent=round(completion * 100, 1),
        mean_airmass=mean_airmass,
        mean_moon_sep_deg=mean_moon_sep,
        plan_score=max(0, min(100, int(plan_score))),
    )


def _time_to_datetime(value: object) -> datetime | None:
    if value is None:
        return None
    if isinstance(value, Time):
        return value.to_datetime()
    if _is_missing(value):
        return None
    return pd.Timestamp(value).to_pydatetime()


def _utc_timestamp(value: object) -> pd.Timestamp | None:
    if value is None or _is_missing(value):
        return None
    if isinstance(value, Time):
        value = value.to_datetime()
    ts = pd.Timestamp(value)
    if ts.tzinfo is None:
        return ts.tz_localize("UTC")
    return ts.tz_convert("UTC")


def _utc_iso(value: object) -> str | None:
    ts = _utc_timestamp(value)
    if ts is None:
        return None
    return ts.isoformat().replace("+00:00", "Z")


def _local_iso(value: object, timezone: str) -> str | None:
    ts = _utc_timestamp(value)
    if ts is None:
        return None
    return ts.tz_convert(ZoneInfo(timezone)).isoformat()


def _clean_float(value: object) -> float | None:
    if value is None or _is_missing(value):
        return None
    parsed = float(value)
    if not math.isfinite(parsed):
        return None
    return round(parsed, 4)


def _is_missing(value: object) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False
