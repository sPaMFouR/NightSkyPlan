from __future__ import annotations

from datetime import date
from typing import Optional

from pydantic import BaseModel, Field


class HealthResponse(BaseModel):
    status: str
    service: str


class ObservatoryOut(BaseModel):
    id: str
    name: str
    latitude_deg: float
    longitude_deg: float
    elevation_m: float
    timezone: str
    horizon_deg: float
    zenith_deg: float


class TargetResolveRequest(BaseModel):
    query: str = Field(min_length=1, max_length=120)


class TargetResolveResponse(BaseModel):
    name: str
    ra_deg: float
    dec_deg: float
    tns_name: str = ""
    prefix: str = ""
    objid: str = ""
    transient_type: str = ""
    redshift: str = ""
    host_name: str = ""
    aliases: list[str] = Field(default_factory=list)


class TargetInput(BaseModel):
    name: str = Field(min_length=1, max_length=120)
    ra_deg: float = Field(ge=0.0, lt=360.0)
    dec_deg: float = Field(ge=-90.0, le=90.0)
    exposure_min: float = Field(default=30.0, gt=0.0, le=1440.0)
    priority: float = Field(default=1.0, ge=0.0, le=10.0)


class ConstraintInput(BaseModel):
    twilight_alt_deg: float = Field(default=-18.0, ge=-24.0, le=-1.0)
    max_airmass: float = Field(default=2.5, ge=1.0, le=5.0)
    min_moon_sep_deg: float = Field(default=30.0, ge=0.0, le=180.0)
    ha_limit_hour: float = Field(default=6.0, ge=0.25, le=12.0)


class ScheduleRequest(BaseModel):
    observatory_id: str = Field(default="HCT", min_length=2, max_length=24)
    date: date
    cadence_min: int = Field(default=5, ge=1, le=60)
    overhead_percent: float = Field(default=20.0, ge=0.0, le=200.0)
    constraints: ConstraintInput = Field(default_factory=ConstraintInput)
    targets: list[TargetInput] = Field(min_length=1, max_length=200)


class TimeWindowOut(BaseModel):
    start_utc: str
    end_utc: str
    start_local: str
    end_local: str
    duration_min: float


class TrackSampleOut(BaseModel):
    time_utc: str
    time_local: str
    alt_deg: Optional[float]
    az_deg: Optional[float]
    airmass: Optional[float]
    ha_hour: Optional[float]
    moon_sep_deg: Optional[float]
    sun_alt_deg: Optional[float]
    valid: bool


class TargetTrackOut(BaseModel):
    target: str
    windows: list[TimeWindowOut]
    samples: list[TrackSampleOut]


class ScheduleBlockOut(BaseModel):
    target: str
    status: str
    start_utc: Optional[str] = None
    end_utc: Optional[str] = None
    start_local: Optional[str] = None
    end_local: Optional[str] = None
    duration_min_with_overhead: float
    reason: str = ""
    mean_alt_deg: Optional[float] = None
    mean_airmass: Optional[float] = None
    mean_moon_sep_deg: Optional[float] = None
    mean_abs_ha_hour: Optional[float] = None


class DiagnosticOut(BaseModel):
    target: str
    observable_slots: int
    window_count: int
    first_window_start_utc: Optional[str] = None
    first_window_start_local: Optional[str] = None
    last_window_end_utc: Optional[str] = None
    last_window_end_local: Optional[str] = None
    longest_window_slots: int
    duration_min_with_overhead: float
    reason: str


class SolarLunarContextOut(BaseModel):
    dusk_utc: Optional[str] = None
    dusk_local: Optional[str] = None
    dawn_utc: Optional[str] = None
    dawn_local: Optional[str] = None
    dark_window_min: Optional[float] = None


class ScoreSummaryOut(BaseModel):
    total_targets: int
    scheduled_targets: int
    completion_percent: float
    mean_airmass: Optional[float] = None
    mean_moon_sep_deg: Optional[float] = None
    plan_score: int


class ScheduleResponse(BaseModel):
    observatory: ObservatoryOut
    date: str
    cadence_min: int
    overhead_percent: float
    constraints: ConstraintInput
    context: SolarLunarContextOut
    tracks: list[TargetTrackOut]
    schedule: list[ScheduleBlockOut]
    diagnostics: list[DiagnosticOut]
    score: ScoreSummaryOut


class ExportRequest(BaseModel):
    schedule: list[ScheduleBlockOut] = Field(default_factory=list)
    diagnostics: list[DiagnosticOut] = Field(default_factory=list)
