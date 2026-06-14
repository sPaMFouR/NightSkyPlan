export interface Observatory {
  id: string;
  name: string;
  latitude_deg: number;
  longitude_deg: number;
  elevation_m: number;
  timezone: string;
  horizon_deg: number;
  zenith_deg: number;
}

export interface PlannerTarget {
  name: string;
  ra_deg: number;
  dec_deg: number;
  exposure_min: number;
  priority: number;
  aliases?: string[];
  transient_type?: string;
  host_name?: string;
  redshift?: string;
}

export interface TargetResolveResponse extends PlannerTarget {
  tns_name: string;
  prefix: string;
  objid: string;
}

export interface ConstraintInput {
  twilight_alt_deg: number;
  max_airmass: number;
  min_moon_sep_deg: number;
  ha_limit_hour: number;
}

export interface ScheduleRequest {
  observatory_id: string;
  date: string;
  cadence_min: number;
  overhead_percent: number;
  constraints: ConstraintInput;
  targets: PlannerTarget[];
}

export interface TrackWindow {
  start_utc: string;
  end_utc: string;
  start_local: string;
  end_local: string;
  duration_min: number;
}

export interface TrackSample {
  time_utc: string;
  time_local: string;
  alt_deg: number | null;
  az_deg: number | null;
  airmass: number | null;
  ha_hour: number | null;
  moon_sep_deg: number | null;
  sun_alt_deg: number | null;
  valid: boolean;
}

export interface TargetTrack {
  target: string;
  windows: TrackWindow[];
  samples: TrackSample[];
}

export interface ScheduleBlock {
  target: string;
  status: string;
  start_utc: string | null;
  end_utc: string | null;
  start_local: string | null;
  end_local: string | null;
  duration_min_with_overhead: number;
  reason: string;
  mean_alt_deg: number | null;
  mean_airmass: number | null;
  mean_moon_sep_deg: number | null;
  mean_abs_ha_hour: number | null;
}

export interface Diagnostic {
  target: string;
  observable_slots: number;
  window_count: number;
  first_window_start_utc: string | null;
  first_window_start_local: string | null;
  last_window_end_utc: string | null;
  last_window_end_local: string | null;
  longest_window_slots: number;
  duration_min_with_overhead: number;
  reason: string;
}

export interface ScheduleResponse {
  observatory: Observatory;
  date: string;
  cadence_min: number;
  overhead_percent: number;
  constraints: ConstraintInput;
  context: {
    dusk_utc: string | null;
    dusk_local: string | null;
    dawn_utc: string | null;
    dawn_local: string | null;
    dark_window_min: number | null;
  };
  tracks: TargetTrack[];
  schedule: ScheduleBlock[];
  diagnostics: Diagnostic[];
  score: {
    total_targets: number;
    scheduled_targets: number;
    completion_percent: number;
    mean_airmass: number | null;
    mean_moon_sep_deg: number | null;
    plan_score: number;
  };
}
