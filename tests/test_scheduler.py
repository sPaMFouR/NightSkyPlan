import pandas as pd

from nightskyplan.scheduling import build_greedy_plan, contiguous_windows, greedy_schedule


def test_greedy_schedule_blocks_do_not_overlap():
    times = pd.date_range("2026-01-01 20:00", periods=8, freq="10min")
    rows = []
    for name, priority in [("A", 2), ("B", 1)]:
        for i, timestamp in enumerate(times):
            rows.append(
                {
                    "name": name,
                    "time": timestamp.to_pydatetime(),
                    "valid": True,
                    "alt_deg": 45.0 + i,
                    "airmass": 1.3,
                    "moon_sep_deg": 80.0,
                    "ha_hour": 0.5,
                    "exposure_min": 20.0,
                    "priority": priority,
                }
            )
    schedule = greedy_schedule(pd.DataFrame(rows), cadence_min=10, overhead_percent=0)
    scheduled = schedule[schedule["status"] == "SCHEDULED"].sort_values("start")

    assert len(scheduled) == 2
    previous_end = None
    for _, row in scheduled.iterrows():
        if previous_end is not None:
            assert row["start"] >= previous_end
        previous_end = row["end"]


def test_greedy_schedule_applies_overhead_to_block_duration():
    times = pd.date_range("2026-01-01 20:00", periods=8, freq="10min")
    tracks = _tracks_for("A", times, valid=True, exposure_min=20.0)

    schedule, diagnostics = build_greedy_plan(tracks, cadence_min=10, overhead_percent=50)

    row = schedule.iloc[0]
    assert row["status"] == "SCHEDULED"
    assert row["duration_min_with_overhead"] == 30.0
    assert (pd.Timestamp(row["end"]) - pd.Timestamp(row["start"])).total_seconds() / 60 == 30.0
    assert diagnostics.iloc[0]["reason"] == "scheduled"


def test_greedy_schedule_reports_impossible_target_reason():
    times = pd.date_range("2026-01-01 20:00", periods=4, freq="10min")
    tracks = _tracks_for("A", times, valid=False, exposure_min=10.0)
    tracks["pass_altitude"] = False
    tracks["pass_airmass"] = True
    tracks["pass_moon_sep"] = True
    tracks["pass_hour_angle"] = True
    tracks["pass_twilight"] = True

    schedule, diagnostics = build_greedy_plan(tracks, cadence_min=10, overhead_percent=0)

    assert schedule.iloc[0]["status"] == "UNSCHEDULED"
    assert schedule.iloc[0]["reason"] == "never passes altitude"
    assert diagnostics.iloc[0]["reason"] == "never passes altitude"


def test_contiguous_windows_handles_window_at_end_of_grid():
    times = pd.date_range("2026-01-01 20:00", periods=5, freq="10min").to_numpy()
    windows = contiguous_windows(times, valid=[False, True, True, False, True])

    assert windows == [(times[1], times[3], 2), (times[4], times[4], 1)]


def _tracks_for(name: str, times: pd.DatetimeIndex, valid: bool, exposure_min: float) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "name": name,
            "time": [timestamp.to_pydatetime() for timestamp in times],
            "valid": valid,
            "alt_deg": 50.0,
            "airmass": 1.2,
            "moon_sep_deg": 90.0,
            "ha_hour": 0.5,
            "exposure_min": exposure_min,
            "priority": 1,
        }
    )
