from __future__ import annotations

import numpy as np
import pandas as pd


def greedy_schedule(tracks: pd.DataFrame, cadence_min: int, overhead_percent: float) -> pd.DataFrame:
    schedule, _ = build_greedy_plan(tracks, cadence_min, overhead_percent)
    return schedule


def build_greedy_plan(tracks: pd.DataFrame, cadence_min: int, overhead_percent: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Simple deterministic scheduler. Not globally optimal; good baseline for interactive triage."""
    times = np.array(sorted(tracks["time"].unique()))
    occupied = pd.Series(False, index=times)
    schedule_rows = []
    diagnostics = []
    targets = tracks.groupby("name").first()[["exposure_min", "priority"]].sort_values("priority", ascending=False)

    for name, meta in targets.iterrows():
        need_min = float(meta.exposure_min) * (1 + overhead_percent / 100.0)
        nslot = max(1, int(np.ceil(need_min / cadence_min)))
        sub = tracks[tracks["name"] == name].set_index("time").reindex(times)
        base_valid = sub["valid"].fillna(False).to_numpy()
        windows = contiguous_windows(times, base_valid)
        valid = base_valid & (~occupied.to_numpy())
        best = _best_window(sub, valid, nslot, float(meta.priority))
        if best is None:
            reason = _unscheduled_reason(sub, base_valid, windows, nslot, cadence_min)
            schedule_rows.append(
                {
                    "target": name,
                    "status": "UNSCHEDULED",
                    "start": None,
                    "end": None,
                    "duration_min_with_overhead": need_min,
                    "reason": reason,
                }
            )
            diagnostics.append(_diagnostic_row(name, sub, windows, need_min, reason))
            continue

        _, start_idx, window = best
        end_idx = start_idx + nslot
        occupied.iloc[start_idx:end_idx] = True
        schedule_rows.append(
            {
                "target": name,
                "status": "SCHEDULED",
                "start": times[start_idx],
                "end": times[end_idx],
                "duration_min_with_overhead": need_min,
                "reason": "",
                "mean_alt_deg": window["alt_deg"].mean(),
                "mean_airmass": window["airmass"].mean(),
                "mean_moon_sep_deg": window["moon_sep_deg"].mean(),
                "mean_abs_ha_hour": window["ha_hour"].abs().mean(),
            }
        )
        diagnostics.append(_diagnostic_row(name, sub, windows, need_min, "scheduled"))
    return pd.DataFrame(schedule_rows), pd.DataFrame(diagnostics)


def contiguous_windows(times: np.ndarray, valid: np.ndarray) -> list[tuple[object, object, int]]:
    windows = []
    start = None
    for idx, is_valid in enumerate(valid):
        if is_valid and start is None:
            start = idx
        if start is not None and (not is_valid or idx == len(valid) - 1):
            stop = idx if not is_valid else idx + 1
            if stop > start:
                end_idx = min(stop, len(times) - 1)
                windows.append((times[start], times[end_idx], stop - start))
            start = None
    return windows


def _best_window(sub: pd.DataFrame, valid: np.ndarray, nslot: int, priority: float) -> tuple[float, int, pd.DataFrame] | None:
    best = None
    for start_idx in range(0, len(sub) - nslot):
        if not valid[start_idx : start_idx + nslot].all():
            continue
        window = sub.iloc[start_idx : start_idx + nslot]
        score = (
            priority * 1000
            + window["alt_deg"].mean()
            - 10 * window["airmass"].mean()
            + 0.05 * window["moon_sep_deg"].mean()
            - 2 * window["ha_hour"].abs().mean()
        )
        if best is None or score > best[0]:
            best = (float(score), start_idx, window)
    return best


def _unscheduled_reason(
    sub: pd.DataFrame,
    base_valid: np.ndarray,
    windows: list[tuple[object, object, int]],
    nslot: int,
    cadence_min: int,
) -> str:
    if not base_valid.any():
        return _constraint_failure_summary(sub)
    longest_window_min = max((slots for _, _, slots in windows), default=0) * cadence_min
    required_min = nslot * cadence_min
    if longest_window_min < required_min:
        return f"observable windows too short; longest {longest_window_min:.0f} min, required {required_min:.0f} min"
    return "observable window exists but was occupied by higher-scoring targets"


def _constraint_failure_summary(sub: pd.DataFrame) -> str:
    checks = [
        ("twilight", "pass_twilight"),
        ("altitude", "pass_altitude"),
        ("airmass", "pass_airmass"),
        ("moon separation", "pass_moon_sep"),
        ("hour angle", "pass_hour_angle"),
    ]
    failed = [label for label, col in checks if col in sub and not sub[col].fillna(False).any()]
    if failed:
        return "never passes " + ", ".join(failed)
    return "never passes combined constraints"


def _diagnostic_row(
    name: str,
    sub: pd.DataFrame,
    windows: list[tuple[object, object, int]],
    duration_min: float,
    reason: str,
) -> dict[str, object]:
    return {
        "target": name,
        "observable_slots": int(sub["valid"].fillna(False).sum()) if "valid" in sub else 0,
        "window_count": len(windows),
        "first_window_start": windows[0][0] if windows else None,
        "last_window_end": windows[-1][1] if windows else None,
        "longest_window_slots": max((slots for _, _, slots in windows), default=0),
        "duration_min_with_overhead": duration_min,
        "reason": reason,
    }
