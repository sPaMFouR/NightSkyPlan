from __future__ import annotations

import csv
import io
from datetime import datetime, timezone

from webapp.backend.models import DiagnosticOut, ExportRequest, ScheduleBlockOut


def schedule_csv(payload: ExportRequest) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(
        [
            "target",
            "status",
            "start_utc",
            "end_utc",
            "start_local",
            "end_local",
            "duration_min_with_overhead",
            "mean_alt_deg",
            "mean_airmass",
            "mean_moon_sep_deg",
            "mean_abs_ha_hour",
            "reason",
        ]
    )
    for block in payload.schedule:
        writer.writerow(_schedule_row(block))
    writer.writerow([])
    writer.writerow(["diagnostics"])
    writer.writerow(
        [
            "target",
            "observable_slots",
            "window_count",
            "first_window_start_utc",
            "last_window_end_utc",
            "longest_window_slots",
            "duration_min_with_overhead",
            "reason",
        ]
    )
    for diagnostic in payload.diagnostics:
        writer.writerow(_diagnostic_row(diagnostic))
    return buffer.getvalue()


def schedule_ics(payload: ExportRequest) -> str:
    lines = [
        "BEGIN:VCALENDAR",
        "VERSION:2.0",
        "PRODID:-//NightSkyPlan//Professional Scheduler//EN",
        "CALSCALE:GREGORIAN",
        "METHOD:PUBLISH",
    ]
    stamp = _ics_now()
    for block in payload.schedule:
        if block.status != "SCHEDULED" or not block.start_utc or not block.end_utc:
            continue
        uid = f"{_ics_time(block.start_utc)}-{_safe_text(block.target)}@nightskyplan"
        description = (
            f"Airmass {block.mean_airmass or 'n/a'}; "
            f"Moon separation {block.mean_moon_sep_deg or 'n/a'} deg; "
            f"HA {block.mean_abs_ha_hour or 'n/a'} h"
        )
        lines.extend(
            [
                "BEGIN:VEVENT",
                f"UID:{uid}",
                f"DTSTAMP:{stamp}",
                f"DTSTART:{_ics_time(block.start_utc)}",
                f"DTEND:{_ics_time(block.end_utc)}",
                f"SUMMARY:Observe {_escape_ics(block.target)}",
                f"DESCRIPTION:{_escape_ics(description)}",
                "END:VEVENT",
            ]
        )
    lines.append("END:VCALENDAR")
    return "\r\n".join(lines) + "\r\n"


def _schedule_row(block: ScheduleBlockOut) -> list[object]:
    return [
        block.target,
        block.status,
        block.start_utc or "",
        block.end_utc or "",
        block.start_local or "",
        block.end_local or "",
        block.duration_min_with_overhead,
        block.mean_alt_deg or "",
        block.mean_airmass or "",
        block.mean_moon_sep_deg or "",
        block.mean_abs_ha_hour or "",
        block.reason,
    ]


def _diagnostic_row(diagnostic: DiagnosticOut) -> list[object]:
    return [
        diagnostic.target,
        diagnostic.observable_slots,
        diagnostic.window_count,
        diagnostic.first_window_start_utc or "",
        diagnostic.last_window_end_utc or "",
        diagnostic.longest_window_slots,
        diagnostic.duration_min_with_overhead,
        diagnostic.reason,
    ]


def _ics_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _ics_time(value: str) -> str:
    return datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _safe_text(value: str) -> str:
    return "".join(char if char.isalnum() else "-" for char in value).strip("-").lower()


def _escape_ics(value: str) -> str:
    return value.replace("\\", "\\\\").replace(";", "\\;").replace(",", "\\,").replace("\n", "\\n")
