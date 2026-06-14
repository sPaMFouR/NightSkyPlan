from fastapi.testclient import TestClient

from webapp.backend.app import app


client = TestClient(app)


def test_health_and_observatories_endpoints():
    health = client.get("/api/health")
    observatories = client.get("/api/observatories")

    assert health.status_code == 200
    assert health.json()["status"] == "ok"
    assert observatories.status_code == 200
    assert any(item["id"] == "HCT" for item in observatories.json())


def test_resolve_target_reports_missing_tns_credentials(monkeypatch):
    monkeypatch.delenv("TNS_API_KEY", raising=False)
    monkeypatch.delenv("TNS_BOT_ID", raising=False)
    monkeypatch.delenv("TNS_BOT_NAME", raising=False)

    response = client.post("/api/targets/resolve", json={"query": "SN 2023ixf"})

    assert response.status_code == 503
    assert "TNS_API_KEY" in response.json()["detail"]


def test_schedule_endpoint_returns_tracks_schedule_and_diagnostics():
    response = client.post(
        "/api/schedule",
        json={
            "observatory_id": "HCT",
            "date": "2026-06-14",
            "cadence_min": 30,
            "overhead_percent": 20,
            "constraints": {
                "twilight_alt_deg": -12,
                "max_airmass": 4.5,
                "min_moon_sep_deg": 0,
                "ha_limit_hour": 12,
            },
            "targets": [
                {"name": "SN2020tlf", "ra_deg": 220.0418, "dec_deg": 42.7777, "exposure_min": 30, "priority": 2},
                {"name": "SN2018zd", "ra_deg": 94.5133, "dec_deg": 78.3669, "exposure_min": 20, "priority": 1},
            ],
        },
    )

    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["observatory"]["id"] == "HCT"
    assert payload["score"]["total_targets"] == 2
    assert len(payload["tracks"]) == 2
    assert len(payload["schedule"]) == 2
    assert len(payload["diagnostics"]) == 2


def test_export_endpoints_return_downloadable_content():
    payload = {
        "schedule": [
            {
                "target": "SN2020tlf",
                "status": "SCHEDULED",
                "start_utc": "2026-06-14T19:00:00Z",
                "end_utc": "2026-06-14T19:30:00Z",
                "start_local": "2026-06-15T00:30:00+05:30",
                "end_local": "2026-06-15T01:00:00+05:30",
                "duration_min_with_overhead": 30,
                "reason": "",
                "mean_alt_deg": 50,
                "mean_airmass": 1.3,
                "mean_moon_sep_deg": 80,
                "mean_abs_ha_hour": 1.2,
            }
        ],
        "diagnostics": [],
    }

    csv_response = client.post("/api/exports/csv", json=payload)
    ics_response = client.post("/api/exports/ics", json=payload)

    assert csv_response.status_code == 200
    assert "SN2020tlf" in csv_response.text
    assert ics_response.status_code == 200
    assert "BEGIN:VCALENDAR" in ics_response.text
    assert "Observe SN2020tlf" in ics_response.text
