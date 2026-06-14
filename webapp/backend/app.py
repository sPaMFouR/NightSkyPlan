from __future__ import annotations

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response

from webapp.backend.exports import schedule_csv, schedule_ics
from webapp.backend.models import (
    ExportRequest,
    HealthResponse,
    ObservatoryOut,
    ScheduleRequest,
    ScheduleResponse,
    TargetResolveRequest,
    TargetResolveResponse,
)
from webapp.backend.observatories import list_observatories
from webapp.backend.scheduler import build_schedule_response
from webapp.backend.tns import TNSClient, TNSCredentialsError, TNSLookupError


app = FastAPI(title="NightSkyPlan Professional API", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://127.0.0.1:5173",
        "http://localhost:5173",
    ],
    allow_credentials=False,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["*"],
)


@app.get("/api/health", response_model=HealthResponse)
def health() -> HealthResponse:
    return HealthResponse(status="ok", service="nightskyplan-webapp")


@app.get("/api/observatories", response_model=list[ObservatoryOut])
def observatories() -> list[dict[str, object]]:
    return list_observatories()


@app.post("/api/targets/resolve", response_model=TargetResolveResponse)
def resolve_target(request: TargetResolveRequest) -> TargetResolveResponse:
    try:
        return TNSClient().lookup(request.query).to_response()
    except TNSCredentialsError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except TNSLookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.post("/api/schedule", response_model=ScheduleResponse)
def schedule(request: ScheduleRequest) -> ScheduleResponse:
    try:
        return build_schedule_response(request)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/exports/csv")
def export_csv(request: ExportRequest) -> Response:
    return Response(
        schedule_csv(request),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="nightskyplan_schedule.csv"'},
    )


@app.post("/api/exports/ics")
def export_ics(request: ExportRequest) -> Response:
    return Response(
        schedule_ics(request),
        media_type="text/calendar",
        headers={"Content-Disposition": 'attachment; filename="nightskyplan_schedule.ics"'},
    )
