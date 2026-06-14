import type {
  Observatory,
  ScheduleRequest,
  ScheduleResponse,
  TargetResolveResponse,
} from "./types";

const API_BASE = import.meta.env.VITE_API_BASE_URL ?? "";

async function requestJson<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${API_BASE}${path}`, {
    ...init,
    headers: {
      "Content-Type": "application/json",
      ...(init?.headers ?? {}),
    },
  });
  const data = await response.json().catch(() => ({}));
  if (!response.ok) {
    throw new Error(data.detail || data.error || `Request failed with ${response.status}`);
  }
  return data as T;
}

export function fetchObservatories(): Promise<Observatory[]> {
  return requestJson<Observatory[]>("/api/observatories");
}

export function resolveTarget(query: string): Promise<TargetResolveResponse> {
  return requestJson<TargetResolveResponse>("/api/targets/resolve", {
    method: "POST",
    body: JSON.stringify({ query }),
  });
}

export function runSchedule(payload: ScheduleRequest): Promise<ScheduleResponse> {
  return requestJson<ScheduleResponse>("/api/schedule", {
    method: "POST",
    body: JSON.stringify(payload),
  });
}

export async function downloadExport(
  format: "csv" | "ics",
  result: Pick<ScheduleResponse, "schedule" | "diagnostics">
): Promise<Blob> {
  const response = await fetch(`${API_BASE}/api/exports/${format}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ schedule: result.schedule, diagnostics: result.diagnostics }),
  });
  if (!response.ok) {
    const text = await response.text();
    throw new Error(text || `Export failed with ${response.status}`);
  }
  return response.blob();
}
