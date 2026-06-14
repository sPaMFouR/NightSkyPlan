import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";
import App from "./App";
import type { ScheduleResponse } from "./types";

const observatories = [
  {
    id: "HCT",
    name: "HCT / Hanle",
    latitude_deg: 32.7794,
    longitude_deg: 78.9642,
    elevation_m: 4486,
    timezone: "Asia/Kolkata",
    horizon_deg: 25,
    zenith_deg: 85,
  },
];

const scheduleResponse: ScheduleResponse = {
  observatory: observatories[0],
  date: "2026-06-14",
  cadence_min: 5,
  overhead_percent: 20,
  constraints: {
    twilight_alt_deg: -18,
    max_airmass: 2.5,
    min_moon_sep_deg: 30,
    ha_limit_hour: 6,
  },
  context: {
    dusk_utc: "2026-06-14T14:00:00Z",
    dusk_local: "2026-06-14T19:30:00+05:30",
    dawn_utc: "2026-06-14T22:00:00Z",
    dawn_local: "2026-06-15T03:30:00+05:30",
    dark_window_min: 480,
  },
  tracks: [],
  schedule: [
    {
      target: "SN2022jli",
      status: "SCHEDULED",
      start_utc: "2026-06-14T19:00:00Z",
      end_utc: "2026-06-14T19:45:00Z",
      start_local: "2026-06-15T00:30:00+05:30",
      end_local: "2026-06-15T01:15:00+05:30",
      duration_min_with_overhead: 45,
      reason: "",
      mean_alt_deg: 56,
      mean_airmass: 1.2,
      mean_moon_sep_deg: 80,
      mean_abs_ha_hour: 1.1,
    },
  ],
  diagnostics: [],
  score: {
    total_targets: 3,
    scheduled_targets: 1,
    completion_percent: 33.3,
    mean_airmass: 1.2,
    mean_moon_sep_deg: 80,
    plan_score: 72,
  },
};

describe("NightSkyPlan webapp", () => {
  beforeEach(() => {
    HTMLCanvasElement.prototype.getContext = vi.fn(() => canvasContext) as unknown as typeof HTMLCanvasElement.prototype.getContext;
    window.matchMedia = vi.fn().mockImplementation(query => ({
      matches: false,
      media: query,
      onchange: null,
      addListener: vi.fn(),
      removeListener: vi.fn(),
      addEventListener: vi.fn(),
      removeEventListener: vi.fn(),
      dispatchEvent: vi.fn(),
    }));
    vi.stubGlobal(
      "fetch",
      vi.fn(async (input: RequestInfo | URL) => {
        const url = String(input);
        if (url.endsWith("/api/observatories")) {
          return jsonResponse(observatories);
        }
        if (url.endsWith("/api/schedule")) {
          return jsonResponse(scheduleResponse);
        }
        if (url.endsWith("/api/targets/resolve")) {
          return jsonResponse({
            name: "SN 2023ixf",
            ra_deg: 210.910675,
            dec_deg: 54.311651,
            tns_name: "2023ixf",
            prefix: "SN",
            objid: "123",
            aliases: ["ZTF23abc"],
            transient_type: "SN II",
            redshift: "0.0008",
            host_name: "M101",
          });
        }
        return new Response("not found", { status: 404 });
      })
    );
  });

  it("adds a manual target to the editable queue", async () => {
    render(<App />);

    fireEvent.change(screen.getByPlaceholderText("Target name"), { target: { value: "Manual SN" } });
    fireEvent.change(screen.getByPlaceholderText("RA deg"), { target: { value: "120" } });
    fireEvent.change(screen.getByPlaceholderText("Dec deg"), { target: { value: "-22" } });
    fireEvent.click(screen.getByRole("button", { name: /Add target/i }));

    expect(await screen.findByDisplayValue("Manual SN")).toBeTruthy();
  });

  it("runs the scheduler and enables exports", async () => {
    render(<App />);

    fireEvent.click(screen.getByRole("button", { name: /Run automatic scheduler/i }));

    await waitFor(() => expect(screen.getByText(/Scheduled 1\/3 targets with score 72/i)).toBeTruthy());
    expect((screen.getByRole("button", { name: /Export CSV/i }) as HTMLButtonElement).disabled).toBe(false);
    expect((screen.getByRole("button", { name: /Export ICS/i }) as HTMLButtonElement).disabled).toBe(false);
  });
});

const gradient = {
  addColorStop: vi.fn(),
};

const canvasContext = {
  setTransform: vi.fn(),
  clearRect: vi.fn(),
  createRadialGradient: vi.fn(() => gradient),
  createLinearGradient: vi.fn(() => gradient),
  fillRect: vi.fn(),
  beginPath: vi.fn(),
  arc: vi.fn(),
  fill: vi.fn(),
  stroke: vi.fn(),
  moveTo: vi.fn(),
  lineTo: vi.fn(),
  save: vi.fn(),
  clip: vi.fn(),
  restore: vi.fn(),
  fillStyle: "",
  strokeStyle: "",
  globalAlpha: 1,
  lineWidth: 1,
};

function jsonResponse(payload: unknown) {
  return new Response(JSON.stringify(payload), {
    status: 200,
    headers: { "Content-Type": "application/json" },
  });
}
