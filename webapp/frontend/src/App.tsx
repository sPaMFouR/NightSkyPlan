import { type CSSProperties, useEffect, useMemo, useState } from "react";
import {
  Activity,
  CalendarClock,
  Download,
  Globe2,
  Loader2,
  Plus,
  Radar,
  Search,
  SlidersHorizontal,
  Trash2,
} from "lucide-react";
import { BackgroundCanvas } from "./BackgroundCanvas";
import { EarthGlobe } from "./EarthGlobe";
import { downloadExport, fetchObservatories, resolveTarget, runSchedule } from "./api";
import type { ConstraintInput, Observatory, PlannerTarget, ScheduleResponse } from "./types";
import "./styles.css";

const today = new Date().toISOString().slice(0, 10);

const fallbackObservatories: Observatory[] = [
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

const initialTargets: PlannerTarget[] = [
  { name: "SN2022jli", ra_deg: 8.6883, dec_deg: -8.3903, exposure_min: 45, priority: 3 },
  { name: "SN2020tlf", ra_deg: 220.0418, dec_deg: 42.7777, exposure_min: 30, priority: 2 },
  { name: "SN2018zd", ra_deg: 94.5133, dec_deg: 78.3669, exposure_min: 20, priority: 1 },
];

const defaultConstraints: ConstraintInput = {
  twilight_alt_deg: -18,
  max_airmass: 2.5,
  min_moon_sep_deg: 30,
  ha_limit_hour: 6,
};

const emptyManualTarget = {
  name: "",
  ra_deg: "",
  dec_deg: "",
  exposure_min: "30",
  priority: "1",
};

function App() {
  const [observatories, setObservatories] = useState<Observatory[]>(fallbackObservatories);
  const [observatoryId, setObservatoryId] = useState("HCT");
  const [date, setDate] = useState(today);
  const [cadenceMin, setCadenceMin] = useState(5);
  const [overheadPercent, setOverheadPercent] = useState(20);
  const [constraints, setConstraints] = useState<ConstraintInput>(defaultConstraints);
  const [targets, setTargets] = useState<PlannerTarget[]>(initialTargets);
  const [tnsQuery, setTnsQuery] = useState("SN 2023ixf");
  const [manualTarget, setManualTarget] = useState(emptyManualTarget);
  const [result, setResult] = useState<ScheduleResponse | null>(null);
  const [status, setStatus] = useState("Ready to run scheduler");
  const [isScheduling, setIsScheduling] = useState(false);
  const [isResolving, setIsResolving] = useState(false);
  const [isExporting, setIsExporting] = useState<"csv" | "ics" | null>(null);

  useEffect(() => {
    fetchObservatories()
      .then(data => {
        setObservatories(data);
        if (!data.some(item => item.id === observatoryId) && data[0]) {
          setObservatoryId(data[0].id);
        }
      })
      .catch(error => {
        setStatus(`Backend unavailable: ${error.message}`);
      });
  }, [observatoryId]);

  const selectedObservatory = useMemo(
    () => observatories.find(item => item.id === observatoryId) ?? observatories[0],
    [observatories, observatoryId]
  );

  const scheduledBlocks = result?.schedule.filter(item => item.status === "SCHEDULED") ?? [];
  const unscheduledBlocks = result?.schedule.filter(item => item.status !== "SCHEDULED") ?? [];

  async function handleResolveTarget() {
    if (!tnsQuery.trim()) {
      setStatus("Enter a TNS, IAU, or ZTF target name before resolving.");
      return;
    }
    setIsResolving(true);
    setStatus("Resolving target through server-side TNS credentials...");
    try {
      const resolved = await resolveTarget(tnsQuery);
      setTargets(current => [
        ...current,
        {
          name: resolved.name,
          ra_deg: resolved.ra_deg,
          dec_deg: resolved.dec_deg,
          exposure_min: 30,
          priority: 2,
          aliases: resolved.aliases,
          transient_type: resolved.transient_type,
          host_name: resolved.host_name,
          redshift: resolved.redshift,
        },
      ]);
      setStatus(`Added ${resolved.name} from TNS.`);
    } catch (error) {
      setStatus(error instanceof Error ? error.message : "TNS lookup failed.");
    } finally {
      setIsResolving(false);
    }
  }

  function handleAddManualTarget() {
    const parsed = {
      name: manualTarget.name.trim(),
      ra_deg: Number(manualTarget.ra_deg),
      dec_deg: Number(manualTarget.dec_deg),
      exposure_min: Number(manualTarget.exposure_min),
      priority: Number(manualTarget.priority),
    };
    if (!parsed.name || !Number.isFinite(parsed.ra_deg) || !Number.isFinite(parsed.dec_deg)) {
      setStatus("Manual targets need a name, RA degrees, and Dec degrees.");
      return;
    }
    if (parsed.ra_deg < 0 || parsed.ra_deg >= 360 || parsed.dec_deg < -90 || parsed.dec_deg > 90) {
      setStatus("Manual target coordinates must be RA 0-360 deg and Dec -90 to +90 deg.");
      return;
    }
    setTargets(current => [...current, parsed]);
    setManualTarget(emptyManualTarget);
    setStatus(`Added ${parsed.name} to the observing queue.`);
  }

  function updateTarget(index: number, key: keyof PlannerTarget, value: string) {
    setTargets(current =>
      current.map((target, targetIndex) => {
        if (targetIndex !== index) return target;
        if (key === "name") return { ...target, name: value };
        return { ...target, [key]: Number(value) };
      })
    );
  }

  function removeTarget(index: number) {
    setTargets(current => current.filter((_, targetIndex) => targetIndex !== index));
  }

  async function handleRunSchedule() {
    if (targets.length === 0) {
      setStatus("Add at least one target before running the scheduler.");
      return;
    }
    setIsScheduling(true);
    setStatus("Computing observability tracks and schedule windows...");
    try {
      const payload = {
        observatory_id: observatoryId,
        date,
        cadence_min: cadenceMin,
        overhead_percent: overheadPercent,
        constraints,
        targets,
      };
      const schedule = await runSchedule(payload);
      setResult(schedule);
      setStatus(
        `Scheduled ${schedule.score.scheduled_targets}/${schedule.score.total_targets} targets with score ${schedule.score.plan_score}.`
      );
    } catch (error) {
      setStatus(error instanceof Error ? error.message : "Scheduling failed.");
    } finally {
      setIsScheduling(false);
    }
  }

  async function handleExport(format: "csv" | "ics") {
    if (!result) {
      setStatus("Run the scheduler before exporting.");
      return;
    }
    setIsExporting(format);
    try {
      const blob = await downloadExport(format, result);
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = `nightskyplan_schedule.${format}`;
      link.click();
      URL.revokeObjectURL(url);
      setStatus(`Exported ${format.toUpperCase()} schedule.`);
    } catch (error) {
      setStatus(error instanceof Error ? error.message : "Export failed.");
    } finally {
      setIsExporting(null);
    }
  }

  return (
    <div className="app-shell">
      <BackgroundCanvas />
      <header className="hero" id="top">
        <nav className="top-nav" aria-label="Primary">
          <a className="brand" href="#top">
            <span className="brand-mark" />
            NightSkyPlan
          </a>
          <div className="nav-links">
            <a href="#scheduler">Scheduler</a>
            <a href="#targets">Targets</a>
            <a href="#outputs">Outputs</a>
          </div>
          <a className="nav-action" href="#scheduler">
            Open planner
          </a>
        </nav>

        <section className="hero-grid">
          <div className="hero-copy">
            <p className="eyebrow">Professional observatory scheduling</p>
            <h1>Turn transient targets into a night-ready observing queue.</h1>
            <p>
              Resolve TNS names, tune observatory constraints, run a deterministic scheduler,
              inspect failed targets, and export the final sequence for operations.
            </p>
            <div className="hero-actions">
              <button className="button primary" type="button" onClick={handleRunSchedule} disabled={isScheduling}>
                {isScheduling ? <Loader2 className="spin" size={17} /> : <Radar size={17} />}
                Run scheduler
              </button>
              <a className="button secondary" href="#targets">
                <Plus size={17} />
                Add targets
              </a>
            </div>
          </div>

          <aside className="mission-card" aria-label="Selected observatory">
            <div className="mission-card__top">
              <span>Observatory lock</span>
              <strong>{selectedObservatory?.name ?? "Loading observatories"}</strong>
            </div>
            <EarthGlobe observatory={selectedObservatory} />
            <dl className="mission-metrics">
              <div>
                <dt>Score</dt>
                <dd>{result?.score.plan_score ?? "--"}</dd>
              </div>
              <div>
                <dt>Scheduled</dt>
                <dd>{result ? `${result.score.scheduled_targets}/${result.score.total_targets}` : "--"}</dd>
              </div>
              <div>
                <dt>Dark</dt>
                <dd>{formatMinutes(result?.context.dark_window_min)}</dd>
              </div>
            </dl>
          </aside>
        </section>
      </header>

      <main>
        <section className="workspace" id="scheduler">
          <div className="section-heading">
            <span>Scheduler workspace</span>
            <h2>Single-night automatic queue builder</h2>
          </div>

          <div className="planner-grid">
            <section className="panel controls-panel" aria-label="Planner controls">
              <div className="panel-title">
                <SlidersHorizontal size={18} />
                <h3>Constraints</h3>
              </div>
              <div className="field-grid two">
                <label>
                  Observatory
                  <select value={observatoryId} onChange={event => setObservatoryId(event.target.value)}>
                    {observatories.map(observatory => (
                      <option key={observatory.id} value={observatory.id}>
                        {observatory.name}
                      </option>
                    ))}
                  </select>
                </label>
                <label>
                  Local date
                  <input type="date" value={date} onChange={event => setDate(event.target.value)} />
                </label>
                <label>
                  Cadence
                  <input
                    type="number"
                    min={1}
                    max={60}
                    value={cadenceMin}
                    onChange={event => setCadenceMin(Number(event.target.value))}
                  />
                </label>
                <label>
                  Overhead %
                  <input
                    type="number"
                    min={0}
                    max={200}
                    value={overheadPercent}
                    onChange={event => setOverheadPercent(Number(event.target.value))}
                  />
                </label>
              </div>

              <div className="slider-stack">
                <RangeField
                  label="Max airmass"
                  value={constraints.max_airmass}
                  min={1}
                  max={5}
                  step={0.1}
                  suffix=""
                  onChange={value => setConstraints(current => ({ ...current, max_airmass: value }))}
                />
                <RangeField
                  label="Moon separation"
                  value={constraints.min_moon_sep_deg}
                  min={0}
                  max={180}
                  step={1}
                  suffix=" deg"
                  onChange={value => setConstraints(current => ({ ...current, min_moon_sep_deg: value }))}
                />
                <RangeField
                  label="Hour angle"
                  value={constraints.ha_limit_hour}
                  min={0.5}
                  max={12}
                  step={0.5}
                  suffix=" h"
                  onChange={value => setConstraints(current => ({ ...current, ha_limit_hour: value }))}
                />
                <RangeField
                  label="Twilight Sun altitude"
                  value={constraints.twilight_alt_deg}
                  min={-24}
                  max={-1}
                  step={1}
                  suffix=" deg"
                  onChange={value => setConstraints(current => ({ ...current, twilight_alt_deg: value }))}
                />
              </div>

              <button className="button primary full" type="button" onClick={handleRunSchedule} disabled={isScheduling}>
                {isScheduling ? <Loader2 className="spin" size={17} /> : <CalendarClock size={17} />}
                {isScheduling ? "Computing plan" : "Run automatic scheduler"}
              </button>
              <p className="status-line">{status}</p>
            </section>

            <section className="panel output-panel" id="outputs" aria-label="Schedule output">
              <div className="panel-title split">
                <div>
                  <Activity size={18} />
                  <h3>Generated queue</h3>
                </div>
                <span>{result ? `${result.score.plan_score}/100` : "Not run"}</span>
              </div>
              <ScoreStrip result={result} />
              <Timeline result={result} />
              <div className="export-row">
                <button className="button secondary" type="button" onClick={() => handleExport("csv")} disabled={!result || !!isExporting}>
                  {isExporting === "csv" ? <Loader2 className="spin" size={16} /> : <Download size={16} />}
                  Export CSV
                </button>
                <button className="button secondary" type="button" onClick={() => handleExport("ics")} disabled={!result || !!isExporting}>
                  {isExporting === "ics" ? <Loader2 className="spin" size={16} /> : <Download size={16} />}
                  Export ICS
                </button>
              </div>
            </section>
          </div>
        </section>

        <section className="target-section" id="targets">
          <div className="section-heading">
            <span>Target workbench</span>
            <h2>TNS lookup and manual coordinates</h2>
          </div>

          <div className="target-grid">
            <section className="panel">
              <div className="panel-title">
                <Search size={18} />
                <h3>Resolve from TNS</h3>
              </div>
              <div className="inline-form">
                <input value={tnsQuery} onChange={event => setTnsQuery(event.target.value)} aria-label="TNS target name" />
                <button className="button primary" type="button" onClick={handleResolveTarget} disabled={isResolving}>
                  {isResolving ? <Loader2 className="spin" size={16} /> : <Search size={16} />}
                  Resolve
                </button>
              </div>
              <p className="hint">Uses server-side TNS credentials only; credentials never leave the backend.</p>
            </section>

            <section className="panel">
              <div className="panel-title">
                <Plus size={18} />
                <h3>Manual target</h3>
              </div>
              <div className="field-grid compact">
                <input
                  placeholder="Target name"
                  value={manualTarget.name}
                  onChange={event => setManualTarget(current => ({ ...current, name: event.target.value }))}
                />
                <input
                  placeholder="RA deg"
                  value={manualTarget.ra_deg}
                  onChange={event => setManualTarget(current => ({ ...current, ra_deg: event.target.value }))}
                />
                <input
                  placeholder="Dec deg"
                  value={manualTarget.dec_deg}
                  onChange={event => setManualTarget(current => ({ ...current, dec_deg: event.target.value }))}
                />
                <input
                  placeholder="Exposure min"
                  value={manualTarget.exposure_min}
                  onChange={event => setManualTarget(current => ({ ...current, exposure_min: event.target.value }))}
                />
                <input
                  placeholder="Priority"
                  value={manualTarget.priority}
                  onChange={event => setManualTarget(current => ({ ...current, priority: event.target.value }))}
                />
                <button className="button secondary" type="button" onClick={handleAddManualTarget}>
                  <Plus size={16} />
                  Add target
                </button>
              </div>
            </section>
          </div>

          <section className="panel queue-panel">
            <div className="panel-title split">
              <div>
                <Globe2 size={18} />
                <h3>Observing queue</h3>
              </div>
              <span>{targets.length} targets</span>
            </div>
            <div className="target-table" role="table" aria-label="Targets">
              <div className="target-row target-row--head" role="row">
                <span>Name</span>
                <span>RA</span>
                <span>Dec</span>
                <span>Exp</span>
                <span>Priority</span>
                <span />
              </div>
              {targets.map((target, index) => (
                <div className="target-row" role="row" key={`${target.name}-${index}`}>
                  <input value={target.name} onChange={event => updateTarget(index, "name", event.target.value)} aria-label="Target name" />
                  <input
                    type="number"
                    value={target.ra_deg}
                    onChange={event => updateTarget(index, "ra_deg", event.target.value)}
                    aria-label="RA degrees"
                  />
                  <input
                    type="number"
                    value={target.dec_deg}
                    onChange={event => updateTarget(index, "dec_deg", event.target.value)}
                    aria-label="Dec degrees"
                  />
                  <input
                    type="number"
                    min={1}
                    value={target.exposure_min}
                    onChange={event => updateTarget(index, "exposure_min", event.target.value)}
                    aria-label="Exposure minutes"
                  />
                  <input
                    type="number"
                    min={0}
                    max={10}
                    value={target.priority}
                    onChange={event => updateTarget(index, "priority", event.target.value)}
                    aria-label="Priority"
                  />
                  <button className="icon-button" type="button" onClick={() => removeTarget(index)} aria-label={`Remove ${target.name}`}>
                    <Trash2 size={16} />
                  </button>
                </div>
              ))}
            </div>
          </section>
        </section>

        <section className="diagnostics-section">
          <div className="section-heading">
            <span>Diagnostics</span>
            <h2>Know exactly why a target did or did not schedule</h2>
          </div>
          <div className="diagnostics-grid">
            <DiagnosticList title="Scheduled" blocks={scheduledBlocks} empty="Run the scheduler to populate scheduled targets." />
            <DiagnosticList title="Unscheduled" blocks={unscheduledBlocks} empty="No blocked targets in the current result." />
          </div>
        </section>
      </main>
    </div>
  );
}

function RangeField({
  label,
  value,
  min,
  max,
  step,
  suffix,
  onChange,
}: {
  label: string;
  value: number;
  min: number;
  max: number;
  step: number;
  suffix: string;
  onChange: (value: number) => void;
}) {
  return (
    <label className="range-field">
      <span>
        {label}
        <strong>
          {value}
          {suffix}
        </strong>
      </span>
      <input type="range" min={min} max={max} step={step} value={value} onChange={event => onChange(Number(event.target.value))} />
    </label>
  );
}

function ScoreStrip({ result }: { result: ScheduleResponse | null }) {
  const items = [
    ["Scheduled", result ? `${result.score.scheduled_targets}/${result.score.total_targets}` : "--"],
    ["Mean X", result?.score.mean_airmass?.toFixed(2) ?? "--"],
    ["Moon sep", result?.score.mean_moon_sep_deg ? `${Math.round(result.score.mean_moon_sep_deg)} deg` : "--"],
    ["Dark window", formatMinutes(result?.context.dark_window_min)],
  ];
  return (
    <dl className="score-strip">
      {items.map(([label, value]) => (
        <div key={label}>
          <dt>{label}</dt>
          <dd>{value}</dd>
        </div>
      ))}
    </dl>
  );
}

function Timeline({ result }: { result: ScheduleResponse | null }) {
  if (!result) {
    return <div className="empty-state">Run the scheduler to render the observing timeline.</div>;
  }
  const scheduled = result.schedule.filter(block => block.status === "SCHEDULED" && block.start_utc && block.end_utc);
  if (scheduled.length === 0) {
    return <div className="empty-state">No target clears the current constraints.</div>;
  }
  const bounds = scheduled.reduce(
    (range, block) => {
      const start = new Date(block.start_utc ?? "").getTime();
      const end = new Date(block.end_utc ?? "").getTime();
      return { start: Math.min(range.start, start), end: Math.max(range.end, end) };
    },
    { start: Number.POSITIVE_INFINITY, end: Number.NEGATIVE_INFINITY }
  );
  const span = Math.max(1, bounds.end - bounds.start);
  return (
    <div className="timeline">
      {scheduled.map(block => {
        const start = new Date(block.start_utc ?? "").getTime();
        const end = new Date(block.end_utc ?? "").getTime();
        const left = ((start - bounds.start) / span) * 100;
        const width = Math.max(8, ((end - start) / span) * 100);
        return (
          <article
            className="timeline-block"
            key={`${block.target}-${block.start_utc}`}
            style={{ "--left": `${left}%`, "--width": `${width}%` } as CSSProperties}
          >
            <strong>{block.target}</strong>
            <span>
              {formatClock(block.start_local)} - {formatClock(block.end_local)}
            </span>
          </article>
        );
      })}
    </div>
  );
}

function DiagnosticList({
  title,
  blocks,
  empty,
}: {
  title: string;
  blocks: Array<{
    target: string;
    status: string;
    start_local: string | null;
    end_local: string | null;
    reason: string;
    mean_airmass: number | null;
    mean_moon_sep_deg: number | null;
  }>;
  empty: string;
}) {
  return (
    <section className="panel">
      <div className="panel-title split">
        <h3>{title}</h3>
        <span>{blocks.length}</span>
      </div>
      <div className="diagnostic-list">
        {blocks.length === 0 ? (
          <p className="empty-inline">{empty}</p>
        ) : (
          blocks.map(block => (
            <article className="diagnostic-item" key={`${block.target}-${block.status}`}>
              <div>
                <strong>{block.target}</strong>
                <span>
                  {block.status === "SCHEDULED"
                    ? `${formatClock(block.start_local)} - ${formatClock(block.end_local)}`
                    : block.reason || "No valid schedule slot"}
                </span>
              </div>
              <small>
                X {block.mean_airmass?.toFixed(2) ?? "--"} / Moon{" "}
                {block.mean_moon_sep_deg ? `${Math.round(block.mean_moon_sep_deg)} deg` : "--"}
              </small>
            </article>
          ))
        )}
      </div>
    </section>
  );
}

function formatMinutes(value: number | null | undefined) {
  if (!value || value <= 0) return "--";
  const hours = Math.floor(value / 60);
  const minutes = Math.round(value % 60);
  return `${hours}h ${minutes.toString().padStart(2, "0")}m`;
}

function formatClock(value: string | null | undefined) {
  if (!value) return "--:--";
  return new Intl.DateTimeFormat(undefined, { hour: "2-digit", minute: "2-digit" }).format(new Date(value));
}

export default App;
