"""
NightSkyPlan 2.0 prototype: interactive visibility + hour-angle + greedy scheduler.

Run:
    pip install streamlit astropy pandas numpy plotly
    streamlit run nightskyplan_streamlit.py

Target CSV columns:
    name,ra,dec,exposure_min,priority
RA examples: 12:34:56.7 or 188.736 deg
Dec examples: -12:34:56 or -12.582 deg
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, time, timedelta
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time


@dataclass(frozen=True)
class Observatory:
    name: str
    latitude_deg: float
    longitude_deg: float
    elevation_m: float
    timezone: str
    horizon_deg: float = 25.0
    zenith_deg: float = 88.0

    @property
    def location(self) -> EarthLocation:
        return EarthLocation(
            lat=self.latitude_deg * u.deg,
            lon=self.longitude_deg * u.deg,
            height=self.elevation_m * u.m,
        )


def parse_angle(value: str, is_ra: bool) -> u.Quantity:
    """Parse decimal degrees or sexagesimal RA/Dec."""
    s = str(value).strip()
    if ":" in s or " " in s:
        return SkyCoord(s, "0d", unit=(u.hourangle, u.deg)).ra if is_ra else SkyCoord("0h", s, unit=(u.hourangle, u.deg)).dec
    return float(s) * u.deg


def parse_targets(df: pd.DataFrame) -> pd.DataFrame:
    required = {"name", "ra", "dec"}
    missing = required - set(df.columns.str.lower())
    df = df.rename(columns={c: c.lower() for c in df.columns})
    if missing:
        raise ValueError(f"Missing target columns: {sorted(missing)}")
    if "exposure_min" not in df:
        df["exposure_min"] = 30.0
    if "priority" not in df:
        df["priority"] = 1.0
    df["coord"] = [
        SkyCoord(parse_angle(ra, True), parse_angle(dec, False), frame="icrs")
        for ra, dec in zip(df["ra"], df["dec"])
    ]
    df["exposure_min"] = pd.to_numeric(df["exposure_min"], errors="coerce").fillna(30.0)
    df["priority"] = pd.to_numeric(df["priority"], errors="coerce").fillna(1.0)
    return df


def make_time_grid(obs_date: date, timezone: str, cadence_min: int = 5) -> tuple[Time, pd.DatetimeIndex]:
    """Grid from local noon to next local noon, covering the whole observing night."""
    tz = ZoneInfo(timezone)
    start_local = datetime.combine(obs_date, time(12, 0), tzinfo=tz)
    end_local = start_local + timedelta(days=1)
    local_index = pd.date_range(start_local, end_local, freq=f"{cadence_min}min", inclusive="both")
    return Time(local_index.to_pydatetime()), local_index


def crossing_time(times: Time, y: np.ndarray, threshold: float, direction: str) -> Time | None:
    z = y - threshold
    if direction == "down":
        idx = np.where((z[:-1] > 0) & (z[1:] <= 0))[0]
    else:
        idx = np.where((z[:-1] < 0) & (z[1:] >= 0))[0]
    if len(idx) == 0:
        return None
    i = idx[0] if direction == "down" else idx[-1]
    frac = abs(z[i]) / (abs(z[i]) + abs(z[i + 1]))
    return times[i] + frac * (times[i + 1] - times[i])


def solar_lunar_context(times: Time, obs: Observatory, twilight_alt_deg: float) -> dict[str, object]:
    frame = AltAz(obstime=times, location=obs.location)
    sun_alt = get_sun(times).transform_to(frame).alt.deg
    moon = get_body("moon", times, obs.location)
    moon_alt = moon.transform_to(frame).alt.deg
    moon_illum_proxy = None  # intentionally not faked; add astroplan or pyephem for true illumination if required
    dusk = crossing_time(times, sun_alt, twilight_alt_deg, "down")
    dawn = crossing_time(times, sun_alt, twilight_alt_deg, "up")
    return {"sun_alt": sun_alt, "moon": moon, "moon_alt": moon_alt, "dusk": dusk, "dawn": dawn, "moon_illum_proxy": moon_illum_proxy}


def target_track(target: pd.Series, times: Time, obs: Observatory, context: dict[str, object]) -> pd.DataFrame:
    coord: SkyCoord = target["coord"]
    frame = AltAz(obstime=times, location=obs.location)
    altaz = coord.transform_to(frame)
    alt = altaz.alt.deg
    az = altaz.az.deg
    # Good enough for screening. For low altitude use a proper airmass model or reject via horizon.
    airmass = np.where(alt > 0, 1.0 / np.cos(np.deg2rad(90.0 - alt)), np.nan)
    lst = times.sidereal_time("apparent", longitude=obs.location.lon)
    ha = (lst - coord.ra).wrap_at(12 * u.hourangle).hour
    moon_sep = coord.separation(context["moon"]).deg
    return pd.DataFrame(
        {
            "name": target["name"],
            "time": times.to_datetime(),
            "alt_deg": alt,
            "az_deg": az,
            "airmass": airmass,
            "ha_hour": ha,
            "moon_sep_deg": moon_sep,
            "priority": float(target["priority"]),
            "exposure_min": float(target["exposure_min"]),
        }
    )


def compute_tracks(targets: pd.DataFrame, times: Time, obs: Observatory, context: dict[str, object]) -> pd.DataFrame:
    return pd.concat([target_track(row, times, obs, context) for _, row in targets.iterrows()], ignore_index=True)


def apply_constraints(
    tracks: pd.DataFrame,
    times: Time,
    context: dict[str, object],
    obs: Observatory,
    twilight_alt_deg: float,
    max_airmass: float,
    min_moon_sep_deg: float,
    ha_limit_hour: float,
) -> pd.DataFrame:
    sun_alt_map = pd.Series(context["sun_alt"], index=times.to_datetime())
    tracks = tracks.copy()
    tracks["sun_alt_deg"] = tracks["time"].map(sun_alt_map)
    tracks["valid"] = (
        (tracks["sun_alt_deg"] < twilight_alt_deg)
        & (tracks["alt_deg"] >= obs.horizon_deg)
        & (tracks["alt_deg"] <= obs.zenith_deg)
        & (tracks["airmass"] <= max_airmass)
        & (tracks["moon_sep_deg"] >= min_moon_sep_deg)
        & (tracks["ha_hour"].abs() <= ha_limit_hour)
    )
    return tracks


def add_twilight_shapes(fig: go.Figure, context: dict[str, object], y0: float, y1: float) -> None:
    dusk, dawn = context.get("dusk"), context.get("dawn")
    if dusk is not None and dawn is not None:
        fig.add_vrect(x0=dusk.datetime, x1=dawn.datetime, fillcolor="gray", opacity=0.12, line_width=0)
        fig.add_vline(x=dusk.datetime, line_dash="dash", annotation_text="twilight start")
        fig.add_vline(x=dawn.datetime, line_dash="dash", annotation_text="twilight end")
    fig.update_yaxes(range=[y0, y1])


def visibility_figure(tracks: pd.DataFrame, context: dict[str, object], selected: list[str], obs: Observatory) -> go.Figure:
    fig = go.Figure()
    for name in selected:
        sub = tracks[tracks["name"] == name]
        fig.add_trace(go.Scatter(
            x=sub["time"], y=sub["alt_deg"], mode="lines", name=name,
            customdata=np.c_[sub["airmass"], sub["ha_hour"], sub["moon_sep_deg"], sub["valid"]],
            hovertemplate="%{x}<br>alt=%{y:.1f} deg<br>X=%{customdata[0]:.2f}<br>HA=%{customdata[1]:.2f} h<br>Moon sep=%{customdata[2]:.0f} deg<br>valid=%{customdata[3]}<extra>%{fullData.name}</extra>",
        ))
    add_twilight_shapes(fig, context, 0, 90)
    fig.add_hline(y=obs.horizon_deg, line_dash="dot", annotation_text="horizon")
    fig.add_hline(y=obs.zenith_deg, line_dash="dot", annotation_text="zenith limit")
    fig.update_layout(title="Visibility", xaxis_title="Time", yaxis_title="Altitude [deg]", hovermode="x unified")
    return fig


def hour_angle_figure(tracks: pd.DataFrame, context: dict[str, object], selected: list[str], ha_limit_hour: float) -> go.Figure:
    fig = go.Figure()
    for name in selected:
        sub = tracks[tracks["name"] == name]
        fig.add_trace(go.Scatter(x=sub["time"], y=sub["ha_hour"], mode="lines", name=name))
    fig.add_hline(y=0, line_dash="solid", annotation_text="meridian")
    fig.add_hline(y=ha_limit_hour, line_dash="dot", annotation_text="+HA limit")
    fig.add_hline(y=-ha_limit_hour, line_dash="dot", annotation_text="-HA limit")
    add_twilight_shapes(fig, context, -12, 12)
    fig.update_layout(title="Hour angle", xaxis_title="Time", yaxis_title="HA [hour]", hovermode="x unified")
    return fig


def greedy_schedule(tracks: pd.DataFrame, cadence_min: int, overhead_percent: float) -> pd.DataFrame:
    """Simple deterministic scheduler. Not globally optimal; good baseline for interactive triage."""
    times = np.array(sorted(tracks["time"].unique()))
    occupied = pd.Series(False, index=times)
    out = []
    targets = tracks.groupby("name").first()[["exposure_min", "priority"]].sort_values("priority", ascending=False)
    for name, meta in targets.iterrows():
        need_min = float(meta.exposure_min) * (1 + overhead_percent / 100.0)
        nslot = max(1, int(np.ceil(need_min / cadence_min)))
        sub = tracks[tracks["name"] == name].set_index("time").reindex(times)
        valid = sub["valid"].fillna(False).to_numpy() & (~occupied.to_numpy())
        best = None
        for i in range(0, len(times) - nslot + 1):
            if valid[i : i + nslot].all():
                window = sub.iloc[i : i + nslot]
                # Science-neutral score: priority dominates, altitude breaks ties, moon avoidance small bonus.
                score = float(meta.priority) * 1000 + window["alt_deg"].mean() + 0.05 * window["moon_sep_deg"].mean()
                if best is None or score > best[0]:
                    best = (score, i, window)
        if best is None:
            out.append({"target": name, "status": "UNSCHEDULED", "start": None, "end": None, "duration_min_with_overhead": need_min})
            continue
        _, i, window = best
        occupied.iloc[i : i + nslot] = True
        out.append({
            "target": name,
            "status": "SCHEDULED",
            "start": times[i],
            "end": times[min(i + nslot, len(times) - 1)],
            "duration_min_with_overhead": need_min,
            "mean_alt_deg": window["alt_deg"].mean(),
            "mean_airmass": window["airmass"].mean(),
            "mean_moon_sep_deg": window["moon_sep_deg"].mean(),
        })
    return pd.DataFrame(out)


def schedule_figure(schedule: pd.DataFrame) -> go.Figure:
    ok = schedule[schedule["status"] == "SCHEDULED"].copy()
    fig = go.Figure()
    for _, row in ok.iterrows():
        fig.add_trace(go.Bar(
            x=[(pd.Timestamp(row["end"]) - pd.Timestamp(row["start"])).total_seconds() / 60],
            y=[row["target"]],
            base=[row["start"]],
            orientation="h",
            name=row["target"],
            hovertemplate="%{y}<br>%{base} + %{x:.0f} min<extra></extra>",
        ))
    fig.update_layout(title="Greedy observing sequence", xaxis_title="Time", yaxis_title="Target", showlegend=False)
    return fig


EXAMPLE_TARGETS = pd.DataFrame(
    {
        "name": ["SN2022jli", "SN2020tlf", "SN2018zd"],
        "ra": ["00:34:45.2", "14:40:10.03", "06:18:03.18"],
        "dec": ["-08:23:25.0", "+42:46:39.6", "+78:22:00.9"],
        "exposure_min": [45, 30, 20],
        "priority": [3, 2, 1],
    }
)


def main() -> None:
    st.set_page_config(page_title="NightSkyPlan 2.0", layout="wide")
    st.title("NightSkyPlan 2.0: visibility, hour angle, and scheduling")

    with st.sidebar:
        st.header("Observatory")
        obs = Observatory(
            name=st.text_input("Name", "HCT / Hanle"),
            latitude_deg=st.number_input("Latitude [deg]", value=32.7794, format="%.6f"),
            longitude_deg=st.number_input("Longitude east [deg]", value=78.9642, format="%.6f"),
            elevation_m=st.number_input("Elevation [m]", value=4486.0),
            timezone=st.text_input("IANA timezone", "Asia/Kolkata"),
            horizon_deg=st.number_input("Telescope horizon [deg]", value=25.0),
            zenith_deg=st.number_input("Zenith avoidance limit [deg]", value=85.0),
        )
        obs_date = st.date_input("Local observing date", value=date.today())
        cadence_min = st.slider("Cadence [min]", 1, 30, 5)
        twilight_alt = st.selectbox("Night definition", [-18.0, -12.0, -6.0], index=0, format_func=lambda x: f"Sun < {x:.0f} deg")
        max_airmass = st.slider("Max airmass", 1.0, 5.0, 2.5, 0.1)
        min_moon_sep = st.slider("Minimum Moon separation [deg]", 0, 180, 30)
        ha_limit = st.slider("Hour-angle limit [h]", 0.5, 12.0, 6.0, 0.5)
        overhead = st.slider("Overheads [% of exposure]", 0, 100, 20)

    upload = st.file_uploader("Upload target CSV", type=["csv"])
    targets_raw = pd.read_csv(upload) if upload else EXAMPLE_TARGETS.copy()
    with st.expander("Targets", expanded=True):
        edited = st.data_editor(targets_raw, num_rows="dynamic", use_container_width=True)

    try:
        targets = parse_targets(edited)
        times, local_index = make_time_grid(obs_date, obs.timezone, cadence_min)
        context = solar_lunar_context(times, obs, twilight_alt)
        tracks = compute_tracks(targets, times, obs, context)
        tracks = apply_constraints(tracks, times, context, obs, twilight_alt, max_airmass, min_moon_sep, ha_limit)
    except Exception as exc:
        st.error(f"Could not compute plan: {exc}")
        st.stop()

    names = targets["name"].tolist()
    selected = st.multiselect("Targets to show", names, default=names[: min(6, len(names))])
    tab1, tab2, tab3, tab4 = st.tabs(["Visibility", "Hour angle", "Scheduler", "Raw table"])

    with tab1:
        st.plotly_chart(visibility_figure(tracks, context, selected, obs), use_container_width=True)
    with tab2:
        st.plotly_chart(hour_angle_figure(tracks, context, selected, ha_limit), use_container_width=True)
    with tab3:
        schedule = greedy_schedule(tracks, cadence_min, overhead)
        st.plotly_chart(schedule_figure(schedule), use_container_width=True)
        st.dataframe(schedule, use_container_width=True)
        st.download_button("Download schedule CSV", schedule.to_csv(index=False), "schedule.csv", "text/csv")
    with tab4:
        st.dataframe(tracks.drop(columns=[]), use_container_width=True)
        st.download_button("Download visibility table CSV", tracks.to_csv(index=False), "visibility_table.csv", "text/csv")


if __name__ == "__main__":
    main()
