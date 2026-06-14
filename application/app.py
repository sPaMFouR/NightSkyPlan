from __future__ import annotations

from datetime import date

import pandas as pd
import streamlit as st

from nightskyplan.constraints import apply_constraints, compute_tracks
from nightskyplan.ephemeris import make_time_grid, solar_lunar_context
from nightskyplan.models import Observatory
from nightskyplan.observatory import DEFAULT_OBSERVATORY
from nightskyplan.plotting import airmass_figure, hour_angle_figure, moon_separation_figure, schedule_figure, visibility_figure
from nightskyplan.scheduling import build_greedy_plan
from nightskyplan.targets import EXAMPLE_TARGETS, parse_targets


def main() -> None:
    st.set_page_config(page_title="NightSkyPlan", layout="wide")
    st.title("NightSkyPlan: visibility and observing schedule")

    with st.sidebar:
        st.header("Observatory")
        obs = Observatory(
            name=st.text_input("Name", DEFAULT_OBSERVATORY.name),
            latitude_deg=st.number_input("Latitude [deg]", value=DEFAULT_OBSERVATORY.latitude_deg, format="%.6f"),
            longitude_deg=st.number_input("Longitude east [deg]", value=DEFAULT_OBSERVATORY.longitude_deg, format="%.6f"),
            elevation_m=st.number_input("Elevation [m]", value=DEFAULT_OBSERVATORY.elevation_m),
            timezone=st.text_input("IANA timezone", DEFAULT_OBSERVATORY.timezone),
            horizon_deg=st.number_input("Telescope horizon [deg]", value=DEFAULT_OBSERVATORY.horizon_deg),
            zenith_deg=st.number_input("Zenith avoidance limit [deg]", value=DEFAULT_OBSERVATORY.zenith_deg),
        )
        obs_date = st.date_input("Local observing date", value=date.today())
        cadence_min = st.slider("Cadence [min]", 1, 30, 5)
        twilight_alt = st.selectbox(
            "Night definition",
            [-18.0, -12.0, -6.0],
            index=0,
            format_func=lambda x: f"Sun < {x:.0f} deg",
        )
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
        times, _ = make_time_grid(obs_date, obs.timezone, cadence_min)
        context = solar_lunar_context(times, obs, twilight_alt)
        tracks = compute_tracks(targets, times, obs, context)
        tracks = apply_constraints(tracks, times, context, obs, twilight_alt, max_airmass, min_moon_sep, ha_limit)
        schedule, diagnostics = build_greedy_plan(tracks, cadence_min, overhead)
    except Exception as exc:
        st.error(f"Could not compute plan: {exc}")
        st.stop()

    names = targets["name"].tolist()
    selected = st.multiselect("Targets to show", names, default=names[: min(6, len(names))])
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(
        ["Visibility", "Airmass", "Hour angle", "Moon separation", "Schedule", "Raw table"]
    )

    with tab1:
        st.plotly_chart(visibility_figure(tracks, context, selected, obs), use_container_width=True)
    with tab2:
        st.plotly_chart(airmass_figure(tracks, context, selected, max_airmass), use_container_width=True)
    with tab3:
        st.plotly_chart(hour_angle_figure(tracks, context, selected, ha_limit), use_container_width=True)
    with tab4:
        st.plotly_chart(moon_separation_figure(tracks, context, selected, min_moon_sep), use_container_width=True)
    with tab5:
        st.plotly_chart(schedule_figure(schedule), use_container_width=True)
        st.dataframe(schedule, use_container_width=True)
        st.download_button("Download schedule CSV", schedule.to_csv(index=False), "schedule.csv", "text/csv")
        st.subheader("Target diagnostics")
        st.dataframe(diagnostics, use_container_width=True)
        st.download_button(
            "Download diagnostics CSV",
            diagnostics.to_csv(index=False),
            "failed_target_diagnostics.csv",
            "text/csv",
        )
    with tab6:
        st.dataframe(tracks, use_container_width=True)
        st.download_button("Download visibility table CSV", tracks.to_csv(index=False), "visibility_table.csv", "text/csv")


if __name__ == "__main__":
    main()
