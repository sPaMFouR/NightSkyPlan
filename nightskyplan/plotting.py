from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from nightskyplan.models import Observatory


def add_twilight_shapes(fig: go.Figure, context: dict[str, object], y0: float, y1: float) -> None:
    dusk, dawn = context.get("dusk"), context.get("dawn")
    if dusk is not None and dawn is not None:
        dusk_dt = dusk.datetime
        dawn_dt = dawn.datetime
        fig.add_vrect(x0=dusk_dt, x1=dawn_dt, fillcolor="gray", opacity=0.12, line_width=0)
        fig.add_vline(x=dusk_dt, line_dash="dash")
        fig.add_vline(x=dawn_dt, line_dash="dash")
        fig.add_annotation(x=dusk_dt, y=y1, text="twilight start", showarrow=False, yanchor="bottom")
        fig.add_annotation(x=dawn_dt, y=y1, text="twilight end", showarrow=False, yanchor="bottom")
    fig.update_yaxes(range=[y0, y1])


def visibility_figure(tracks: pd.DataFrame, context: dict[str, object], selected: list[str], obs: Observatory) -> go.Figure:
    fig = go.Figure()
    for name in selected:
        sub = tracks[tracks["name"] == name]
        fig.add_trace(
            go.Scatter(
                x=sub["time"],
                y=sub["alt_deg"],
                mode="lines",
                name=name,
                customdata=np.c_[sub["airmass"], sub["ha_hour"], sub["moon_sep_deg"], sub["valid"]],
                hovertemplate="%{x}<br>alt=%{y:.1f} deg<br>X=%{customdata[0]:.2f}<br>HA=%{customdata[1]:.2f} h<br>Moon sep=%{customdata[2]:.0f} deg<br>valid=%{customdata[3]}<extra>%{fullData.name}</extra>",
            )
        )
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


def moon_separation_figure(
    tracks: pd.DataFrame,
    context: dict[str, object],
    selected: list[str],
    min_moon_sep_deg: float,
) -> go.Figure:
    fig = go.Figure()
    for name in selected:
        sub = tracks[tracks["name"] == name]
        fig.add_trace(go.Scatter(x=sub["time"], y=sub["moon_sep_deg"], mode="lines", name=name))
    fig.add_hline(y=min_moon_sep_deg, line_dash="dot", annotation_text="minimum separation")
    add_twilight_shapes(fig, context, 0, 180)
    fig.update_layout(title="Moon separation", xaxis_title="Time", yaxis_title="Separation [deg]", hovermode="x unified")
    return fig


def airmass_figure(tracks: pd.DataFrame, context: dict[str, object], selected: list[str], max_airmass: float) -> go.Figure:
    fig = go.Figure()
    for name in selected:
        sub = tracks[tracks["name"] == name]
        fig.add_trace(go.Scatter(x=sub["time"], y=sub["airmass"], mode="lines", name=name))
    fig.add_hline(y=max_airmass, line_dash="dot", annotation_text="max airmass")
    add_twilight_shapes(fig, context, 1, max(5, max_airmass + 0.5))
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(title="Airmass", xaxis_title="Time", yaxis_title="Airmass", hovermode="x unified")
    return fig


def schedule_figure(schedule: pd.DataFrame) -> go.Figure:
    ok = schedule[schedule["status"] == "SCHEDULED"].copy()
    fig = go.Figure()
    for _, row in ok.iterrows():
        fig.add_trace(
            go.Bar(
                x=[(pd.Timestamp(row["end"]) - pd.Timestamp(row["start"])).total_seconds() / 60],
                y=[row["target"]],
                base=[row["start"]],
                orientation="h",
                name=row["target"],
                hovertemplate="%{y}<br>%{base} + %{x:.0f} min<extra></extra>",
            )
        )
    fig.update_layout(title="Greedy observing sequence", xaxis_title="Time", yaxis_title="Target", showlegend=False)
    return fig
