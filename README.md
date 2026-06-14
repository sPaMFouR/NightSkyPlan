# NightSkyPlan

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Tests](https://github.com/sPaMFouR/NightSkyPlan/actions/workflows/tests.yml/badge.svg?branch=web-dev)](https://github.com/sPaMFouR/NightSkyPlan/actions/workflows/tests.yml)
[![Streamlit app](https://img.shields.io/badge/app-Streamlit-FF4B4B?logo=streamlit&logoColor=white)](https://streamlit.io/)
[![Plotly charts](https://img.shields.io/badge/charts-Plotly-3F4F75?logo=plotly&logoColor=white)](https://plotly.com/python/)
[![Branch](https://img.shields.io/badge/branch-web--dev-69E1FF)](https://github.com/sPaMFouR/NightSkyPlan/tree/web-dev)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-2F81F7)](LICENSE.txt)

NightSkyPlan helps plan ground-based night sky observations. It combines the
original observing scripts with a reusable Python planning package, a Streamlit
planner application, and a portfolio-ready webpage for the web experience.

## Repository Layout

```text
nightskyplan/   Core planning package: targets, ephemerides, constraints, plots, scheduling
application/    Streamlit observing planner
webpage/        Static portfolio landing page for the web experience
tests/          Package tests
*.py            Legacy command-line and GUI planning scripts
*.dat           Example observatory, target, and date inputs
```

Generated files such as `*.egg-info/`, `__pycache__/`, `*.pyc`, build outputs,
and `.DS_Store` are ignored and should not be committed.

## Streamlit Planner

The interactive planner lives in `application/app.py`. It uses the
`nightskyplan` package for target parsing, ephemerides, twilight context,
constraints, Plotly charts, and baseline scheduling.

Install the planner dependencies:

```bash
python3 -m pip install -r application/requirements.txt
```

Run the planner from the repository root:

```bash
streamlit run application/app.py
```

Target uploads should be CSV files with these columns:

```text
name,ra,dec,exposure_min,priority
```

RA and Dec can be sexagesimal values such as `12:34:56.7`, `-12:34:56`, or
decimal degrees. The planner shows visibility, airmass, hour angle, Moon
separation, a greedy observing schedule, diagnostics, and downloadable CSV
tables.

## Webpage

The modern static webpage lives in `webpage/index.html`. It is an independent
browser-native planner with a project-local hero image, responsive CSS, native
canvas motion, an illuminated observatory globe, and JavaScript schedule
calculation.

Open it directly in a browser:

```bash
open webpage/index.html
```

The webpage runs directly in the browser and does not need the Streamlit app.

## Legacy Scripts

The original scripts remain available for PDF and table workflows:

| Script | Purpose | Outputs |
| --- | --- | --- |
| `NightSkyPlan.py` | Nightly target visibility from a selected observatory | `NightSkyPlan_DATE.pdf` |
| `YearlyPlan.py` | Long-range target observability across a proposal cycle | `YearlyPlan_StartDATEToEndDATE.pdf` |
| `CalcTwilightTime.py` | Sunset, sunrise, twilight, night duration, and Moon phase | `NightDuration_StartDATEToEndDATE.pdf`, `TwilightTimes_StartDATEToEndDATE.asc` |
| `CalcMoonAnglePhase.py` | Moon phase and Moon separation for listed dates | `MoonPhaseAngle.asc` |

Input examples are provided in `TelescopeList.dat`, `TargetList.dat`, and
`DateList.dat`.

## Development

Run the package tests from the repository root:

```bash
python3 -m pytest
```

The active web branch is `web-dev`. Keep the Streamlit planner code in
`application/`, the static landing page in `webpage/`, and generated Python
metadata out of version control.
