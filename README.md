# NightSkyPlan

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Tests](https://github.com/sPaMFouR/NightSkyPlan/actions/workflows/tests.yml/badge.svg?branch=web-dev)](https://github.com/sPaMFouR/NightSkyPlan/actions/workflows/tests.yml)
[![FastAPI](https://img.shields.io/badge/API-FastAPI-009688?logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![React](https://img.shields.io/badge/web-React%20%2B%20TypeScript-61DAFB?logo=react&logoColor=0B1020)](https://react.dev/)
[![Branch](https://img.shields.io/badge/branch-web--dev-69E1FF)](https://github.com/sPaMFouR/NightSkyPlan/tree/web-dev)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-2F81F7)](LICENSE.txt)

NightSkyPlan helps plan ground-based night sky observations. It combines the
original observing scripts with a reusable Python planning package and a
professional FastAPI + React webapp for observatory scheduling.

## Repository Layout

```text
nightskyplan/   Core planning package: targets, ephemerides, constraints, plots, scheduling
webapp/         FastAPI backend and React/TypeScript professional scheduler
application/    Streamlit observing planner
webpage/        Static standalone page kept separate from the webapp
tests/          Package tests
*.py            Legacy command-line and GUI planning scripts
*.dat           Example observatory, target, and date inputs
```

Generated files such as `*.egg-info/`, `__pycache__/`, `*.pyc`, build outputs,
and `.DS_Store` are ignored and should not be committed.

## Professional Webapp

The active web version lives in `webapp/`. It is built for professional
observatories that need target resolution, constraints, single-night automatic
scheduling, diagnostics, and CSV/ICS exports.

Install Python dependencies from the repository root:

```bash
python3 -m pip install --upgrade pip
python3 -m pip install -e .
```

Set TNS credentials for target-name lookup:

```bash
export TNS_API_KEY="..."
export TNS_BOT_ID="..."
export TNS_BOT_NAME="..."
```

Run the FastAPI backend:

```bash
uvicorn webapp.backend.app:app --reload --host 127.0.0.1 --port 8000
```

Run the React frontend:

```bash
cd webapp/frontend
npm install
npm run dev
```

Open `http://127.0.0.1:5173`. The Vite dev server proxies `/api` requests to
the FastAPI backend.

The V1 webapp includes:

- TNS target lookup with server-side credentials only
- manual RA/Dec target entry
- observatory presets for Hanle, Devasthal, Kanata/Hiroshima, La Palma, Mauna
  Kea, Paranal, La Silla, and Palomar
- airmass, Moon separation, hour-angle, twilight, cadence, and overhead controls
- deterministic single-night scheduling with diagnostics
- CSV and ICS export endpoints

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

## Static Webpage

The static standalone page lives in `webpage/index.html`. It remains separate
from the FastAPI + React webapp and does not use the Streamlit application.

Open it directly in a browser:

```bash
open webpage/index.html
```

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

Run frontend checks:

```bash
cd webapp/frontend
npm test
npm run build
```

The active web branch is `web-dev`. Keep the professional webapp in `webapp/`,
the Streamlit planner code in `application/`, and generated metadata or build
outputs out of version control.
