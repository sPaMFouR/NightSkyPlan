# NightSkyPlan Professional Webapp

This directory contains the active web application:

- `backend/`: FastAPI service that wraps the `nightskyplan` astronomy package.
- `frontend/`: React + TypeScript planner UI built with Vite.

## Run

From the repository root:

```bash
python3 -m pip install --upgrade pip
python3 -m pip install -e .
uvicorn webapp.backend.app:app --reload --host 127.0.0.1 --port 8000
```

In a second terminal:

```bash
cd webapp/frontend
npm install
npm run dev
```

Open `http://127.0.0.1:5173`.

## TNS

Target lookup uses server-side TNS credentials only:

```bash
export TNS_API_KEY="..."
export TNS_BOT_ID="..."
export TNS_BOT_NAME="..."
```

If credentials are missing, `/api/targets/resolve` returns a setup error and
manual RA/Dec target entry remains available.
