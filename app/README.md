# Live demo

An interactive FastAPI + React app for exploring MAPA: view the resulting
calibration curve, step through the pooling process band by band, tune the
parameters, and try it on your own data. The backend is a thin wrapper around
[`reference/python/mapa.py`](../reference/python/mapa.py) - no part of the
algorithm is reimplemented here.

This runs locally only; there is no hosted deployment.

## Backend

```bash
cd app/backend
pip install -r requirements.txt
uvicorn main:app --reload
```

Runs on `http://localhost:8000`. Routes:

- `GET /api/example` - the bundled example observations (the same dataset
  used by [`examples/`](../examples)).
- `POST /api/calibrate` - runs `run_pipeline` and returns the calibrated
  bands plus a smoothed score-to-PD curve.
- `POST /api/pooling-steps` - runs the initial pooling pass step by step,
  for the pooling-process visualization.

## Frontend

```bash
cd app/frontend
npm install
npm run dev
```

Runs on `http://localhost:5173` and talks to the backend at
`http://localhost:8000` (CORS is enabled for this origin in `main.py`).

## Usage

- The demo loads with the bundled example dataset (2846 observations) and
  the same default parameters as the reference fixtures (`min_obs=50`,
  `min_bads=10`, `k=10`).
- Adjust `min_obs`, `min_bads`, `k`, the monotonicity direction, or enable
  confidence-based pooling - the calibration chart and pooling stepper
  update accordingly.
- Upload a CSV with `score,bad` columns to run the pipeline on your own
  data instead.
