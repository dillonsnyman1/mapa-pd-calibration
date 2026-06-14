"""Plot the unsmoothed (step) and smoothed (log-odds interpolated) PD
mappings produced by the Python reference's run_pipeline, using the shared
reference fixtures.

Requires matplotlib (see requirements.txt) - the reference implementations
themselves have no external dependencies; this is a demo-only script.

Usage:
    pip install -r requirements.txt
    python3 example_plot_python.py
"""

import csv
import sys
from pathlib import Path

from plotting import plot_calibration

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "reference" / "python"))
from mapa import run_pipeline, interpolate_pd

FIXTURES_DIR = Path(__file__).resolve().parent.parent / "reference" / "fixtures"


def load_observations():
    with open(FIXTURES_DIR / "raw_observations.csv", newline="") as f:
        reader = csv.DictReader(f)
        return [(float(row["score"]), int(row["bad"])) for row in reader]


def main():
    observations = load_observations()
    result = run_pipeline(observations, k=10, min_obs=50, min_bads=10)

    bands = [(b.score_min, b.score_max, b.pd) for b in result.bands]

    score_min = result.bands[0].score_min
    score_max = result.bands[-1].score_max
    smoothed = [
        (s, interpolate_pd(result.bands, s))
        for s in (score_min + i * (score_max - score_min) / 500 for i in range(501))
    ]

    out_path = Path(__file__).resolve().parent / "output" / "python_calibration.png"
    plot_calibration(bands, smoothed, "MAPA score-to-PD calibration (Python)", out_path)


if __name__ == "__main__":
    main()
