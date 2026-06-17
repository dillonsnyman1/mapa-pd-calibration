"""Thin wrapper around the validated reference/python/mapa.py implementation."""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "reference" / "python"))

from mapa import CalibratedBin, run_pipeline, interpolate_pd  # noqa: E402

FIXTURES_DIR = Path(__file__).resolve().parent.parent.parent / "reference" / "fixtures"

Observation = Tuple[float, int, float]


def load_example_observations() -> List[Observation]:
    with open(FIXTURES_DIR / "raw_observations.csv", newline="") as f:
        reader = csv.DictReader(f)
        return [(float(row["score"]), int(row["bad"]), 1.0) for row in reader]


def load_example_observations_weighted() -> List[Observation]:
    with open(FIXTURES_DIR / "raw_observations_weighted.csv", newline="") as f:
        reader = csv.DictReader(f)
        return [(float(row["score"]), int(row["bad"]), float(row["weight"])) for row in reader]


def compute_smoothed(bands: List[CalibratedBin], num_points: int = 200) -> list[tuple[float, float]]:
    score_min = bands[0].score_min
    score_max = bands[-1].score_max
    if score_max > score_min:
        step = (score_max - score_min) / (num_points - 1)
        scores = [score_min + i * step for i in range(num_points)]
    else:
        scores = [score_min]
    return [(s, interpolate_pd(bands, s)) for s in scores]


def run_calibration(
    observations: List[Observation],
    min_obs: float,
    min_bads: float,
    k: float,
    min_confidence: float | None,
    increasing: bool,
    use_counts: bool = True,
    num_smoothed_points: int = 200,
) -> tuple[list[CalibratedBin], list[tuple[float, float]]]:
    result = run_pipeline(
        observations,
        k=k,
        min_obs=min_obs,
        min_bads=min_bads,
        increasing=increasing,
        min_confidence=min_confidence,
        use_counts=use_counts,
    )

    bands = result.bands
    smoothed = compute_smoothed(bands, num_smoothed_points)

    return bands, smoothed
