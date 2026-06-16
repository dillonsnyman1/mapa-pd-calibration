"""FastAPI backend for the MAPA interactive demo."""

from __future__ import annotations

import sys
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "reference" / "python"))

from mapa import (  # noqa: E402
    Bin,
    apply_bayesian_adjustment,
    bins_from_observations,
    interpolate_pd,
    mapa,
    repool_calibrated_bins,
)

from pipeline import compute_smoothed, load_example_observations, run_calibration
from pooling_steps import minimum_size_steps, pooling_steps, repool_pd_steps
from schemas import (
    Band,
    BayesianBand,
    CalibrationRequest,
    CalibrationResponse,
    Metrics,
    Observation,
    PipelineRequest,
    PipelineResponse,
    ScorePd,
    SmoothingStage,
    Step,
)

app = FastAPI(title="MAPA demo API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/api/example")
def get_example() -> list[Observation]:
    return load_example_observations()


@app.post("/api/calibrate")
def calibrate(request: CalibrationRequest) -> CalibrationResponse:
    observations = request.observations if request.observations is not None else load_example_observations()
    params = request.params

    bands, smoothed = run_calibration(
        observations,
        min_obs=params.min_obs,
        min_bads=params.min_bads,
        k=params.k,
        min_confidence=params.min_confidence,
        increasing=params.increasing,
    )

    actual_bins = bins_from_observations(observations)

    deviations = [abs(b.n_bads / b.n_obs - interpolate_pd(bands, b.score_min)) for b in actual_bins]
    mad = sum(deviations) / len(deviations)

    return CalibrationResponse(
        bands=[
            Band(score_min=b.score_min, score_max=b.score_max, n_obs=b.n_obs, n_bads=b.n_bads, pd=b.pd)
            for b in bands
        ],
        smoothed=[ScorePd(score=s, pd=pd) for s, pd in smoothed],
        actual=[ScorePd(score=b.score_min, pd=b.n_bads / b.n_obs) for b in actual_bins],
        metrics=Metrics(mad=mad),
    )


@app.post("/api/pipeline")
def pipeline_endpoint(request: PipelineRequest) -> PipelineResponse:
    """Run the full MAPA pipeline and return a step-by-step trace of every
    stage, for the pipeline visualization."""
    observations = request.observations if request.observations is not None else load_example_observations()
    params = request.params

    initial_bins = bins_from_observations(observations)
    pooled = mapa(initial_bins, params.increasing, params.min_confidence)
    pooling_stage = pooling_steps(initial_bins, params.increasing, params.min_confidence)

    size_steps = minimum_size_steps(pooled, params.min_obs, params.min_bads)
    sized_intermediate = _bins_from_step(size_steps[-1])
    resize_steps = pooling_steps(sized_intermediate, params.increasing, params.min_confidence)
    sized = mapa(sized_intermediate, params.increasing, params.min_confidence)
    minimum_size_stage = size_steps + resize_steps

    calibrated = apply_bayesian_adjustment(sized, params.k)
    bayesian_stage = [
        BayesianBand(
            score_min=b.score_min,
            score_max=b.score_max,
            n_obs=b.n_obs,
            n_bads=b.n_bads,
            bad_rate=b.n_bads / b.n_obs,
            pd=b.pd,
        )
        for b in calibrated
    ]

    repooling_stage = repool_pd_steps(calibrated, params.increasing)
    repooled = repool_calibrated_bins(calibrated, params.increasing)
    smoothed = compute_smoothed(repooled)

    return PipelineResponse(
        pooling=pooling_stage,
        minimum_size=minimum_size_stage,
        bayesian=bayesian_stage,
        repooling=repooling_stage,
        smoothing=SmoothingStage(
            bands=[
                Band(score_min=b.score_min, score_max=b.score_max, n_obs=b.n_obs, n_bads=b.n_bads, pd=b.pd)
                for b in repooled
            ],
            smoothed=[ScorePd(score=s, pd=pd) for s, pd in smoothed],
        ),
    )


def _bins_from_step(step: Step) -> list[Bin]:
    return [Bin(score_min=b.score_min, score_max=b.score_max, n_obs=b.n_obs, n_bads=b.n_bads) for b in step.stack]
