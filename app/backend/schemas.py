"""Pydantic request/response models for the demo API."""

from __future__ import annotations

from typing import List, Literal, Optional, Tuple

from pydantic import BaseModel

Observation = Tuple[float, int]


class CalibrationParams(BaseModel):
    min_obs: int = 50
    min_bads: int = 10
    k: float = 10
    min_confidence: Optional[float] = None
    increasing: bool = False


class CalibrationRequest(BaseModel):
    observations: Optional[List[Observation]] = None
    params: CalibrationParams = CalibrationParams()


class Band(BaseModel):
    score_min: float
    score_max: float
    n_obs: int
    n_bads: int
    pd: float


class ScorePd(BaseModel):
    score: float
    pd: float


class Metrics(BaseModel):
    mad: float


class CalibrationResponse(BaseModel):
    bands: List[Band]
    smoothed: List[ScorePd]
    actual: List[ScorePd]
    metrics: Metrics


class StepBin(BaseModel):
    score_min: float
    score_max: float
    n_obs: int
    n_bads: int
    bad_rate: float


class Step(BaseModel):
    action: Literal["push", "merge"]
    stack: List[StepBin]
    reason: Optional[str] = None


class PdStepBin(BaseModel):
    score_min: float
    score_max: float
    n_obs: int
    n_bads: int
    pd: float


class PdStep(BaseModel):
    action: Literal["push", "merge"]
    stack: List[PdStepBin]
    reason: Optional[str] = None


class BayesianBand(BaseModel):
    score_min: float
    score_max: float
    n_obs: int
    n_bads: int
    bad_rate: float
    pd: float


class SmoothingStage(BaseModel):
    bands: List[Band]
    smoothed: List[ScorePd]


class PipelineRequest(BaseModel):
    observations: Optional[List[Observation]] = None
    params: CalibrationParams = CalibrationParams()


class PipelineResponse(BaseModel):
    pooling: List[Step]
    minimum_size: List[Step]
    bayesian: List[BayesianBand]
    repooling: List[PdStep]
    smoothing: SmoothingStage
