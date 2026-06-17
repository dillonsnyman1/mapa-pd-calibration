export type Observation = [score: number, bad: number, weight: number];

export type WeightingMode = "number" | "value";

export interface CalibrationParams {
  min_obs: number;
  min_bads: number;
  k: number;
  min_confidence: number | null;
  increasing: boolean;
  use_counts_for_thresholds: boolean;
}

export interface Band {
  score_min: number;
  score_max: number;
  n_obs: number;
  n_bads: number;
  pd: number;
}

export interface ScorePd {
  score: number;
  pd: number;
}

export interface Metrics {
  mad: number;
}

export interface CalibrationResponse {
  bands: Band[];
  smoothed: ScorePd[];
  actual: ScorePd[];
  metrics: Metrics;
}

export interface StepBin {
  score_min: number;
  score_max: number;
  n_obs: number;
  n_bads: number;
  bad_rate: number;
  count: number;
  count_bads: number;
}

export interface Step {
  action: "push" | "merge";
  stack: StepBin[];
  reason: string | null;
}

export interface PdStepBin {
  score_min: number;
  score_max: number;
  n_obs: number;
  n_bads: number;
  pd: number;
  count: number;
  count_bads: number;
}

export interface PdStep {
  action: "push" | "merge";
  stack: PdStepBin[];
  reason: string | null;
}

export interface BayesianBand {
  score_min: number;
  score_max: number;
  n_obs: number;
  n_bads: number;
  bad_rate: number;
  pd: number;
}

export interface SmoothingStage {
  bands: Band[];
  smoothed: ScorePd[];
}

export interface PipelineResponse {
  pooling: Step[];
  minimum_size: Step[];
  bayesian: BayesianBand[];
  repooling: PdStep[];
  smoothing: SmoothingStage;
}
