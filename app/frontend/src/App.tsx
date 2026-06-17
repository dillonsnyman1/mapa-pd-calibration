import { useEffect, useState } from "react";
import { fetchCalibration, fetchExample, fetchExampleWeighted, fetchPipeline } from "./api";
import { BacktestChart } from "./components/BacktestChart";
import { CalibrationChart } from "./components/CalibrationChart";
import { ControlsPanel } from "./components/ControlsPanel";
import { PipelineView } from "./components/PipelineView";
import type { CalibrationParams, CalibrationResponse, Observation, PipelineResponse, WeightingMode } from "./types";

const DEFAULT_PARAMS: CalibrationParams = {
  min_obs: 50,
  min_bads: 10,
  k: 10,
  min_confidence: null,
  increasing: false,
  use_counts_for_thresholds: true,
};

function useDebounced<T>(value: T, delayMs: number): T {
  const [debounced, setDebounced] = useState(value);
  useEffect(() => {
    const timer = setTimeout(() => setDebounced(value), delayMs);
    return () => clearTimeout(timer);
  }, [value, delayMs]);
  return debounced;
}

export default function App() {
  const [params, setParams] = useState<CalibrationParams>(DEFAULT_PARAMS);
  const [weightingMode, setWeightingMode] = useState<WeightingMode>("number");
  const [customObservations, setCustomObservations] = useState<Observation[] | null>(null);
  const [exampleObs, setExampleObs] = useState<Observation[] | null>(null);
  const [exampleWeightedObs, setExampleWeightedObs] = useState<Observation[] | null>(null);
  const [calibration, setCalibration] = useState<CalibrationResponse | null>(null);
  const [pipeline, setPipeline] = useState<PipelineResponse | null>(null);
  const [error, setError] = useState<string | null>(null);

  const isCustom = customObservations !== null;
  const observations: Observation[] | null = isCustom
    ? (weightingMode === "number" ? customObservations.map(([s, b]) => [s, b, 1]) : customObservations)
    : weightingMode === "value"
      ? exampleWeightedObs
      : exampleObs;

  const debouncedParams = useDebounced(params, 300);
  const debouncedObservations = useDebounced(observations, 300);

  useEffect(() => {
    fetchExample().then(setExampleObs).catch((e) => setError(String(e)));
    fetchExampleWeighted().then(setExampleWeightedObs).catch((e) => setError(String(e)));
  }, []);

  useEffect(() => {
    fetchCalibration(debouncedObservations, debouncedParams)
      .then(setCalibration)
      .catch((e) => setError(String(e)));
  }, [debouncedObservations, debouncedParams]);

  useEffect(() => {
    fetchPipeline(debouncedObservations, debouncedParams)
      .then(setPipeline)
      .catch((e) => setError(String(e)));
  }, [debouncedObservations, debouncedParams]);

  const exampleCount = isCustom ? null : observations?.length ?? null;
  const datasetLabel = isCustom
    ? `Custom upload (${customObservations.length} observations)`
    : exampleCount !== null
      ? `Bundled example (${exampleCount} observations)`
      : "Loading...";

  return (
    <div>
      <header>
        <h1>MAPA score-to-PD calibration demo</h1>
        <p>
          Interactive demo of the Monotone Adjacent Pooling Algorithm (MAPA). See{" "}
          <a href="https://github.com/dillonsnyman1/mapa-pd-calibration" target="_blank" rel="noreferrer">
            the repository
          </a>{" "}
          for the reference implementations and methodology.
        </p>
      </header>

      <div className="card">
        <ControlsPanel
          params={params}
          onParamsChange={setParams}
          onObservationsChange={setCustomObservations}
          onError={setError}
          datasetLabel={datasetLabel}
          error={error}
          weightingMode={weightingMode}
          onWeightingModeChange={setWeightingMode}
        />
      </div>

      <div className="card">
        <h2>Calibration curve</h2>
        <p>Unsmoothed (pooled-band) and smoothed (log-odds interpolated) score-to-PD mappings.</p>
        {calibration ? (
          <CalibrationChart bands={calibration.bands} smoothed={calibration.smoothed} />
        ) : (
          <p>Loading...</p>
        )}
      </div>

      <div className="card">
        <h2>Backtest: predicted vs observed</h2>
        <p>Calibrated (smoothed) PD curve against the actual observed bad rate at each individual score.</p>
        {calibration ? (
          <>
            <BacktestChart smoothed={calibration.smoothed} actual={calibration.actual} />
            <p className="muted-note">
              Weighted mean absolute deviation (predicted vs. observed bad rate): {calibration.metrics.mad.toFixed(4)}
            </p>
          </>
        ) : (
          <p>Loading...</p>
        )}
      </div>

      <div className="card">
        <h2>Pipeline walkthrough</h2>
        <p>
          Steps through the full calibration pipeline: initial pooling, minimum-size enforcement, Bayesian
          shrinkage, re-pooling on the adjusted PDs, and final smoothing.
        </p>
        {pipeline ? <PipelineView pipeline={pipeline} /> : <p>Loading...</p>}
      </div>
    </div>
  );
}
