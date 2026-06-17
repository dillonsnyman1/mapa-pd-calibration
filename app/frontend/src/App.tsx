import { useEffect, useState } from "react";
import { fetchCalibration, fetchExample, fetchPipeline } from "./api";
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
  const [rawObservations, setRawObservations] = useState<Observation[] | null>(null);
  const [observations, setObservations] = useState<Observation[] | null>(null);
  const [exampleCount, setExampleCount] = useState<number | null>(null);
  const [calibration, setCalibration] = useState<CalibrationResponse | null>(null);
  const [pipeline, setPipeline] = useState<PipelineResponse | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (rawObservations === null) {
      setObservations(null);
      return;
    }
    if (weightingMode === "number") {
      setObservations(rawObservations.map(([s, b]) => [s, b, 1]));
    } else {
      setObservations(rawObservations);
    }
  }, [rawObservations, weightingMode]);

  const debouncedParams = useDebounced(params, 300);
  const debouncedObservations = useDebounced(observations, 300);

  useEffect(() => {
    fetchExample().then((obs) => setExampleCount(obs.length)).catch((e) => setError(String(e)));
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

  const datasetLabel = observations
    ? `Custom upload (${observations.length} observations)`
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
          onObservationsChange={setRawObservations}
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
              Mean absolute deviation (predicted vs. actual, per score, averaged): {calibration.metrics.mad.toFixed(4)}
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
