import type { CalibrationParams, Observation, WeightingMode } from "../types";

interface Props {
  params: CalibrationParams;
  onParamsChange: (params: CalibrationParams) => void;
  onObservationsChange: (observations: Observation[] | null) => void;
  onError: (message: string | null) => void;
  datasetLabel: string;
  error: string | null;
  weightingMode: WeightingMode;
  onWeightingModeChange: (mode: WeightingMode) => void;
}

function parseCsv(text: string): Observation[] {
  const lines = text.trim().split(/\r?\n/);
  const observations: Observation[] = [];

  for (const [i, line] of lines.entries()) {
    if (!line.trim()) continue;
    const parts = line.split(",").map((p) => p.trim());
    if (i === 0 && Number.isNaN(Number(parts[0]))) continue; // header row

    if (parts.length < 2) {
      throw new Error(`Line ${i + 1}: expected at least "score,bad", got "${line}"`);
    }
    const score = Number(parts[0]);
    const bad = Number(parts[1]);
    if (Number.isNaN(score) || (bad !== 0 && bad !== 1)) {
      throw new Error(`Line ${i + 1}: invalid row "${line}" (expected numeric score and bad of 0 or 1)`);
    }

    let weight = 1;
    if (parts.length >= 3) {
      weight = Number(parts[2]);
      if (Number.isNaN(weight) || weight < 0) {
        throw new Error(`Line ${i + 1}: invalid weight "${parts[2]}" (expected non-negative number)`);
      }
    }

    observations.push([score, bad, weight]);
  }

  if (observations.length === 0) {
    throw new Error("No observations found in file");
  }
  return observations;
}

export function ControlsPanel({
  params, onParamsChange, onObservationsChange, onError, datasetLabel, error,
  weightingMode, onWeightingModeChange,
}: Props) {
  const update = <K extends keyof CalibrationParams>(key: K, value: CalibrationParams[K]) => {
    onParamsChange({ ...params, [key]: value });
  };

  const handleFile = async (file: File) => {
    try {
      const observations = parseCsv(await file.text());
      onError(null);
      onObservationsChange(observations);
    } catch (e) {
      onError(e instanceof Error ? e.message : String(e));
    }
  };

  return (
    <fieldset>
      <legend>Parameters</legend>

      <div className="control-group">
        <div className="control-field">
          <label>Dataset</label>
          <p className="dataset-label">{datasetLabel}</p>
          <div style={{ marginTop: "0.5rem", display: "flex", gap: "0.5rem", alignItems: "center", flexWrap: "wrap" }}>
            <input
              type="file"
              accept=".csv"
              onChange={(e) => {
                const file = e.target.files?.[0];
                if (file) handleFile(file);
              }}
            />
            {!datasetLabel.startsWith("Bundled example") && (
              <button
                type="button"
                onClick={() => {
                  onError(null);
                  onObservationsChange(null);
                }}
              >
                Reset to bundled example
              </button>
            )}
          </div>
          {error && <p className="error">{error}</p>}
          <p className="muted-note" style={{ marginTop: "0.5rem" }}>
            Upload a CSV with <code>score,bad</code>
            {weightingMode === "value" ? <code>,weight</code> : null} columns.
          </p>
        </div>
      </div>

      <div className="control-group">
        <div className="controls-grid">
          <div className="control-field">
            <label>Weighting</label>
            <div style={{ display: "flex", gap: "1rem", marginTop: "0.25rem" }}>
              <label className="checkbox-label">
                <input
                  type="radio"
                  name="weighting"
                  checked={weightingMode === "number"}
                  onChange={() => onWeightingModeChange("number")}
                />
                Number weighted
              </label>
              <label className="checkbox-label">
                <input
                  type="radio"
                  name="weighting"
                  checked={weightingMode === "value"}
                  onChange={() => onWeightingModeChange("value")}
                />
                Value weighted
              </label>
            </div>
          </div>
        </div>
      </div>

      <div className="control-group">
        <div className="controls-grid">
          <div className="control-field">
            <label>
              Minimum observations per band <span className="value">({params.min_obs})</span>
            </label>
            <input
              type="range"
              min={0}
              max={200}
              value={params.min_obs}
              onChange={(e) => update("min_obs", Number(e.target.value))}
            />
          </div>

          <div className="control-field">
            <label>
              Minimum bads per band <span className="value">({params.min_bads})</span>
            </label>
            <input
              type="range"
              min={0}
              max={50}
              value={params.min_bads}
              onChange={(e) => update("min_bads", Number(e.target.value))}
            />
          </div>

          <div className="control-field">
            <label>
              Bayesian credibility, k <span className="value">({params.k})</span>
            </label>
            <input
              type="range"
              min={0}
              max={100}
              value={params.k}
              onChange={(e) => update("k", Number(e.target.value))}
            />
          </div>

          {params.min_confidence !== null && (
            <div className="control-field">
              <label>
                Confidence level <span className="value">({params.min_confidence.toFixed(3)})</span>
              </label>
              <input
                type="range"
                min={0.5}
                max={0.999}
                step={0.001}
                value={params.min_confidence}
                onChange={(e) => update("min_confidence", Number(e.target.value))}
              />
            </div>
          )}
        </div>
      </div>

      <div className="control-group">
        <div className="controls-grid">
          <label className="checkbox-label">
            <input
              type="checkbox"
              checked={params.increasing}
              onChange={(e) => update("increasing", e.target.checked)}
            />
            Bad rate increases with score
          </label>

          <label className="checkbox-label">
            <input
              type="checkbox"
              checked={params.min_confidence !== null}
              onChange={(e) => update("min_confidence", e.target.checked ? 0.95 : null)}
            />
            Confidence-based pooling
          </label>

          {weightingMode === "value" && (
            <label className="checkbox-label">
              <input
                type="checkbox"
                checked={!params.use_counts_for_thresholds}
                onChange={(e) => update("use_counts_for_thresholds", !e.target.checked)}
              />
              Apply minimum thresholds to weighted sums
            </label>
          )}
        </div>
      </div>
    </fieldset>
  );
}
