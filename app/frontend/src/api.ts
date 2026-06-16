import type {
  CalibrationParams,
  CalibrationResponse,
  Observation,
  PipelineResponse,
} from "./types";

const BASE_URL = import.meta.env.VITE_API_BASE_URL ?? "http://localhost:8000";

export async function fetchExample(): Promise<Observation[]> {
  const res = await fetch(`${BASE_URL}/api/example`);
  if (!res.ok) throw new Error(`/api/example failed: ${res.status}`);
  return res.json();
}

export async function fetchCalibration(
  observations: Observation[] | null,
  params: CalibrationParams,
): Promise<CalibrationResponse> {
  const res = await fetch(`${BASE_URL}/api/calibrate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ observations, params }),
  });
  if (!res.ok) throw new Error(`/api/calibrate failed: ${res.status}`);
  return res.json();
}

export async function fetchPipeline(
  observations: Observation[] | null,
  params: CalibrationParams,
): Promise<PipelineResponse> {
  const res = await fetch(`${BASE_URL}/api/pipeline`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ observations, params }),
  });
  if (!res.ok) throw new Error(`/api/pipeline failed: ${res.status}`);
  return res.json();
}
