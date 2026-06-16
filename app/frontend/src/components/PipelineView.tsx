import { useEffect, useMemo, useState } from "react";
import { StackStepper } from "./StackStepper";
import { BayesianTransition } from "./BayesianTransition";
import { SmoothingView } from "./SmoothingView";
import type { PipelineResponse } from "../types";

interface Props {
  pipeline: PipelineResponse;
}

type Frame =
  | { stage: 0; step: number }
  | { stage: 1; step: number }
  | { stage: 2; current: number; shrunk: boolean }
  | { stage: 3; step: number }
  | { stage: 4; phase: "step" | "anchors" | "smoothed" };

const STAGE_LABELS = [
  "1. Pooling",
  "2. Minimum size enforcement",
  "3. Bayesian adjustment",
  "4. Re-pooling",
  "5. Smoothing",
];

const STAGE_DESCRIPTIONS: Record<0 | 1 | 2 | 3 | 4, string> = {
  0: "Pooling: every score starts as its own bin. Working from the lowest score upward, each bin is pushed onto a stack. Whenever the newly pushed bin's bad rate would make the stack non-monotonic with score (or the bin can't be confidently distinguished from the one below it), the top two bins are merged into one. This repeats until the stack is fully monotonic before the next bin is pushed.",
  1: "Minimum-size enforcement: any band with fewer observations or defaults than the configured minimums is merged into whichever neighboring band has the closer bad rate. Because this can change bad rates, the bands are then walked through again with the same push/merge logic as stage 1 to restore monotonicity.",
  2: "Bayesian adjustment: each band's empirical bad rate is shrunk toward the overall portfolio bad rate (the prior). The amount of shrinkage depends on the band's size relative to the credibility weight k - smaller bands (less data) shrink more toward the prior, larger bands shrink less and stay closer to their observed rate.",
  3: "Re-pooling: shrinking each band's rate independently can break the monotonic ordering that pooling established. The bands are walked through again, left to right, merging any adjacent pair whose Bayesian-adjusted PDs violate monotonicity - the same push/merge logic as stage 1, applied to the adjusted PDs.",
  4: "Smoothing: the re-pooled bands form a step function, where every score within a band maps to that band's PD and jumps abruptly at each band boundary. Log-odds interpolation replaces those jumps with a smooth curve, so scores near a boundary get gradually changing PDs instead of a sudden step.",
};

function describeStep(frame: Frame, pipeline: PipelineResponse): string | null {
  switch (frame.stage) {
    case 2: {
      const band = pipeline.bayesian[frame.current];
      const range = `${band.score_min}-${band.score_max}`;
      return frame.shrunk
        ? `Band ${range}: its empirical bad rate of ${band.bad_rate.toFixed(2)} (blue) is shrunk toward the green dashed portfolio prior line, producing a Bayesian-adjusted PD of ${band.pd.toFixed(2)} (orange).`
        : `Band ${range}: starting point is its empirical bad rate of ${band.bad_rate.toFixed(2)} (blue), before being shrunk toward the green dashed portfolio prior line - the overall bad rate across all bands.`;
    }
    case 4:
      switch (frame.phase) {
        case "step":
          return "This is the step function from re-pooling: a flat PD across each band, jumping abruptly at the boundaries.";
        case "anchors":
          return "Each band is reduced to a single anchor point (green dot): its midpoint score, (score_min + score_max) / 2, paired with the band's PD. These anchors are what the interpolation connects.";
        case "smoothed":
          return "Interpolating in log-odds space between consecutive anchor points produces the continuous orange curve - scores between two midpoints get a PD that blends linearly (in log-odds) between the neighboring anchors.";
      }
    default:
      return null;
  }
}

function durationFor(frame: Frame): number {
  switch (frame.stage) {
    case 0:
    case 1:
    case 3:
      return 1100;
    case 2:
      return frame.shrunk ? 1000 : 800;
    case 4:
      return frame.phase === "smoothed" ? 1400 : 1100;
  }
}

export function PipelineView({ pipeline }: Props) {
  const frames = useMemo<Frame[]>(() => {
    const f: Frame[] = [];
    pipeline.pooling.forEach((_, i) => f.push({ stage: 0, step: i }));
    pipeline.minimum_size.forEach((_, i) => f.push({ stage: 1, step: i }));
    pipeline.bayesian.forEach((_, i) => {
      f.push({ stage: 2, current: i, shrunk: false });
      f.push({ stage: 2, current: i, shrunk: true });
    });
    pipeline.repooling.forEach((_, i) => f.push({ stage: 3, step: i }));
    f.push({ stage: 4, phase: "step" });
    f.push({ stage: 4, phase: "anchors" });
    f.push({ stage: 4, phase: "smoothed" });
    return f;
  }, [pipeline]);

  const [frameIndex, setFrameIndex] = useState(0);
  const [playing, setPlaying] = useState(false);

  useEffect(() => {
    setFrameIndex(0);
    setPlaying(false);
  }, [frames]);

  useEffect(() => {
    if (!playing) return;
    const duration = durationFor(frames[frameIndex]);
    const timer = window.setTimeout(() => {
      setFrameIndex((i) => {
        if (i >= frames.length - 1) {
          setPlaying(false);
          return i;
        }
        return i + 1;
      });
    }, duration);
    return () => clearTimeout(timer);
  }, [playing, frameIndex, frames]);

  const frame = frames[frameIndex];
  const atEnd = frameIndex === frames.length - 1;

  return (
    <div>
      <div className="stepper-controls">
        <button
          type="button"
          onClick={() => {
            setFrameIndex(0);
            setPlaying(false);
          }}
          disabled={frameIndex === 0 && !playing}
        >
          Reset
        </button>
        <button
          type="button"
          onClick={() => setFrameIndex((i) => Math.max(0, i - 1))}
          disabled={frameIndex === 0}
        >
          Prev
        </button>
        <button
          type="button"
          onClick={() => setFrameIndex((i) => Math.min(frames.length - 1, i + 1))}
          disabled={atEnd}
        >
          Next
        </button>
        <button
          type="button"
          onClick={() => setPlaying((p) => !p)}
          disabled={atEnd && !playing}
        >
          {playing ? "Pause" : "Play"}
        </button>
        <span className="step-count">
          {STAGE_LABELS[frame.stage]} — Step {frameIndex + 1} / {frames.length}
        </span>
      </div>

      <p className="stepper-reason">{STAGE_DESCRIPTIONS[frame.stage]}</p>
      {describeStep(frame, pipeline) && <p className="stepper-reason">{describeStep(frame, pipeline)}</p>}

      {frame.stage === 0 && (
        <StackStepper
          steps={pipeline.pooling}
          index={frame.step}
          valueOf={(b) => b.bad_rate}
          valueLabel="bad rate"
          pushMessage={(step) => {
            const b = step.stack[step.stack.length - 1];
            return `Pushed the ${b.score_min}-${b.score_max} bin (bad rate ${b.bad_rate.toFixed(3)}) onto the stack.`;
          }}
        />
      )}

      {frame.stage === 1 && (
        <StackStepper
          steps={pipeline.minimum_size}
          index={frame.step}
          valueOf={(b) => b.bad_rate}
          valueLabel="bad rate"
          pushMessage={(step) => {
            if (step.stack.length > 1) {
              return "Starting point: the bands produced by pooling, before checking minimum size and default-count requirements.";
            }
            const b = step.stack[0];
            return `Re-checking monotonicity after a merge: pushed the ${b.score_min}-${b.score_max} band (bad rate ${b.bad_rate.toFixed(3)}) onto the stack.`;
          }}
          referenceBands={pipeline.pooling[pipeline.pooling.length - 1].stack.map((b) => ({
            score_min: b.score_min,
            score_max: b.score_max,
            value: b.bad_rate,
          }))}
          referenceLabel="Post-pooling bad rate (before minimum-size enforcement)"
        />
      )}

      {frame.stage === 2 && (
        <BayesianTransition bands={pipeline.bayesian} current={frame.current} shrunk={frame.shrunk} />
      )}

      {frame.stage === 3 && (
        <StackStepper
          steps={pipeline.repooling}
          index={frame.step}
          valueOf={(b) => b.pd}
          valueLabel="pd"
          pushMessage={(step) => {
            const b = step.stack[step.stack.length - 1];
            return `Pushed the ${b.score_min}-${b.score_max} band (PD ${b.pd.toFixed(3)}) onto the stack.`;
          }}
          referenceBands={pipeline.bayesian.map((b) => ({
            score_min: b.score_min,
            score_max: b.score_max,
            value: b.pd,
          }))}
          referenceLabel="Post-Bayesian PD (before re-pooling)"
        />
      )}

      {frame.stage === 4 && <SmoothingView smoothing={pipeline.smoothing} phase={frame.phase} />}
    </div>
  );
}
