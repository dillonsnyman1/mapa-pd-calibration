import { CalibrationChart } from "./CalibrationChart";
import type { SmoothingStage } from "../types";

interface Props {
  smoothing: SmoothingStage;
  phase: "step" | "anchors" | "smoothed";
  disableAllAnimation?: boolean;
}

export function SmoothingView({ smoothing, phase, disableAllAnimation = false }: Props) {
  const anchors =
    phase === "step"
      ? undefined
      : smoothing.bands.map((b) => ({ score: (b.score_min + b.score_max) / 2, pd: b.pd }));

  return (
    <div>
      <CalibrationChart
        bands={smoothing.bands}
        smoothed={phase === "smoothed" ? smoothing.smoothed : []}
        animateSmoothed={phase === "smoothed"}
        anchors={anchors}
        showBandRanges
        disableAllAnimation={disableAllAnimation}
      />
    </div>
  );
}
