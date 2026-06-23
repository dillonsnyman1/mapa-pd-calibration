import {
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ReferenceArea,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import type { Band, ScorePd } from "../types";

interface Props {
  bands: Band[];
  smoothed: ScorePd[];
  animateSmoothed?: boolean;
  anchors?: ScorePd[];
  showBandRanges?: boolean;
  disableAllAnimation?: boolean;
}

interface Point {
  score: number;
  unsmoothed?: number;
  smoothed?: number;
  anchor?: number;
}

export function CalibrationChart({ bands, smoothed, animateSmoothed = false, anchors, showBandRanges = false, disableAllAnimation = false }: Props) {
  const stepPoints: Point[] = [];
  for (const b of bands) {
    stepPoints.push({ score: b.score_min, unsmoothed: b.pd });
    stepPoints.push({ score: b.score_max, unsmoothed: b.pd });
  }

  const smoothedPoints: Point[] = smoothed.map((p) => ({ score: p.score, smoothed: p.pd }));
  const anchorPoints: Point[] = (anchors ?? []).map((p) => ({ score: p.score, anchor: p.pd }));

  const points = [...stepPoints, ...smoothedPoints, ...anchorPoints].sort((a, b) => a.score - b.score);

  return (
    <div className="chart-wrap">
    <ResponsiveContainer width="100%" height="100%" debounce={1}>
      <LineChart data={points} margin={{ top: 10, right: 20, bottom: 20, left: 10 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
        {showBandRanges &&
          bands.map((b, i) => (
            <ReferenceArea
              key={`${b.score_min}-${b.score_max}-${i}`}
              x1={b.score_min}
              x2={b.score_max}
              fill="#3461eb"
              fillOpacity={i % 2 === 0 ? 0.07 : 0.03}
              stroke="#3461eb"
              strokeOpacity={0.15}
            />
          ))}
        <XAxis
          dataKey="score"
          type="number"
          domain={["dataMin", "dataMax"]}
          tick={{ fontSize: 12, fill: "#6b7280" }}
          label={{ value: "Score", position: "insideBottom", offset: -10, fill: "#374151" }}
        />
        <YAxis
          domain={[0, 1]}
          tick={{ fontSize: 12, fill: "#6b7280" }}
          label={{ value: "PD", angle: -90, position: "insideLeft", fill: "#374151" }}
        />
        <Tooltip />
        <Legend wrapperStyle={{ fontSize: "0.85rem", paddingTop: "0.5rem" }} />
        <Line
          type="linear"
          dataKey="unsmoothed"
          name="Unsmoothed (pooled bands)"
          stroke="#3461eb"
          strokeWidth={2}
          dot={false}
          isAnimationActive={false}
          connectNulls
        />
        <Line
          type="monotone"
          dataKey="smoothed"
          name="Smoothed (log-odds interpolation)"
          stroke="#f97316"
          strokeWidth={2}
          dot={false}
          isAnimationActive={animateSmoothed && !disableAllAnimation}
          connectNulls
        />
        {anchors && (
          <Line
            dataKey="anchor"
            name="Band midpoints (interpolation anchors)"
            stroke="none"
            dot={{ r: 5, fill: "#16a34a", stroke: "#fff", strokeWidth: 1 }}
            isAnimationActive={false}
            connectNulls
          />
        )}
      </LineChart>
    </ResponsiveContainer>
    </div>
  );
}
