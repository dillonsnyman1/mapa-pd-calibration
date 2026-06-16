import {
  CartesianGrid,
  ComposedChart,
  Legend,
  Line,
  ResponsiveContainer,
  Scatter,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import type { ScorePd } from "../types";

interface Props {
  smoothed: ScorePd[];
  actual: ScorePd[];
}

export function BacktestChart({ smoothed, actual }: Props) {
  const scoreMin = Math.min(smoothed[0]?.score ?? 0, actual[0]?.score ?? 0);
  const scoreMax = Math.max(
    smoothed[smoothed.length - 1]?.score ?? 0,
    actual[actual.length - 1]?.score ?? 0,
  );

  return (
    <div className="chart-wrap">
    <ResponsiveContainer width="100%" height="100%">
      <ComposedChart margin={{ top: 10, right: 20, bottom: 20, left: 10 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
        <XAxis
          dataKey="score"
          type="number"
          domain={[scoreMin, scoreMax]}
          tick={{ fontSize: 12, fill: "#6b7280" }}
          label={{ value: "Score", position: "insideBottom", offset: -10, fill: "#374151" }}
        />
        <YAxis
          dataKey="pd"
          domain={[0, 1]}
          tick={{ fontSize: 12, fill: "#6b7280" }}
          label={{ value: "PD / bad rate", angle: -90, position: "insideLeft", fill: "#374151" }}
        />
        <Tooltip />
        <Legend wrapperStyle={{ fontSize: "0.85rem", paddingTop: "0.5rem" }} />
        <Scatter
          data={actual}
          dataKey="pd"
          name="Actual bad rate (per score)"
          fill="#3461eb"
          fillOpacity={0.5}
          shape={(props: { cx?: number; cy?: number }) => (
            <circle cx={props.cx} cy={props.cy} r={2} fill="#3461eb" fillOpacity={0.5} />
          )}
          isAnimationActive={false}
        />
        <Line
          data={smoothed}
          type="monotone"
          dataKey="pd"
          name="Calibrated PD (smoothed)"
          stroke="#f97316"
          strokeWidth={2}
          dot={false}
          isAnimationActive={false}
        />
      </ComposedChart>
    </ResponsiveContainer>
    </div>
  );
}
