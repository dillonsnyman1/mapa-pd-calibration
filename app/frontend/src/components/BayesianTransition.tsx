import {
  CartesianGrid,
  Legend,
  Line,
  LineChart,
  ReferenceArea,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import type { BayesianBand } from "../types";

interface Props {
  bands: BayesianBand[];
  current: number;
  shrunk: boolean;
  disableAllAnimation?: boolean;
}

interface Point {
  score: number;
  bad_rate?: number;
  pd?: number;
  highlight?: number;
}

export function BayesianTransition({ bands, current, shrunk, disableAllAnimation = false }: Props) {
  const scoreMin = bands[0].score_min;
  const scoreMax = bands[bands.length - 1].score_max;

  const totalObs = bands.reduce((sum, b) => sum + b.n_obs, 0);
  const totalBads = bands.reduce((sum, b) => sum + b.n_bads, 0);
  const prior = totalBads / totalObs;

  const points: Point[] = [];
  bands.forEach((b, i) => {
    const isCurrent = i === current;
    const pdValue = i < current || (isCurrent && shrunk) ? b.pd : b.bad_rate;
    points.push({ score: b.score_min, bad_rate: b.bad_rate, pd: pdValue, highlight: isCurrent ? pdValue : undefined });
    points.push({ score: b.score_max, bad_rate: b.bad_rate, pd: pdValue, highlight: isCurrent ? pdValue : undefined });
  });

  return (
    <div>
      <div className="chart-wrap">
      <ResponsiveContainer width="100%" height="100%">
        <LineChart data={points} margin={{ top: 10, right: 20, bottom: 20, left: 10 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
          {bands.map((b, i) => {
            const isCurrent = i === current;
            return (
              <ReferenceArea
                key={`${b.score_min}-${b.score_max}-${i}`}
                x1={b.score_min}
                x2={b.score_max}
                fill={isCurrent ? "#f97316" : "#3461eb"}
                fillOpacity={isCurrent ? 0.18 : i % 2 === 0 ? 0.07 : 0.03}
                stroke={isCurrent ? "#f97316" : "#3461eb"}
                strokeOpacity={isCurrent ? 0.4 : 0.15}
              />
            );
          })}
          <XAxis
            dataKey="score"
            type="number"
            domain={[scoreMin, scoreMax]}
            tick={{ fontSize: 12, fill: "#6b7280" }}
            label={{ value: "Score", position: "insideBottom", offset: -10, fill: "#374151" }}
          />
          <YAxis
            domain={[0, 1]}
            tick={{ fontSize: 12, fill: "#6b7280" }}
            label={{ value: "Rate", angle: -90, position: "insideLeft", fill: "#374151" }}
          />
          <Tooltip />
          <Legend wrapperStyle={{ fontSize: "0.85rem", paddingTop: "0.5rem" }} />
          <ReferenceLine
            y={prior}
            stroke="#10b981"
            strokeWidth={1.5}
            strokeDasharray="6 3"
            label={{
              value: `Portfolio prior: ${prior.toFixed(3)}`,
              position: "insideTopRight",
              fill: "#10b981",
              fontSize: 11,
            }}
          />
          <Line
            type="linear"
            dataKey="bad_rate"
            name="Empirical bad rate"
            stroke="#3461eb"
            strokeWidth={2}
            dot={false}
            isAnimationActive={false}
            connectNulls
          />
          <Line
            type="linear"
            dataKey="pd"
            name="Bayesian-adjusted PD"
            stroke="#f97316"
            strokeWidth={2}
            dot={false}
            isAnimationActive={!disableAllAnimation}
            animationDuration={250}
            animationEasing="ease-in-out"
            connectNulls
          />
          <Line
            type="linear"
            dataKey="highlight"
            name="Current band"
            stroke="#f97316"
            strokeWidth={4}
            dot={false}
            isAnimationActive={false}
            connectNulls={false}
          />
        </LineChart>
      </ResponsiveContainer>
      </div>
    </div>
  );
}
