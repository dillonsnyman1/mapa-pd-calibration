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

interface BinLike {
  score_min: number;
  score_max: number;
  n_obs: number;
  n_bads: number;
}

interface StepLike<T extends BinLike> {
  action: "push" | "merge";
  stack: T[];
  reason: string | null;
}

interface ReferenceBand {
  score_min: number;
  score_max: number;
  value: number;
}

interface Props<T extends BinLike> {
  steps: StepLike<T>[];
  index: number;
  valueOf: (b: T) => number;
  valueLabel: string;
  pushMessage: (step: StepLike<T>) => string;
  referenceBands?: ReferenceBand[];
  referenceLabel?: string;
}

interface Point {
  score: number;
  value?: number;
  highlight?: number;
  reference?: number;
}

export function StackStepper<T extends BinLike>({
  steps,
  index,
  valueOf,
  valueLabel,
  pushMessage,
  referenceBands,
  referenceLabel,
}: Props<T>) {
  if (steps.length === 0) {
    return <p>No steps to show.</p>;
  }

  const clampedIndex = Math.min(index, steps.length - 1);
  const step = steps[clampedIndex];
  const scoreMin = Math.min(
    ...steps[0].stack.map((b) => b.score_min),
    ...(referenceBands ?? []).map((b) => b.score_min),
  );
  const scoreMax = Math.max(
    ...steps[steps.length - 1].stack.map((b) => b.score_max),
    ...(referenceBands ?? []).map((b) => b.score_max),
  );

  const points: Point[] = [];
  step.stack.forEach((b, i) => {
    const value = valueOf(b);
    const isLast = i === step.stack.length - 1;
    points.push({ score: b.score_min, value, highlight: isLast ? value : undefined });
    points.push({ score: b.score_max, value, highlight: isLast ? value : undefined });
  });

  const referencePoints: Point[] = [];
  (referenceBands ?? []).forEach((b) => {
    referencePoints.push({ score: b.score_min, reference: b.value });
    referencePoints.push({ score: b.score_max, reference: b.value });
  });

  const mergedPoints = [...points, ...referencePoints].sort((a, b) => a.score - b.score);

  return (
    <div>
      <p className="stepper-reason">{step.action === "push" ? pushMessage(step) : step.reason}</p>

      <div className="chart-wrap">
      <ResponsiveContainer width="100%" height="100%" debounce={1}>
        <LineChart data={mergedPoints} margin={{ top: 10, right: 20, bottom: 20, left: 10 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#e5e7eb" />
          {step.stack.map((b, i) => {
            const isLast = i === step.stack.length - 1;
            return (
              <ReferenceArea
                key={`${b.score_min}-${b.score_max}-${i}`}
                x1={b.score_min}
                x2={b.score_max}
                fill={isLast ? "#f97316" : "#3461eb"}
                fillOpacity={isLast ? 0.18 : i % 2 === 0 ? 0.07 : 0.03}
                stroke={isLast ? "#f97316" : "#3461eb"}
                strokeOpacity={isLast ? 0.4 : 0.15}
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
            label={{ value: valueLabel, angle: -90, position: "insideLeft", fill: "#374151" }}
          />
          <Tooltip />
          {referenceBands && (
            <Legend wrapperStyle={{ fontSize: "0.85rem", paddingTop: "0.5rem" }} />
          )}
          {referenceBands && (
            <Line
              type="linear"
              dataKey="reference"
              name={referenceLabel ?? "Reference"}
              stroke="#9ca3af"
              strokeWidth={1.5}
              strokeDasharray="4 4"
              dot={false}
              isAnimationActive={false}
              connectNulls
            />
          )}
          <Line
            type="linear"
            dataKey="value"
            name="Pooled bands"
            stroke="#3461eb"
            strokeWidth={2}
            dot={false}
            isAnimationActive={false}
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
