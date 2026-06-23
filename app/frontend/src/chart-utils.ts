export function niceYMax(dataMax: number): number {
  if (dataMax >= 0.9) return 1;
  const padded = dataMax * 1.15;
  const step = padded > 0.2 ? 0.1 : padded > 0.05 ? 0.05 : 0.02;
  return Math.min(1, Math.ceil(padded / step) * step);
}
