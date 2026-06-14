"""Monotone Adjacent Pooling Algorithm (MAPA).

MAPA turns raw (score, bad) observations into a score-to-PD calibration
curve that is guaranteed to be monotone: as the score improves, the
calibrated PD never gets worse.

The full pipeline:

1. `bins_from_observations` groups the raw observations into one bin per
   unique score - the finest possible starting point.
2. `mapa` repeatedly merges ("pools") adjacent bins whose bad rates
   violate the required monotonicity, until none remain. This is the same
   idea as the Pool Adjacent Violators Algorithm (PAVA) used in isotonic
   regression. An optional `min_confidence` also merges adjacent bins
   whose bad rates aren't statistically distinguishable, even if they
   don't violate monotonicity.
3. `enforce_minimum_size` further pools any bins that fall below given
   `min_obs`/`min_bads` thresholds into their closer-rate neighbour, then
   re-runs `mapa` to restore monotonicity.
4. `apply_bayesian_adjustment` shrinks each bin's bad rate toward a prior
   (typically the overall bad rate) using a credibility weight `k`,
   producing a `CalibratedBin.pd` per bin.
5. `repool_calibrated_bins` re-pools any bins whose `pd` became
   non-monotone after shrinkage, restoring a monotone `pd` sequence.
6. `interpolate_pd` smooths the resulting step function into a
   continuous PD curve via log-odds interpolation.

`calibrate` runs steps 1-2. `run_pipeline` runs steps 1-6 and returns a
`CalibrationResult` bundling the final bands and a `pd_for_score` lookup.

See ../../docs/mapa-methodology.md for background and attribution.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple


@dataclass
class Bin:
    """A score band spanning [score_min, score_max], with its observation
    and bad counts."""

    score_min: float
    score_max: float
    n_obs: int
    n_bads: int

    @property
    def bad_rate(self) -> float:
        return self.n_bads / self.n_obs


@dataclass
class CalibratedBin:
    """A pooled bin together with its Bayesian-adjusted PD, as produced by
    `apply_bayesian_adjustment`."""

    score_min: float
    score_max: float
    n_obs: int
    n_bads: int
    pd: float


@dataclass
class CalibrationResult:
    """The output of `run_pipeline`: a calibrated band table, plus a smoothed
    per-score PD curve derived from it.

    `bands` is the step-function calibration table - the typical deliverable
    for reporting and governance. `pd_for_score` gives a smoothed,
    continuous PD for an individual score, via `interpolate_pd`. Both are
    views of the same underlying bands; which one to use depends on the
    consumer.
    """

    bands: List[CalibratedBin]

    def pd_for_score(self, score: float) -> float:
        return interpolate_pd(self.bands, score)


def bins_from_observations(observations: Iterable[Tuple[float, int]]) -> List[Bin]:
    """Group raw (score, bad) observations into one bin per unique score,
    ordered by score ascending. This is the finest possible starting point
    for `mapa`.

    Args:
        observations: An iterable of (score, bad) pairs, where `bad` is 1
            for a default and 0 otherwise.
    """
    counts: dict[float, list[int]] = defaultdict(lambda: [0, 0])
    for score, bad in observations:
        entry = counts[score]
        entry[0] += 1
        entry[1] += int(bad)

    return [
        Bin(score_min=score, score_max=score, n_obs=n_obs, n_bads=n_bads)
        for score, (n_obs, n_bads) in sorted(counts.items())
    ]


def mapa(
    bins: List[Bin],
    increasing: bool = False,
    min_confidence: Optional[float] = None,
) -> List[Bin]:
    """Run the Monotone Adjacent Pooling Algorithm.

    Args:
        bins: Bins ordered by score, ascending (e.g. from
            `bins_from_observations`).
        increasing: If False (the default), the bad rate is required to be
            non-increasing as score increases - the standard credit scoring
            convention (higher score = lower risk). If True, the bad rate
            is required to be non-decreasing instead.
        min_confidence: If given (e.g. 0.95 for 95%), adjacent bins whose
            bad rates do not differ at this confidence level (two-proportion
            z-test) are merged as well, even if they don't violate
            monotonicity. This produces fewer, larger bins whose bad rates
            are more reliably distinguishable from their neighbours. If not
            given (the default), only monotonicity violations are merged.

    Returns:
        Pooled bins, in score order, whose bad rates satisfy the requested
        monotonicity (and, if `min_confidence` is given, are pairwise
        distinguishable from their neighbours at that confidence level).
        Together they partition the full score range and population of the
        input bins.
    """
    stack: List[Bin] = []

    for b in bins:
        stack.append(b)
        while len(stack) >= 2 and (
            _violates(stack[-2], stack[-1], increasing)
            or (min_confidence is not None and _not_significant(stack[-2], stack[-1], min_confidence))
        ):
            top = stack.pop()
            below = stack.pop()
            stack.append(_merge(below, top))

    return stack


def calibrate(
    observations: Iterable[Tuple[float, int]],
    increasing: bool = False,
    min_confidence: Optional[float] = None,
) -> List[Bin]:
    """Convenience wrapper: group raw observations into per-score bins, then
    pool them with `mapa`."""
    return mapa(bins_from_observations(observations), increasing, min_confidence)


def enforce_minimum_size(
    bins: List[Bin],
    min_obs: int = 0,
    min_bads: int = 0,
    increasing: bool = False,
    min_confidence: Optional[float] = None,
) -> List[Bin]:
    """Further pool bins that don't meet minimum size thresholds, even if
    they don't violate monotonicity.

    A bin "violates" if `n_obs < min_obs` or `n_bads < min_bads`. Each
    violating bin is repeatedly merged into whichever adjacent bin has the
    closer bad rate (minimizing the distortion to the calibration curve),
    until every remaining bin meets both thresholds or only one bin
    remains.

    Merging toward the closer-rate neighbour is not guaranteed to preserve
    the monotonicity established by `mapa`, so the result is passed back
    through `mapa` before being returned.

    Args:
        bins: Pooled bins, in score order, typically the output of `mapa`.
        min_obs: Minimum number of observations required per bin.
        min_bads: Minimum number of bads required per bin.
        increasing: Passed through to the final `mapa` pass; see `mapa`.
        min_confidence: Passed through to the final `mapa` pass; see `mapa`.
    """
    bins = list(bins)

    while len(bins) > 1:
        violator = next(
            (i for i, b in enumerate(bins) if b.n_obs < min_obs or b.n_bads < min_bads),
            None,
        )
        if violator is None:
            break

        if violator == 0:
            neighbour = 1
        elif violator == len(bins) - 1:
            neighbour = violator - 1
        else:
            rate = bins[violator].bad_rate
            left_diff = abs(rate - bins[violator - 1].bad_rate)
            right_diff = abs(rate - bins[violator + 1].bad_rate)
            neighbour = violator - 1 if left_diff <= right_diff else violator + 1

        i, j = sorted((violator, neighbour))
        bins = bins[:i] + [_merge(bins[i], bins[j])] + bins[j + 1 :]

    return mapa(bins, increasing, min_confidence)


def apply_bayesian_adjustment(
    bins: List[Bin], k: float, prior: Optional[float] = None
) -> List[CalibratedBin]:
    """Shrink each bin's empirical bad rate toward a prior using Bayesian
    (credibility) weighting.

    Each bin's adjusted PD is:

        pd = (n_bads + k * prior) / (n_obs + k)

    `k` is the credibility weight, expressed as a number of "equivalent
    observations" of the prior: a bin with `n_obs == k` is shrunk halfway
    between its own bad rate and the prior, bins with `n_obs >> k` are
    barely adjusted, and bins with `n_obs << k` end up close to the prior.

    Args:
        bins: Pooled bins, typically the output of `mapa`.
        k: Credibility weight (in equivalent observations). Larger values
            apply stronger shrinkage.
        prior: The PD to shrink toward. If not given, defaults to the
            overall bad rate across all bins (sum(n_bads) / sum(n_obs)).

    Returns:
        One `CalibratedBin` per input bin, in the same order, with `pd` set
        to the shrunk estimate.

    Note:
        Shrinking each bin independently toward a single global prior is
        not guaranteed to preserve the monotonicity established by `mapa`:
        a small bin whose empirical rate is close to a much larger
        neighbour's can be pulled past it toward the prior. See
        ../../docs/mapa-methodology.md for discussion.
    """
    if prior is None:
        total_obs = sum(b.n_obs for b in bins)
        total_bads = sum(b.n_bads for b in bins)
        prior = total_bads / total_obs

    return [
        CalibratedBin(
            score_min=b.score_min,
            score_max=b.score_max,
            n_obs=b.n_obs,
            n_bads=b.n_bads,
            pd=(b.n_bads + k * prior) / (b.n_obs + k),
        )
        for b in bins
    ]


def repool_calibrated_bins(
    bins: List[CalibratedBin], increasing: bool = False
) -> List[CalibratedBin]:
    """Re-apply pooling to Bayesian-adjusted bins, restoring monotonicity of
    `pd`.

    `apply_bayesian_adjustment` shrinks each bin's bad rate toward a shared
    prior independently, which is not guaranteed to preserve the
    monotonicity established by `mapa` (see its docstring). This function
    runs the same adjacent-pooling algorithm again, but on `pd` instead of
    `bad_rate`, merging violating bins by taking the `n_obs`-weighted
    average of their `pd` values.

    Args:
        bins: Calibrated bins, in score order, typically the output of
            `apply_bayesian_adjustment`.
        increasing: Same meaning as in `mapa`, applied to `pd` instead of
            `bad_rate`.

    Returns:
        Pooled calibrated bins, in score order, whose `pd` values satisfy
        the requested monotonicity.
    """
    stack: List[CalibratedBin] = []

    for b in bins:
        stack.append(b)
        while len(stack) >= 2 and _violates_pd(stack[-2], stack[-1], increasing):
            top = stack.pop()
            below = stack.pop()
            stack.append(_merge_calibrated(below, top))

    return stack


def run_pipeline(
    observations: Iterable[Tuple[float, int]],
    k: float,
    min_obs: int = 0,
    min_bads: int = 0,
    prior: Optional[float] = None,
    increasing: bool = False,
    min_confidence: Optional[float] = None,
) -> CalibrationResult:
    """Run the full MAPA pipeline: bin, pool, enforce minimum size, apply
    Bayesian adjustment, and re-pool.

    This chains `calibrate`, `enforce_minimum_size`,
    `apply_bayesian_adjustment` and `repool_calibrated_bins`. The result
    bundles the resulting band table (`bands`) together with a smoothed,
    continuous PD curve derived from it (`pd_for_score`, via
    `interpolate_pd`) - use whichever representation suits the consumer.

    Args:
        observations: Raw (score, bad) observations; see
            `bins_from_observations`.
        k: Bayesian credibility weight; see `apply_bayesian_adjustment`.
        min_obs: Minimum observations per bin; see `enforce_minimum_size`.
        min_bads: Minimum bads per bin; see `enforce_minimum_size`.
        prior: PD to shrink toward; see `apply_bayesian_adjustment`.
        increasing: Direction of monotonicity; see `mapa`.
        min_confidence: Confidence-based pooling threshold; see `mapa`.
    """
    pooled = calibrate(observations, increasing, min_confidence)
    sized = enforce_minimum_size(pooled, min_obs, min_bads, increasing, min_confidence)
    calibrated = apply_bayesian_adjustment(sized, k, prior)
    repooled = repool_calibrated_bins(calibrated, increasing)
    return CalibrationResult(bands=repooled)


def interpolate_pd(bins: List[CalibratedBin], score: float) -> float:
    """Return a smoothed PD for an individual score via log-odds
    interpolation between pools.

    The pooled, calibrated PD curve is a step function: every score within
    a pool gets that pool's single PD, with a discontinuous jump at each
    pool boundary. This produces a smooth, continuous PD instead.

    Each pool is represented by a single anchor point: its midpoint score,
    `(score_min + score_max) / 2`, and the log-odds of its `pd`,

        log_odds = ln((1 - pd) / pd)

    For the given `score`:

    - if it is at or before the first pool's midpoint, or at or after the
      last pool's midpoint, the nearest pool's `pd` is returned unchanged
      (flat extrapolation beyond the anchor points).
    - otherwise, `log_odds` is linearly interpolated between the midpoints
      of the two pools bracketing `score`, and converted back to a PD via
      `pd = 1 / (1 + exp(log_odds))`.

    Because log-odds is a monotonic transform of `pd`, this preserves the
    monotonicity of a monotone input sequence of pool PDs.

    Args:
        bins: Calibrated bins, in score order, typically the output of
            `repool_calibrated_bins`.
        score: The individual score to compute a smoothed PD for.
    """
    mids = [(b.score_min + b.score_max) / 2 for b in bins]
    log_odds = [math.log((1 - b.pd) / b.pd) for b in bins]

    if score <= mids[0]:
        return bins[0].pd
    if score >= mids[-1]:
        return bins[-1].pd

    for i in range(len(mids) - 1):
        if mids[i] <= score <= mids[i + 1]:
            t = (score - mids[i]) / (mids[i + 1] - mids[i])
            interpolated = log_odds[i] + t * (log_odds[i + 1] - log_odds[i])
            return 1 / (1 + math.exp(interpolated))

    raise AssertionError("unreachable: score must lie between mids[0] and mids[-1] here")


def _violates(lower: Bin, upper: Bin, increasing: bool) -> bool:
    """Whether `upper` (the higher-scoring bin) violates monotonicity relative to `lower`."""
    if increasing:
        return upper.bad_rate < lower.bad_rate
    return upper.bad_rate > lower.bad_rate


def _not_significant(a: Bin, b: Bin, confidence: float) -> bool:
    """Whether the bad rates of `a` and `b` are statistically indistinguishable
    at the given confidence level, via a two-proportion z-test.

    Returns True (i.e. "merge these") when the observed difference in bad
    rates is small relative to its standard error - either because the
    rates are genuinely close, or because both bins are too small to tell
    them apart.
    """
    pooled_rate = (a.n_bads + b.n_bads) / (a.n_obs + b.n_obs)
    if pooled_rate <= 0 or pooled_rate >= 1:
        return True

    se = math.sqrt(pooled_rate * (1 - pooled_rate) * (1 / a.n_obs + 1 / b.n_obs))
    z = abs(a.bad_rate - b.bad_rate) / se
    z_critical = _inverse_normal_cdf((1 + confidence) / 2)
    return z < z_critical


def _inverse_normal_cdf(p: float) -> float:
    """Approximate the inverse of the standard normal CDF (the quantile
    function / probit), via Acklam's rational approximation.

    Used to convert a confidence level (e.g. 0.95) into a critical z-value
    for `_not_significant`. Accurate to about 1.15e-9.
    """
    a = (
        -3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
        1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00,
    )
    b = (
        -5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
        6.680131188771972e01, -1.328068155288572e01,
    )
    c = (
        -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
        -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00,
    )
    d = (
        7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
        3.754408661907416e00,
    )

    p_low = 0.02425
    p_high = 1 - p_low

    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1
        )
    if p <= p_high:
        q = p - 0.5
        r = q * q
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (
            ((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1
        )
    q = math.sqrt(-2 * math.log(1 - p))
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
        (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1
    )


def _merge(a: Bin, b: Bin) -> Bin:
    return Bin(
        score_min=a.score_min,
        score_max=b.score_max,
        n_obs=a.n_obs + b.n_obs,
        n_bads=a.n_bads + b.n_bads,
    )


def _violates_pd(lower: CalibratedBin, upper: CalibratedBin, increasing: bool) -> bool:
    """Whether `upper` (the higher-scoring bin) violates monotonicity of `pd` relative to `lower`."""
    if increasing:
        return upper.pd < lower.pd
    return upper.pd > lower.pd


def _merge_calibrated(a: CalibratedBin, b: CalibratedBin) -> CalibratedBin:
    n_obs = a.n_obs + b.n_obs
    n_bads = a.n_bads + b.n_bads
    return CalibratedBin(
        score_min=a.score_min,
        score_max=b.score_max,
        n_obs=n_obs,
        n_bads=n_bads,
        pd=(a.pd * a.n_obs + b.pd * b.n_obs) / n_obs,
    )
