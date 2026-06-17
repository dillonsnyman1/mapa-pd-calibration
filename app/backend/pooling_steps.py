"""Step-by-step trace of the MAPA pooling loop, for the interactive demo.

Reuses the validated merge-decision helpers from reference/python/mapa.py
(`_violates`, `_not_significant`, `_merge`) so the recorded steps follow
exactly the same logic as `mapa()` - this module only adds instrumentation
around that loop, it does not reimplement the algorithm.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "reference" / "python"))

from mapa import Bin, CalibratedBin, _merge, _merge_calibrated, _not_significant, _violates, _violates_pd  # noqa: E402

from schemas import PdStep, PdStepBin, Step, StepBin  # noqa: E402


def _snapshot(stack: List[Bin]) -> List[StepBin]:
    return [
        StepBin(
            score_min=b.score_min, score_max=b.score_max,
            n_obs=b.n_obs, n_bads=b.n_bads, bad_rate=b.bad_rate,
            count=b.count, count_bads=b.count_bads,
        )
        for b in stack
    ]


def _pd_snapshot(stack: List[CalibratedBin]) -> List[PdStepBin]:
    return [
        PdStepBin(
            score_min=b.score_min, score_max=b.score_max,
            n_obs=b.n_obs, n_bads=b.n_bads, pd=b.pd,
            count=b.count, count_bads=b.count_bads,
        )
        for b in stack
    ]


def pooling_steps(bins: List[Bin], increasing: bool, min_confidence: Optional[float]) -> List[Step]:
    stack: List[Bin] = []
    steps: List[Step] = []

    for b in bins:
        stack.append(b)
        steps.append(Step(action="push", stack=_snapshot(stack)))

        while len(stack) >= 2:
            violates = _violates(stack[-2], stack[-1], increasing)
            not_significant = min_confidence is not None and _not_significant(stack[-2], stack[-1], min_confidence)
            if not (violates or not_significant):
                break

            top = stack.pop()
            below = stack.pop()

            if violates:
                reason = (
                    f"Merged bands at {below.score_min:g}-{below.score_max:g} and "
                    f"{top.score_min:g}-{top.score_max:g}: bad rate "
                    f"{'decreased' if increasing else 'increased'} from "
                    f"{below.bad_rate:.3f} to {top.bad_rate:.3f}, violating monotonicity"
                )
            else:
                reason = (
                    f"Merged bands at {below.score_min:g}-{below.score_max:g} and "
                    f"{top.score_min:g}-{top.score_max:g}: bad rates "
                    f"{below.bad_rate:.3f} and {top.bad_rate:.3f} not significantly "
                    f"different at the requested confidence level"
                )

            stack.append(_merge(below, top))
            steps.append(Step(action="merge", stack=_snapshot(stack), reason=reason))

    return steps


def minimum_size_steps(bins: List[Bin], min_obs: float, min_bads: float, use_counts: bool = True) -> List[Step]:
    """Trace the size-violation merges from `enforce_minimum_size`'s while loop.

    Mirrors the loop in reference/python/mapa.py's `enforce_minimum_size`,
    using `_merge` so the merge result matches exactly. Does not include the
    final `mapa` re-pooling pass - call `pooling_steps` on the returned bins
    for that.
    """
    bins = list(bins)
    steps: List[Step] = [Step(action="push", stack=_snapshot(bins))]

    def _obs(b: Bin) -> float:
        return b.count if use_counts else b.n_obs

    def _bads(b: Bin) -> float:
        return b.count_bads if use_counts else b.n_bads

    while len(bins) > 1:
        violator = next(
            (i for i, b in enumerate(bins) if _obs(b) < min_obs or _bads(b) < min_bads),
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

        obs_label = "count" if use_counts else "n_obs"
        bads_label = "count_bads" if use_counts else "n_bads"
        if _obs(bins[violator]) < min_obs:
            shortfall = obs_label
            shortfall_value = _obs(bins[violator])
            threshold = min_obs
        else:
            shortfall = bads_label
            shortfall_value = _bads(bins[violator])
            threshold = min_bads
        reason = (
            f"Band {bins[violator].score_min:g}-{bins[violator].score_max:g} has "
            f"{shortfall}={shortfall_value:g} (below the minimum of {threshold:g}); "
            f"merged into its closer-rate neighbour at "
            f"{bins[neighbour].score_min:g}-{bins[neighbour].score_max:g}"
        )

        i, j = sorted((violator, neighbour))
        bins = bins[:i] + [_merge(bins[i], bins[j])] + bins[j + 1 :]
        steps.append(Step(action="merge", stack=_snapshot(bins), reason=reason))

    return steps


def repool_pd_steps(bins: List[CalibratedBin], increasing: bool) -> List[PdStep]:
    """Trace `repool_calibrated_bins`'s re-pooling of Bayesian-adjusted `pd`
    values, using `_violates_pd` and `_merge_calibrated` so the merge
    decisions match exactly.
    """
    stack: List[CalibratedBin] = []
    steps: List[PdStep] = []

    for b in bins:
        stack.append(b)
        steps.append(PdStep(action="push", stack=_pd_snapshot(stack)))

        while len(stack) >= 2 and _violates_pd(stack[-2], stack[-1], increasing):
            top = stack.pop()
            below = stack.pop()

            reason = (
                f"Re-pooled bands at {below.score_min:g}-{below.score_max:g} and "
                f"{top.score_min:g}-{top.score_max:g}: adjusted PD "
                f"{'decreased' if increasing else 'increased'} from "
                f"{below.pd:.3f} to {top.pd:.3f} after Bayesian shrinkage, "
                f"violating monotonicity"
            )

            stack.append(_merge_calibrated(below, top))
            steps.append(PdStep(action="merge", stack=_pd_snapshot(stack), reason=reason))

    return steps
