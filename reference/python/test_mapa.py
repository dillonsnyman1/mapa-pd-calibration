import csv
import math
from pathlib import Path

from mapa import (
    Bin,
    CalibratedBin,
    apply_bayesian_adjustment,
    bins_from_observations,
    calibrate,
    enforce_minimum_size,
    interpolate_pd,
    mapa,
    repool_calibrated_bins,
    run_pipeline,
)

FIXTURES_DIR = Path(__file__).resolve().parent.parent / "fixtures"
BAYESIAN_K = 10
MIN_OBS = 50
MIN_BADS = 10
MIN_CONFIDENCE = 0.95


def _load_observations() -> list[tuple[float, int]]:
    with open(FIXTURES_DIR / "raw_observations.csv", newline="") as f:
        reader = csv.DictReader(f)
        return [(float(row["score"]), int(row["bad"])) for row in reader]


def _load_weighted_observations() -> list[tuple[float, int, float]]:
    with open(FIXTURES_DIR / "raw_observations_weighted.csv", newline="") as f:
        reader = csv.DictReader(f)
        return [(float(row["score"]), int(row["bad"]), float(row["weight"])) for row in reader]


def _load_bins(filename: str) -> list[Bin]:
    with open(FIXTURES_DIR / filename, newline="") as f:
        reader = csv.DictReader(f)
        return [
            Bin(float(row["score_min"]), float(row["score_max"]), int(row["n_obs"]), int(row["n_bads"]))
            for row in reader
        ]


def _load_weighted_bins(filename: str) -> list[Bin]:
    with open(FIXTURES_DIR / filename, newline="") as f:
        reader = csv.DictReader(f)
        return [
            Bin(
                float(row["score_min"]), float(row["score_max"]),
                float(row["n_obs"]), float(row["n_bads"]),
                int(row["count"]), int(row["count_bads"]),
            )
            for row in reader
        ]


def test_bins_from_observations_matches_expected():
    result = bins_from_observations(_load_observations())
    expected = _load_bins("expected_initial_bins.csv")

    assert result == expected


def test_calibrate_matches_expected_pooled_bins():
    result = calibrate(_load_observations())
    expected = _load_bins("expected_pooled_bins.csv")

    assert result == expected


def test_result_is_monotone_non_increasing():
    result = calibrate(_load_observations())

    rates = [b.bad_rate for b in result]
    assert rates == sorted(rates, reverse=True)


def test_pooling_preserves_totals():
    observations = _load_observations()
    initial = bins_from_observations(observations)
    result = mapa(initial)

    assert sum(b.n_obs for b in result) == sum(b.n_obs for b in initial)
    assert sum(b.n_bads for b in result) == sum(b.n_bads for b in initial)
    assert sum(b.n_obs for b in result) == len(observations)


def test_enforce_minimum_size_matches_expected():
    pooled = calibrate(_load_observations())

    result = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    expected = _load_bins("expected_min_size_bins.csv")

    assert result == expected


def test_enforce_minimum_size_satisfies_thresholds():
    pooled = calibrate(_load_observations())

    result = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)

    # A single remaining bin may still fall short if the whole population
    # can't meet the thresholds.
    if len(result) > 1:
        for b in result:
            assert b.n_obs >= MIN_OBS
            assert b.n_bads >= MIN_BADS


def test_enforce_minimum_size_preserves_totals_and_monotonicity():
    pooled = calibrate(_load_observations())

    result = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)

    assert sum(b.n_obs for b in result) == sum(b.n_obs for b in pooled)
    assert sum(b.n_bads for b in result) == sum(b.n_bads for b in pooled)

    rates = [b.bad_rate for b in result]
    assert rates == sorted(rates, reverse=True)


def test_enforce_minimum_size_is_noop_with_default_thresholds():
    pooled = calibrate(_load_observations())

    assert enforce_minimum_size(pooled) == pooled


def test_bayesian_adjustment_matches_expected():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)

    result = apply_bayesian_adjustment(sized, k=BAYESIAN_K)

    with open(FIXTURES_DIR / "expected_calibrated_bins.csv", newline="") as f:
        expected = list(csv.DictReader(f))

    assert len(result) == len(expected)
    for calibrated, row in zip(result, expected):
        assert calibrated.score_min == float(row["score_min"])
        assert calibrated.score_max == float(row["score_max"])
        assert calibrated.n_obs == int(row["n_obs"])
        assert calibrated.n_bads == int(row["n_bads"])
        assert math.isclose(calibrated.pd, float(row["pd"]), rel_tol=1e-9)


def test_bayesian_adjustment_shrinks_toward_prior():
    pooled = calibrate(_load_observations())
    prior = sum(b.n_bads for b in pooled) / sum(b.n_obs for b in pooled)

    result = apply_bayesian_adjustment(pooled, k=BAYESIAN_K, prior=prior)

    # Every adjusted PD must lie strictly between the bin's own bad rate
    # and the prior (or equal one of them, for k == 0).
    for original, adjusted in zip(pooled, result):
        lo, hi = sorted((original.bad_rate, prior))
        assert lo <= adjusted.pd <= hi


def _load_calibrated_bins(filename: str) -> list[CalibratedBin]:
    with open(FIXTURES_DIR / filename, newline="") as f:
        reader = csv.DictReader(f)
        return [
            CalibratedBin(
                float(row["score_min"]),
                float(row["score_max"]),
                int(row["n_obs"]),
                int(row["n_bads"]),
                float(row["pd"]),
            )
            for row in reader
        ]


def test_repool_calibrated_bins_matches_expected():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    calibrated = apply_bayesian_adjustment(sized, k=BAYESIAN_K)

    result = repool_calibrated_bins(calibrated)
    expected = _load_calibrated_bins("expected_repooled_calibrated_bins.csv")

    assert len(result) == len(expected)
    for r, e in zip(result, expected):
        assert r.score_min == e.score_min
        assert r.score_max == e.score_max
        assert r.n_obs == e.n_obs
        assert r.n_bads == e.n_bads
        assert math.isclose(r.pd, e.pd, rel_tol=1e-9)


def test_repool_calibrated_bins_restores_monotonicity():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    calibrated = apply_bayesian_adjustment(sized, k=BAYESIAN_K)

    # The bundled fixture deliberately crosses after Bayesian adjustment.
    pds = [b.pd for b in calibrated]
    assert pds != sorted(pds, reverse=True)

    result = repool_calibrated_bins(calibrated)
    pds = [b.pd for b in result]
    assert pds == sorted(pds, reverse=True)


def test_repool_calibrated_bins_preserves_totals():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    calibrated = apply_bayesian_adjustment(sized, k=BAYESIAN_K)

    result = repool_calibrated_bins(calibrated)

    assert sum(b.n_obs for b in result) == sum(b.n_obs for b in calibrated)
    assert sum(b.n_bads for b in result) == sum(b.n_bads for b in calibrated)


def test_interpolate_pd_matches_expected():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    calibrated = apply_bayesian_adjustment(sized, k=BAYESIAN_K)
    repooled = repool_calibrated_bins(calibrated)

    with open(FIXTURES_DIR / "expected_smoothed_pds.csv", newline="") as f:
        expected = list(csv.DictReader(f))

    for row in expected:
        score = float(row["score"])
        assert math.isclose(interpolate_pd(repooled, score), float(row["pd"]), rel_tol=1e-9)


def test_interpolate_pd_is_monotone_non_increasing():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    calibrated = apply_bayesian_adjustment(sized, k=BAYESIAN_K)
    repooled = repool_calibrated_bins(calibrated)

    scores = [b.score_min for b in repooled] + [b.score_max for b in repooled]
    pds = [interpolate_pd(repooled, s) for s in sorted(set(scores))]

    assert pds == sorted(pds, reverse=True)


def test_interpolate_pd_at_pool_midpoint_matches_pool_pd():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    calibrated = apply_bayesian_adjustment(sized, k=BAYESIAN_K)
    repooled = repool_calibrated_bins(calibrated)

    for b in repooled:
        midpoint = (b.score_min + b.score_max) / 2
        assert math.isclose(interpolate_pd(repooled, midpoint), b.pd, rel_tol=1e-9)


def test_run_pipeline_bands_match_repooled_calibrated_bins():
    pooled = calibrate(_load_observations())
    sized = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS)
    calibrated = apply_bayesian_adjustment(sized, k=BAYESIAN_K)
    repooled = repool_calibrated_bins(calibrated)

    result = run_pipeline(_load_observations(), k=BAYESIAN_K, min_obs=MIN_OBS, min_bads=MIN_BADS)

    assert result.bands == repooled


def test_run_pipeline_pd_for_score_matches_interpolate_pd():
    result = run_pipeline(_load_observations(), k=BAYESIAN_K, min_obs=MIN_OBS, min_bads=MIN_BADS)

    with open(FIXTURES_DIR / "expected_smoothed_pds.csv", newline="") as f:
        expected = list(csv.DictReader(f))

    for row in expected:
        score = float(row["score"])
        assert math.isclose(result.pd_for_score(score), interpolate_pd(result.bands, score), rel_tol=1e-9)
        assert math.isclose(result.pd_for_score(score), float(row["pd"]), rel_tol=1e-9)


def test_mapa_min_confidence_matches_expected():
    initial = bins_from_observations(_load_observations())

    result = mapa(initial, min_confidence=MIN_CONFIDENCE)
    expected = _load_bins("expected_pooled_bins_confidence.csv")

    assert result == expected


def test_mapa_min_confidence_preserves_totals_and_monotonicity():
    initial = bins_from_observations(_load_observations())

    result = mapa(initial, min_confidence=MIN_CONFIDENCE)

    assert sum(b.n_obs for b in result) == sum(b.n_obs for b in initial)
    assert sum(b.n_bads for b in result) == sum(b.n_bads for b in initial)

    rates = [b.bad_rate for b in result]
    assert rates == sorted(rates, reverse=True)


def test_mapa_min_confidence_merges_at_least_as_much_as_plain_mapa():
    initial = bins_from_observations(_load_observations())

    plain = mapa(initial)
    confident = mapa(initial, min_confidence=MIN_CONFIDENCE)

    assert len(confident) <= len(plain)


def test_increasing_direction():
    # Flip each score's sign so that "higher score = lower risk" becomes
    # "higher score = higher risk", and check the non-decreasing variant.
    observations = [(-score, bad) for score, bad in _load_observations()]

    result = calibrate(observations, increasing=True)

    rates = [b.bad_rate for b in result]
    assert rates == sorted(rates)
    assert sum(b.n_obs for b in result) == len(observations)


# ---------------------------------------------------------------------------
# Value-weighted tests
# ---------------------------------------------------------------------------


def test_weighted_bins_from_observations_matches_expected():
    result = bins_from_observations(_load_weighted_observations())
    expected = _load_weighted_bins("expected_initial_bins_weighted.csv")

    assert len(result) == len(expected)
    for r, e in zip(result, expected):
        assert r.score_min == e.score_min
        assert r.score_max == e.score_max
        assert math.isclose(r.n_obs, e.n_obs, rel_tol=1e-9)
        assert math.isclose(r.n_bads, e.n_bads, rel_tol=1e-9)
        assert r.count == e.count
        assert r.count_bads == e.count_bads


def test_weighted_bins_differ_from_unweighted():
    weighted = _load_weighted_observations()
    unweighted = [(s, b) for s, b, _ in weighted]

    w_bins = bins_from_observations(weighted)
    u_bins = bins_from_observations(unweighted)

    assert len(w_bins) == len(u_bins)
    diffs = sum(1 for w, u in zip(w_bins, u_bins) if not math.isclose(w.bad_rate, u.bad_rate, rel_tol=1e-9))
    assert diffs > 0


def test_weighted_pooling_preserves_totals():
    observations = _load_weighted_observations()
    initial = bins_from_observations(observations)
    result = mapa(initial)

    assert math.isclose(sum(b.n_obs for b in result), sum(b.n_obs for b in initial), rel_tol=1e-9)
    assert math.isclose(sum(b.n_bads for b in result), sum(b.n_bads for b in initial), rel_tol=1e-9)
    assert sum(b.count for b in result) == sum(b.count for b in initial)
    assert sum(b.count_bads for b in result) == sum(b.count_bads for b in initial)


def test_weighted_pooling_is_monotone():
    initial = bins_from_observations(_load_weighted_observations())
    result = mapa(initial)

    rates = [b.bad_rate for b in result]
    assert rates == sorted(rates, reverse=True)


def test_weighted_pooling_matches_expected():
    initial = bins_from_observations(_load_weighted_observations())
    result = mapa(initial)
    expected = _load_weighted_bins("expected_pooled_bins_weighted.csv")

    assert len(result) == len(expected)
    for r, e in zip(result, expected):
        assert r.score_min == e.score_min
        assert r.score_max == e.score_max
        assert math.isclose(r.n_obs, e.n_obs, rel_tol=1e-9)
        assert math.isclose(r.n_bads, e.n_bads, rel_tol=1e-9)
        assert r.count == e.count
        assert r.count_bads == e.count_bads


def test_weighted_enforce_minimum_size_uses_counts():
    initial = bins_from_observations(_load_weighted_observations())
    pooled = mapa(initial)

    result = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS, use_counts=True)

    if len(result) > 1:
        for b in result:
            assert b.count >= MIN_OBS
            assert b.count_bads >= MIN_BADS


def test_weighted_enforce_minimum_size_matches_expected():
    initial = bins_from_observations(_load_weighted_observations())
    pooled = mapa(initial)
    result = enforce_minimum_size(pooled, min_obs=MIN_OBS, min_bads=MIN_BADS, use_counts=True)
    expected = _load_weighted_bins("expected_min_size_bins_weighted.csv")

    assert len(result) == len(expected)
    for r, e in zip(result, expected):
        assert r.score_min == e.score_min
        assert r.score_max == e.score_max
        assert math.isclose(r.n_obs, e.n_obs, rel_tol=1e-9)
        assert math.isclose(r.n_bads, e.n_bads, rel_tol=1e-9)
        assert r.count == e.count
        assert r.count_bads == e.count_bads


def test_weighted_run_pipeline_matches_expected():
    observations = _load_weighted_observations()
    result = run_pipeline(observations, k=BAYESIAN_K, min_obs=MIN_OBS, min_bads=MIN_BADS, use_counts=True)

    with open(FIXTURES_DIR / "expected_repooled_calibrated_bins_weighted.csv", newline="") as f:
        expected = list(csv.DictReader(f))

    assert len(result.bands) == len(expected)
    for band, row in zip(result.bands, expected):
        assert band.score_min == float(row["score_min"])
        assert band.score_max == float(row["score_max"])
        assert math.isclose(band.n_obs, float(row["n_obs"]), rel_tol=1e-9)
        assert math.isclose(band.n_bads, float(row["n_bads"]), rel_tol=1e-9)
        assert band.count == int(row["count"])
        assert band.count_bads == int(row["count_bads"])
        assert math.isclose(band.pd, float(row["pd"]), rel_tol=1e-9)


def test_weighted_smoothed_pds_match_expected():
    result = run_pipeline(
        _load_weighted_observations(), k=BAYESIAN_K, min_obs=MIN_OBS, min_bads=MIN_BADS, use_counts=True
    )

    with open(FIXTURES_DIR / "expected_smoothed_pds_weighted.csv", newline="") as f:
        expected = list(csv.DictReader(f))

    for row in expected:
        score = float(row["score"])
        assert math.isclose(result.pd_for_score(score), float(row["pd"]), rel_tol=1e-9)
