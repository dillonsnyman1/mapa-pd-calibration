#pragma once

#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace mapa {

// A score band spanning [score_min, score_max], with its observation and
// bad counts.
//
// n_obs / n_bads are weighted sums (equal to count / count_bads when every
// observation has weight 1).  count / count_bads are raw observation counts,
// used for statistical tests (z-test sample sizes) and minimum-size checks
// when use_counts is true.
struct Bin {
    double score_min;
    double score_max;
    double n_obs;
    double n_bads;
    long count;
    long count_bads;

    double bad_rate() const { return n_bads / n_obs; }
};

// A pooled bin together with its Bayesian-adjusted PD, as produced by
// apply_bayesian_adjustment().
struct CalibratedBin {
    double score_min;
    double score_max;
    double n_obs;
    double n_bads;
    long count;
    long count_bads;
    double pd;
};

// The output of run_pipeline(): a calibrated band table, plus a smoothed
// per-score PD curve derived from it.
//
// `bands` is the step-function calibration table - the typical deliverable
// for reporting and governance. `pd_for_score()` gives a smoothed,
// continuous PD for an individual score, via interpolate_pd(). Both are
// views of the same underlying bands; which one to use depends on the
// consumer.
struct CalibrationResult {
    std::vector<CalibratedBin> bands;

    double pd_for_score(double score) const;
};

// Groups raw (score, bad, weight) observations into one bin per unique score,
// ordered by score ascending. This is the finest possible starting point
// for mapa(). `bad` should be 1 for a default and 0 otherwise. `weight` is
// the value weight of the observation (e.g. exposure at default).
std::vector<Bin> bins_from_observations(const std::vector<std::tuple<double, int, double>>& observations);

// Convenience overload for unweighted observations (weight = 1 for all).
std::vector<Bin> bins_from_observations(const std::vector<std::pair<double, int>>& observations);

// Monotone Adjacent Pooling Algorithm.
//
// `bins` must be ordered by score, ascending (e.g. from
// bins_from_observations()). By default the bad rate (PD) is required to
// be non-increasing as score increases - the standard credit scoring
// convention (higher score = lower risk). Pass `increasing = true` to
// require a non-decreasing bad rate instead.
//
// Repeatedly merges adjacent bins that violate the required monotonicity
// until the resulting sequence of bad rates is monotone. The result
// partitions the full score range and population of the input bins.
//
// If `min_confidence` is given (e.g. 0.95 for 95%), adjacent bins whose bad
// rates do not differ at this confidence level (two-proportion z-test) are
// merged as well, even if they don't violate monotonicity. This produces
// fewer, larger bins whose bad rates are more reliably distinguishable from
// their neighbours. If not given (the default), only monotonicity
// violations are merged.
//
// See ../../docs/mapa-methodology.md for background and attribution.
std::vector<Bin> mapa(const std::vector<Bin>& bins, bool increasing = false,
                       std::optional<double> min_confidence = std::nullopt);

// Convenience wrapper: group raw weighted observations into per-score bins,
// then pool them with mapa().
std::vector<Bin> calibrate(const std::vector<std::tuple<double, int, double>>& observations,
                            bool increasing = false,
                            std::optional<double> min_confidence = std::nullopt);

// Convenience overload for unweighted observations (weight = 1 for all).
std::vector<Bin> calibrate(const std::vector<std::pair<double, int>>& observations,
                            bool increasing = false,
                            std::optional<double> min_confidence = std::nullopt);

// Further pool bins that don't meet minimum size thresholds, even if they
// don't violate monotonicity.
//
// When `use_counts` is true (the default), a bin "violates" if
// count < min_obs or count_bads < min_bads (raw observation counts).
// When false, the weighted sums n_obs and n_bads are compared instead.
//
// Each violating bin is repeatedly merged into whichever adjacent bin has
// the closer bad rate (minimizing the distortion to the calibration curve),
// until every remaining bin meets both thresholds or only one bin remains.
//
// Merging toward the closer-rate neighbour is not guaranteed to preserve
// the monotonicity established by mapa(), so the result is passed back
// through mapa() before being returned. `increasing` is forwarded to that
// final pass; see mapa().
//
// `min_confidence` is forwarded to that final pass; see mapa().
//
// `bins` is typically the output of mapa().
std::vector<Bin> enforce_minimum_size(const std::vector<Bin>& bins, double min_obs = 0,
                                       double min_bads = 0, bool increasing = false,
                                       std::optional<double> min_confidence = std::nullopt,
                                       bool use_counts = true);

// Shrink each bin's empirical bad rate toward a prior using Bayesian
// (credibility) weighting:
//
//     pd = (n_bads + k * prior) / (n_obs + k)
//
// `k` is the credibility weight, expressed as a number of "equivalent
// observations" of the prior: a bin with n_obs == k is shrunk halfway
// between its own bad rate and the prior, bins with n_obs >> k are barely
// adjusted, and bins with n_obs << k end up close to the prior.
//
// `bins` is typically the output of mapa(). If `prior` is not given, it
// defaults to the overall bad rate across all bins
// (sum(n_bads) / sum(n_obs)).
//
// Note: shrinking each bin independently toward a single global prior is
// not guaranteed to preserve the monotonicity established by mapa(): a
// small bin whose empirical rate is close to a much larger neighbour's can
// be pulled past it toward the prior. See ../../docs/mapa-methodology.md
// for discussion.
std::vector<CalibratedBin> apply_bayesian_adjustment(const std::vector<Bin>& bins, double k,
                                                      std::optional<double> prior = std::nullopt);

// Re-applies pooling to Bayesian-adjusted bins, restoring monotonicity of
// `pd`.
//
// apply_bayesian_adjustment() shrinks each bin's bad rate toward a shared
// prior independently, which is not guaranteed to preserve the
// monotonicity established by mapa() (see its comment). This runs the same
// adjacent-pooling algorithm again, but on `pd` instead of `bad_rate`,
// merging violating bins by taking the n_obs-weighted average of their `pd`
// values. `increasing` has the same meaning as in mapa(), applied to `pd`.
//
// `bins` is typically the output of apply_bayesian_adjustment().
std::vector<CalibratedBin> repool_calibrated_bins(const std::vector<CalibratedBin>& bins,
                                                   bool increasing = false);

// Returns a smoothed PD for an individual score via log-odds interpolation
// between pools.
//
// The pooled, calibrated PD curve is a step function: every score within a
// pool gets that pool's single PD, with a discontinuous jump at each pool
// boundary. This produces a smooth, continuous PD instead.
//
// Each pool is represented by a single anchor point: its midpoint score,
// (score_min + score_max) / 2, and the log-odds of its pd,
// log_odds = ln((1 - pd) / pd). For the given `score`:
//
// - if it is at or before the first pool's midpoint, or at or after the
//   last pool's midpoint, the nearest pool's pd is returned unchanged (flat
//   extrapolation beyond the anchor points).
// - otherwise, log_odds is linearly interpolated between the midpoints of
//   the two pools bracketing `score`, and converted back to a PD via
//   pd = 1 / (1 + exp(log_odds)).
//
// Because log-odds is a monotonic transform of pd, this preserves the
// monotonicity of a monotone input sequence of pool PDs.
//
// `bins` is typically the output of repool_calibrated_bins().
double interpolate_pd(const std::vector<CalibratedBin>& bins, double score);

// Runs the full MAPA pipeline: bin, pool, enforce minimum size, apply
// Bayesian adjustment, and re-pool.
//
// This chains calibrate(), enforce_minimum_size(),
// apply_bayesian_adjustment() and repool_calibrated_bins(). The result
// bundles the resulting band table (`bands`) together with a smoothed,
// continuous PD curve derived from it (`pd_for_score()`, via
// interpolate_pd()) - use whichever representation suits the consumer.
CalibrationResult run_pipeline(const std::vector<std::tuple<double, int, double>>& observations,
                                double k, double min_obs = 0, double min_bads = 0,
                                std::optional<double> prior = std::nullopt,
                                bool increasing = false,
                                std::optional<double> min_confidence = std::nullopt,
                                bool use_counts = true);

// Convenience overload for unweighted observations (weight = 1 for all).
CalibrationResult run_pipeline(const std::vector<std::pair<double, int>>& observations, double k,
                                double min_obs = 0, double min_bads = 0,
                                std::optional<double> prior = std::nullopt,
                                bool increasing = false,
                                std::optional<double> min_confidence = std::nullopt,
                                bool use_counts = true);

}  // namespace mapa
