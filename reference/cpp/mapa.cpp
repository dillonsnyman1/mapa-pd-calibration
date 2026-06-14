#include "mapa.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace mapa {

namespace {

bool violates(const Bin& lower, const Bin& upper, bool increasing) {
    if (increasing) {
        return upper.bad_rate() < lower.bad_rate();
    }
    return upper.bad_rate() > lower.bad_rate();
}

Bin merge(const Bin& a, const Bin& b) {
    return Bin{a.score_min, b.score_max, a.n_obs + b.n_obs, a.n_bads + b.n_bads};
}

bool violates_pd(const CalibratedBin& lower, const CalibratedBin& upper, bool increasing) {
    if (increasing) {
        return upper.pd < lower.pd;
    }
    return upper.pd > lower.pd;
}

CalibratedBin merge_calibrated(const CalibratedBin& a, const CalibratedBin& b) {
    long n_obs = a.n_obs + b.n_obs;
    long n_bads = a.n_bads + b.n_bads;
    double pd = (a.pd * static_cast<double>(a.n_obs) + b.pd * static_cast<double>(b.n_obs)) /
                static_cast<double>(n_obs);
    return CalibratedBin{a.score_min, b.score_max, n_obs, n_bads, pd};
}

// Approximates the inverse of the standard normal CDF (the quantile
// function / probit), via Acklam's rational approximation. Used to convert
// a confidence level (e.g. 0.95) into a critical z-value for
// not_significant(). Accurate to about 1.15e-9.
double inverse_normal_cdf(double p) {
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                                1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                                6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                                -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                                3.754408661907416e+00};

    const double p_low = 0.02425;
    const double p_high = 1.0 - p_low;

    if (p < p_low) {
        double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p <= p_high) {
        double q = p - 0.5;
        double r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
    double q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
}

// Whether the bad rates of `a` and `b` are statistically indistinguishable
// at the given confidence level, via a two-proportion z-test. Returns true
// (i.e. "merge these") when the observed difference in bad rates is small
// relative to its standard error - either because the rates are genuinely
// close, or because both bins are too small to tell them apart.
bool not_significant(const Bin& a, const Bin& b, double confidence) {
    double pooled_rate = static_cast<double>(a.n_bads + b.n_bads) / static_cast<double>(a.n_obs + b.n_obs);
    if (pooled_rate <= 0.0 || pooled_rate >= 1.0) {
        return true;
    }

    double se = std::sqrt(pooled_rate * (1.0 - pooled_rate) *
                           (1.0 / static_cast<double>(a.n_obs) + 1.0 / static_cast<double>(b.n_obs)));
    double z = std::fabs(a.bad_rate() - b.bad_rate()) / se;
    double z_critical = inverse_normal_cdf((1.0 + confidence) / 2.0);
    return z < z_critical;
}

}  // namespace

std::vector<Bin> bins_from_observations(const std::vector<std::pair<double, int>>& observations) {
    // std::map keeps entries sorted by score, giving us ascending order for free.
    std::map<double, std::pair<long, long>> counts;
    for (const auto& [score, bad] : observations) {
        auto& entry = counts[score];
        entry.first += 1;
        entry.second += bad;
    }

    std::vector<Bin> bins;
    bins.reserve(counts.size());
    for (const auto& [score, count] : counts) {
        bins.push_back(Bin{score, score, count.first, count.second});
    }
    return bins;
}

std::vector<Bin> mapa(const std::vector<Bin>& bins, bool increasing,
                       std::optional<double> min_confidence) {
    std::vector<Bin> stack;
    stack.reserve(bins.size());

    for (const auto& b : bins) {
        stack.push_back(b);
        while (stack.size() >= 2 &&
               (violates(stack[stack.size() - 2], stack.back(), increasing) ||
                (min_confidence.has_value() &&
                 not_significant(stack[stack.size() - 2], stack.back(), *min_confidence)))) {
            Bin top = stack.back();
            stack.pop_back();
            stack.back() = merge(stack.back(), top);
        }
    }

    return stack;
}

std::vector<Bin> calibrate(const std::vector<std::pair<double, int>>& observations,
                            bool increasing, std::optional<double> min_confidence) {
    return mapa(bins_from_observations(observations), increasing, min_confidence);
}

std::vector<Bin> enforce_minimum_size(const std::vector<Bin>& input, long min_obs, long min_bads,
                                       bool increasing, std::optional<double> min_confidence) {
    std::vector<Bin> bins = input;

    while (bins.size() > 1) {
        size_t violator = bins.size();
        for (size_t i = 0; i < bins.size(); ++i) {
            if (bins[i].n_obs < min_obs || bins[i].n_bads < min_bads) {
                violator = i;
                break;
            }
        }
        if (violator == bins.size()) {
            break;
        }

        size_t neighbour;
        if (violator == 0) {
            neighbour = 1;
        } else if (violator == bins.size() - 1) {
            neighbour = violator - 1;
        } else {
            double rate = bins[violator].bad_rate();
            double left_diff = std::fabs(rate - bins[violator - 1].bad_rate());
            double right_diff = std::fabs(rate - bins[violator + 1].bad_rate());
            neighbour = (left_diff <= right_diff) ? violator - 1 : violator + 1;
        }

        size_t i = std::min(violator, neighbour);
        size_t j = std::max(violator, neighbour);
        bins[i] = merge(bins[i], bins[j]);
        bins.erase(bins.begin() + static_cast<long>(j));
    }

    return mapa(bins, increasing, min_confidence);
}

std::vector<CalibratedBin> apply_bayesian_adjustment(const std::vector<Bin>& bins, double k,
                                                      std::optional<double> prior) {
    double p0;
    if (prior.has_value()) {
        p0 = *prior;
    } else {
        long total_obs = 0, total_bads = 0;
        for (const auto& b : bins) {
            total_obs += b.n_obs;
            total_bads += b.n_bads;
        }
        p0 = static_cast<double>(total_bads) / static_cast<double>(total_obs);
    }

    std::vector<CalibratedBin> result;
    result.reserve(bins.size());
    for (const auto& b : bins) {
        double pd = (static_cast<double>(b.n_bads) + k * p0) / (static_cast<double>(b.n_obs) + k);
        result.push_back(CalibratedBin{b.score_min, b.score_max, b.n_obs, b.n_bads, pd});
    }
    return result;
}

std::vector<CalibratedBin> repool_calibrated_bins(const std::vector<CalibratedBin>& bins,
                                                   bool increasing) {
    std::vector<CalibratedBin> stack;
    stack.reserve(bins.size());

    for (const auto& b : bins) {
        stack.push_back(b);
        while (stack.size() >= 2 &&
               violates_pd(stack[stack.size() - 2], stack.back(), increasing)) {
            CalibratedBin top = stack.back();
            stack.pop_back();
            stack.back() = merge_calibrated(stack.back(), top);
        }
    }

    return stack;
}

double interpolate_pd(const std::vector<CalibratedBin>& bins, double score) {
    std::vector<double> mids;
    std::vector<double> log_odds;
    mids.reserve(bins.size());
    log_odds.reserve(bins.size());
    for (const auto& b : bins) {
        mids.push_back((b.score_min + b.score_max) / 2.0);
        log_odds.push_back(std::log((1.0 - b.pd) / b.pd));
    }

    if (score <= mids.front()) {
        return bins.front().pd;
    }
    if (score >= mids.back()) {
        return bins.back().pd;
    }

    for (size_t i = 0; i + 1 < mids.size(); ++i) {
        if (mids[i] <= score && score <= mids[i + 1]) {
            double t = (score - mids[i]) / (mids[i + 1] - mids[i]);
            double interpolated = log_odds[i] + t * (log_odds[i + 1] - log_odds[i]);
            return 1.0 / (1.0 + std::exp(interpolated));
        }
    }

    // Unreachable: score must lie between mids.front() and mids.back() here.
    return bins.back().pd;
}

double CalibrationResult::pd_for_score(double score) const { return interpolate_pd(bands, score); }

CalibrationResult run_pipeline(const std::vector<std::pair<double, int>>& observations, double k,
                                long min_obs, long min_bads, std::optional<double> prior,
                                bool increasing, std::optional<double> min_confidence) {
    std::vector<Bin> pooled = calibrate(observations, increasing, min_confidence);
    std::vector<Bin> sized = enforce_minimum_size(pooled, min_obs, min_bads, increasing, min_confidence);
    std::vector<CalibratedBin> calibrated = apply_bayesian_adjustment(sized, k, prior);
    std::vector<CalibratedBin> repooled = repool_calibrated_bins(calibrated, increasing);
    return CalibrationResult{std::move(repooled)};
}

}  // namespace mapa
