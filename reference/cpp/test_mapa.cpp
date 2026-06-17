#include "mapa.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

// Minimal CSV reader for the fixed-schema fixture files used by this test.
// Skips the header row; does not handle quoting or embedded commas, which
// the fixtures never contain.
std::vector<std::vector<std::string>> read_csv(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Could not open fixture file: " + path);
    }

    std::vector<std::vector<std::string>> rows;
    std::string line;
    bool first = true;
    while (std::getline(file, line)) {
        if (first) {
            first = false;  // header row
            continue;
        }
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ',')) {
            fields.push_back(field);
        }
        rows.push_back(fields);
    }
    return rows;
}

std::vector<std::pair<double, int>> load_observations(const std::string& path) {
    std::vector<std::pair<double, int>> observations;
    for (const auto& row : read_csv(path)) {
        observations.emplace_back(std::stod(row[0]), std::stoi(row[1]));
    }
    return observations;
}

std::vector<mapa::Bin> load_bins(const std::string& path) {
    std::vector<mapa::Bin> bins;
    for (const auto& row : read_csv(path)) {
        bins.push_back(mapa::Bin{std::stod(row[0]), std::stod(row[1]), std::stol(row[2]),
                                  std::stol(row[3])});
    }
    return bins;
}

struct ExpectedCalibratedBin {
    double score_min;
    double score_max;
    long n_obs;
    long n_bads;
    double pd;
};

std::vector<ExpectedCalibratedBin> load_calibrated_bins(const std::string& path) {
    std::vector<ExpectedCalibratedBin> bins;
    for (const auto& row : read_csv(path)) {
        bins.push_back(ExpectedCalibratedBin{std::stod(row[0]), std::stod(row[1]),
                                               std::stol(row[2]), std::stol(row[3]),
                                               std::stod(row[4])});
    }
    return bins;
}

struct ExpectedScorePd {
    double score;
    double pd;
};

std::vector<ExpectedScorePd> load_score_pds(const std::string& path) {
    std::vector<ExpectedScorePd> rows;
    for (const auto& row : read_csv(path)) {
        rows.push_back(ExpectedScorePd{std::stod(row[0]), std::stod(row[1])});
    }
    return rows;
}

bool bins_equal(const mapa::Bin& a, const mapa::Bin& b) {
    return a.score_min == b.score_min && a.score_max == b.score_max && a.n_obs == b.n_obs &&
           a.n_bads == b.n_bads;
}

std::vector<std::tuple<double, int, double>> load_weighted_observations(const std::string& path) {
    std::vector<std::tuple<double, int, double>> observations;
    for (const auto& row : read_csv(path)) {
        observations.emplace_back(std::stod(row[0]), std::stoi(row[1]), std::stod(row[2]));
    }
    return observations;
}

struct ExpectedWeightedBin {
    double score_min, score_max, n_obs, n_bads;
    long count, count_bads;
};

std::vector<ExpectedWeightedBin> load_weighted_bins(const std::string& path) {
    std::vector<ExpectedWeightedBin> bins;
    for (const auto& row : read_csv(path)) {
        bins.push_back(ExpectedWeightedBin{
            std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]),
            std::stol(row[4]), std::stol(row[5])
        });
    }
    return bins;
}

struct ExpectedWeightedCalibratedBin {
    double score_min, score_max, n_obs, n_bads;
    long count, count_bads;
    double pd;
};

std::vector<ExpectedWeightedCalibratedBin> load_weighted_calibrated_bins(const std::string& path) {
    std::vector<ExpectedWeightedCalibratedBin> bins;
    for (const auto& row : read_csv(path)) {
        bins.push_back(ExpectedWeightedCalibratedBin{
            std::stod(row[0]), std::stod(row[1]), std::stod(row[2]), std::stod(row[3]),
            std::stol(row[4]), std::stol(row[5]), std::stod(row[6])
        });
    }
    return bins;
}

}  // namespace

int main() {
    const std::string fixtures_dir = FIXTURES_DIR;

    const auto observations = load_observations(fixtures_dir + "/raw_observations.csv");
    const auto expected_initial = load_bins(fixtures_dir + "/expected_initial_bins.csv");
    const auto expected_pooled = load_bins(fixtures_dir + "/expected_pooled_bins.csv");

    // bins_from_observations matches expected per-score bins.
    const auto initial = mapa::bins_from_observations(observations);
    assert(initial.size() == expected_initial.size());
    for (size_t i = 0; i < initial.size(); ++i) {
        assert(bins_equal(initial[i], expected_initial[i]));
    }

    // calibrate() (binning + pooling) matches expected pooled bins.
    const auto result = mapa::calibrate(observations);
    assert(result.size() == expected_pooled.size());
    for (size_t i = 0; i < result.size(); ++i) {
        assert(bins_equal(result[i], expected_pooled[i]));
    }

    // The resulting bad rate sequence must be non-increasing.
    for (size_t i = 1; i < result.size(); ++i) {
        assert(result[i].bad_rate() <= result[i - 1].bad_rate());
    }

    // Pooling must preserve totals.
    long total_obs = 0, total_bads = 0;
    for (const auto& b : result) {
        total_obs += b.n_obs;
        total_bads += b.n_bads;
    }
    assert(total_obs == static_cast<long>(observations.size()));

    long expected_bads = 0;
    for (const auto& [score, bad] : observations) {
        expected_bads += bad;
    }
    assert(total_bads == expected_bads);

    // enforce_minimum_size matches expected min-size bins.
    constexpr long kMinObs = 50;
    constexpr long kMinBads = 10;
    const auto expected_min_size = load_bins(fixtures_dir + "/expected_min_size_bins.csv");
    const auto min_size_result = mapa::enforce_minimum_size(result, kMinObs, kMinBads);

    assert(min_size_result.size() == expected_min_size.size());
    for (size_t i = 0; i < min_size_result.size(); ++i) {
        assert(bins_equal(min_size_result[i], expected_min_size[i]));
    }

    // Thresholds are satisfied (more than one bin remains) and totals/monotonicity preserved.
    if (min_size_result.size() > 1) {
        for (const auto& b : min_size_result) {
            assert(b.n_obs >= kMinObs);
            assert(b.n_bads >= kMinBads);
        }
    }
    long min_size_obs = 0, min_size_bads = 0;
    for (const auto& b : min_size_result) {
        min_size_obs += b.n_obs;
        min_size_bads += b.n_bads;
    }
    assert(min_size_obs == total_obs);
    assert(min_size_bads == total_bads);
    for (size_t i = 1; i < min_size_result.size(); ++i) {
        assert(min_size_result[i].bad_rate() <= min_size_result[i - 1].bad_rate());
    }

    // Default thresholds are a no-op.
    const auto noop_result = mapa::enforce_minimum_size(result);
    assert(noop_result.size() == result.size());
    for (size_t i = 0; i < result.size(); ++i) {
        assert(bins_equal(noop_result[i], result[i]));
    }

    // Bayesian adjustment matches expected calibrated PDs.
    constexpr double kBayesianK = 10.0;
    const auto expected_calibrated = load_calibrated_bins(fixtures_dir + "/expected_calibrated_bins.csv");
    const auto calibrated = mapa::apply_bayesian_adjustment(min_size_result, kBayesianK);

    assert(calibrated.size() == expected_calibrated.size());
    for (size_t i = 0; i < calibrated.size(); ++i) {
        const auto& c = calibrated[i];
        const auto& e = expected_calibrated[i];
        assert(c.score_min == e.score_min);
        assert(c.score_max == e.score_max);
        assert(c.n_obs == e.n_obs);
        assert(c.n_bads == e.n_bads);
        assert(std::fabs(c.pd - e.pd) < 1e-9);
    }

    // Each adjusted PD lies between the bin's own bad rate and the prior.
    long prior_total_obs = 0, prior_total_bads = 0;
    for (const auto& b : min_size_result) {
        prior_total_obs += b.n_obs;
        prior_total_bads += b.n_bads;
    }
    const double prior = static_cast<double>(prior_total_bads) / static_cast<double>(prior_total_obs);
    for (size_t i = 0; i < min_size_result.size(); ++i) {
        const double rate = min_size_result[i].bad_rate();
        const double lo = std::min(rate, prior);
        const double hi = std::max(rate, prior);
        assert(calibrated[i].pd >= lo - 1e-12 && calibrated[i].pd <= hi + 1e-12);
    }

    // The bundled fixture's calibrated PDs deliberately cross between two
    // adjacent bins.
    bool has_crossing = false;
    for (size_t i = 1; i < calibrated.size(); ++i) {
        if (calibrated[i].pd > calibrated[i - 1].pd) {
            has_crossing = true;
            break;
        }
    }
    assert(has_crossing);

    // repool_calibrated_bins matches the expected re-pooled output and
    // restores monotonicity.
    const auto expected_repooled =
        load_calibrated_bins(fixtures_dir + "/expected_repooled_calibrated_bins.csv");
    const auto repooled = mapa::repool_calibrated_bins(calibrated);

    assert(repooled.size() == expected_repooled.size());
    for (size_t i = 0; i < repooled.size(); ++i) {
        const auto& c = repooled[i];
        const auto& e = expected_repooled[i];
        assert(c.score_min == e.score_min);
        assert(c.score_max == e.score_max);
        assert(c.n_obs == e.n_obs);
        assert(c.n_bads == e.n_bads);
        assert(std::fabs(c.pd - e.pd) < 1e-9);
    }

    for (size_t i = 1; i < repooled.size(); ++i) {
        assert(repooled[i].pd <= repooled[i - 1].pd);
    }

    long repooled_obs = 0, repooled_bads = 0;
    for (const auto& b : repooled) {
        repooled_obs += b.n_obs;
        repooled_bads += b.n_bads;
    }
    assert(repooled_obs == total_obs);
    assert(repooled_bads == total_bads);

    // interpolate_pd matches the expected smoothed PDs.
    const auto expected_smoothed = load_score_pds(fixtures_dir + "/expected_smoothed_pds.csv");
    for (const auto& row : expected_smoothed) {
        double pd = mapa::interpolate_pd(repooled, row.score);
        assert(std::fabs(pd - row.pd) < 1e-9);
    }

    // interpolate_pd at each pool's midpoint returns that pool's own pd.
    for (const auto& b : repooled) {
        double midpoint = (b.score_min + b.score_max) / 2.0;
        assert(std::fabs(mapa::interpolate_pd(repooled, midpoint) - b.pd) < 1e-9);
    }

    // interpolate_pd is monotone non-increasing across the pool boundaries.
    {
        std::vector<double> scores;
        for (const auto& b : repooled) {
            scores.push_back(b.score_min);
            scores.push_back(b.score_max);
        }
        std::sort(scores.begin(), scores.end());
        scores.erase(std::unique(scores.begin(), scores.end()), scores.end());

        double prev = mapa::interpolate_pd(repooled, scores.front());
        for (size_t i = 1; i < scores.size(); ++i) {
            double pd = mapa::interpolate_pd(repooled, scores[i]);
            assert(pd <= prev);
            prev = pd;
        }
    }

    // run_pipeline() bundles both representations: bands matches
    // repool_calibrated_bins() output, and pd_for_score() matches
    // interpolate_pd().
    {
        const auto pipeline_result = mapa::run_pipeline(observations, kBayesianK, kMinObs, kMinBads);
        assert(pipeline_result.bands.size() == repooled.size());
        for (size_t i = 0; i < repooled.size(); ++i) {
            const auto& a = pipeline_result.bands[i];
            const auto& b = repooled[i];
            assert(a.score_min == b.score_min);
            assert(a.score_max == b.score_max);
            assert(a.n_obs == b.n_obs);
            assert(a.n_bads == b.n_bads);
            assert(std::fabs(a.pd - b.pd) < 1e-9);
        }
        for (const auto& row : expected_smoothed) {
            assert(std::fabs(pipeline_result.pd_for_score(row.score) - row.pd) < 1e-9);
        }
    }

    // mapa() with min_confidence merges adjacent bins whose bad rates are
    // not statistically distinguishable, in addition to monotonicity
    // violations.
    {
        const auto expected_confidence = load_bins(fixtures_dir + "/expected_pooled_bins_confidence.csv");
        const auto confident = mapa::mapa(initial, /*increasing=*/false, /*min_confidence=*/0.95);

        assert(confident.size() == expected_confidence.size());
        for (size_t i = 0; i < confident.size(); ++i) {
            assert(bins_equal(confident[i], expected_confidence[i]));
        }

        long confident_obs = 0, confident_bads = 0;
        for (const auto& b : confident) {
            confident_obs += b.n_obs;
            confident_bads += b.n_bads;
        }
        assert(confident_obs == total_obs);
        assert(confident_bads == total_bads);

        for (size_t i = 1; i < confident.size(); ++i) {
            assert(confident[i].bad_rate() <= confident[i - 1].bad_rate());
        }

        assert(confident.size() <= result.size());
    }

    // ---------------------------------------------------------------
    // Weighted test cases
    // ---------------------------------------------------------------

    const auto weighted_obs = load_weighted_observations(fixtures_dir + "/raw_observations_weighted.csv");
    const auto expected_initial_weighted = load_weighted_bins(fixtures_dir + "/expected_initial_bins_weighted.csv");

    // Weighted bins_from_observations matches expected per-score bins.
    {
        const auto w_initial = mapa::bins_from_observations(weighted_obs);
        assert(w_initial.size() == expected_initial_weighted.size());
        for (size_t i = 0; i < w_initial.size(); ++i) {
            const auto& a = w_initial[i];
            const auto& e = expected_initial_weighted[i];
            assert(a.score_min == e.score_min);
            assert(a.score_max == e.score_max);
            assert(std::fabs(a.n_obs - e.n_obs) < 1e-6);
            assert(std::fabs(a.n_bads - e.n_bads) < 1e-6);
            assert(a.count == e.count);
            assert(a.count_bads == e.count_bads);
        }
    }

    // Weighted bins differ from unweighted (at least one bad rate differs).
    {
        std::vector<std::pair<double, int>> unweighted_pairs;
        for (const auto& [score, bad, weight] : weighted_obs) {
            unweighted_pairs.emplace_back(score, bad);
        }
        const auto w_bins = mapa::bins_from_observations(weighted_obs);
        const auto uw_bins = mapa::bins_from_observations(unweighted_pairs);
        assert(w_bins.size() == uw_bins.size());
        bool any_differ = false;
        for (size_t i = 0; i < w_bins.size(); ++i) {
            if (std::fabs(w_bins[i].bad_rate() - uw_bins[i].bad_rate()) > 1e-12) {
                any_differ = true;
                break;
            }
        }
        assert(any_differ);
    }

    // Weighted pooling preserves totals.
    {
        const auto w_initial = mapa::bins_from_observations(weighted_obs);
        const auto w_pooled = mapa::calibrate(weighted_obs);

        double init_n_obs = 0, init_n_bads = 0;
        long init_count = 0, init_count_bads = 0;
        for (const auto& b : w_initial) {
            init_n_obs += b.n_obs;
            init_n_bads += b.n_bads;
            init_count += b.count;
            init_count_bads += b.count_bads;
        }

        double pooled_n_obs = 0, pooled_n_bads = 0;
        long pooled_count = 0, pooled_count_bads = 0;
        for (const auto& b : w_pooled) {
            pooled_n_obs += b.n_obs;
            pooled_n_bads += b.n_bads;
            pooled_count += b.count;
            pooled_count_bads += b.count_bads;
        }

        assert(std::fabs(pooled_n_obs - init_n_obs) < 1e-6);
        assert(std::fabs(pooled_n_bads - init_n_bads) < 1e-6);
        assert(pooled_count == init_count);
        assert(pooled_count_bads == init_count_bads);
    }

    // Weighted pooling is monotone (non-increasing bad rate).
    {
        const auto w_pooled = mapa::calibrate(weighted_obs);
        for (size_t i = 1; i < w_pooled.size(); ++i) {
            assert(w_pooled[i].bad_rate() <= w_pooled[i - 1].bad_rate());
        }
    }

    // Weighted pooling matches expected.
    {
        const auto expected_pooled_weighted = load_weighted_bins(fixtures_dir + "/expected_pooled_bins_weighted.csv");
        const auto w_pooled = mapa::calibrate(weighted_obs);
        assert(w_pooled.size() == expected_pooled_weighted.size());
        for (size_t i = 0; i < w_pooled.size(); ++i) {
            const auto& a = w_pooled[i];
            const auto& e = expected_pooled_weighted[i];
            assert(a.score_min == e.score_min);
            assert(a.score_max == e.score_max);
            assert(std::fabs(a.n_obs - e.n_obs) < 1e-6);
            assert(std::fabs(a.n_bads - e.n_bads) < 1e-6);
            assert(a.count == e.count);
            assert(a.count_bads == e.count_bads);
        }
    }

    // Weighted enforce_minimum_size uses counts (use_counts=true).
    {
        const auto w_pooled = mapa::calibrate(weighted_obs);
        const auto w_min_size = mapa::enforce_minimum_size(w_pooled, kMinObs, kMinBads,
                                                            /*increasing=*/false,
                                                            /*min_confidence=*/std::nullopt,
                                                            /*use_counts=*/true);
        if (w_min_size.size() > 1) {
            for (const auto& b : w_min_size) {
                assert(b.count >= kMinObs);
                assert(b.count_bads >= kMinBads);
            }
        }
    }

    // Weighted enforce_minimum_size matches expected.
    {
        const auto expected_min_size_weighted = load_weighted_bins(fixtures_dir + "/expected_min_size_bins_weighted.csv");
        const auto w_pooled = mapa::calibrate(weighted_obs);
        const auto w_min_size = mapa::enforce_minimum_size(w_pooled, kMinObs, kMinBads,
                                                            /*increasing=*/false,
                                                            /*min_confidence=*/std::nullopt,
                                                            /*use_counts=*/true);
        assert(w_min_size.size() == expected_min_size_weighted.size());
        for (size_t i = 0; i < w_min_size.size(); ++i) {
            const auto& a = w_min_size[i];
            const auto& e = expected_min_size_weighted[i];
            assert(a.score_min == e.score_min);
            assert(a.score_max == e.score_max);
            assert(std::fabs(a.n_obs - e.n_obs) < 1e-6);
            assert(std::fabs(a.n_bads - e.n_bads) < 1e-6);
            assert(a.count == e.count);
            assert(a.count_bads == e.count_bads);
        }
    }

    // Weighted run_pipeline matches expected repooled calibrated bins.
    {
        const auto expected_repooled_weighted =
            load_weighted_calibrated_bins(fixtures_dir + "/expected_repooled_calibrated_bins_weighted.csv");
        const auto w_pipeline = mapa::run_pipeline(weighted_obs, kBayesianK, kMinObs, kMinBads,
                                                    /*prior=*/std::nullopt,
                                                    /*increasing=*/false,
                                                    /*min_confidence=*/std::nullopt,
                                                    /*use_counts=*/true);
        assert(w_pipeline.bands.size() == expected_repooled_weighted.size());
        for (size_t i = 0; i < w_pipeline.bands.size(); ++i) {
            const auto& a = w_pipeline.bands[i];
            const auto& e = expected_repooled_weighted[i];
            assert(a.score_min == e.score_min);
            assert(a.score_max == e.score_max);
            assert(std::fabs(a.n_obs - e.n_obs) < 1e-6);
            assert(std::fabs(a.n_bads - e.n_bads) < 1e-6);
            assert(a.count == e.count);
            assert(a.count_bads == e.count_bads);
            assert(std::fabs(a.pd - e.pd) < 1e-9);
        }
    }

    // Weighted smoothed PDs match expected.
    {
        const auto expected_smoothed_weighted = load_score_pds(fixtures_dir + "/expected_smoothed_pds_weighted.csv");
        const auto w_pipeline = mapa::run_pipeline(weighted_obs, kBayesianK, kMinObs, kMinBads,
                                                    /*prior=*/std::nullopt,
                                                    /*increasing=*/false,
                                                    /*min_confidence=*/std::nullopt,
                                                    /*use_counts=*/true);
        for (const auto& row : expected_smoothed_weighted) {
            double pd = w_pipeline.pd_for_score(row.score);
            assert(std::fabs(pd - row.pd) < 1e-9);
        }
    }

    std::cout << "All tests passed (unweighted: " << result.size()
              << " pooled bins; weighted tests included).\n";
    return 0;
}
