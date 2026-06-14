// Runs the C++ reference's run_pipeline() on the shared fixtures and writes
// the resulting band table and a densely-sampled smoothed PD curve to CSV,
// for plotting by example_plot_cpp.py.
//
// Build (from this directory):
//   g++ -std=c++17 -I../../reference/cpp generate_output.cpp \
//       ../../reference/cpp/mapa.cpp -o generate_output
//
// Run:
//   ./generate_output

#include "mapa.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::pair<double, int>> load_observations(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Could not open fixture file: " + path);
    }

    std::vector<std::pair<double, int>> observations;
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

        std::stringstream ss(line);
        std::string score_str, bad_str;
        std::getline(ss, score_str, ',');
        std::getline(ss, bad_str, ',');
        observations.emplace_back(std::stod(score_str), std::stoi(bad_str));
    }
    return observations;
}

}  // namespace

int main() {
    const std::string fixtures_dir = "../../reference/fixtures";
    const std::string output_dir = "../output";

    const auto observations = load_observations(fixtures_dir + "/raw_observations.csv");
    const auto result = mapa::run_pipeline(observations, /*k=*/10.0, /*min_obs=*/50, /*min_bads=*/10);

    std::ofstream bands_file(output_dir + "/cpp_bands.csv");
    bands_file << std::setprecision(15);
    bands_file << "score_min,score_max,pd\n";
    for (const auto& b : result.bands) {
        bands_file << b.score_min << "," << b.score_max << "," << b.pd << "\n";
    }
    std::cout << "Wrote " << output_dir << "/cpp_bands.csv\n";

    std::ofstream smoothed_file(output_dir + "/cpp_smoothed.csv");
    smoothed_file << std::setprecision(15);
    smoothed_file << "score,pd\n";
    const double score_min = result.bands.front().score_min;
    const double score_max = result.bands.back().score_max;
    constexpr int kNumPoints = 501;
    for (int i = 0; i < kNumPoints; ++i) {
        const double score = score_min + i * (score_max - score_min) / (kNumPoints - 1);
        smoothed_file << score << "," << result.pd_for_score(score) << "\n";
    }
    std::cout << "Wrote " << output_dir << "/cpp_smoothed.csv\n";

    return 0;
}
