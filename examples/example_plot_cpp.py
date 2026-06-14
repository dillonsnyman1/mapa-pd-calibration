"""Plot the unsmoothed (step) and smoothed (log-odds interpolated) PD
mappings produced by the C++ reference's run_pipeline().

Requires the C++ generator to have been run first (see cpp/generate_output.cpp)
to produce output/cpp_bands.csv and output/cpp_smoothed.csv.

Usage:
    pip install -r requirements.txt
    cd cpp && g++ -std=c++17 -I../../reference/cpp generate_output.cpp \
        ../../reference/cpp/mapa.cpp -o generate_output && ./generate_output && cd ..
    python3 example_plot_cpp.py
"""

import csv
from pathlib import Path

from plotting import plot_calibration

OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def main():
    with open(OUTPUT_DIR / "cpp_bands.csv", newline="") as f:
        bands = [
            (float(row["score_min"]), float(row["score_max"]), float(row["pd"]))
            for row in csv.DictReader(f)
        ]

    with open(OUTPUT_DIR / "cpp_smoothed.csv", newline="") as f:
        smoothed = [(float(row["score"]), float(row["pd"])) for row in csv.DictReader(f)]

    out_path = OUTPUT_DIR / "cpp_calibration.png"
    plot_calibration(bands, smoothed, "MAPA score-to-PD calibration (C++)", out_path)


if __name__ == "__main__":
    main()
