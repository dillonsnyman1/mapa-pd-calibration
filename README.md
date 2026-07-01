# MAPA: Monotone Adjacent Pooling Algorithm for PD Calibration

[![CI/CD](https://github.com/dillonsnyman1/mapa-pd-calibration/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/dillonsnyman1/mapa-pd-calibration/actions/workflows/ci-cd.yml)

A reference implementation and interactive demo of the **Monotone Adjacent
Pooling Algorithm (MAPA)** - a method for turning empirical bad rates per
score band into a monotone score-to-PD calibration curve. Supports both
number-weighted (standard) and value-weighted observations (e.g. for
IFRS 9 exposure-weighted PD calibration).

> **Attribution**: MAPA, as implemented here, is based on the method
> described by **Raymond Anderson** in *The Credit Scoring Toolkit* (Oxford
> University Press, 2007). All credit for the underlying methodology belongs
> to Anderson; this repository is an independent, from-scratch
> implementation for educational and portfolio purposes, using synthetic
> data, and is not derived from or representative of any proprietary model
> or implementation. See [`docs/mapa-methodology.md`](docs/mapa-methodology.md)
> for full attribution.

---

## What's here

- [`docs/mapa-methodology.md`](docs/mapa-methodology.md) - explanation of
  the algorithm, its relationship to PAVA / isotonic regression, and
  attribution.
- [`reference/`](reference/) - clean, dependency-light, side-by-side
  implementations of the core algorithm in **Python**, **C++**, **R**,
  **MATLAB / GNU Octave** and **SAS**, all validated against the same shared
  fixture data. Intended as a readable starting point for anyone wanting to
  implement or port MAPA themselves.
- [`examples/`](examples/) - scripts that plot the resulting score-to-PD
  calibration curves (unsmoothed and smoothed) for the Python and C++
  implementations.
- [`app/`](app/) - a full-stack interactive demo (FastAPI + React)
  visualizing the pooling process and the resulting calibration curve on a
  synthetic scored population, with adjustable parameters and CSV upload.

## Methodology

For background and attribution, see
[`docs/mapa-methodology.md`](docs/mapa-methodology.md). In short: MAPA is
the Pool Adjacent Violators Algorithm applied to score-band bad rates,
merging adjacent bands until the resulting bad rate sequence is monotone.
It is based on the method described by Raymond Anderson in *The Credit
Scoring Toolkit* (Oxford University Press, 2007).

## Quickstart

### Python

```bash
cd reference/python
pip install pytest
pytest
```

### C++

```bash
cd reference/cpp
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

### R

```bash
cd reference/r
Rscript -e "install.packages('testthat', repos='https://cloud.r-project.org')"
Rscript test_mapa.R
```

### MATLAB / GNU Octave

**MATLAB:** Add `reference/matlab/` to your path, then run `test_mapa` from the command window.

**GNU Octave** (free, open source — [octave.org](https://octave.org)):

```bash
pkg install -forge datatypes   # one-time setup
cd reference/matlab
octave --no-gui test_mapa.m
```

### SAS

See [`reference/sas/README.md`](reference/sas/README.md).

### Live demo

> **Live demo**: [dcg14fdv56g8g.cloudfront.net](https://dcg14fdv56g8g.cloudfront.net)
>
> The backend is fully stateless - each request is processed and returned in one go, with no data written to disk or stored anywhere.

See [`app/README.md`](app/README.md) to run the interactive FastAPI + React
demo locally.
