# MAPA: Monotone Adjacent Pooling Algorithm for PD Calibration

A reference implementation and interactive demo of the **Monotone Adjacent
Pooling Algorithm (MAPA)** - a method for turning empirical bad rates per
score band into a monotone score-to-PD calibration curve.

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
  implementations of the core algorithm in **Python**, **C++** and **SAS**,
  validated against shared fixture data. Intended as a readable starting
  point for anyone wanting to implement or port MAPA themselves.
- [`examples/`](examples/) - scripts that plot the resulting score-to-PD
  calibration curves (unsmoothed and smoothed) for the Python and C++
  implementations.
- `app/` - *(planned)* a full-stack interactive demo (FastAPI + React)
  visualizing the pooling process and the resulting calibration curve on a
  synthetic scored population.

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

### SAS

See [`reference/sas/README.md`](reference/sas/README.md).
