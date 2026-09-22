Code and data archive, version 4.1.2 (22 September 2026), for the article

  Complete-positivity regions of Runge-Kutta discretizations of
  phase-covariant qubit dynamics
  G. Blake Pierpoint, Olivier Bernard, and Yichen Liu
  (prepared for BIT Numerical Mathematics)

Reproduce in a clean virtual environment, from the folder that contains this
file (PowerShell: use .venv\Scripts\python.exe in place of .venv/bin/python):
  python -m venv .venv
  .venv/bin/python -m pip install -e ".[test]"
  .venv/bin/python reproduce_bit.py all

Expected result: 4 stages pass and the manuscript-contract stage is skipped
(it needs the article's LaTeX sources, which are not included); 384 tests
pass and 4 are skipped for the same reason.

Single stages:
  certificates | tests | audits | figures | manuscript-contract

The historical v3.5 chain remains available: bash reproduce.sh all

The scientific results and certified data of versions 3.4/3.5 are unchanged
in version 4, which adds an independent claim verifier, the Richardson-
extrapolation corollary, figures and tables in the Springer format with
accessibility audits, contract tests, and the cross-platform entry point.
Evidence classes (exact symbolic, certified rational/interval, high precision,
binary64, regression test, visual) are explained in README.md and mapped to
the article in BIT_MANUSCRIPT_CONTRACT.md.

Versioned public release:
https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.1.2
Zenodo concept DOI (all versions): https://doi.org/10.5281/zenodo.21522057
