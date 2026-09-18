Code and data archive, version 4.0.1 (BIT Numerical Mathematics revision, 2026-09-18).

Article: Complete-positivity regions of Runge-Kutta discretizations of
phase-covariant qubit dynamics (not submitted).

Install the pinned environment with:
  python -m pip install -e ".[test]"

Run the full reproduction with (PowerShell or POSIX):
  python reproduce_bit.py all

Or execute modular stages:
  certificates | tests | audits | figures | manuscript-contract

The historical v3.5 chain remains available: bash reproduce.sh all

The scientific results and certified data of versions 3.4/3.5 are unchanged in
4.0; the revision adds an independent claim verifier, the Richardson-
extrapolation corollary, BIT-format figures and tables with accessibility
audits, contract tests, and the cross-platform entry point. Evidence classes
(exact symbolic, certified rational/interval, high precision, binary64,
regression test, visual) are explained in README.md and mapped to the article
in BIT_MANUSCRIPT_CONTRACT.md.

Versioned public release:
https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.0.1
Zenodo concept DOI (latest version): https://doi.org/10.5281/zenodo.21522057

The authors approved the manuscript and declarations and confirmed that the
work is neither under consideration nor accepted elsewhere. Publication of
this reproducibility archive is separate from journal submission.
