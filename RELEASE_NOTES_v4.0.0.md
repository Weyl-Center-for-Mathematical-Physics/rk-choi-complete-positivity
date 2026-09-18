# v4.0.0 — BIT reproducibility release

This release supports *Complete-positivity regions of Runge-Kutta discretizations of phase-covariant
qubit dynamics*, by G. Blake Pierpoint, Olivier Bernard, and Yichen Liu.

It adds the reviewed BIT reproducibility tree under `bit/`: exact symbolic claim checks, Richardson
extrapolation analysis, candidate-map and endpoint certificates, retained benchmark records,
deterministic figures and tables, and regression tests. The root historical source/data/figure tree and
the v1.0.0 tag/release are preserved. The earlier release is archived at
[10.5281/zenodo.21522058](https://doi.org/10.5281/zenodo.21522058); that DOI is not presented as a DOI for v4.0.0.

The attached `BIT_Code_and_Data.zip` contains the exact reviewed archive. The same 408 files appear under
`bit/`, preserved byte-for-byte. Historical manuscript-contract and flattened article-source fixtures
are intentionally included for reproducibility tests. Private revision reports, development histories,
virtual environments, caches and journal-submission ZIPs are not included.

Use `python reproduce_bit.py all` from the extracted archive or `bit/`, after installing the pinned
environment described in its README. Workspace-only manuscript comparisons are explicitly skipped when
the separate manuscript directory is absent. This is a public reproducibility release, not a journal
submission.

Release URL:
https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.0.0

License: MIT. The exact asset SHA-256 and final verification totals are included in the GitHub release body.
