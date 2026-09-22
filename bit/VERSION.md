# Version

## Versioning rule

The archive follows semantic versioning of its *reproduction interface and manuscript identity*:
a major version changes when the supported manuscript, the entry point, or the artifact set changes;
a minor version adds results without changing existing ones; a patch version changes prose, packaging,
or tests only. Historical version directories (`results/v26` … `results/v34`, `results/jcp_figures_v34`,
`figures_jcp`) and the historical entry point `reproduce.sh` are preserved unchanged across versions.

## 4.1.1 — corrections release (2026-09-22)

Corrections to the mathematical text of the article and of Online Resource 1, fixes to the figures, and
clearer documentation. Certified values and data records are unchanged. The documentation now states how to
reproduce the results from a clean environment and what a stand-alone run reports. The cover letter is no
longer part of the reproduction manifest, so its date can change without a new archive. The package builder
writes the plain-text cover letter from its LaTeX source with paragraph breaks, sets the line numbers of the
review copy further from the text, no longer writes file-name comments into the flattened article source, and
dates every packaged file 22 September 2026. The flattened article source in `manuscript_source/` is
regenerated from the corrected article. The archive is cited as
https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.1.1
and by the Zenodo concept DOI 10.5281/zenodo.21522057.

## 4.1.0 — figure restyle and prose release (2026-09-19)

The seven figures were restyled for the typeset article (Computer Modern lettering, solid fills instead of
hatching, consistent line weights and panel labels); the plotted data, series and panels are unchanged, and
the figure contract tests were updated for the new export settings. No certified value, data record,
verification script or table changed; the figure script changed with the restyle. The article and Online
Resource prose was edited; the flattened article source shipped in `manuscript_source/` was regenerated. The
archive is cited as https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.1.0
and by the Zenodo concept DOI 10.5281/zenodo.21522057; Zenodo archives this version as 10.5281/zenodo.22846581.

## 4.0.1 — documentation patch (2026-09-18)

No change to code, tests, certified data, figures or tables. The article's Code availability statement,
this archive's citation metadata and documentation now cite the versioned release
https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.0.1
and the Zenodo concept DOI 10.5281/zenodo.21522057, which resolves to the latest archived version
(version 4.0.0 is archived as 10.5281/zenodo.22837831). Zenodo did not archive this documentation-only
release. The flattened article source shipped in `manuscript_source/` was regenerated from the updated article.

## 4.0.0 — version for the BIT Numerical Mathematics article (2026-09-18)

Supports the article *Complete-positivity regions of Runge–Kutta discretizations of phase-covariant qubit
dynamics* and its Online Resource 1. The mathematical results, certified data, and controller
implementation of 3.4/3.5 are unchanged. Added: `scripts/verify_bit_claims.py` (independent first-principles
verification of the article's central claims, 69 checks), `src/rk_choi_margin/extrapolation.py` and
`scripts/verify_bit_extension.py` (n-fold Richardson-extrapolation corollary of the signed first-defect law),
`scripts/generate_bit_figures.py` and `scripts/generate_bit_tables.py` (Springer-format figures and tables
with overlap, grayscale, and colour-vision audits), `scripts/sync_article_numbers.py`,
`scripts/audit_bit_source.py`, `scripts/build_bit_package.py`, the cross-platform entry point
`reproduce_bit.py`, the contract tests `tests/test_bit_*.py`, and `BIT_MANUSCRIPT_CONTRACT.md`.
Corrected in the manuscript (not in the code): the attached-branch/ceiling crossover of the RK4
classification is stated explicitly (q_c ≈ 3.5282, ϖ_c ≈ 0.7950).

This version also adds corrected rational proof brackets, the complete Richardson sign
chart, and a whole-interval slack hypothesis for the persistence result. It also strengthens input
provenance and packaging checks and edits the article-facing prose. Certified numerical values and
historical scientific inputs are unchanged. The public release is available at
https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.0.0.

## 3.5.0 — editorial-polish release (2026-08-20)

This patch release retains the v3.4 mathematical results, certified data, and controller implementation.
It revises the manuscript-facing prose and graphical abstract of the earlier manuscript version,
synchronizes the source contract, and updates the source audit without changing any theorem, numerical
value, or scientific scope. The buffered-component panel reserves a non-overlapping legend band and carries
a render-level regression test.
