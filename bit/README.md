# Code and data archive

**Article:** Complete-positivity regions of Runge–Kutta discretizations of phase-covariant qubit dynamics
**Article prepared for:** BIT Numerical Mathematics
**Archive version:** 4.1.1, 22 September 2026 (see `VERSION.md`)
**Authors:** G. Blake Pierpoint, Olivier Bernard, and Yichen Liu
**Corresponding author:** G. Blake Pierpoint, pierpogb@odu.edu
**License:** MIT (`LICENSE`); citation metadata in `CITATION.cff`

**Versioned release:** https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.1.1
**Archived:** Zenodo, concept DOI https://doi.org/10.5281/zenodo.21522057 (all versions; resolves to the latest)

This archive contains the exact symbolic calculations, certified numerical records, deterministic figure
and table generators, and regression tests supporting the article and its Online Resource 1. Version 4 of
the archive accompanies the BIT Numerical Mathematics article: the mathematical results and certified data
of versions 3.4/3.5 are unchanged, an independent first-principles verifier and the Richardson-extrapolation
corollary were added, the figures and tables were regenerated for the Springer format, and a cross-platform
reproduction entry point (`reproduce_bit.py`) was added. Historical version directories (`results/v26` …
`results/v34`, `results/jcp_figures_v34`, `figures_jcp`) are retained as provenance and are never edited.

The Choi-block criterion applies to the autonomous phase-covariant qubit family. The signed first-defect
law and its extrapolation corollary apply to real-coefficient Runge–Kutta methods on the boundary rays.
The six-regime classification concerns classical RK4 on the laboratory-frame one-way, zero-dephasing
ray. Four transverse-field couplings, three of them nonzero, have separate pointwise certificates. The
adaptive experiments measure physicality, error and operation counts; they do not measure runtime efficiency.

## Reproducing the results

Requirements: Python 3.11 or newer (the pinned NumPy and SciPy wheels require 3.11); the validated
environment is CPython 3.13.5 with the pins in `pyproject.toml` (SymPy 1.14.0, NumPy 2.3.5, SciPy 1.16.1,
mpmath 1.3.0, Matplotlib 3.10.8; tests: pytest 9.0.2, pypdf 5.9.0). No LaTeX installation is needed.

Work in the folder that contains this README: the folder `BIT_Code_and_Data/` of the unpacked release asset
`BIT_Code_and_Data.zip`, or `bit/` in the repository. Create a clean virtual environment, install the
archive with its pinned dependencies, and run every stage.

PowerShell (Windows):

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe .\reproduce_bit.py all
```

POSIX shell (Linux/macOS):

```bash
python3 -m venv .venv
./.venv/bin/python -m pip install -e ".[test]"
./.venv/bin/python reproduce_bit.py all
```

**Expected result.** Four stages pass (`certificates`, `tests`, `audits`, `figures`) and the
`manuscript-contract` stage is skipped, because it audits the article's LaTeX sources, which this archive does
not contain. The test stage reports 384 tests passed and 4 skipped; the four skipped tests also compare against
the article's LaTeX sources. A full run takes a few minutes, and its last line is

```text
==> stages: 4 passed, 0 failed, 1 skipped; tests passed: 384
```

`reproduce_bit.py` runs the stages below in this order, prints every command it executes, stops at the
first failure with that command's exit code, writes a timestamped log per stage to
`results/bit_revision/logs/`, and (for `all`) writes `results/bit_revision/manifest.json` with pass/fail
counts, the number of tests collected and passed, and the SHA-256 of every generated artifact and frozen
input. Run a single stage with `reproduce_bit.py <stage>`; `reproduce_bit.py all --list` prints the commands
without running them.

| Stage | What it runs | Evidence class |
|---|---|---|
| `certificates` | `scripts/verify_bit_claims.py` (69 SymPy checks of the central claims, rebuilt from Butcher tableaux and the generator without importing the package) and `scripts/verify_bit_extension.py` (extrapolation corollary for eight methods, n = 2, 3, 4) | symbolic, certified |
| `tests` | every pytest module in its own process (16 historical modules + the BIT contract modules), after a collection step that counts the tests | regression |
| `audits` | the retained v2.x–v3.3 audit scripts: theorem identities, phase-covariant identities, method audit, 15,000-sample Choi audit, executed-candidate contract, endpoint-domain adversarial audit | symbolic, floating point |
| `figures` | `scripts/generate_bit_figures.py` (Fig1–Fig3, FigS1–FigS4 with overlap, grayscale, and colour-vision audits) and `scripts/generate_bit_tables.py` (Online Resource tables) | visual, deterministic regeneration |
| `manuscript-contract` | `scripts/audit_bit_source.py` on the article's LaTeX sources, passed with `--manuscript` or found in `../manuscript` (skipped otherwise) | mechanical source audit |

Direct commands (PowerShell shown; use `./.venv/bin/python` on POSIX):

```powershell
.\.venv\Scripts\python.exe -m pytest -q                                  # whole suite in one process
.\.venv\Scripts\python.exe .\scripts\generate_bit_figures.py            # regenerate figures_bit/ and the figure manifest
.\.venv\Scripts\python.exe .\scripts\generate_bit_tables.py             # regenerate results/bit_revision/tables/
.\.venv\Scripts\python.exe .\scripts\verify_bit_claims.py               # independent claim verification
.\.venv\Scripts\python.exe .\reproduce_bit.py manuscript-contract --manuscript <folder with the article's main.tex>
```

The four skipped tests and the `manuscript-contract` stage run when the article's LaTeX sources are
available: `reproduce_bit.py` uses the folder given with `--manuscript`, or a folder `manuscript/` next to
this archive folder, where the tests also look. Three of those tests compare the Online Resource tables with
their typeset copies, and one audits the article source; pytest lists them as skipped with `-rs`. The five
article-table tests of `tests/test_bit_table_contract.py` read the typeset Tables 1 and 2 from
`manuscript_source/main.tex` (the flattened article source) and therefore run without the LaTeX sources.

The historical v3.5 chain is retained as `bash reproduce.sh all`, with the known issue described below.
Building the journal package (`scripts/build_bit_package.py`) additionally needs `latexmk` with pdfTeX and
the Poppler tools `pdffonts`/`pdftotext`; an audit of supplied PDFs fails if font inspection is unavailable.
The stages of `reproduce_bit.py` need none of these tools.

## What kind of evidence each result is

- **Exact symbolic mathematics** (SymPy over the rationals or algebraic numbers): stability functions from
  Butcher tableaux, the Choi block form, trace-preservation identities, the signed first-defect law and its
  frequency factors, the RK4 cubic factorization and its discriminant, the extrapolation-corollary
  coefficients. Outputs: `results/bit_revision/claim_verification.{json,md}`,
  `results/bit_revision/extrapolation_corollary.{json,md}`, `results/symbolic_verification.txt`.
- **Certified interval/rational computation**: real-root isolation by Sturm sequences with exact rational
  brackets, exact rational evaluation of the constraint polynomials at midpoints and at executed binary64
  steps, exact algebraic root identity. Outputs: `results/closeout/root_sign_certificate_v2_7.json`,
  `results/v27/v27_verification.json`, `results/v34/v34_verification.json` (`candidate_regions_varpi_2`,
  `richardson`), the transverse-field endpoint records.
- **High-precision floating point** (mpmath, 60 digits): the fixed-horizon negative-mass curves of Fig. S3.
- **Binary64 floating point**: the adaptive-controller runs (`results/jcp_figures_v34/figure5_*.csv`,
  `results/v34/adaptive_benchmark_v34.csv`); every acceptance in the guarded runs is nevertheless decided by
  exact rational certification of the executed step.
- **Regression tests**: `tests/` (see `BIT_MANUSCRIPT_CONTRACT.md` for the statement-to-test map).
- **Visual inspection**: figure renders and the grayscale/deuteranopia proofs in `figures_bit/audit/`.

## Principal outputs

- `results/bit_revision/claim_verification.json`: 69 named checks with their exact values.
- `results/bit_revision/extrapolation_corollary.json`: per-method extrapolation defects (Online Resource Table S2).
- `results/bit_revision/figures/manifest.json`: source hashes, output hashes, overlap-audit results.
- `results/bit_revision/tables/*.tex` and `manifest.json`: Online Resource tables and their input hashes.
- `figures_bit/Fig1–Fig3, FigS1–FigS4` (`.eps`, `.pdf`, `.png`) and `figures_bit/audit/` proofs.
- `results/bit_revision/manifest.json`: written by `reproduce_bit.py all`.
- `BIT_MANUSCRIPT_CONTRACT.md`: every theorem, figure, and table of the article mapped to its generator,
  data, tests, and certificates.

Frozen v3.4 records consumed by the BIT generators: `results/jcp_figures_v34/*.csv`,
`results/v27/v27_verification.json`, `results/v34/v34_verification.json` (SHA-256 pinned in
`tests/test_bit_figure_contract.py` and in the manifests). Re-running the historical
`scripts/verify_v34_extensions.py` regenerates these records identically except for `wall_time_seconds`
fields.

## Historical v3.5 workflow (retained)

The v3.5 archive (2026-08-20) supported an earlier manuscript version prepared for a different journal.
Its reproduction script is kept operational:

```bash
bash reproduce.sh all          # or: certificates | tests | audits | legacy-heavy | figures | source-audit
```

`reproduce.sh tests` collects the historical 240-test suite (the current tree collects more tests because
of the BIT modules; the historical stage counts only its own 16 modules). `reproduce.sh figures`
regenerates the historical figure set in `figures_jcp/` from `scripts/generate_jcp_figures_v34.py`;
`reproduce.sh source-audit` audits the sources of the earlier manuscript version in `manuscript_contract/`,
which the historical regression tests also read. Known issue: the historical check
`independent_validation/remediation_v34.py` in `reproduce.sh audits` fails on a source-text heuristic about
the historical figure script (it expects three marker-only series; the script has four). Its scientific
checks pass, and `reproduce_bit.py` does not run it.

The files named `V3_*`, `TEST_*`, `VERIFY_*` and `*.log` at the top level of the archive are the unchanged
records of the earlier 3.x versions.

## Closed-interval root contract

The general-qubit, phase-covariant, and candidate locators establish exact algebraic membership in the
requested closed interval before restricting any Sturm bracket. Roots strictly outside the interval are
discarded. A genuine upper-endpoint root is retained and jointly certified. Strict interior roots are
re-isolated inside the requested domain. Unresolved membership, root identity, ordering, overlap, or sign
returns `UNCERTAIN` with no usable component set. The regression suite permanently covers two failure modes
found in earlier audits (an out-of-domain root cannot leave a phantom clipped bracket; an exact root at the
upper boundary is retained) and verifies fail-closed behaviour under unresolved root separation.

## Controller telemetry

Every component-locator invocation increments `locator_calls`, including empty and uncertain outcomes:

```text
locator_calls = successful_projections + empty_component_outcomes + uncertain_outcomes
```

The benchmark also reports `root_isolation_ops`, `certification_ops`, `binary_recertifications`,
`precision_escalations`, and `fallback_events`. Across the four standard boundary tolerances the frozen data
contain 15 locator calls, 15 empty-component outcomes, 45 root-isolation operations, 117 certification
operations, and 15 fallback events.

## Executed-step contract

A positive-width exact component guarantees an interior rational step specification. It does not guarantee
that a finite search will discover a representable IEEE-754 binary64 value accepted by the implementation.
Projection is accepted only after the chosen float is interpreted as its exact binary rational and the
corresponding candidate is recertified. Exhaustion or failure returns `UNCERTAIN` with diagnostic
`EXECUTED_BINARY_RECERTIFICATION_FAILED`.

## Transverse-field scope

The transverse-field study provides independent pointwise certificates only at
`Omega_x / Gamma = 0, 0.05, 0.10, 0.20`; phase covariance is broken at the three nonzero values, and
`Omega_x = 0` is the phase-covariant reference. Those four points each retain two finite-step
complete-positivity components. The archive makes no uniform claim over the intervening interval and no generic claim for
arbitrary noncommuting Lindbladians.

## Provenance

No external dataset is required. All symbolic expressions, rational brackets, exact interval samples,
benchmark records, and plotted data are included. Figures and tables are deterministic outputs of versioned
scripts from certified data (`SOURCE_DATE_EPOCH=1787184000`). `CITATION.cff` identifies the versioned release.
The reproduction manifest records 30 generated artifacts, 12 frozen inputs, the computation inputs and,
when available, the article and Online Resource sources (not the cover letter). The root
`CODE_ARCHIVE_ALLOWLIST.txt` lists the 408 delivered files.

## License

The software and associated documentation are distributed under the MIT License; see `LICENSE`.
