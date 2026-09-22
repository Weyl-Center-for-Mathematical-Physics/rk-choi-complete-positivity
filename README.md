# Runge–Kutta complete-positivity: code and data

## Current release: v4.1.2

[Release v4.1.2](https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.1.2) contains the code and data for *Complete-positivity regions of Runge–Kutta
discretizations of phase-covariant qubit dynamics* by G. Blake Pierpoint, Olivier Bernard and Yichen Liu.
The code is in [`bit/`](bit/); the release asset `BIT_Code_and_Data.zip` holds the same 408 files. The
folders `bit/manuscript_contract/` (sources of an earlier manuscript version) and `bit/manuscript_source/`
(the flattened article source) are included because the regression tests read them.

```bash
cd bit
python -m venv .venv
.venv/bin/python -m pip install -e ".[test]"   # on Windows: .venv\Scripts\python.exe
.venv/bin/python reproduce_bit.py all
```

Run from `bit/` or from the extracted release asset, `reproduce_bit.py all` passes four stages and skips the
audit of the article's LaTeX sources, with 384 tests passed and 4 skipped (the skipped tests compare against
those sources). See [`bit/README.md`](bit/README.md) and [`bit/VERSION.md`](bit/VERSION.md).

License: MIT. Versions 1.0.0, 4.0.0, 4.1.0, 4.1.1 and 4.1.2 are archived on Zenodo under the concept DOI
[10.5281/zenodo.21522057](https://doi.org/10.5281/zenodo.21522057), which resolves to the latest version (v4.1.2: [10.5281/zenodo.22906473](https://doi.org/10.5281/zenodo.22906473)); the
documentation-only v4.0.1 was not archived. Version 4.1.2 revises the wording of the article and its supplement; version 4.1.1
corrected statements in both, refined the figures and cleaned up the documentation. Certified values, data
records and tables are unchanged since v4.0.0.

## Version 1.0.0: zero-precession archive

The [v1.0.0 tag and release](https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v1.0.0) and its
[archived DOI](https://doi.org/10.5281/zenodo.21522058) are unchanged. The root `src/`, `data/`, `figures/`,
`supplement/`, `AUTHORS.md` and the historical reproduction workflow belong to that release, and the
documentation below describes it. The current code in `bit/` is reproduced separately.

### Version 1.0.0 documentation

This repository contains the exact symbolic calculations, certified boundary data, plotting data, and figures accompanying our zero-precession research note (Weyl Center for Mathematical Physics, Research Note 001).

**Authors:** G. Blake Pierpoint, Olivier Bernard, and Yichen Liu  
**Archive version:** 1.0.0  
**Release date:** 23 July 2026

Author affiliations and ORCIDs are listed in `AUTHORS.md`.

## Scope

The repository evaluates eight Runge-Kutta formulas for the pure-amplitude-damping boundary case. For each method, the reproduction script:

1. derives the stability function from exact rational Butcher coefficients;
2. computes the first nonzero defect between the stability function and the exponential;
3. constructs the population factor \(a(x)=R(-x)\) and complete-positivity margin
   \[
   M_R(x)=R(-x)-R(-x/2)^2;
   \]
4. derives primitive boundary polynomials for \(a=0\), \(a=1\), and \(M_R=0\);
5. certifies every listed positive-root bracket with exact Sturm counts;
6. evaluates the sign chart at exact rational sample points;
7. derives the CPTP-admissible components; and
8. regenerates the article figures and fixed-horizon defect data.

No external datasets, web access, or proprietary software are required.

## Repository contents

- `src/reproduce_pla_results.py` performs the exact verification and regenerates all tracked outputs.
- `data/certified_boundaries.json` contains tableau coefficients, exact symbolic expressions, primitive boundary polynomials and scalar factors, root certificates, multiplicities, rational sign samples, and interval truth values.
- `data/method_audit.csv` is the compact eight-method summary table.
- `data/fixed_horizon_defect.csv` contains values computed at 80-digit working precision and serialized to 30 significant digits for the fixed-horizon structural-defect figure.
- `supplement/root_sign_certificate.md` presents the root-isolation and sign-chart results in human-readable form.
- `figures/` contains publication-ready PDF figures and 600 dpi PNG copies.
- `CITATION.cff` and `.zenodo.json` provide citation and archival metadata.
- `.github/workflows/reproduce.yml` checks the exact calculations and tracked outputs on GitHub.
- `LICENSE` states the reuse terms.

## Reproduce the archive

A clean virtual environment is recommended:

```bash
python -m venv .venv
```

Activate it using the command appropriate to your platform, then install the pinned dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Regenerate every data file and figure:

```bash
python src/reproduce_pla_results.py
```

A successful run ends with:

```text
Reproduced 8 tableau-derived method audits, 14 positive-root certificates, and 3 figures.
Exact symbolic, Sturm, sign-chart, and fixed-horizon checks passed.
```

To verify the mathematics and confirm that the committed outputs are current without modifying the repository, run:

```bash
python src/reproduce_pla_results.py --check
```

The check is independent of Python's `assert` mechanism and therefore remains active when Python is invoked with optimization flags.

## Tested environment

The released outputs were generated with:

- CPython 3.13.5
- SymPy 1.14.0
- NumPy 2.3.5
- mpmath 1.3.0
- Matplotlib 3.10.8
- Pillow 12.3.0

The direct and rendering dependency versions are pinned in `requirements.txt`. PDF metadata use a fixed release timestamp, and JSON output uses fixed line endings, so that all tracked outputs are byte-reproducible in the tested environment.

## Interpretation of the certificate

The exact proof obligations are the boundary-polynomial identities, Sturm root counts, and rational interval signs. Decimal endpoint values are provided for readability; the corresponding decimal brackets are interpreted as exact rational intervals during verification.

The Markdown supplement summarizes the verified results but is not generated automatically. The JSON file is the complete machine-readable record produced by the script.

## License

Unless otherwise noted, the code, data, documentation, and figures in this repository are released under the MIT License; see `LICENSE`.

## Citation and persistent deposit

The public repository is:

`https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity`

Version 1.0.0 is archived at https://doi.org/10.5281/zenodo.21522058; later versions are listed at the top of this page.
