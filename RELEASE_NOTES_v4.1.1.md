# v4.1.1 — corrections to the article text, figures and documentation

Code and data for *Complete-positivity regions of Runge–Kutta discretizations of phase-covariant qubit
dynamics*, by G. Blake Pierpoint, Olivier Bernard, and Yichen Liu.

- **Article and Online Resource 1.** Corrected statements: the remark on the halving threshold now
  describes both reentrant-band roots of the resultant correctly; the substepping statements carry the
  condition that the substep population multiplier is nonnegative (with a forward Euler counterexample
  otherwise); the transverse-field evidence is stated as three certificates at nonzero field; the scope of
  the abstract, the acceptance test and several summary sentences now matches the theorems. Terms used
  informally before (Choi margin, boundary rays, buffer, small-step orientation) are defined. Algorithm 1 is
  restated as a single loop with the same logic.
- **Online Resource 1, Section 9.** The certificate appendix is now generated from the archived records
  (`scripts/generate_root_certificate.py`), with consistent method names, true minus signs and aligned
  polynomials. Tables S2–S4 use consistent number formats.
- **Figures.** Computer Modern ϖ in every figure, no guide line through the Fig. 1 legend, cleaner inset
  ticks, both detached components of Fig. 2a drawn in the detached style, a visible square-root reference
  in Fig. S2a, and consistent labels.
- **Packaging and documentation.** The cover letter is no longer bound into the reproduction manifest; the
  archive README, VERSION.md and CITATION.cff are updated and self-contained.

Certified values, data records and table values are unchanged since v4.0.0.

The attached `BIT_Code_and_Data.zip` holds the same 408 files as `bit/` at this tag, byte for byte.
Run from either, `python reproduce_bit.py all` passes four stages and skips the audit of the article's
LaTeX sources, with 384 tests passed and 4 skipped. Versions 1.0.0, 4.0.0, 4.1.0 and 4.1.1
are archived on Zenodo under the concept DOI
[10.5281/zenodo.21522057](https://doi.org/10.5281/zenodo.21522057).

Release URL: https://github.com/Weyl-Center-for-Mathematical-Physics/rk-choi-complete-positivity/releases/tag/v4.1.1

License: MIT.
