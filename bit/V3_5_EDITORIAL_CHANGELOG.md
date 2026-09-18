# JCP v3.5 editorial change report

**Date:** 2026-08-18; Figure 1 layout correction 2026-08-20  
**Baseline:** immutable v3.4 submission package dated 2026-08-15  
**Disposition:** **TECHNICAL PASS — editorial revision complete; journal upload remains a human decision**

## Scope and approach

Version 3.5 is a publication-facing editorial patch over the scientifically
validated v3.4 baseline. It applies every HIGH and MODERATE recommendation in
`JCP_v3_4_Hostile_Human_Prose_Audit.md`, while preserving equations, theorem
scope, numerical values, definitions, formal declarations, and required
limitations. The source package and the manuscript contract embedded in the
code archive were kept byte-identical for every manifest-listed file.

## Applied prose changes

| Severity | Location | Resolution |
|---|---|---|
| MODERATE | Abstract | Replaced the vague candidate-failure transition and closing validation inventory with the concrete disjoint-component result and the finite-step consequence. The abstract is now 180 words. |
| HIGH | Introduction | Compressed the contribution matrix and seven-sentence roadmap into three short paragraphs while preserving the full-family/ray distinction. |
| HIGH | Executed-map semantics | Replaced the internal schema and gate language with the certified mathematical object, executed binary value, and component responsibilities. |
| HIGH | Root/component isolation | Reduced the remediation-style checklist to endpoint comparison, interior re-isolation, equality testing, and one fail-closed rule. |
| HIGH | Benchmark discussion | Kept the work-proxy definition and decisive results; moved the counter taxonomy back to the telemetry table. |
| MODERATE | Verification paragraph | Retained exact-data provenance and the independent Choi reconstruction while removing a three-list software inventory. |
| HIGH | Conclusion | Removed the repeated telemetry catalogue and the generic “broader lesson / principled bridge” ending; the paper now closes on candidate-map certification. |
| MODERATE | Supplement opening | Replaced the section-by-section itinerary with purpose, grouped navigation, and reproduction instructions. |
| MODERATE | Supplement certification | Kept operational detail in the supplement but shortened duplicated root/failure prose. |
| MODERATE | Figure 5 caption | Limited the caption to panel descriptions and the decisive safety result; a short clarification now identifies the intentional equal-work pair. |
| HIGH | Cover letter | Compressed controller evidence and archive/disclosure prose into editor-facing result statements. |
| HIGH | Graphical abstract | Replaced process-assurance language with three parameterized findings and the fallback rule. |

## Additional quality fixes

- Added `lmodern` to the line-numbered review source so MiKTeX/pdflatex uses
  scalable fonts with `microtype`; this resolves the prior font-expansion
  failure without changing the submission manuscript.
- Selected Latin Modern Mono for e-mail addresses and execution diagnostics,
  replacing two embedded bitmap EC fallbacks while retaining the Elsevier
  Times body. Every final PDF font is now embedded and vector-based.
- Declared and pinned `scipy==1.16.1`, which three bundled independent-audit
  scripts import but v3.4 did not list as a dependency.
- Updated version/date metadata to 3.5.0 / 2026-08-18 and introduced the v3.5
  source-audit script and report names.
- Regenerated tables, figures, graphical abstract, audit reports, and test logs
  from the pinned environment. Figure 4 reflects only machine-roundoff-level
  recomputation; the graphical abstract carries the approved new copy.
- Reauthored Figures 1--5 at a 7.45-inch native width, close to Elsevier's
  190-mm double-column target. Main figure text is at least 7.6 pt before the
  small final-placement reduction; the resulting normal text is approximately
  7 pt in both manuscript layouts.
- Standardized the figures on embedded Arial text with STIX sans-serif math,
  an Okabe--Ito colorblind-safe palette, stronger line/marker weights, and
  color-independent marker, line-style, and hatch cues.
- Replaced transparency in Figures 1--5 with opaque light tints and grids, and
  emitted EPS text as portable embedded Type 3 glyphs. Ghostscript/MiKTeX
  round-trip renders retain every label and the intended visual hierarchy.
- Relocated the Figure 1a reentrant-band inset into the empty central region,
  moved its zoom label inside the opaque inset, and reduced the callout to one
  short local leader. A geometry regression rejects inset collisions with the
  component regions, legend, and annotations; tick-label overlap with the
  highlighted source band; and overlong or annotation-crossing leaders.
- Moved the Figure 3 legends into data-free regions and added a render-level
  regression that fails if the panel-a legend covers a component bar.
- Labeled the two intentional equal-work error-only observations in Figure 5
  ($10^{-2}$ and $3\times10^{-3}$ at work 60), and separated the panel-b
  telemetry labels so adjacent bars do not collide.
- Recentered Table 1 at its natural width with balanced numeric-column spacing,
  eliminating the stretched final column without changing any values.
- Focused Figure 5a on the nonzero benchmark curves while identifying the
  commuting Strang result as a roundoff-floor reference omitted from the log
  scale. Figure 5b now suppresses the all-zero exact-fallback stack. Both facts
  are stated in the caption instead of occupying the plotted data regions.
- Repositioned the lower Figure 5a tolerance label and added a rendered-geometry
  regression that rejects annotation intersections with the panel-a legend,
  the error-only curve, or panel-b bars.
- Added a review-only centered figure-width override (80% of paper width),
  preserving the standard two-column placement in the clean manuscript while
  making every figure substantially easier to read in the line-numbered copy.
- Replaced stale v3.4 branding in upload-facing PDF and PNG metadata with
  publication-facing, version-neutral generator metadata dated 2026-08-18;
  the v3.4-named source script remains unchanged in name to preserve the
  validated scientific provenance.
- Added the validated PDF graphical abstract to the upload carrier alongside
  the high-resolution PNG and normalized the official title punctuation across
  the title page, portal draft, package guide, and reproducibility README.

## Scientific invariance check

The regenerated `results/v34/v34_verification.json` has no nonnumeric changes
after excluding wall-clock fields. Eleven floating-point leaves differ from the
v3.4 artifact, with maximum absolute difference `4.43e-16`; every difference is
within `atol=1e-14, rtol=1e-12`. The single largest relative report is caused by
`0` versus `-3.23e-30`, both numerical zero at the stated scale. Tests, exact
certificates, scope guards, and independent audits all pass.

## Validation completed

- 234/234 pytest tests passed in the package-prescribed isolated invocations.
- Exact certificate, theorem, phase-covariant, method, endpoint, executed-map,
  v3.4 remediation, and legacy-heavy noncommuting stages passed.
- Independent audit: 20/20 checks; 15,000 randomized scalar-versus-Choi
  comparisons; zero mismatches.
- Source audit: 20 TeX files, 99 unique labels, 42 reference occurrences,
  64/64 bibliography entries cited, no missing assets, and five Highlights at
  or below 85 characters.
- Final PDFs: 16-page article, 36-page line-numbered review copy, 14-page
  supplement, one-page cover letter, and one-page graphical abstract (68 core
  pages), plus five one-page separate figure PDFs. Figure 1 and its clean page
  6 and review page 12 were freshly rasterized after the inset correction; the
  prior Figure 5 correction remains verified on its standalone PDF, clean page
  12, and review page 25. No overlap, clipping, or page-edge collision remains
  in the affected outputs.
- LaTeX logs contain no fatal errors, undefined references/citations, duplicate
  destinations, missing characters, or overfull boxes. All PDF fonts are
  embedded.

## File-level ledger

### Author-facing source modified

- `jcp_article.tex`
- `jcp_preamble.tex`
- `sections/01_introduction.tex`
- `sections/04_rk4_topology.tex`
- `sections/05_robustness.tex`
- `sections/06_candidate_maps.tex`
- `sections/06_step_control.tex`
- `sections/07_comparison.tex`
- `sections/08_conclusion.tex`
- `Supplemental_Material_JCP.tex`
- `Cover_Letter_JCP.tex`
- `Figure_Captions.txt`
- `Manuscript_JCP_review.tex`
- `README.txt`, `BUILD_INSTRUCTIONS.txt`, and `Title_Page.txt`

The same files were modified under `manuscript_contract/` in the code archive.

### Reproducibility and metadata modified

- `CITATION.cff`, `pyproject.toml`, `README.md`, `README.txt`, `VERSION.md`,
  `ENVIRONMENT.txt`, and `reproduce.sh`
- `scripts/generate_jcp_figures_v34.py` and
  `scripts/verify_v34_extensions.py`
- `tests/test_v35_editorial_layout.py` now guards publication-facing metadata,
  Figure 3a legend placement, finished-size typography, Figure 5 series
  distinguishability and equal-work labels, suppression of scale-distorting
  roundoff references and all-zero action categories, rendered annotation
  clearance from legends, curves, and bars, non-color component encodings,
  portable EPS glyph embedding, transparency-free main figures, and the
  explicit Figure 1 inset zoom cue, including source-band/tick clearance,
  component and legend clearance, and a short noncrossing leader.
- `scripts_source_audit_v35.py` replaced the v3.4-named source-audit entrypoint;
  `SOURCE_AUDIT_v3_5.json` replaced the v3.4-named report.
- Generated figures under `figures_jcp/` and `manuscript_contract/figures/`,
  generated tables/macros, verification JSON/CSV/Markdown records, independent
  audit reports, and stage/test logs were refreshed by the validated workflow.

### Preserved unchanged

The v3.4 ZIP archives, v3.4 download checksum record, v3.4 remediation report,
formal CRediT/interest DOCX files, Highlights text, equations, theorem/proof
content outside the listed prose locations, and all human-signoff requirements
remain unchanged.

## Remaining human gate

No upload, repository deposition, or journal submission was performed. Authors
must still confirm metadata, order and contributions, affiliations/ORCIDs,
funding and conflicts, AI-use wording, graphical-abstract provenance,
prior/concurrent submission status, repository identifier if any, reviewers and
conflicts, and the Editorial Manager-generated proof.
