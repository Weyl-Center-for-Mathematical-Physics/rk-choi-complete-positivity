# Claim-to-evidence map

Article: *Complete-positivity regions of Runge–Kutta discretizations of phase-covariant qubit dynamics*
(manuscript for BIT Numerical Mathematics, version of 22 September 2026). Every numbered statement, figure,
and table of the article and of Online Resource 1 maps below to the script that produces or verifies it,
the data it reads, the tests that guard it, and the certificate files that record it. Evidence classes:
**symbolic** = exact SymPy identity; **certified** = exact rational/algebraic root isolation (Sturm) or exact
rational evaluation; **high precision** = mpmath at stated precision; **floating point** = binary64 runs;
**visual** = inspected renders.

The independent checker `scripts/verify_bit_claims.py` (69 checks, output
`results/bit_revision/claim_verification.{json,md}`) rebuilds every object from first principles without
importing the package; check names below refer to its JSON keys. `scripts/verify_bit_extension.py` (output
`results/bit_revision/extrapolation_corollary.{json,md}`) verifies the extrapolation corollary.

## Article statements

| Identifier (article) | Evidence class | Producer / verifier | Source data | Tests |
|---|---|---|---|---|
| Eq. (1)–(3) generator, modes, spectrum | symbolic | `verify_bit_claims.py`: `mode_equation_population`, `mode_equation_coherence`, `generator_spectrum`; `scripts/verify_phase_covariant.py` | — | `tests/test_core.py` |
| Definition 1 (candidate maps) and Eq. (4)–(8): stability function, multipliers, direct map, extrapolated candidate | symbolic | `rk4_direct_map_equals_claimed_block_form`, `stage_domain_is_spectral_condition`; `src/rk_choi_margin` map constructors | Butcher tableaux in `src/rk_choi_margin` | `tests/test_core.py`, `tests/test_v28_candidates.py` |
| Theorem 1 (CPTP criterion), Eq. (9)–(14): Choi matrix, admissible set, criterion, one-way margin, exact margin | symbolic | `rk4_choi_matrix_equals_claimed_form`, `trace_preservation_identities`, `choi_block_determinant_identity`, `criterion_reduction_structure`, `exact_semigroup_margin`; `scripts/verify_theorems.py`, `scripts/verify_phase_covariant.py` | — | `tests/test_core.py`, `tests/test_v29_theory.py` |
| Eq. (15) buffers | symbolic | `dephasing_buffer_first_order`, `excitation_buffer_coefficient` | — | `tests/test_v26.py` |
| Theorem 2 (signed first-defect law), Eq. (16)–(19); Proposition (degenerate frequency, ESM) | symbolic | `signed_first_defect_law_generic_m2_to_m8`, `G4_G5_G6_G7_closed_forms`, `unit_circle_rule_q4`, `rk4_threshold_factorization`, `ssp3_threshold_factorization`; `src/rk_choi_margin/extrapolation.py::frequency_factor` | — | `tests/test_bit_scientific_contract.py` (closed forms, `G_m(0)`), `tests/test_v26.py` |
| Corollary 3 (extrapolated candidates), Eq. (21); Table 1 | symbolic | `scripts/verify_bit_extension.py`; `src/rk_choi_margin/extrapolation.py::richardson_extrapolate`, `predicted_extrapolation_defect`, `local_orientation_coefficient`; Table 1 values = `results/bit_revision/extrapolation_corollary.json` | tableaux | `tests/test_bit_scientific_contract.py` (8 methods × n = 2, 3, 4, parity rule, degenerate midpoint), `tests/test_bit_table_contract.py::test_article_methods_table_matches_corollary_record` |
| Eq. (22)–(23) population multiplier, α₄ | certified | `rk4_population_multiplier_positive`, `alpha4_population_ceiling`; `scripts/generate_root_certificate.py` → `results/closeout/root_sign_certificate_v2_7.json` | — | `tests/test_closeout.py` |
| Proposition 4 (margin factorization), Eq. (24)–(26): cubic 𝒞_q, discriminant, q±, ϖ±, critical factorizations | symbolic | `rk4_margin_cubic_factorization`, `cubic_discriminant`, `q_pm_roots`, `varpi_pm_closed_forms`, `critical_cubic_factorizations` | — | `tests/test_v27_certification.py`, `tests/test_bit_scientific_contract.py` |
| Theorem 5 (six regimes), crossover q_c / ϖ_c, isolated points, Eq. (27) window asymptotics | certified + symbolic | `cubic_root_counts_by_regime`, `empty_positive_step_regime_4_to_qplus`, `attached_branch_ceiling_crossover` (exact certificate: reduced discriminant in ℚ(α₄), resultant with the minimal polynomial Sturm-isolated, sign change), `detached_endpoints_below_alpha4` (y* inequality and sign samples), `isolated_points_q4_and_qplus`, `window_edge_asymptotics`, `rk4_imaginary_axis_boundary`; archive certificate `results/closeout/root_sign_certificate_v2_7.json`, `results/v27/v27_verification.json` (`varpi2_detached_interval`) | `results/jcp_figures_v34/figure1_exact_branches.csv` (plot only) | `tests/test_closeout.py`, `tests/test_topology_operational.py`, `tests/test_v27_certification.py`, `tests/test_bit_scientific_contract.py::test_rk4_regime_one_crossover_frequency` |
| Proposition 6 (representation dependence), Eq. (28) | symbolic | `hamiltonian_and_dissipative_parts_commute`, `generator_splits`, `unitary_channel_equals_exp_hamiltonian_part`, `lab_and_rotating_rk4_maps_differ`, `rotating_candidate_admissibility_is_omega_free`, `frame_candidates_differ_prop6ii` | — | `tests/test_v28_candidates.py` |
| Proposition 7 (composition and dilation), Eq. (29) | symbolic | `composition_calculus` | — | `tests/test_v28_candidates.py` |
| Proposition 8 (substep ladder), Example (ϖ = 2), Eq. (30), halving threshold ϖ₁/₂ | certified | `varpi2_direct_interval`, `varpi2_two_half_step_disjoint`, `substep_ladder_overlap_for_n_ge_2`, `halving_threshold_resultant`, `halving_threshold_value`; `results/v34/v34_verification.json::candidate_regions_varpi_2` | `results/jcp_figures_v34/figure2_candidate_regions.csv` | `tests/test_v28_candidates.py`, `tests/test_v34_remediation.py`, `tests/test_bit_scientific_contract.py::test_varpi2_full_and_two_half_step_positive_components_are_disjoint` |
| Theorem 9 (loss of CP under extrapolation), Eq. (31)–(33), P₁₀ | symbolic + certified | `richardson_stability_polynomial`, `richardson_first_defect`, `richardson_margin_factorization`, `richardson_first_positive_root_above_alpha4`, `richardson_remote_interval`, `constituents_cptp_on_common_interval`, `rk4_nonrotating_margin_root`, `richardson_orientation_bands`, `dp5_orientation_bands_reversed`; `results/v34/v34_verification.json::richardson` | `results/jcp_figures_v34/figure2_richardson_margin.csv` | `tests/test_v34_remediation.py`, `tests/test_v28_candidates.py` |
| Remark (convergence from outside; fixed horizon) | high precision (mpmath, 60 digits) | `scripts/generate_bit_figures.py::figS3`; prior note (cited) | `results/jcp_figures_v34/figureS2_operational.csv` | `tests/test_topology_operational.py` |
| Theorem 10 (persistence under buffers); buffered roots at ϖ = 2 | symbolic + certified | `table1_buffered_roots_varpi2`, `persistence_sign_samples`; `scripts/verify_v27_extensions.py` | `results/v27/v27_verification.json` | `tests/test_v27_certification.py` |
| Eq. (34)–(35), Lemma 11 (noncommuting, pointwise) | certified (exact characteristic polynomials of the 4×4 Choi matrix) | `noncommuting_eigenvalue_signs_omega_x_0p1`, `noncommuting_components_omega_x_0p1`, `transverse_table_endpoints`; `scripts/verify_v27_extensions.py` | `results/v27/v27_verification.json`, `results/jcp_figures_v34/figure3_noncommuting_endpoints.csv` | `tests/test_v27_noncommuting_adaptive.py` |
| Eq. (36) tri-state guard, Algorithm 1, the projection paragraph of §6.3 | certified + floating point | `src/rk_choi_margin` locator/guard; `independent_validation/executed_candidate_contract_v31.py` (82 checks), `endpoint_membership_v33.py`; `printed_decimal_vs_binary64_margin` | — | `tests/test_v29_executed_candidate.py`, `tests/test_v30_execution_contract.py`, `tests/test_v33_endpoint_membership.py` |
| Table 2 (adaptive tests) and §6.4 numbers (117 certification operations, 15 locator calls) | floating point (binary64 runs) with exact certification | `scripts/verify_v34_extensions.py` (regenerates the records; wall-clock fields differ between runs) | `results/jcp_figures_v34/figure5_tolerance_sweep.csv`, `results/v34/v34_verification.json` (`scientifically_useful_candidate_guard_rows`, `guard_rotating_fallback_action_totals`) | `tests/test_bit_table_contract.py::test_article_adaptive_table_matches_frozen_records`, `tests/test_v34_remediation.py` |
| Nine-candidate admissible sets on nonrotating one-way damping (ESM §9: eight formulas from the prior note, plus the Richardson candidate) | certified | `nine_candidate_admissible_sets`, `richardson_global_sign_chart`; certificate `nine_candidate_certificate.tex` in the Online Resource source archive | — | `tests/test_closeout.py`, `tests/test_bit_scientific_contract.py::test_richardson_global_chart_accounts_for_all_positive_events_and_tail` |

## Figures and tables

| Item | Generator | Source data (SHA-256 pinned in `results/bit_revision/figures/manifest.json`) | Tests |
|---|---|---|---|
| Fig. 1 (RK4 classification, inset) | `scripts/generate_bit_figures.py::fig1` | `figure1_exact_branches.csv`, `v27_verification.json`, exact critical values | `tests/test_bit_figure_contract.py` (sources, size, typography, overlap audit, series counts, determinism) |
| Fig. 2 (candidate maps at ϖ = 2; Richardson margin) | `::fig2` | `figure2_candidate_regions.csv`, `figure2_richardson_margin.csv` | same |
| Fig. 3 (buffered components; transverse certificates) | `::fig3` | `v27_verification.json`, `figure3_noncommuting_endpoints.csv` | same |
| Fig. S4 of the Online Resource (accuracy–work; accepted candidate types; moved from the article to the Online Resource) | `::fig4` (output `FigS4`) | `figure5_tolerance_sweep.csv`, `figure5_action_composition.csv`, `v34_verification.json` | same |
| Fig. S1 (orientation bands) | `::figS1` | exact `G_m` polynomials (symbolic) | same |
| Fig. S2 (conditioning; signed polynomial; exact vs binary64) | `::figS2` | `figure1_conditioning.csv`, `figure4_conditioning.csv`, `v27_verification.json` | same |
| Fig. S3 (witness fraction; fixed horizon) | `::figS3` | `figureS2_operational.csv`, mpmath recomputation | same |
| Table 1 (article) | hand-typeset from `extrapolation_corollary.json` | `results/bit_revision/extrapolation_corollary.json` | `tests/test_bit_table_contract.py` |
| Table 2 (article) | hand-typeset from frozen records | `figure5_tolerance_sweep.csv`, `v34_verification.json` | `tests/test_bit_table_contract.py` |
| Tables S2–S4 of the Online Resource (extrapolation per method; tolerance sweep; telemetry); Table S1 (regime samples) is typeset from the `cubic_root_counts_by_regime` record | `scripts/generate_bit_tables.py` → `results/bit_revision/tables/*.tex` | as above (`results/bit_revision/tables/manifest.json`) | `tests/test_bit_table_contract.py` |
| Grayscale / deuteranopia proofs | `scripts/generate_bit_figures.py::proofs` → `figures_bit/audit/` | rendered PNGs | visual (proofs in `figures_bit/audit/`) |

## Reproduction

`python reproduce_bit.py all` runs, in order, `certificates` (the two verifiers above), `tests` (every
module in isolated processes), `audits` (the retained v3.x audit scripts), `figures` (figures and tables),
and `manuscript-contract` (`scripts/audit_bit_source.py`, only when the article's LaTeX sources are
supplied with `--manuscript` or found in a folder `manuscript/` next to the archive), writing timestamped logs to `results/bit_revision/logs/` and the summary
`results/bit_revision/manifest.json`. The historical v3.5 chain (`bash reproduce.sh all`) is unchanged;
its `verify_v34_extensions.py` step regenerates the frozen v3.4 records identically except for
`wall_time_seconds` fields, and its historical `remediation_v34.py` check fails on a stale source-text
heuristic about the historical figure script (see `README.md`); `reproduce_bit.py` does not run that check.
