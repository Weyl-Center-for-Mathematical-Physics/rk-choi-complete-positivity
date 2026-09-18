# v3.4 hostile-audit remediation changelog

## Closed-interval exact root semantics

- Replaced the general-qubit locator's early bracket clipping with exact
  algebraic membership testing before interval restriction.
- Roots strictly outside the requested closed interval are discarded.
- Exact upper-boundary roots are retained and jointly certified.
- Unresolved separation, membership, ordering, overlap, or sign fails closed as
  `UNCERTAIN` and yields no usable component decomposition.
- Added permanent regression tests for the phantom-clipped-root defect, the
  exact upper-endpoint defect, and unresolved-root fail-closed behavior.

## Controller telemetry

- `locator_calls` now increments on every component-locator invocation.
- Added `successful_projections`, `empty_component_outcomes`,
  `uncertain_outcomes`, `root_isolation_ops`, `certification_ops`,
  `binary_recertifications`, `precision_escalations`, and `fallback_events`.
- Retained corrected compatibility aliases for legacy downstream readers.
- Regenerated the standard tolerance sweep and all manuscript macros, tables,
  CSV records, and performance figures.

## Scope and execution contract

- Restricted the transverse-field evidence to four independent pointwise
  certificates at `Omega_x/Gamma = 0, 0.05, 0.10, 0.20`.
- Removed interpolating curves and filled interval bands from the robustness
  figure.
- Distinguished existence of an interior rational step from finite discovery
  and exact recertification of an IEEE-754 binary64 step.
- Made `EXECUTED_BINARY_RECERTIFICATION_FAILED` an explicit fail-closed outcome.

## Manuscript and release engineering

- Standardized the classified set as the laboratory-frame one-way,
  zero-dephasing ray; reserved rotating-frame terminology for the analytically
  transformed scheme.
- Reframed the candidate-aware controller as a validated safety proof of
  concept, not an optimized runtime engine.
- Stated that the plotted work proxy excludes algebraic root-isolation and
  certification costs, which are reported separately.
- Installed the 186-word abstract, five compliant Highlights, Sturm/root-
  isolation and software references, AI-use disclosure, and deterministic
  graphical-abstract provenance statement.
- Synchronized the complete v3.4 manuscript source into `manuscript_contract/`.
