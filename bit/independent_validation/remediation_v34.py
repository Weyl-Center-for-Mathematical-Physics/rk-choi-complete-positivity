#!/usr/bin/env python3
"""Independent release gate for the v3.4 hostile-audit remediation."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import sympy as sp

from rk_choi_margin import liouvillian
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import certified_candidate_guard_step

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "independent_validation" / "remediation_v34_report.json"

checks: dict[str, object] = {}
h = sp.Symbol("h", nonnegative=True, real=True)
original = liouvillian.rk4_gksl_principal_minor_polynomials
try:
    upper = sp.Rational(1414213562, 10**9)

    def outside_polynomials(**_kwargs):
        return h, sp.Matrix([[1]]), {(0,): sp.Poly(h**2 - 2, h, domain=sp.QQ)}

    liouvillian.rk4_gksl_principal_minor_polynomials = outside_polynomials
    outside = liouvillian.locate_rk4_gksl_cptp_components_exact(upper=upper, root_digits=8)
    checks["outside_root_discarded_before_restriction"] = (
        outside["status"] is CPStatus.PASS
        and outside["roots"] == tuple()
        and outside["components"] == ((sp.Rational(0), upper, upper / 2),)
    )

    def boundary_polynomials(**_kwargs):
        return h, sp.Matrix([[1 - h]]), {(0,): sp.Poly(1 - h, h, domain=sp.QQ)}

    liouvillian.rk4_gksl_principal_minor_polynomials = boundary_polynomials
    boundary = liouvillian.locate_rk4_gksl_cptp_components_exact(upper=1, root_digits=8)
    root = boundary["roots"][0] if boundary["roots"] else None
    checks["exact_upper_boundary_retained"] = bool(
        boundary["status"] is CPStatus.PASS
        and root is not None
        and root.lower == root.upper == root.exact_root == 1
        and root.point_status is CPStatus.PASS
    )
finally:
    liouvillian.rk4_gksl_principal_minor_polynomials = original

decision = certified_candidate_guard_step(
    method="rk4",
    candidate_kind="equal-substeps",
    substeps=2,
    proposed_x=1.4,
    theta=0.0,
    kappa=0.0,
    varpi=2.0,
)
telemetry = decision.telemetry
checks["locator_outcome_partition"] = telemetry["locator_calls"] == (
    telemetry["successful_projections"]
    + telemetry["empty_component_outcomes"]
    + telemetry["uncertain_outcomes"]
)
checks["granular_telemetry_keys"] = all(
    key in telemetry
    for key in (
        "locator_calls",
        "successful_projections",
        "empty_component_outcomes",
        "uncertain_outcomes",
        "root_isolation_ops",
        "certification_ops",
        "binary_recertifications",
        "precision_escalations",
    )
)

v34 = json.loads((ROOT / "results" / "v34" / "v34_verification.json").read_text())
boundary_totals = v34["guard_rotating_fallback_action_totals"]
checks["boundary_telemetry_totals"] = {
    key: int(boundary_totals[key])
    for key in (
        "locator_calls",
        "successful_projections",
        "empty_component_outcomes",
        "uncertain_outcomes",
        "root_isolation_ops",
        "certification_ops",
        "binary_recertifications",
        "precision_escalations",
        "fallback_events",
    )
}
checks["boundary_expected_counts"] = (
    int(boundary_totals["locator_calls"]) == 15
    and int(boundary_totals["empty_component_outcomes"]) == 15
    and int(boundary_totals["root_isolation_ops"]) == 45
    and int(boundary_totals["certification_ops"]) == 117
    and int(boundary_totals["fallback_events"]) == 15
)

with (ROOT / "results" / "jcp_figures_v34" / "figure3_noncommuting_endpoints.csv").open(newline="") as handle:
    transverse_rows = list(csv.DictReader(handle))
transverse_values = [float(row["omega_x"]) for row in transverse_rows]
checks["transverse_pointwise_values"] = transverse_values
checks["exactly_four_transverse_certificates"] = transverse_values == [0.0, 0.05, 0.1, 0.2]

figure_script = (ROOT / "scripts" / "generate_jcp_figures_v34.py").read_text(errors="replace")
# The robustness panel must use marker-only scatter calls and must not construct
# a fill_between confidence band for the transverse data.
fig3_block = figure_script.split("def fig2_robustness_noncommuting():", 1)[1].split("def fig4_certified_guard_conditioning():", 1)[0]
checks["marker_only_transverse_plot"] = (
    fig3_block.count('linestyle="None"') == 3
    and "fill_between" not in fig3_block
    and "interpolation" in fig3_block
)

checks["pass"] = all(value is True for key, value in checks.items() if key not in {
    "boundary_telemetry_totals", "transverse_pointwise_values"
})
REPORT.write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks, indent=2))
raise SystemExit(0 if checks["pass"] else 1)
