#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from rk_choi_margin.candidates import (
    CandidateSpec,
    equal_substep_candidate,
    locate_candidate_components,
)
from rk_choi_margin.certification import CPStatus, locate_cptp_components
from rk_choi_margin.methods import RKMethod

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "independent_validation" / "external_review_regressions_v33_report.json"

checks: list[dict[str, object]] = []

def record(name: str, passed: bool, detail: object) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")

# Reviewer's exact counterexample: R(z)=1+z+2z^2+z^3.
method = RKMethod(
    key="external-review-counterexample",
    name="external-review-counterexample",
    A=((0,0,0),(1,0,0),(0,1,0)),
    b=(-1,1,1),
    order=1, family="explicit",
)
loc = locate_cptp_components(
    method, theta=0, kappa=0, varpi=1, upper=2, root_digits=30
)
record("bad isolated root not admitted", not any(r.exact_root == 1 and r.point_status is CPStatus.PASS for r in loc.isolated_points),
       {"status": loc.status.value, "isolated": [str(r.exact_root) for r in loc.isolated_points]})
root = next(r for r in loc.roots if r.lower <= 1 <= r.upper)
record("bad root exact point status is FAIL", root.point_status is CPStatus.FAIL,
       {"point_status": root.point_status.value, "sources": root.source})

cloc = locate_candidate_components(
    "direct", method=method, theta=0, kappa=0, varpi=1, upper=2, root_digits=30
)
record("candidate locator also rejects false isolated point", not any(r.exact_root == 1 and r.point_status is CPStatus.PASS for r in cloc.isolated_points),
       {"status": cloc.status.value, "isolated": [str(r.exact_root) for r in cloc.isolated_points]})

# Positive control with an exact isolated CPTP point at x=1.
genuine = RKMethod(
    key="genuine-isolated-control", name="genuine isolated control",
    A=((0,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,1,0)),
    b=(sp.Rational(59,36), sp.Rational(37,36), -sp.Rational(8,9), -sp.Rational(7,9)),
    order=1, family="explicit",
)
gloc = locate_candidate_components(
    "direct", method=genuine, theta=0, kappa=0, varpi=0,
    upper=2, root_digits=30
)
record("genuine isolated CPTP equality retained",
       any(root.exact_root == 1 and root.point_status is CPStatus.PASS
           for root in gloc.isolated_points),
       {"status": gloc.status.value,
        "isolated": [str(root.exact_root) for root in gloc.isolated_points]})

# Exact-integer public boundary.
for bad in (True, 2.0, 2.9, 0, -1):
    try:
        equal_substep_candidate(method="rk4", total_x=1.0, substeps=bad,
                                theta=0, kappa=0, varpi=0)
    except (TypeError, ValueError):
        record(f"invalid substeps {bad!r} rejected early", True, type(bad).__name__)
    else:
        record(f"invalid substeps {bad!r} rejected early", False, "accepted")

# A valid integer remains accepted.
valid = equal_substep_candidate(method="rk4", total_x=1.0, substeps=2,
                                theta=0, kappa=0, varpi=0)
record("valid exact integer substep accepted", valid.a is not None, {"a": str(valid.a)})

payload = {
    "version": "3.3.0",
    "checks": checks,
    "passed": sum(1 for row in checks if row["passed"]),
    "failed": sum(1 for row in checks if not row["passed"]),
}
OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
