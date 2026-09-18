#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from rk_choi_margin.certification import (
    CPStatus,
    certify_cptp_step,
    isolate_constraint_roots_exact,
    locate_cptp_components,
)
from rk_choi_margin.candidates import locate_candidate_components
from rk_choi_margin.methods import RKMethod

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "independent_validation" / "endpoint_membership_v33_report.json"
Q = sp.Rational
X = sp.Symbol("x", real=True)
EPS = Q(1, 10**35)

checks: list[dict[str, object]] = []

def record(name: str, passed: bool, detail: object) -> None:
    row = {"name": name, "passed": bool(passed), "detail": detail}
    checks.append(row)
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def method_at(root: sp.Rational) -> RKMethod:
    c2 = 2 + 1 / root
    c3 = 2 / root
    return RKMethod(
        key=f"endpoint-{root}",
        name="Endpoint-domain adversary",
        A=((Q(0), Q(0), Q(0)), (Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0))),
        b=(1 - c2, c2 - c3, c3),
        order=1,
        family="explicit adversarial",
    )


def locate(root: sp.Rational):
    method = method_at(root)
    params = dict(theta=Q(0), kappa=Q(3, 4), varpi=Q(0), upper=Q(1), root_digits=30)
    return (
        certify_cptp_step(method, x=Q(1), theta=0, kappa=Q(3, 4), varpi=0),
        locate_cptp_components(method, **params),
        locate_candidate_components("direct", method=method, substeps=2, **params),
    )

for label, root, expected_fixed in (
    ("inside", Q(1) - EPS, CPStatus.PASS),
    ("equal", Q(1), CPStatus.PASS),
    ("outside", Q(1) + EPS, CPStatus.FAIL),
):
    fixed, generic, candidate = locate(root)
    record(f"upper {label}: fixed certificate", fixed.status is expected_fixed, fixed.status.value)
    for kind, location in (("generic", generic), ("candidate", candidate)):
        isolated_upper = any(
            point.lower == Q(1) == point.upper and point.point_status is CPStatus.PASS
            for point in location.isolated_points
        )
        terminal = any(
            component.lower_value == root and component.upper_value == Q(1)
            for component in location.components
        )
        if label == "inside":
            passed = terminal and not isolated_upper
        elif label == "equal":
            passed = isolated_upper
        else:
            passed = not isolated_upper and all(r.exact_root != root for r in location.roots)
        record(
            f"upper {label}: {kind} locator",
            passed,
            {
                "status": location.status.value,
                "components": [[str(c.lower_value), str(c.upper_value)] for c in location.components],
                "isolated": [str(p.exact_root) for p in location.isolated_points],
            },
        )

# Symmetric lower-domain membership.
for label, poly in (
    ("below", sp.Poly(X + EPS, X, domain=sp.QQ)),
    ("equal", sp.Poly(X, X, domain=sp.QQ)),
    ("inside", sp.Poly(X - EPS, X, domain=sp.QQ)),
):
    roots, complete = isolate_constraint_roots_exact({"p": poly}, Q(1), digits=30)
    if label in {"below", "equal"}:
        passed = complete and roots == []
    else:
        passed = complete and len(roots) == 1 and roots[0].exact_root == EPS and roots[0].lower == EPS == roots[0].upper
    record(f"lower {label}: exact domain membership", passed, [str(r.exact_root) for r in roots])

# Near-coincident roots must be separated in exact order, independent of input order.
polys = {
    "right": sp.Poly(X - (Q(1) + EPS), X, domain=sp.QQ),
    "left": sp.Poly(X - (Q(1) - EPS), X, domain=sp.QQ),
}
roots, complete = isolate_constraint_roots_exact(polys, Q(2), digits=30, max_digits=100)
record(
    "near-coincident roots separated and ordered",
    complete
    and [r.exact_root for r in roots] == [Q(1) - EPS, Q(1) + EPS]
    and roots[0].upper < roots[1].lower,
    [[str(r.lower), str(r.upper), str(r.exact_root)] for r in roots],
)

payload = {
    "version": "3.3.0",
    "passed": sum(row["passed"] for row in checks),
    "failed": sum(not row["passed"] for row in checks),
    "checks": checks,
}
OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(json.dumps(payload, indent=2))
