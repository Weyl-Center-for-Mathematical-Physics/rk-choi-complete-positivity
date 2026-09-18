from __future__ import annotations

import sympy as sp
import pytest

from rk_choi_margin.candidates import (
    equal_substep_candidate,
    locate_candidate_components,
)
from rk_choi_margin.certification import (
    CPStatus,
    isolate_constraint_roots_exact,
    locate_cptp_components,
)
from rk_choi_margin.methods import RKMethod

Q = sp.Rational
X = sp.Symbol("x", real=True)


def _false_isolated_method() -> RKMethod:
    # Stability function R(z)=1+z+2z^2+z^3.  At theta=kappa=0,
    # varpi=1 and x=1, 1-a has an even root but F=-45/64.
    return RKMethod(
        key="false-isolated",
        name="false isolated-point adversary",
        A=((Q(0), Q(0), Q(0)), (Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0))),
        b=(-Q(1), Q(1), Q(1)),
        order=1,
        family="explicit adversarial",
    )


def _genuine_isolated_method() -> RKMethod:
    # R(z)=1+z-23 z^2/36-5 z^3/3-7 z^4/9 and
    # M(x)=-x^2(x-1)^2 Q_4(x)/20736, with Q_4(x)>0 at x=1.
    return RKMethod(
        key="genuine-isolated",
        name="genuine isolated-point control",
        A=(
            (Q(0), Q(0), Q(0), Q(0)),
            (Q(1), Q(0), Q(0), Q(0)),
            (Q(0), Q(1), Q(0), Q(0)),
            (Q(0), Q(0), Q(1), Q(0)),
        ),
        b=(Q(59, 36), Q(37, 36), -Q(8, 9), -Q(7, 9)),
        order=1,
        family="explicit isolated control",
    )


def test_false_even_root_is_not_reported_as_cptp_isolated_point() -> None:
    method = _false_isolated_method()
    location = locate_cptp_components(
        method, theta=0, kappa=0, varpi=1, upper=2, root_digits=20
    )
    assert location.status is CPStatus.PASS
    assert not any(root.exact_root == 1 for root in location.isolated_points)
    root = next(root for root in location.roots if root.lower <= 1 <= root.upper)
    assert root.point_status is CPStatus.FAIL
    assert "one_minus_a" in root.source


def test_candidate_locator_uses_same_joint_isolated_point_semantics() -> None:
    method = _false_isolated_method()
    location = locate_candidate_components(
        "direct", method=method, theta=0, kappa=0, varpi=1, upper=2,
        root_digits=20,
    )
    assert location.status is CPStatus.PASS
    assert not any(root.exact_root == 1 for root in location.isolated_points)
    root = next(root for root in location.roots if root.lower <= 1 <= root.upper)
    assert root.point_status is CPStatus.FAIL


def test_genuine_isolated_cptp_point_is_retained() -> None:
    method = _genuine_isolated_method()
    generic = locate_cptp_components(
        method, theta=0, kappa=0, varpi=0, upper=2, root_digits=20
    )
    candidate = locate_candidate_components(
        "direct", method=method, theta=0, kappa=0, varpi=0, upper=2,
        root_digits=20,
    )
    assert generic.status is CPStatus.PASS
    assert candidate.status is CPStatus.PASS
    assert any(root.exact_root == 1 and root.point_status is CPStatus.PASS
               for root in generic.isolated_points)
    assert any(root.exact_root == 1 and root.point_status is CPStatus.PASS
               for root in candidate.isolated_points)


@pytest.mark.parametrize("value", [True, 2.0, 2.9, 0, -1])
def test_equal_substep_public_api_rejects_non_exact_positive_integers(value) -> None:
    with pytest.raises((TypeError, ValueError), match="substeps"):
        equal_substep_candidate(
            "rk4", total_x=1, substeps=value, theta=0, kappa=0, varpi=0
        )


def test_adaptively_separates_nearby_algebraic_roots() -> None:
    eps = Q(1, 10**35)
    polynomials = {
        "left": sp.Poly(X**2 - 2, X, domain=sp.QQ),
        "right": sp.Poly((X - eps) ** 2 - 2, X, domain=sp.QQ),
    }
    roots, separated = isolate_constraint_roots_exact(
        polynomials, Q(2), digits=8, max_digits=100
    )
    positive = [root for root in roots if root.lower > 1]
    assert separated
    assert len(positive) == 2
    assert positive[0].upper < positive[1].lower
