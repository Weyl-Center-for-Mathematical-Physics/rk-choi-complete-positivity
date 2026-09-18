from __future__ import annotations

import sympy as sp
import pytest

from rk_choi_margin.certification import (
    CPStatus,
    certify_cptp_step,
    isolate_constraint_roots_exact,
    locate_cptp_components,
)
from rk_choi_margin.candidates import locate_candidate_components
from rk_choi_margin.methods import RKMethod

Q = sp.Rational
X = sp.Symbol("x", real=True)
EPS = Q(1, 10**35)
PARAMS = dict(theta=Q(0), kappa=Q(3, 4), varpi=Q(0), upper=Q(1), root_digits=30)


def _method(root: sp.Rational) -> RKMethod:
    """Consistent explicit method with a controllable root of 1-R(-x).

    Its stability polynomial satisfies
        1-R(-x) = x(2x-1)(x-root)/root.
    """

    c2 = 2 + 1 / root
    c3 = 2 / root
    return RKMethod(
        key=f"endpoint-{root}",
        name="Endpoint membership adversary",
        A=((Q(0), Q(0), Q(0)), (Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0))),
        b=(1 - c2, c2 - c3, c3),
        order=1,
        family="explicit adversarial",
    )


def _locations(root: sp.Rational):
    method = _method(root)
    fixed = certify_cptp_step(method, x=Q(1), theta=0, kappa=Q(3, 4), varpi=0)
    generic = locate_cptp_components(method, **PARAMS)
    candidate = locate_candidate_components("direct", method=method, substeps=2, **PARAMS)
    return fixed, generic, candidate


def _has_upper_isolated(location) -> bool:
    return any(
        item.lower == Q(1) == item.upper and item.point_status is CPStatus.PASS
        for item in location.isolated_points
    )


def _has_terminal_component(location, lower: sp.Rational) -> bool:
    return any(item.lower_value == lower and item.upper_value == Q(1) for item in location.components)


@pytest.mark.parametrize("which", ["generic", "candidate"])
def test_root_just_outside_upper_is_excluded(which: str):
    fixed, generic, candidate = _locations(Q(1) + EPS)
    location = generic if which == "generic" else candidate
    assert fixed.status is CPStatus.FAIL
    assert not _has_upper_isolated(location)
    assert all(root.exact_root != Q(1) + EPS for root in location.roots)


@pytest.mark.parametrize("which", ["generic", "candidate"])
def test_root_exactly_at_upper_is_retained_as_isolated_equality(which: str):
    fixed, generic, candidate = _locations(Q(1))
    location = generic if which == "generic" else candidate
    assert fixed.status is CPStatus.PASS
    assert _has_upper_isolated(location)


@pytest.mark.parametrize("which", ["generic", "candidate"])
def test_root_just_inside_upper_exposes_terminal_component(which: str):
    root = Q(1) - EPS
    fixed, generic, candidate = _locations(root)
    location = generic if which == "generic" else candidate
    assert fixed.status is CPStatus.PASS
    assert _has_terminal_component(location, root)


def test_lower_root_just_below_domain_is_excluded_before_clipping():
    roots, complete = isolate_constraint_roots_exact(
        {"p": sp.Poly(X + EPS, X, domain=sp.QQ)}, Q(1), digits=30
    )
    assert complete
    assert roots == []


def test_lower_root_exactly_at_zero_is_not_reported_as_positive_root():
    roots, complete = isolate_constraint_roots_exact(
        {"p": sp.Poly(X, X, domain=sp.QQ)}, Q(1), digits=30
    )
    assert complete
    assert roots == []


def test_lower_root_just_inside_domain_is_retained_exactly():
    roots, complete = isolate_constraint_roots_exact(
        {"p": sp.Poly(X - EPS, X, domain=sp.QQ)}, Q(1), digits=30
    )
    assert complete
    assert len(roots) == 1
    root = roots[0]
    assert root.exact_root == EPS
    assert root.lower == EPS == root.upper
    assert root.lower > 0


@pytest.mark.parametrize("which", ["generic", "candidate"])
def test_public_locator_does_not_clip_strict_lower_interior_root_to_zero(which: str):
    method = _method(EPS)
    generic = locate_cptp_components(method, theta=0, kappa=Q(3, 4), varpi=0, upper=Q(1, 10), root_digits=30)
    candidate = locate_candidate_components(
        "direct", method=method, substeps=2, theta=0, kappa=Q(3, 4), varpi=0,
        upper=Q(1, 10), root_digits=30,
    )
    location = generic if which == "generic" else candidate
    exact = [root for root in location.roots if root.exact_root == EPS]
    assert len(exact) == 1
    assert exact[0].lower == EPS == exact[0].upper
    assert not any(point.lower == 0 == point.upper for point in location.isolated_points)


@pytest.mark.parametrize("which", ["generic", "candidate"])
def test_public_locator_excludes_controllable_root_just_below_zero(which: str):
    method = _method(-EPS)
    generic = locate_cptp_components(method, theta=0, kappa=Q(3, 4), varpi=0, upper=Q(1, 10), root_digits=30)
    candidate = locate_candidate_components(
        "direct", method=method, substeps=2, theta=0, kappa=Q(3, 4), varpi=0,
        upper=Q(1, 10), root_digits=30,
    )
    location = generic if which == "generic" else candidate
    assert all(root.exact_root != -EPS for root in location.roots)
    assert not any(point.lower == 0 == point.upper for point in location.isolated_points)


def test_near_coincident_roots_are_adaptively_separated():
    polynomials = {
        "left": sp.Poly(X - (Q(1) - EPS), X, domain=sp.QQ),
        "right": sp.Poly(X - (Q(1) + EPS), X, domain=sp.QQ),
    }
    roots, complete = isolate_constraint_roots_exact(polynomials, Q(2), digits=30, max_digits=100)
    assert complete
    assert [root.exact_root for root in roots] == [Q(1) - EPS, Q(1) + EPS]
    assert roots[0].upper < roots[1].lower


def test_boundary_straddling_near_coincident_roots_keep_only_in_domain_root():
    polynomials = {
        "inside": sp.Poly(X - (Q(1) - EPS), X, domain=sp.QQ),
        "outside": sp.Poly(X - (Q(1) + EPS), X, domain=sp.QQ),
    }
    roots, complete = isolate_constraint_roots_exact(polynomials, Q(1), digits=30, max_digits=100)
    assert complete
    assert len(roots) == 1
    assert roots[0].exact_root == Q(1) - EPS
    assert roots[0].upper < Q(1)


def test_near_coincident_root_order_is_independent_of_mapping_insertion():
    polynomials = {
        "right": sp.Poly(X - (Q(1) + EPS), X, domain=sp.QQ),
        "left": sp.Poly(X - (Q(1) - EPS), X, domain=sp.QQ),
    }
    roots, complete = isolate_constraint_roots_exact(
        polynomials, Q(2), digits=30, max_digits=100
    )
    assert complete
    assert [root.exact_root for root in roots] == [Q(1) - EPS, Q(1) + EPS]
    assert roots[0].upper < roots[1].lower
