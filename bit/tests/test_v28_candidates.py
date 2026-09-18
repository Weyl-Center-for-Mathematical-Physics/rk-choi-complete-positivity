from __future__ import annotations

import numpy as np
import sympy as sp

from rk_choi_margin.candidates import (
    affine_candidate,
    certify_candidate,
    direct_candidate,
    equal_substep_candidate,
    locate_candidate_components,
    phase_covariant_superoperator_numeric,
    richardson_rk4_margin_nonrotating,
    rk4_step_doubling_bundle,
)
from rk_choi_margin.certification import CPStatus, as_rational
from rk_choi_margin.controllers import run_adaptive_channel_benchmark
from rk_choi_margin.liouvillian import min_choi_eigenvalue_numeric


def test_float_inputs_are_certified_as_exact_ieee_binary_values() -> None:
    value = 1.270334626973889
    numerator, denominator = value.as_integer_ratio()
    assert as_rational(value) == sp.Rational(numerator, denominator)
    assert as_rational(str(value)) == sp.Rational(str(value))
    assert as_rational(value) != as_rational(str(value))


def test_binary_boundary_case_is_classified_for_executed_value() -> None:
    value = 1.270334626973889
    binary = direct_candidate("rk4", total_x=value, theta=0, kappa=0, varpi=2)
    decimal = direct_candidate("rk4", total_x=str(value), theta=0, kappa=0, varpi=2)
    assert certify_candidate(binary).status is CPStatus.FAIL
    assert certify_candidate(decimal).status is CPStatus.PASS


def test_full_and_two_half_step_candidate_regions_are_non_nested() -> None:
    at_one = rk4_step_doubling_bundle(total_x=1, theta=0, kappa=0, varpi=2)
    assert certify_candidate(at_one.coarse).status is CPStatus.PASS
    assert certify_candidate(at_one.fine).status is CPStatus.FAIL

    at_two = rk4_step_doubling_bundle(total_x=2, theta=0, kappa=0, varpi=2)
    assert certify_candidate(at_two.coarse).status is CPStatus.FAIL
    assert certify_candidate(at_two.fine).status is CPStatus.PASS


def test_two_half_step_total_region_is_dilation_of_direct_region() -> None:
    direct = locate_candidate_components(
        "direct", method="rk4", theta=0, kappa=0, varpi=2, upper=3
    )
    fine = locate_candidate_components(
        "equal-substeps", method="rk4", theta=0, kappa=0, varpi=2, upper=3, substeps=2
    )
    assert len(direct.components) == 1
    assert len(fine.components) == 1
    d = direct.components[0]
    f = fine.components[0]
    assert abs(float(f.lower_value) - 2 * float(d.lower_value)) < 1e-11
    assert abs(float(f.upper_value) - 2 * float(d.upper_value)) < 1e-11


def test_convex_mixture_of_cptp_candidates_is_cptp() -> None:
    c1 = direct_candidate("rk4", total_x="0.2", theta=0, kappa=0, varpi=0)
    c2 = equal_substep_candidate("rk4", total_x="0.2", substeps=2, theta=0, kappa=0, varpi=0)
    assert certify_candidate(c1).status is CPStatus.PASS
    assert certify_candidate(c2).status is CPStatus.PASS
    mixture = affine_candidate((c1, c2), (sp.Rational(1, 3), sp.Rational(2, 3)), kind="convex")
    assert certify_candidate(mixture).status is CPStatus.PASS


def test_richardson_extrapolate_is_non_cptp_while_constituents_are_cptp() -> None:
    for value in (sp.Rational(1, 100), sp.Rational(1, 10), sp.Integer(1), sp.Rational(5, 2)):
        bundle = rk4_step_doubling_bundle(total_x=value, theta=0, kappa=0, varpi=0)
        assert certify_candidate(bundle.coarse).status is CPStatus.PASS
        assert certify_candidate(bundle.fine).status is CPStatus.PASS
        assert certify_candidate(bundle.extrapolated).status is CPStatus.FAIL


def test_richardson_margin_factorization_and_leading_coefficient() -> None:
    x = sp.Symbol("x", nonnegative=True, real=True)
    margin = sp.factor(richardson_rk4_margin_nonrotating())
    polynomial = (
        x**10 - 64*x**9 + 2304*x**8 - 59392*x**7 + 1183744*x**6
        - 19169280*x**5 + 258932736*x**4 - 2960916480*x**3
        + 19888865280*x**2 - 97391738880*x + 280850595840
    )
    expected = -x**6 * polynomial / sp.Integer(1252412463513600)
    assert sp.factor(margin - expected) == 0
    assert sp.expand(margin).coeff(x, 6) == -sp.Rational(31, 138240)
    roots = sp.Poly(polynomial, x, domain=sp.QQ).intervals()
    assert all(not (lo < sp.Rational("2.786") and hi > 0) for (lo, hi), _ in roots)


def test_candidate_matrix_is_the_map_whose_certificate_is_checked() -> None:
    bundle = rk4_step_doubling_bundle(total_x=1, theta=0, kappa=0, varpi=2)
    fine_matrix = phase_covariant_superoperator_numeric(bundle.fine)
    assert min_choi_eigenvalue_numeric(fine_matrix) < 0
    assert certify_candidate(bundle.fine).status is CPStatus.FAIL


def test_adaptive_candidate_guard_never_accumulates_uncertified_fine_map() -> None:
    guarded = run_adaptive_channel_benchmark(
        final_time=2.0,
        initial_h=1.0,
        tolerance=0.01,
        theta=0,
        kappa=0,
        varpi=2,
        policy="candidate_guard",
        fallback="rotating_frame",
    )
    assert guarded.min_step_choi_eigenvalue > -1e-11
    assert all(row.candidate_kind == row.certified_kind or row.action.startswith("fallback") for row in guarded.records)
