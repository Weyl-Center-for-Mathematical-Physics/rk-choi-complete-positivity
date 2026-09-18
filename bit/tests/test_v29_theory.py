from __future__ import annotations

import math
import sympy as sp

from rk_choi_margin.candidates import (
    certify_candidate,
    locate_candidate_components,
    richardson_rk4_margin_nonrotating,
    richardson_rk4_stability_function,
    rk4_step_doubling_bundle,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import halving_threshold_certificate
from rk_choi_margin.symbolic import first_exponential_defect
from rk_choi_margin.topology import frequency_factor_q


def test_richardson_is_a_stability_function_with_expected_defect() -> None:
    s = sp.Symbol("s")
    Rext = sp.expand(richardson_rk4_stability_function())
    expected = (
        1 + s + s**2/sp.Integer(2) + s**3/sp.Integer(6)
        + s**4/sp.Integer(24) + s**5/sp.Integer(120)
        + s**6/sp.Integer(864) + s**7/sp.Integer(8640)
        + s**8/sp.Integer(138240)
    )
    assert sp.expand(Rext - expected) == 0
    defect = first_exponential_defect(Rext, x=s)
    assert defect.index == 6
    assert defect.coefficient == -sp.Rational(1, 4320)


def test_richardson_frequency_factor_gives_outward_inward_outward() -> None:
    q = sp.Symbol("q", real=True)
    G6 = sp.factor(frequency_factor_q(6, q))
    assert G6 == q * (q**2 - 18*q + 48) / 32
    eta6 = -sp.Rational(1, 4320)
    for varpi, expected_sign in [(0, -1), (1, 1), (sp.Rational(3,2), 1), (2, -1), (5, -1)]:
        qv = 1 + 4*varpi**2
        coeff = sp.factor(eta6 * G6.subs(q, qv))
        assert sp.sign(coeff) == expected_sign


def test_richardson_global_nonrotating_set_and_common_interval_failure() -> None:
    loc = locate_candidate_components(
        "richardson-4", method="rk4", theta=0, kappa=0, varpi=0, upper=8
    )
    assert len(loc.components) == 1
    component = loc.components[0]
    assert abs(float(component.lower_value) - 6.034523958649) < 1e-10
    assert abs(float(component.upper_value) - 6.459127767826) < 1e-10
    for x in (sp.Rational(1,100), sp.Rational(1,10), sp.Integer(1), sp.Rational(5,2)):
        bundle = rk4_step_doubling_bundle(total_x=x, theta=0, kappa=0, varpi=0)
        assert certify_candidate(bundle.coarse).status is CPStatus.PASS
        assert certify_candidate(bundle.fine).status is CPStatus.PASS
        assert certify_candidate(bundle.extrapolated).status is CPStatus.FAIL


def test_richardson_exact_margin_factorization() -> None:
    x = sp.Symbol("x", nonnegative=True, real=True)
    p10 = (
        x**10 - 64*x**9 + 2304*x**8 - 59392*x**7 + 1183744*x**6
        - 19169280*x**5 + 258932736*x**4 - 2960916480*x**3
        + 19888865280*x**2 - 97391738880*x + 280850595840
    )
    margin = sp.factor(richardson_rk4_margin_nonrotating())
    assert sp.factor(margin + x**6*p10/sp.Integer(1252412463513600)) == 0
    assert sp.expand(margin).coeff(x,6) == -sp.Rational(31,138240)


def test_disjoint_full_and_fine_threshold_matches_halving_threshold() -> None:
    cert = halving_threshold_certificate(42)
    threshold = float(cert["varpi_mid"])
    assert abs(threshold - 2.2482546045149805) < 1e-14
    direct = locate_candidate_components(
        "direct", method="rk4", theta=0, kappa=0, varpi=2, upper=3
    ).components[0]
    ratio = float(direct.upper_value/direct.lower_value)
    assert ratio < 2
    # Consecutive n- and n+1-substep candidate intervals overlap iff
    # x_+/x_- >= (n+1)/n.  At varpi=2 only the 1-to-2 pair is disjoint.
    assert ratio < 2/1
    assert ratio >= 3/2


def test_margin_leading_coefficient_from_first_defect_theorem() -> None:
    
    Rext = richardson_rk4_stability_function()
    s = next(iter(Rext.free_symbols))
    defect = first_exponential_defect(Rext, x=s)
    G6_0 = sp.Rational(31,32)
    coefficient = (-1)**defect.index * defect.coefficient * G6_0
    assert coefficient == -sp.Rational(31,138240)


def test_richardson_remote_component_active_boundaries() -> None:
    loc = locate_candidate_components(
        "richardson-4", method="rk4", theta=0, kappa=0, varpi=0, upper=8
    )
    component = loc.components[0]
    assert component.left is not None and "F" in component.left.source
    assert component.right is not None and "1-a" in component.right.source


def test_substep_ladder_overlap_condition_at_varpi_two() -> None:
    component = locate_candidate_components(
        "direct", method="rk4", theta=0, kappa=0, varpi=2, upper=3
    ).components[0]
    ratio = sp.N(component.upper_value / component.lower_value, 30)
    for n in range(1, 8):
        predicted = ratio >= sp.Rational(n + 1, n)
        interval_n = (n * float(component.lower_value), n * float(component.upper_value))
        interval_np1 = ((n + 1) * float(component.lower_value), (n + 1) * float(component.upper_value))
        actual = max(interval_n[0], interval_np1[0]) <= min(interval_n[1], interval_np1[1])
        assert bool(predicted) == actual
