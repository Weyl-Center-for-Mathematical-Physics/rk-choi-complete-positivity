from __future__ import annotations

import math

import mpmath as mp
import sympy as sp

from rk_choi_margin.methods import all_methods
from rk_choi_margin.topology import (
    degenerate_next_coefficient,
    frequency_factor_q,
    frequency_sigma,
    normalized_choi_defect_and_error,
    rk4_general_margin,
    rk4_rotating_cubic,
    rk4_rotating_margin,
)

x, q, y, v = sp.symbols("x q y v", real=True)


def test_frequency_recurrence_and_correct_g6() -> None:
    assert frequency_sigma(0, q) == 2
    assert frequency_sigma(1, q) == 2
    for m in range(2, 10):
        assert sp.expand(frequency_sigma(m, q) - 2 * frequency_sigma(m - 1, q) + q * frequency_sigma(m - 2, q)) == 0
    assert sp.factor(frequency_factor_q(4, q)) == -q * (q - 8) / 8
    assert sp.factor(frequency_factor_q(5, q)) == -5 * q * (q - 4) / 16
    assert sp.factor(frequency_factor_q(6, q)) == q * (q**2 - 18 * q + 48) / 32
    assert sp.factor(frequency_factor_q(7, q)) == 7 * q * (q - 4) ** 2 / 64


def test_unit_circle_mod_six_rule() -> None:
    for m in range(2, 31):
        value = sp.simplify(frequency_factor_q(m, sp.Integer(4)))
        assert (value == 0) == (m % 6 in (1, 5))


def test_degenerate_next_coefficients() -> None:
    coeff_ssp = degenerate_next_coefficient(4, -sp.Rational(1, 24), -sp.Rational(1, 120), sp.sqrt(7) / 2)
    coeff_rk4 = degenerate_next_coefficient(5, -sp.Rational(1, 120), -sp.Rational(1, 720), sp.sqrt(3) / 2)
    assert sp.simplify(coeff_ssp - sp.Rational(1, 3)) == 0
    assert sp.simplify(coeff_rk4 + sp.Rational(1, 144)) == 0


def test_rk4_rotating_closed_form() -> None:
    R = all_methods()["rk4"].stability_function(x=sp.Symbol("s"))
    s = sp.Symbol("s")
    direct = sp.factor(
        R.subs(s, -x)
        - R.subs(s, -x * (sp.Rational(1, 2) - sp.I * v))
        * R.subs(s, -x * (sp.Rational(1, 2) + sp.I * v))
    )
    assert sp.simplify(direct - rk4_rotating_margin(x, v)) == 0
    assert sp.simplify(sp.discriminant(rk4_rotating_cubic(q, y), y) - 16384 * (q - 4) * (8 * q**2 - 123 * q + 348)) == 0


def test_rk4_critical_factorizations() -> None:
    qm = (123 - 11 * sp.sqrt(33)) / 16
    qp = (123 + 11 * sp.sqrt(33)) / 16
    Cm = sp.factor(rk4_rotating_cubic(qm, y), extension=sp.sqrt(33))
    Cp = sp.factor(rk4_rotating_cubic(qp, y), extension=sp.sqrt(33))
    assert sp.simplify(Cm - (y - 9 + sp.sqrt(33)) ** 2 * (y - 2 * sp.sqrt(33) + 2)) == 0
    assert sp.simplify(Cp - (y - 9 - sp.sqrt(33)) ** 2 * (y + 2 + 2 * sp.sqrt(33))) == 0
    assert sp.factor(rk4_rotating_cubic(4, y)) == y * (y - 8) ** 2


def _positive_roots(expr: sp.Expr) -> list[float]:
    poly = sp.Poly(sp.together(expr).as_numer_denom()[0], x)
    roots = sp.nroots(poly, n=60, maxsteps=2000)
    return sorted(float(sp.re(r)) for r in roots if abs(float(sp.im(r))) < 1e-20 and float(sp.re(r)) > 1e-10)


def test_exact_open_set_examples_at_varpi_two() -> None:
    cases = [
        (sp.Rational(0), sp.Rational(1, 1000), [0.2533488691158, 0.7406959996519, 1.2695427111396]),
        (sp.Rational(1, 1000), sp.Rational(0), [0.1222689034698, 0.7425583125542, 1.2705443701360]),
        (sp.Rational(1, 1000), sp.Rational(1, 1000), [0.2620333331609, 0.7383646272135, 1.2697535708265]),
        (sp.Rational(0), sp.Rational(1, 100), [0.5193082188198, 0.6816216100681, 1.2622801549800]),
        (sp.Rational(1, 100), sp.Rational(0), [0.2769302758106, 0.7208968050590, 1.2723913949353]),
    ]
    for theta, kappa, expected in cases:
        roots = _positive_roots(rk4_general_margin(x, theta, kappa, sp.Integer(2)))
        assert len(roots) == 3
        for found, target in zip(roots, expected):
            assert abs(found - target) < 3e-12


def test_bell_defect_is_bounded_by_total_choi_error() -> None:
    mp.mp.dps = 80
    for varpi in (1.0, 2.0, 5.0):
        for step in (0.1, 0.01, 0.001):
            defect, error = normalized_choi_defect_and_error(step, varpi)
            assert defect >= 0
            assert error >= defect - mp.mpf("1e-60")
    defect, error = normalized_choi_defect_and_error(0.01, 2.0)
    ratio = defect / error
    assert mp.mpf("0.90") < ratio < mp.mpf("0.92")


def test_high_frequency_asymptotic_window() -> None:
    # Numerical roots of the exact cubic approach the stated scaled limits.
    for varpi in (25.0, 100.0):
        qv = 1 + 4 * varpi**2
        roots = sp.nroots(sp.Poly(rk4_rotating_cubic(sp.Float(qv, 50), y), y), n=50, maxsteps=1000)
        positive = sorted(float(sp.re(r)) / qv for r in roots if abs(float(sp.im(r))) < 1e-15 and float(sp.re(r)) > 0)
        assert len(positive) == 2
        lo, hi = positive
        assert abs(lo * varpi**2 - 3.0) < 0.01
        assert abs(hi * varpi - 2 * math.sqrt(2)) < 0.03


def _real_roots_y(q_value: sp.Expr) -> list[float]:
    roots = sp.nroots(sp.Poly(rk4_rotating_cubic(q_value, y), y), n=70, maxsteps=2000)
    return sorted(float(sp.re(r)) for r in roots if abs(float(sp.im(r))) < 1e-18)


def test_rk4_complete_topology_sign_regimes() -> None:
    alpha4 = 2.785293563405282
    qm = (123 - 11 * sp.sqrt(33)) / 16
    qp = (123 + 11 * sp.sqrt(33)) / 16

    # Connected regime before the first discriminant transition.
    roots = _real_roots_y(sp.Rational(2))
    positive = [r for r in roots if r > 0]
    assert len(positive) == 1
    assert positive[0] / 2 > alpha4  # population ceiling is active at q=2

    # Three positive roots create attached and detached components.
    qv = sp.Rational(19, 5)  # 3.8 lies between q_- and 4
    positive = [r for r in _real_roots_y(qv) if r > 0]
    assert len(positive) == 3
    samples = [positive[0] / 2, (positive[0] + positive[1]) / 2,
               (positive[1] + positive[2]) / 2, positive[2] + 1]
    signs = [sp.sign(sp.N(rk4_rotating_cubic(qv, sp.Float(s, 50)), 50)) for s in samples]
    assert signs == [-1, 1, -1, 1]
    assert positive[2] / float(qv) < alpha4

    # At q=4 the two components collapse to isolated admissible points.
    assert sp.factor(rk4_rotating_cubic(4, y)) == y * (y - 8) ** 2

    # Between 4 and q_+ there is no positive margin interval.
    roots = _real_roots_y(sp.Rational(5))
    assert all(r < 0 for r in roots)

    # At q_+ a double positive root gives an isolated remote point.
    qp_factor = sp.factor(rk4_rotating_cubic(qp, y), extension=sp.sqrt(33))
    assert qp_factor == (y - 9 - sp.sqrt(33)) ** 2 * (y + 2 + 2 * sp.sqrt(33))
    isolated = float(sp.N((9 + sp.sqrt(33)) / qp, 40))
    assert 0 < isolated < alpha4

    # Beyond q_+ the remote point opens into a detached interval.
    qv = sp.Integer(16)
    roots = _real_roots_y(qv)
    positive = [r for r in roots if r > 0]
    negative = [r for r in roots if r < 0]
    assert len(negative) == 1 and len(positive) == 2
    assert positive[1] / 16 < alpha4
