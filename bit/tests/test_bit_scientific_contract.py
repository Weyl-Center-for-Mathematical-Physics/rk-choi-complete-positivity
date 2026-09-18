"""BIT-revision scientific contract tests (Task 5/6).

These tests pin the corrected statements, explicit constants, and the one targeted extension adopted for
the BIT manuscript.  They are written before the implementation (test-driven) and must stay fast.
"""
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest
import sympy as sp

from rk_choi_margin.methods import all_methods
from rk_choi_margin.symbolic import first_exponential_defect
from rk_choi_margin.extrapolation import (
    frequency_factor,
    local_orientation_coefficient,
    predicted_extrapolation_defect,
    richardson_extrapolate,
)

Z = sp.Symbol("x")
Q = sp.Rational


def _defect_pair(R, p):
    """Return (eta_{p+1}, eta_{p+2}) of R(z)-exp(z)."""
    ser = sp.series(R - sp.exp(Z), Z, 0, p + 3).removeO()
    return sp.simplify(ser.coeff(Z, p + 1)), sp.simplify(ser.coeff(Z, p + 2))


# ---------------------------------------------------------------------------
# Extension: n-fold Richardson extrapolation of an order-p RK method
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("key", ["euler", "rk2", "rk3", "rk4", "dp5", "dp4", "be", "im"])
@pytest.mark.parametrize("n", [2, 3, 4])
def test_extrapolation_defect_formula(key, n):
    method = all_methods()[key]
    R = method.stability_function()
    p = method.order
    eta1, eta2 = _defect_pair(R, p)
    E = richardson_extrapolate(R, p, n)
    defect = first_exponential_defect(E)
    predicted = predicted_extrapolation_defect(eta1, eta2, p, n)
    if eta1 == eta2:
        # degenerate case: the (p+2) coefficient cancels; the defect order is at least p+3
        assert defect.index >= p + 3
    else:
        assert defect.index == p + 2
        assert sp.simplify(defect.coefficient - predicted) == 0


def test_rk4_two_fold_extrapolate_recovers_archive_certificate():
    R4 = all_methods()["rk4"].stability_function()
    E = sp.expand(richardson_extrapolate(R4, 4, 2))
    assert E.coeff(Z, 6) == Q(1, 864) and E.coeff(Z, 7) == Q(1, 8640) and E.coeff(Z, 8) == Q(1, 138240)
    defect = first_exponential_defect(E)
    assert (defect.index, defect.coefficient) == (6, -Q(1, 4320))
    assert predicted_extrapolation_defect(-Q(1, 120), -Q(1, 720), 4, 2) == -Q(1, 4320)


def test_euler_two_fold_extrapolate_is_heun_polynomial():
    euler = all_methods()["euler"].stability_function()
    heun = all_methods()["rk2"].stability_function()
    assert sp.expand(richardson_extrapolate(euler, 1, 2) - heun) == 0


def test_truncated_exponential_orientation_by_parity():
    """For truncated-exponential methods of order p, the 2-fold extrapolate is locally non-CP at
    varpi=0 when p is even and locally CP when p is odd."""
    for key, p in (("euler", 1), ("rk2", 2), ("rk3", 3), ("rk4", 4)):
        R = all_methods()[key].stability_function()
        E = richardson_extrapolate(R, p, 2)
        defect = first_exponential_defect(E)
        coeff = local_orientation_coefficient(defect.index, defect.coefficient, 0)
        assert (coeff < 0) == (p % 2 == 0), (key, coeff)


def test_frequency_factor_matches_closed_forms():
    v = sp.Symbol("varpi", real=True)
    q = 1 + 4 * v**2
    assert sp.simplify(frequency_factor(6, v) - q * (q**2 - 18 * q + 48) / 32) == 0
    assert sp.simplify(frequency_factor(5, v) + 5 * q * (q - 4) / 16) == 0
    assert sp.simplify(frequency_factor(2, v) - q / 2) == 0
    assert sp.simplify(frequency_factor(3, v) - 3 * q / 4) == 0
    assert sp.simplify(frequency_factor(4, v) + q * (q - 8) / 8) == 0
    assert sp.simplify(frequency_factor(7, v) - 7 * q * (q - 4) ** 2 / 64) == 0
    # at zero precession G_m(0) = 1 - 2^{1-m}
    assert all(frequency_factor(m, 0) == 1 - Q(1, 2 ** (m - 1)) for m in range(2, 9))


def test_dp5_two_fold_extrapolate_is_locally_outward_at_zero_frequency():
    dp5 = all_methods()["dp5"].stability_function()
    E = richardson_extrapolate(dp5, 5, 2)
    defect = first_exponential_defect(E)
    assert defect.index == 7
    assert local_orientation_coefficient(7, defect.coefficient, 0) < 0


# ---------------------------------------------------------------------------
# Explicit constants adopted in the BIT statements
# ---------------------------------------------------------------------------
def test_rk4_regime_one_crossover_frequency():
    """q_c: exact certificate (reduced discriminant sign, resultant isolation, sign change) and the y* inequality."""
    y, q = sp.symbols("y q", real=True)
    t = sp.Symbol("t")
    Cq = y**3 - 16 * y**2 - 32 * (q - 6) * y + 384 * (q - 4)
    minpoly = t**3 - 4 * t**2 + 12 * t - 24
    alpha4 = sp.CRootOf(minpoly, 0)
    g = sp.expand(Cq.subs(y, q * t))
    disc = sp.rem(sp.expand(sp.discriminant(g, q)), minpoly, t)
    assert sp.expand(disc - (-1774190592 * t**2 + 3472883712 * t + 3170893824)) == 0
    lo, hi = sp.Poly(minpoly, t).refine_root(2, 3, eps=Q(1, 10**30))
    assert disc.subs(t, lo) < 0 and disc.subs(t, hi) < 0
    assert sp.diff(disc, t).subs(t, lo) < 0 and sp.diff(disc, t).subs(t, hi) < 0  # monotone on the bracket
    res = sp.Poly(sp.expand(sp.resultant(g, minpoly, t)), q)
    assert res.degree() == 9
    brackets = res.intervals()
    assert len(brackets) == 1 and brackets[0][1] == 1
    a, b = res.refine_root(*brackets[0][0], eps=Q(1, 10**20))
    qm = (123 - 11 * sp.sqrt(33)) / 16
    assert 1 < a and b < qm

    def sign_g(qv):
        val = sp.rem(sp.expand(g.subs(q, qv)), minpoly, t)
        s0, s1 = sp.sign(val.subs(t, lo)), sp.sign(val.subs(t, hi))
        assert s0 == s1 != 0
        return int(s0)

    assert sign_g(a) == -1 and sign_g(b) == 1
    qc = (a + b) / 2
    assert abs(qc - 3.5282058982943983) < 1e-12
    varpi_c = sp.sqrt((qc - 1) / 4)
    assert abs(varpi_c - 0.7950166505008556) < 1e-12
    # q*alpha4 exceeds the larger critical point of C_q for q >= 10/3: negative discriminant of the quadratic in q
    quad = sp.expand((3 * t * q - 16) ** 2 - (96 * q - 320))
    dq = sp.discriminant(quad, q)
    assert sp.simplify(dq + 2304 * (t - 2) * (5 * t + 2)) == 0
    assert dq.subs(t, lo) < 0 and dq.subs(t, hi) < 0
    # sign of g at rational samples: negative below q_c, positive above
    assert sign_g(Q(2)) == -1 and all(sign_g(v) == 1 for v in (Q(37, 10), Q(39, 10), Q(17), Q(100), Q(10**6)))
    # below q_c the cubic root y3/q exceeds alpha4 (ceiling active); above, it is smaller (margin active)
    for qv, expect_ceiling in ((Q(2), True), (Q(37, 10), False)):
        P = sp.Poly(sp.expand(Cq.subs(q, qv)), y)
        y3 = max(float((lo + hi) / 2) for (lo, hi), m in P.intervals(eps=Q(1, 10**20)) if hi > 0)
        assert (y3 / float(qv) > float(alpha4)) == expect_ceiling


def test_frame_candidates_differ_prop6ii():
    """Proposition 6(ii): laboratory- and rotating-frame coherence multipliers differ except at isolated steps."""
    kap, vp, nu, xv = sp.symbols("kappa varpi nu x", real=True)
    s = sp.Symbol("s", positive=True)
    z = sp.Symbol("z")
    R4 = 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24
    w = -(Q(1, 2) + kap) * xv + sp.I * vp * xv
    lab = sp.Poly(sp.expand(R4.subs(z, w) * R4.subs(z, sp.conjugate(w))), xv).coeff_monomial(xv**8)
    rot = sp.Poly(sp.expand(R4.subs(z, -(Q(1, 2) + kap) * xv) ** 2), xv).coeff_monomial(xv**8)
    assert sp.simplify(lab - rot - (((Q(1, 2) + kap) ** 2 + vp**2) ** 4 - (Q(1, 2) + kap) ** 8) / 576) == 0
    be = sp.simplify(sp.Abs(1 / (1 - (-s + sp.I * nu))) ** 2 - (1 / (1 + s)) ** 2)
    assert sp.simplify(be + nu**2 / ((1 + s) ** 2 * ((1 + s) ** 2 + nu**2))) == 0
    mid = lambda zz: (1 + zz / 2) / (1 - zz / 2)
    mp_ = sp.simplify(sp.Abs(mid(-s + sp.I * nu)) ** 2 - mid(-s) ** 2)
    assert sp.simplify(mp_ - 8 * nu**2 * s / ((2 + s) ** 2 * ((2 + s) ** 2 + nu**2))) == 0


def test_printed_decimal_versus_binary64_margin_at_varpi2_endpoint():
    R4 = all_methods()["rk4"].stability_function()
    x = sp.Symbol("xx", real=True)
    a = R4.subs(Z, -x)
    c = R4.subs(Z, -x * (Q(1, 2) - 2 * sp.I))
    M = sp.expand(a - c * sp.conjugate(c))
    dec = M.subs(x, Q("1.270334626973889"))
    flt = M.subs(x, Q(*Fraction(1.270334626973889).as_integer_ratio()))
    assert dec > 0 > flt


def test_varpi2_full_and_two_half_step_positive_components_are_disjoint():
    R4 = all_methods()["rk4"].stability_function()
    x = sp.Symbol("xx", real=True)
    M = sp.Poly(sp.expand(R4.subs(Z, -x) - R4.subs(Z, -x * (Q(1, 2) - 2 * sp.I)) * R4.subs(Z, -x * (Q(1, 2) + 2 * sp.I))), x)
    roots = [(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in M.intervals(eps=Q(1, 10**18)) if hi > 0]
    lo1, hi1 = roots[0]
    lo2, hi2 = roots[1]
    assert hi2 < 2 * lo1  # x_+ < 2 x_-  <=>  disjoint positive components


def _claim_record(name):
    """Read computed proof data; tests recheck its mathematics, not its pass flag."""
    path = Path(__file__).resolve().parents[1] / "results/bit_revision/claim_verification.json"
    return json.loads(path.read_text(encoding="utf-8"))["checks"][name]


def test_printable_crossover_brackets_certify_the_claimed_signs():
    record = _claim_record("printed_crossover_brackets")
    t, q = sp.symbols("t q")
    f = t**3 - 4*t**2 + 12*t - 24
    D = -1774190592*t**2 + 3472883712*t + 3170893824
    lo, hi = map(Q, record["alpha_bracket"])
    assert sp.Poly(f, t).count_roots(lo, hi) == 1
    assert D.subs(t, lo) < 0 and D.subs(t, hi) < 0
    assert sp.diff(D, t).subs(t, lo) < 0 and sp.diff(D, t, 2) < 0
    g = (q*t)**3 - 16*(q*t)**2 - 32*(q-6)*q*t + 384*(q-4)
    res = sp.Poly(sp.resultant(g, f, t), q)
    qlo, qhi = map(Q, record["q_c_bracket"])
    assert res.count_roots(-sp.oo, sp.oo) == res.count_roots(qlo, qhi) == 1
    assert 1 < qlo < qhi < (123-11*sp.sqrt(33))/16


def test_richardson_global_chart_accounts_for_all_positive_events_and_tail():
    record = _claim_record("richardson_global_sign_chart")
    E = richardson_extrapolate(all_methods()["rk4"].stability_function(), 4, 2)
    a = sp.Poly(E.subs(Z, -Z), Z)
    margin = sp.Poly(sp.expand(a.as_expr()-E.subs(Z, -Z/2)**2), Z)
    polys = {"a": a, "a-1": a-1, "M": margin}
    events = record["events"]
    assert [e["source"] for e in events] == ["a", "a", "M", "a-1", "M"]
    previous = Q(0)
    for event in events:
        lo, hi = map(Q, event["bracket"])
        assert previous < lo < hi
        assert polys[event["source"]].count_roots(lo, hi) == 1
        assert polys[event["source"]].gcd(polys[event["source"]].diff()).count_roots(lo, hi) == 0
        previous = hi
    for source, poly in polys.items():
        without_zero = sp.Poly(poly.as_expr()/Z**min(m[0] for m, c in poly.terms()), Z)
        assert without_zero.count_roots(0, sp.oo) == sum(e["source"] == source for e in events)
    want = [(1,-1,-1), (-1,-1,-1), (1,-1,-1), (1,-1,1), (1,1,1), (1,1,-1)]
    for i, row in enumerate(record["intervals"]):
        sample = Q(row["sample"])
        assert (i == 0 or Q(events[i-1]["bracket"][1]) < sample)
        assert (i == len(events) or sample < Q(events[i]["bracket"][0]))
        assert tuple(int(sp.sign(p.eval(sample))) for p in polys.values()) == want[i]
        assert tuple(row["signs"]) == want[i]
    assert len(record["intervals"]) == 6


def test_rk4_persistence_slack_holds_on_whole_detached_neighbourhood():
    lo, hi = map(Q, _claim_record("persistence_inactive_slack_interval")["interval"])
    R = all_methods()["rk4"].stability_function()
    a = sp.Poly(R.subs(Z, -Z), Z)
    M = sp.Poly(sp.expand(a.as_expr()-R.subs(Z, -(Q(1,2)-2*sp.I)*Z)*R.subs(Z, -(Q(1,2)+2*sp.I)*Z)), Z)
    assert M.count_roots(lo, hi) == 2  # both detached endpoints lie inside the neighbourhood
    for inactive in (a, 1-a):
        assert inactive.count_roots(lo, hi) == 0
        assert inactive.eval(lo) > 0 and inactive.eval(hi) > 0
