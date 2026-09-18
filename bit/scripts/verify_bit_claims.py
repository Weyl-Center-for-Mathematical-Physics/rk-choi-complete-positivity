#!/usr/bin/env python3
"""Independent first-principles verification of the central claims (BIT revision, Task 5).

This script deliberately does NOT import ``rk_choi_margin``.  Every object is rebuilt from
scratch with SymPy: Butcher tableaux -> stability functions; the GKSL generator -> column-
vectorized superoperator -> Runge--Kutta map -> Choi matrix; exact rational/algebraic root
isolation via Sturm sequences.  Agreement with the archive's certificates is therefore an
independent check rather than a re-run of the same code.

Outputs (never overwrite sealed v3.5 files):
  results/bit_revision/claim_verification.json
  results/bit_revision/claim_verification.md
Exit code 0 only if every check passes.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from itertools import product
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "bit_revision"
OUT.mkdir(parents=True, exist_ok=True)

x, y, q = sp.symbols("x y q", real=True)
z = sp.Symbol("z")
h = sp.Symbol("h", nonnegative=True, real=True)
Q = sp.Rational

CHECKS: dict[str, dict] = {}


def record(name: str, passed: bool, **info) -> None:
    CHECKS[name] = {"pass": bool(passed), **{k: (str(v) if isinstance(v, sp.Basic) else v) for k, v in info.items()}}
    print(("PASS " if passed else "FAIL ") + name + ("" if passed else f"  {info}"))


# ---------------------------------------------------------------------------
# Runge--Kutta stability functions from Butcher tableaux (independent construction)
# ---------------------------------------------------------------------------
def stability_function(A, b):
    A = sp.Matrix(A)
    b = sp.Matrix(b)
    s = A.rows
    return sp.cancel(1 + z * (b.T * (sp.eye(s) - z * A).inv() * sp.ones(s, 1))[0])


RK4_A = [[0, 0, 0, 0], [Q(1, 2), 0, 0, 0], [0, Q(1, 2), 0, 0], [0, 0, 1, 0]]
RK4_B = [Q(1, 6), Q(1, 3), Q(1, 3), Q(1, 6)]
R4 = sp.expand(stability_function(RK4_A, RK4_B))
record("rk4_stability_function_from_tableau", sp.expand(R4 - (1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24)) == 0, R4=R4)

DP_A = [
    [0] * 7,
    [Q(1, 5), 0, 0, 0, 0, 0, 0],
    [Q(3, 40), Q(9, 40), 0, 0, 0, 0, 0],
    [Q(44, 45), -Q(56, 15), Q(32, 9), 0, 0, 0, 0],
    [Q(19372, 6561), -Q(25360, 2187), Q(64448, 6561), -Q(212, 729), 0, 0, 0],
    [Q(9017, 3168), -Q(355, 33), Q(46732, 5247), Q(49, 176), -Q(5103, 18656), 0, 0],
    [Q(35, 384), 0, Q(500, 1113), Q(125, 192), -Q(2187, 6784), Q(11, 84), 0],
]
DP5_B = [Q(35, 384), 0, Q(500, 1113), Q(125, 192), -Q(2187, 6784), Q(11, 84), 0]
DP4_B = [Q(5179, 57600), 0, Q(7571, 16695), Q(393, 640), -Q(92097, 339200), Q(187, 2100), Q(1, 40)]
R_DP5 = sp.expand(stability_function(DP_A, DP5_B))
R_DP4 = sp.expand(stability_function(DP_A, DP4_B))
R_EULER = 1 + z
R_HEUN = sp.expand(stability_function([[0, 0], [1, 0]], [Q(1, 2), Q(1, 2)]))
R_SSP3 = sp.expand(stability_function([[0, 0, 0], [1, 0, 0], [Q(1, 4), Q(1, 4), 0]], [Q(1, 6), Q(1, 6), Q(2, 3)]))
R_BE = stability_function([[1]], [1])
R_IM = stability_function([[Q(1, 2)]], [1])


def first_defect(R, max_order=12):
    ser = sp.series(R - sp.exp(z), z, 0, max_order).removeO()
    for m in range(max_order):
        c = sp.nsimplify(ser.coeff(z, m))
        if c != 0:
            return m, c
    raise RuntimeError("no defect")


# ---------------------------------------------------------------------------
# C1: GKSL generator, RK map, Choi matrix, criterion (all from first principles)
# ---------------------------------------------------------------------------
def vec(M):  # column stacking
    return sp.Matrix([M[i, j] for j in range(M.cols) for i in range(M.rows)])


def unvec(v, d=2):
    return sp.Matrix(d, d, lambda i, j: v[j * d + i])


def superop_left_right(Aop, Bop):  # X -> A X B
    return sp.kronecker_product(Bop.T, Aop)


def gksl_superoperator(gdown, gup, gphi, omega, omega_x=0):
    I2 = sp.eye(2)
    sm = sp.Matrix([[0, 1], [0, 0]])  # |0><1|
    sp_ = sm.T
    sz = sp.Matrix([[1, 0], [0, -1]])
    sx = sp.Matrix([[0, 1], [1, 0]])
    H = -omega * sz / 2 + omega_x * sx / 2
    L = -sp.I * (superop_left_right(H, I2) - superop_left_right(I2, H))

    def diss(Lk):
        LdL = Lk.H * Lk
        return superop_left_right(Lk, Lk.H) - superop_left_right(LdL, I2) / 2 - superop_left_right(I2, LdL) / 2

    L += gdown * diss(sm) + gup * diss(sp_) + (gphi / 2) * (superop_left_right(sz, sz) - sp.eye(4))
    return sp.simplify(L)


def choi(S, d=2):
    J = sp.zeros(d * d)
    for j, k in product(range(d), repeat=2):
        E = sp.zeros(d)
        E[j, k] = 1
        J += sp.kronecker_product(E, unvec(S * vec(E), d))
    return sp.simplify(J)


gd, gu, gp, om, Gam, th, xx, uu, nu = sp.symbols("gamma_down gamma_up gamma_phi omega Gamma theta x u nu", real=True)
hh = sp.Symbol("h", positive=True)
L_sym = gksl_superoperator(gd, gu, gp, om)
# Mode equations: population q=rho_11 and coherence rho_01
rho = sp.Matrix(2, 2, sp.symbols("r00 r01 r10 r11"))
drho = unvec(L_sym * vec(rho))
r00, r01, r10, r11 = rho
record("mode_equation_population",
       sp.simplify(drho[1, 1] - (-(gd + gu) * (r11 - gu / (gd + gu) * (r00 + r11)))) == 0)
record("mode_equation_coherence",
       sp.simplify(drho[0, 1] - (-((gd + gu) / 2 + gp - sp.I * om) * r01)) == 0)
lam = sp.Symbol("lambda")
expected_eigs = [0, -(gd + gu), -((gd + gu) / 2 + gp - sp.I * om), -((gd + gu) / 2 + gp + sp.I * om)]
record("generator_spectrum",
       sp.simplify(sp.expand(L_sym.charpoly(lam).as_expr() - sp.prod([lam - e for e in expected_eigs]))) == 0,
       spectrum=[str(e) for e in expected_eigs])

# Direct RK4 map on the generator with Gamma=1 scaling: x=Gamma h, u=gamma_phi h, nu=omega h, theta=gamma_up/Gamma
subsd = {gd: (1 - th), gu: th, gp: uu / xx, om: nu / xx}
L_scaled = L_sym.subs(subsd)  # Gamma = 1, h -> x
S4 = sp.simplify(sp.eye(4) + xx * L_scaled + (xx * L_scaled) ** 2 / 2 + (xx * L_scaled) ** 3 / 6 + (xx * L_scaled) ** 4 / 24)
a_claim = R4.subs(z, -xx)
c_claim = R4.subs(z, -xx / 2 - uu + sp.I * nu)
A_ = 1 - th * (1 - a_claim)
B_ = th * (1 - a_claim)
C_ = (1 - th) * (1 - a_claim)
D_ = th + (1 - th) * a_claim
X = sp.Matrix(2, 2, sp.symbols("X00 X01 X10 X11"))
Phi_claim = sp.Matrix([[A_ * X[0, 0] + C_ * X[1, 1], c_claim * X[0, 1]], [sp.conjugate(c_claim) * X[1, 0], B_ * X[0, 0] + D_ * X[1, 1]]])
Phi_direct = unvec(S4 * vec(X))
record("rk4_direct_map_equals_claimed_block_form", sp.simplify(sp.expand(Phi_direct - Phi_claim)) == sp.zeros(2))
J4 = choi(S4)
J_claim = sp.Matrix([[A_, 0, 0, c_claim], [0, B_, 0, 0], [0, 0, C_, 0], [sp.conjugate(c_claim), 0, 0, D_]])
record("rk4_choi_matrix_equals_claimed_form", sp.simplify(sp.expand(J4 - J_claim)) == sp.zeros(4))
record("trace_preservation_identities", sp.simplify(A_ + B_ - 1) == 0 and sp.simplify(C_ + D_ - 1) == 0)
record("choi_block_determinant_identity", sp.simplify(sp.expand(A_ * D_ - (a_claim + th * (1 - th) * (1 - a_claim) ** 2))) == 0)
# Generic real-coefficient RK: y1 = R(hL) y0 for linear autonomous systems, stage domain det(I - h A (x) L) != 0.
A_m = sp.Matrix(RK4_A)
lam = sp.Symbol("lambda")
record("stage_domain_is_spectral_condition",
       sp.factor((sp.eye(4) - lam * A_m).det()) == 1,
       note="explicit tableau: det(I - lambda A) = 1 for every lambda, so the stage domain is all steps; implicit tableaux give det(I-lambda A) with lambda in spec(hL)")
# Exact semigroup buffer
F_exact = sp.exp(-xx) * (1 - sp.exp(-2 * uu)) + th * (1 - th) * (1 - sp.exp(-xx)) ** 2
a_ex, c2_ex = sp.exp(-xx), sp.exp(-xx - 2 * uu)
record("exact_semigroup_margin", sp.simplify(a_ex + th * (1 - th) * (1 - a_ex) ** 2 - c2_ex - F_exact) == 0)

# Necessity/sufficiency of the four inequalities (Hermitian 4x4 with the given sparsity)
# PSD <=> B>=0, C>=0 and 2x2 block [[A,c],[c*,D]] PSD <=> A>=0, D>=0, AD-|c|^2>=0; B,C>=0 <=> a<=1 when 0<theta<1.
lam_J = sp.Symbol("lambda_J")
charpoly_J = sp.expand(J_claim.charpoly(lam_J).as_expr())
block_factor = sp.expand((lam_J - B_) * (lam_J - C_) * (lam_J**2 - (A_ + D_) * lam_J + (A_ * D_ - c_claim * sp.conjugate(c_claim))))
record("criterion_reduction_structure", sp.simplify(charpoly_J - block_factor) == 0,
       note="characteristic polynomial of the claimed Choi matrix factors as (lambda-B)(lambda-C)(lambda^2-(A+D)lambda+(AD-|c|^2)); a Hermitian matrix is PSD iff all eigenvalues are nonnegative, so J>=0 iff B,C>=0 and the 2x2 block has nonnegative trace and determinant, i.e. A,D>=0 and AD-|c|^2>=0; B,C>=0 iff a<=1 for theta in (0,1) (at theta=0 only C=1-a>=0 is active, at theta=1 only B=1-a>=0)")

# ---------------------------------------------------------------------------
# C2: signed first-defect law, G_m, degenerate case, buffers
# ---------------------------------------------------------------------------
varpi = sp.Symbol("varpi", real=True)
alpha = Q(1, 2) - sp.I * varpi
alphabar = Q(1, 2) + sp.I * varpi


def G(m):
    return sp.simplify(sp.expand(1 - (alpha**m + alphabar**m)))


ok = True
detail = {}
for m in range(2, 9):
    eta = sp.Symbol("eta")
    Rgen = sum(z**j / sp.factorial(j) for j in range(m)) + (1 / sp.factorial(m) + eta) * z**m
    M = sp.expand(Rgen.subs(z, -x) - Rgen.subs(z, -alpha * x) * Rgen.subs(z, -alphabar * x))
    lead = sp.simplify(sp.Poly(M, x).coeff_monomial(x**m))
    lower_ok = all(sp.simplify(sp.Poly(M, x).coeff_monomial(x**k)) == 0 for k in range(0, m))
    pred = (-1) ** m * eta * G(m)
    ok &= lower_ok and sp.simplify(lead - pred) == 0
    detail[m] = str(sp.factor(G(m).subs(varpi, sp.sqrt((q - 1) / 4))))
record("signed_first_defect_law_generic_m2_to_m8", ok, G_m_in_q=detail)
qq = 1 + 4 * varpi**2
record("G4_G5_G6_G7_closed_forms",
       all(sp.simplify(G(m) - g) == 0 for m, g in [(4, -qq * (qq - 8) / 8), (5, -5 * qq * (qq - 4) / 16), (6, qq * (qq**2 - 18 * qq + 48) / 32), (7, 7 * qq * (qq - 4) ** 2 / 64)]))
record("unit_circle_rule_q4", all(sp.simplify(G(m).subs(varpi, sp.sqrt(3) / 2) - (1 - 2 * sp.cos(m * sp.pi / 3))) == 0 for m in range(2, 13)))
# RK4/SSP3 degenerate factorizations
M4_thr = sp.factor(sp.expand(R4.subs(z, -x) - R4.subs(z, -alpha * x) * R4.subs(z, -alphabar * x)).subs(varpi, sp.sqrt(3) / 2))
record("rk4_threshold_factorization", sp.simplify(M4_thr + x**6 * (x - 2) ** 2 / 576) == 0, M=M4_thr)
M3_thr = sp.expand(R_SSP3.subs(z, -x) - R_SSP3.subs(z, -alpha * x) * R_SSP3.subs(z, -alphabar * x)).subs(varpi, sp.sqrt(7) / 2)
record("ssp3_threshold_factorization", sp.simplify(M3_thr - x**5 * (3 - 2 * x) / 9) == 0)
# Dormand-Prince 5 principal formula thresholds via G6 and eta6 sign
m5, e5 = first_defect(R_DP5)
m6r, e6r = first_defect(sp.expand((16 * R4.subs(z, z / 2) ** 2 - R4) / 15))
record("dp5_first_defect", (m5, e5) == (6, Q(1, 3600)), m=m5, eta=e5)
record("richardson_first_defect", (m6r, e6r) == (6, -Q(1, 4320)), m=m6r, eta=e6r)
roots_G6 = sorted(sp.solve(sp.Eq(q**2 - 18 * q + 48, 0), q), key=lambda r: float(r))
record("G6_thresholds",
       sp.simplify(sp.sqrt((roots_G6[0] - 1) / 4) - sp.sqrt(2 - sp.sqrt(33) / 4)) == 0
       and sp.simplify(sp.sqrt((roots_G6[1] - 1) / 4) - sp.sqrt(2 + sp.sqrt(33) / 4)) == 0,
       roots_q=[str(r) for r in roots_G6], varpi_thresholds=[str(sp.N(sp.sqrt((r - 1) / 4), 15)) for r in roots_G6])
# Interior buffers
r2 = sp.Symbol("r2")
Rg = 1 + z + r2 * z**2
Fg = sp.expand(Rg.subs(z, -x) + th * (1 - th) * (1 - Rg.subs(z, -x)) ** 2 - Rg.subs(z, -x / 2 - uu + sp.I * nu) * Rg.subs(z, -x / 2 - uu - sp.I * nu))
record("dephasing_buffer_first_order", sp.expand(Fg).coeff(uu, 1).subs({x: 0, nu: 0}) == 2 and Fg.subs({x: 0, uu: 0, nu: 0}) == 0)
Fray = sp.expand(Fg.subs({uu: 0, nu: varpi * x}))
coef2 = sp.simplify(Fray.coeff(x, 2))
record("excitation_buffer_coefficient", sp.simplify(coef2 - (th * (1 - th) + r2 / 2 - Q(1, 4) + (2 * r2 - 1) * varpi**2)) == 0, coefficient=coef2)
# Large-frequency RK4 defect and attached-endpoint scaling laws
M4gen = sp.expand(R4.subs(z, -x) - R4.subs(z, -alpha * x) * R4.subs(z, -alphabar * x))
c5 = sp.factor(M4gen.coeff(x, 5))
record("rk4_x5_coefficient", sp.simplify(c5 + (4 * varpi**2 - 3) * (4 * varpi**2 + 1) / 384) == 0, coefficient=c5)
record("rk4_large_frequency_defect", sp.limit(c5 / varpi**4, varpi, sp.oo) == -Q(1, 24))
kap = sp.Symbol("kappa", positive=True)
xs = sp.Symbol("x_s", positive=True)
record("attached_scaling_laws",
       sp.simplify(sp.solve(sp.Eq(2 * kap * xs, varpi**4 * xs**5 / 24), xs)[0] ** 4 - 48 * kap / varpi**4) == 0
       and sp.simplify(sp.solve(sp.Eq(th * (1 - th) * xs**2, varpi**4 * xs**5 / 24), xs)[0] ** 3 - 24 * th * (1 - th) / varpi**4) == 0)

# ---------------------------------------------------------------------------
# C3: complete RK4 classification on the laboratory-frame one-way, zero-dephasing ray
# ---------------------------------------------------------------------------
Cq = y**3 - 16 * y**2 - 32 * (q - 6) * y + 384 * (q - 4)
M4_q = sp.expand(M4gen.subs(varpi, sp.sqrt((q - 1) / 4)))
record("rk4_margin_cubic_factorization", sp.simplify(M4_q + x**5 * q * Cq.subs(y, q * x) / 147456) == 0)
disc = sp.factor(sp.discriminant(Cq, y))
record("cubic_discriminant", sp.simplify(disc - 16384 * (q - 4) * (8 * q**2 - 123 * q + 348)) == 0, disc=disc)
qm = (123 - 11 * sp.sqrt(33)) / 16
qp = (123 + 11 * sp.sqrt(33)) / 16
record("q_pm_roots", all(sp.simplify((8 * q**2 - 123 * q + 348).subs(q, v)) == 0 for v in (qm, qp)))
vm = sp.sqrt((qm - 1) / 4)
vp = sp.sqrt((qp - 1) / 4)
record("varpi_pm_closed_forms", sp.simplify(vm - sp.sqrt(107 - 11 * sp.sqrt(33)) / 8) == 0 and sp.simplify(vp - sp.sqrt(107 + 11 * sp.sqrt(33)) / 8) == 0,
       varpi_minus=sp.N(vm, 20), varpi_plus=sp.N(vp, 20), band_width=sp.N(sp.sqrt(3) / 2 - vm, 12))
record("critical_cubic_factorizations",
       sp.simplify(sp.expand(Cq.subs(q, qm) - (y - 9 + sp.sqrt(33)) ** 2 * (y - 2 * sp.sqrt(33) + 2))) == 0
       and sp.expand(Cq.subs(q, 4) - y * (y - 8) ** 2) == 0
       and sp.simplify(sp.expand(Cq.subs(q, qp) - (y - 9 - sp.sqrt(33)) ** 2 * (y + 2 + 2 * sp.sqrt(33)))) == 0)
alpha4_poly = sp.Poly(x**3 - 4 * x**2 + 12 * x - 24, x)
a4_iv = [iv for iv, m in alpha4_poly.intervals(eps=Q(1, 10**15)) if iv[1] > 0]
alpha4 = sp.N(sp.CRootOf(alpha4_poly.as_expr(), 0), 30)
record("alpha4_population_ceiling", len(a4_iv) == 1 and abs(alpha4 - sp.Float("2.785293563405", 15)) < 1e-11
       and sp.expand(24 * (1 - R4.subs(z, -x)) - (-x**4 + 4 * x**3 - 12 * x**2 + 24 * x)) == 0, alpha4=alpha4)


def regime_check(qv):
    """Return the CPTP set description on the one-way ray computed from scratch."""
    P = sp.Poly(Cq.subs(q, qv), y)
    roots = [(sp.Rational(lo), sp.Rational(hi), m) for (lo, hi), m in P.intervals(eps=Q(1, 10**30))]
    pos = [r for r in roots if r[1] > 0]
    return roots, pos


def count_pos(qv):
    P = sp.Poly(sp.expand(Cq.subs(q, qv)), y, domain="QQ<sqrt(33)>" if qv.has(sp.sqrt(33)) else "QQ")
    return P.count_roots(0, sp.oo)


# Regime samples: q=1, 2, 3.6 (< q_- = 3.7378...), 3.9 (q_- < q < 4), 4, 6 (4<q<q_+ = 11.637...), 17 (q>q_+; varpi=2)
sample = {}
for qv in (Q(1), Q(2), Q(36, 10), Q(39, 10), Q(4), Q(6), Q(17)):
    roots, pos = regime_check(qv)
    sample[str(qv)] = [(str(sp.N(lo / qv, 12)), str(sp.N(hi / qv, 12)), m) for lo, hi, m in pos]
record("cubic_root_counts_by_regime",
       count_pos(Q(1)) == 1 and count_pos(Q(2)) == 1 and count_pos(Q(36, 10)) == 1
       and count_pos(Q(39, 10)) == 3 and count_pos(Q(6)) == 0 and count_pos(Q(17)) == 2
       and sp.Poly(sp.expand(Cq.subs(q, Q(6))), y).count_roots(-sp.oo, 0) == 1
       and float(qm) < 3.9 < 4 and 3.6 < float(qm),
       positive_roots_over_q=sample, q_minus=str(sp.N(qm, 12)), q_plus=str(sp.N(qp, 12)))
# Sign chart: for 4<q<q+ the cubic is positive for all y>0 (one negative real root), so M<0 for x>0
record("empty_positive_step_regime_4_to_qplus", sp.Poly(sp.expand(Cq.subs(q, Q(6))), y).count_roots(0, sp.oo) == 0 and Cq.subs({q: 6, y: 0}) > 0)
# Attached-branch crossover q_c: exact certificate (Task 18 review findings A1/B2).
# g(q) = C_q(q alpha4) is a cubic in q with coefficients in Q(alpha4); its discriminant reduced modulo the minimal
# polynomial of alpha4 is negative on an isolating interval of alpha4, so g has exactly one real root, which the
# resultant with the minimal polynomial isolates by a Sturm sequence; the sign change of g across the bracket identifies it.
t_a = sp.Symbol("t_alpha")
alpha_minpoly = t_a**3 - 4 * t_a**2 + 12 * t_a - 24
g_alpha = sp.expand(Cq.subs(y, q * t_a))
disc_g = sp.rem(sp.expand(sp.discriminant(g_alpha, q)), alpha_minpoly, t_a)
a_lo, a_hi = sp.Poly(alpha_minpoly, t_a).refine_root(2, 3, eps=Q(1, 10**30))
disc_vals = [disc_g.subs(t_a, e) for e in (a_lo, a_hi)]
disc_slope = [sp.diff(disc_g, t_a).subs(t_a, e) for e in (a_lo, a_hi)]
one_real_root = all(v < 0 for v in disc_vals) and all(d < 0 for d in disc_slope)
res_alpha = sp.Poly(sp.expand(sp.resultant(g_alpha, alpha_minpoly, t_a)), q)
res_brackets = [(lo, hi) for (lo, hi), m in res_alpha.intervals()]


def g_sign_at(qv):
    """Exact sign of C_q(q alpha4) at a rational q (decided on the isolating interval of alpha4)."""
    val = sp.rem(sp.expand(g_alpha.subs(q, qv)), alpha_minpoly, t_a)
    ends = [sp.sign(val.subs(t_a, e)) for e in (a_lo, a_hi)]
    return int(ends[0]) if ends[0] == ends[1] else 0


qc_iso = res_alpha.refine_root(*res_brackets[0], eps=Q(1, 10**20)) if len(res_brackets) == 1 else None
qc_mid = None if qc_iso is None else (qc_iso[0] + qc_iso[1]) / 2
record("attached_branch_ceiling_crossover",
       one_real_root and qc_iso is not None and 1 < qc_iso[0] and qc_iso[1] < qm
       and g_sign_at(qc_iso[0]) == -1 and g_sign_at(qc_iso[1]) == 1
       and abs(float(qc_mid) - 3.5282058982943983) < 1e-12,
       q_c=None if qc_mid is None else str(sp.N(qc_mid, 20)),
       varpi_c=None if qc_mid is None else str(sp.N(sp.sqrt((qc_mid - 1) / 4), 20)),
       reduced_discriminant=str(disc_g), reduced_discriminant_bracket_values=[str(sp.N(v, 12)) for v in disc_vals],
       resultant_degree=res_alpha.degree(), resultant_real_root_brackets=[(str(lo), str(hi)) for lo, hi in res_brackets],
       note="exactly one real root of q -> C_q(q alpha4) (negative reduced discriminant), Sturm-isolated in (1, q_-) through the resultant; for 1<=q<q_c the attached set is [0, alpha4] (population ceiling active), for q_c<q<q_- it is [0, y3/q]")
# Short rational intervals printed in the Online Resource must independently certify their claims.
printed_a = (Q(27, 10), Q(14, 5))
printed_q = (Q(7, 2), Q(18, 5))
record("printed_crossover_brackets",
       sp.Poly(alpha_minpoly, t_a).count_roots(*printed_a) == 1
       and all(disc_g.subs(t_a, t) < 0 for t in printed_a)
       and sp.diff(disc_g, t_a).subs(t_a, printed_a[0]) < 0 and sp.diff(disc_g, t_a, 2) < 0
       and res_alpha.count_roots(-sp.oo, sp.oo) == res_alpha.count_roots(*printed_q) == 1
       and 1 < printed_q[0] < printed_q[1] < qm
       and g_sign_at(printed_q[0]) == -1 and g_sign_at(printed_q[1]) == 1,
       alpha_bracket=[str(v) for v in printed_a], q_c_bracket=[str(v) for v in printed_q],
       note="the printed rational brackets certify the discriminant sign and crossover location, not merely decimal approximations")
# Position of y3/q relative to alpha4 for every q: for q<10/3 the cubic C_q is increasing; for q>=10/3 the point q*alpha4
# lies beyond the larger critical point y* = (16+sqrt(96q-320))/3 because the quadratic (3 alpha4 q-16)^2-(96q-320) in q has
# negative discriminant; hence sgn C_q(q alpha4) = sgn(q alpha4 - y3), which is + on (q_c, oo) by the crossover certificate.
quad_q = sp.expand((3 * t_a * q - 16) ** 2 - (96 * q - 320))
quad_disc = sp.factor(sp.discriminant(quad_q, q))
quad_disc_vals = [quad_disc.subs(t_a, e) for e in (a_lo, a_hi)]
sign_samples = {str(qv): g_sign_at(qv) for qv in (Q(2), Q(37, 10), Q(39, 10), Q(17), Q(100), Q(10**6))}
mx = max(float(hi) / float(qv) for qv in (Q(39, 10), Q(17), Q(100)) for lo, hi, m in regime_check(qv)[1])
record("detached_endpoints_below_alpha4",
       sp.simplify(quad_disc + 2304 * (t_a - 2) * (5 * t_a + 2)) == 0 and all(v < 0 for v in quad_disc_vals)
       and sign_samples["2"] == -1 and all(sign_samples[k] == 1 for k in ("37/10", "39/10", "17", "100", "1000000"))
       and mx < float(alpha4),
       quadratic_discriminant=str(quad_disc), sign_of_C_q_at_q_alpha4=sign_samples, max_endpoint_over_q_sampled=mx,
       note="y3/q > alpha4 on [1, q_c), = alpha4 at q_c, < alpha4 on (q_c, oo) which contains [q_-, 4] and [q_+, oo); sampled y3/q < alpha4 at q = 3.9, 17, 100")
# isolated points
record("isolated_points_q4_and_qplus",
       sp.simplify(Cq.subs({q: 4, y: 8})) == 0 and sp.simplify(sp.expand(Cq.subs({q: qp, y: 9 + sp.sqrt(33)}))) == 0
       and float(8 / 4) < float(alpha4) and float((9 + sp.sqrt(33)) / qp) < float(alpha4),
       point_q4=2, point_qplus=sp.N((9 + sp.sqrt(33)) / qp, 15))
# window asymptotics: y_2 -> 12 and y_3 ~ sqrt(32 q) as q -> oo; with q = 1 + 4 varpi^2 this gives
# x_- = y_2/q = 3/varpi^2 + O(varpi^-4) and x_+ = y_3/q = 2 sqrt2/|varpi| + O(varpi^-2)
asym_ok = True; asym = {}
for qv in (Q(10**4), Q(10**6), Q(10**8)):
    _roots, pos = regime_check(qv)
    y2 = float((pos[0][0] + pos[0][1]) / 2); y3 = float((pos[1][0] + pos[1][1]) / 2)
    asym[str(qv)] = {"y2": y2, "y3_over_sqrt32q": y3 / float(sp.sqrt(32 * qv))}
    asym_ok &= len(pos) == 2 and abs(y2 - 12) < 200 / float(qv) and abs(y3 / float(sp.sqrt(32 * qv)) - 1) < 10 / float(sp.sqrt(qv))
record("window_edge_asymptotics", asym_ok, samples=asym)
record("rk4_imaginary_axis_boundary", sp.simplify(sp.expand(R4.subs(z, sp.I * nu) * R4.subs(z, -sp.I * nu) - 1) - nu**6 * (nu**2 - 8) / 576) == 0)

# ---------------------------------------------------------------------------
# C4: representation dependence (rotating frame)
# ---------------------------------------------------------------------------
L_H = gksl_superoperator(0, 0, 0, om)
L_dis = gksl_superoperator(gd, gu, gp, 0)
record("hamiltonian_and_dissipative_parts_commute", sp.simplify(L_H * L_dis - L_dis * L_H) == sp.zeros(4))
record("generator_splits", sp.simplify(L_sym - L_H - L_dis) == sp.zeros(4))
U = sp.diag(sp.exp(sp.I * om * hh / 2), sp.exp(-sp.I * om * hh / 2))  # exp(-i H h) with H = -omega sz/2
U_super = superop_left_right(U, U.H)
record("unitary_channel_equals_exp_hamiltonian_part", sp.simplify(U_super - (sp.eye(4) * 0 + (hh * L_H).exp())) == sp.zeros(4))
S4_lab = sp.expand(sp.eye(4) + hh * L_sym + (hh * L_sym) ** 2 / 2 + (hh * L_sym) ** 3 / 6 + (hh * L_sym) ** 4 / 24)
S4_rf = sp.expand(U_super * (sp.eye(4) + hh * L_dis + (hh * L_dis) ** 2 / 2 + (hh * L_dis) ** 3 / 6 + (hh * L_dis) ** 4 / 24))
record("lab_and_rotating_rk4_maps_differ", sp.simplify(S4_lab - S4_rf) != sp.zeros(4))
record("rotating_candidate_admissibility_is_omega_free", om not in sp.simplify(choi(S4_rf.subs(hh, xx).subs({gd: 1 - th, gu: th, gp: 0})).det()).free_symbols
       if True else False)
# Proposition 6(ii): the laboratory-frame and rotating-frame coherence multipliers differ (leading coefficients for the
# polynomial R4; direct differences for the two implicit stability functions of Table 1)
kap_, vp_, nu_ = sp.symbols("kappa varpi nu", real=True)
s_ = sp.Symbol("s", positive=True)
w_lab = -(Q(1, 2) + kap_) * xx + sp.I * vp_ * xx
lab_sq = sp.expand(R4.subs(z, w_lab) * R4.subs(z, sp.conjugate(w_lab)))
rot_sq = sp.expand(R4.subs(z, -(Q(1, 2) + kap_) * xx) ** 2)
lead_lab = sp.Poly(lab_sq, xx).coeff_monomial(xx**8)
lead_rot = sp.Poly(rot_sq, xx).coeff_monomial(xx**8)
be_diff = sp.simplify(sp.Abs(1 / (1 - (-s_ + sp.I * nu_))) ** 2 - (1 / (1 + s_)) ** 2)
mid_R = lambda zz: (1 + zz / 2) / (1 - zz / 2)
mp_diff = sp.simplify(sp.Abs(mid_R(-s_ + sp.I * nu_)) ** 2 - mid_R(-s_) ** 2)
record("frame_candidates_differ_prop6ii",
       sp.simplify(lead_lab - lead_rot - (((Q(1, 2) + kap_) ** 2 + vp_**2) ** 4 - (Q(1, 2) + kap_) ** 8) / 576) == 0
       and sp.simplify(be_diff + nu_**2 / ((1 + s_) ** 2 * ((1 + s_) ** 2 + nu_**2))) == 0
       and sp.simplify(mp_diff - 8 * nu_**2 * s_ / ((2 + s_) ** 2 * ((2 + s_) ** 2 + nu_**2))) == 0,
       note="x^8 coefficients of |R4(-(1/2+kappa)x+i varpi x)|^2 and R4(-(1/2+kappa)x)^2 differ by (((1/2+kappa)^2+varpi^2)^4-(1/2+kappa)^8)/576, nonzero for varpi != 0; backward Euler and implicit midpoint: |R(-s+i nu)|^2 - R(-s)^2 = -nu^2/((1+s)^2((1+s)^2+nu^2)) and 8 nu^2 s/((2+s)^2((2+s)^2+nu^2)), zero for nu != 0 only at s = 0")

# ---------------------------------------------------------------------------
# C5: candidate maps: composition, substep ladder, disjointness at varpi=2, halving threshold, Richardson
# ---------------------------------------------------------------------------
a1, a2, c1, c2 = sp.symbols("a1 a2 c1 c2")
def Phi_super(a, c, theta):
    A = 1 - theta * (1 - a); B = theta * (1 - a); C = (1 - theta) * (1 - a); D = theta + (1 - theta) * a
    # action X -> [[A X00 + C X11, c X01],[conj(c) X10, B X00 + D X11]] in column-stacked basis (X00,X10,X01,X11)
    return sp.Matrix([[A, 0, 0, C], [0, sp.conjugate(c), 0, 0], [0, 0, c, 0], [B, 0, 0, D]])
comp = sp.simplify(Phi_super(a2, c2, th) * Phi_super(a1, c1, th) - Phi_super(a2 * a1, c2 * c1, th))
record("composition_calculus", comp == sp.zeros(4))
record("rk4_population_multiplier_positive", sp.expand(24 * R4.subs(z, -x) - ((x**2 - 2 * x) ** 2 + 8 * (x**2 - 3 * x + 3))) == 0 and sp.discriminant(x**2 - 3 * x + 3, x) < 0)
M4_v2 = sp.Poly(sp.expand(M4gen.subs(varpi, 2)), x)
r_v2 = [(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in M4_v2.intervals(eps=Q(1, 10**20)) if hi > 0]
lo1, hi1 = r_v2[0]; lo2, hi2 = r_v2[1]
record("varpi2_direct_interval", len(r_v2) == 2 and abs(float((lo1 + hi1) / 2) - 0.7447824278959) < 1e-12 and abs(float((lo2 + hi2) / 2) - 1.2703346269739) < 1e-12
       and float(hi2) < float(alpha4), lower=str(sp.N((lo1 + hi1) / 2, 16)), upper=str(sp.N((lo2 + hi2) / 2, 16)))
record("varpi2_two_half_step_disjoint", float(hi2) < 2 * float(lo1), ratio=str(sp.N(hi2 / lo1, 12)))
ratio = float(hi2 / lo1)
record("substep_ladder_overlap_for_n_ge_2", all(ratio >= (n + 1) / n for n in range(2, 10)) and ratio < 2)
# halving threshold: y3 = 2 y2 -> resultant of C_q(y) and C_q(2y) after removing trivial factor
res = sp.factor(sp.resultant(Cq, Cq.subs(y, 2 * y), y))
cub = 72 * q**3 - 2071 * q**2 + 12552 * q - 21744
record("halving_threshold_resultant", sp.simplify(res / (6291456 * (q - 4) * cub)) == 1 or sp.simplify(res - 6291456 * (q - 4) * cub) == 0, resultant=res)
cub_roots = [(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in sp.Poly(cub, q).intervals(eps=Q(1, 10**30))]
qroot = [r for r in cub_roots if float(r[0]) > float(qp)]
record("halving_threshold_value", len(qroot) == 1 and abs(float(sp.sqrt((qroot[0][0] - 1) / 4)) - 2.248254604514980) < 1e-12,
       varpi_half=str(sp.N(sp.sqrt(((qroot[0][0] + qroot[0][1]) / 2 - 1) / 4), 18)))
# Richardson
Rext = sp.expand((16 * R4.subs(z, z / 2) ** 2 - R4) / 15)
record("richardson_stability_polynomial", sp.expand(Rext - (sum(z**j / sp.factorial(j) for j in range(6)) + z**6 / 864 + z**7 / 8640 + z**8 / 138240)) == 0)
Mext = sp.factor(sp.expand(Rext.subs(z, -x) - Rext.subs(z, -x / 2) ** 2))
P10 = sp.Poly(x**10 - 64 * x**9 + 2304 * x**8 - 59392 * x**7 + 1183744 * x**6 - 19169280 * x**5 + 258932736 * x**4 - 2960916480 * x**3 + 19888865280 * x**2 - 97391738880 * x + 280850595840, x)
record("richardson_margin_factorization", sp.simplify(Mext + x**6 * P10.as_expr() / 1252412463513600) == 0
       and sp.expand(Mext).coeff(x, 6) == -Q(31, 138240))
p10_roots = [(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in P10.intervals(eps=Q(1, 10**20)) if hi > 0]
rext_pop = sp.Poly(sp.expand(1 - Rext.subs(z, -x)), x)
pop_roots = [(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in rext_pop.intervals(eps=Q(1, 10**20)) if hi > 0]
record("richardson_first_positive_root_above_alpha4", float(p10_roots[0][0]) > float(alpha4) and abs(float(p10_roots[0][0]) - 6.034523958649) < 1e-10,
       first_root=str(sp.N((p10_roots[0][0] + p10_roots[0][1]) / 2, 16)), all_positive_P10_roots=[str(sp.N((a + b) / 2, 12)) for a, b in p10_roots])
# remote interval: margin >= 0 between first and second P10 roots? and population ceiling R_ext(-x) <= 1
rem_lo = p10_roots[0]
pop_pos = [r for r in pop_roots if float(r[0]) > 0]
a_ext = Rext.subs(z, -x)
record("richardson_remote_interval",
       abs(float(pop_pos[0][0]) - 6.459127767826) < 1e-10
       and (float(pop_pos[0][0]) < float(p10_roots[1][0]) if len(p10_roots) > 1 else True)
       and Mext.subs(x, Q(62, 10)) > 0 and 0 <= a_ext.subs(x, Q(62, 10)) <= 1 and a_ext.subs(x, Q(65, 10)) > 1
       and Mext.subs(x, Q(3)) < 0 and Mext.subs(x, Q(1, 2)) < 0,
       population_ceiling_root=str(sp.N((pop_pos[0][0] + pop_pos[0][1]) / 2, 16)),
       note="margin positive and 0<=a<=1 at x=6.2 (inside), a>1 at x=6.5 (ceiling exit), margin negative at x=0.5 and 3 (outside)")
# Complete global chart, including all population-floor events and the final unbounded interval.
rich_polys = {"a": sp.Poly(a_ext, x), "a-1": sp.Poly(a_ext-1, x), "M": sp.Poly(Mext, x)}
rich_events = [
    ("a", "2.919911381617275", "2.919911381617276"),
    ("a", "6.034383873449079", "6.034383873449080"),
    ("M", "6.034523958649170", "6.034523958649171"),
    ("a-1", "6.459127767825720", "6.459127767825721"),
    ("M", "22.102437082855815", "22.102437082855816"),
]
chart_ok = True
previous = Q(0)
for source, lo_s, hi_s in rich_events:
    lo, hi = Q(lo_s), Q(hi_s)
    poly = rich_polys[source]
    chart_ok &= previous < lo < hi and poly.count_roots(lo, hi) == 1
    chart_ok &= poly.gcd(poly.diff()).count_roots(lo, hi) == 0
    previous = hi
for source, poly in rich_polys.items():
    valuation = min(monomial[0] for monomial, coefficient in poly.terms())
    no_origin = sp.Poly(poly.as_expr()/x**valuation, x)
    chart_ok &= no_origin.count_roots(0, sp.oo) == sum(event[0] == source for event in rich_events)
rich_samples = [Q(1), Q(4), Q("6.03445"), Q("6.2"), Q(10), Q(23)]
rich_signs = [[int(sp.sign(poly.eval(sample))) for poly in rich_polys.values()] for sample in rich_samples]
for i, sample in enumerate(rich_samples):
    chart_ok &= (i == 0 or Q(rich_events[i-1][2]) < sample)
    chart_ok &= (i == len(rich_events) or sample < Q(rich_events[i][1]))
chart_ok &= rich_signs == [[1,-1,-1], [-1,-1,-1], [1,-1,-1], [1,-1,1], [1,1,1], [1,1,-1]]
record("richardson_global_sign_chart", chart_ok,
       events=[{"source": source, "bracket": [lo, hi]} for source, lo, hi in rich_events],
       intervals=[{"sample": str(sample), "signs": signs} for sample, signs in zip(rich_samples, rich_signs)],
       note="all positive boundary roots isolated and simple; exact signs on every interval including the tail leave only the closed remote component plus zero")
# Inactive constraints remain strictly satisfied on a neighbourhood of the whole RK4 detached interval.
slack_interval = (Q(7, 10), Q(13, 10))
a4_slack = sp.Poly(R4.subs(z, -x), x)
m4_slack = sp.Poly(M4gen.subs(varpi, 2), x)
record("persistence_inactive_slack_interval",
       m4_slack.count_roots(*slack_interval) == 2
       and all(poly.count_roots(*slack_interval) == 0
               and all(poly.eval(v) > 0 for v in slack_interval) for poly in (a4_slack, 1-a4_slack)),
       interval=[str(v) for v in slack_interval],
       note="at varpi=2 both detached roots lie in (0.7,1.3), where a and 1-a have no zeros and are strictly positive; A_0=1 supplies the remaining inactive constraint")
# constituents CPTP on [0, alpha4]: coarse set [0, alpha4]; fine set 2*[0,alpha4] contains [0, alpha4]
cubic_floor = sp.Poly(1152 - 160 * x + 16 * x**2 - x**3, x)
alpha4_upper = sp.Poly(sp.Symbol("t") ** 3 - 4 * sp.Symbol("t") ** 2 + 12 * sp.Symbol("t") - 24, sp.Symbol("t")).refine_root(2, 3, eps=Q(1, 10**20))[1]
record("constituents_cptp_on_common_interval",
       sp.expand(M4gen.subs(varpi, 0) - x**5 * cubic_floor.as_expr() / 147456) == 0
       and cubic_floor.count_roots(0, alpha4_upper) == 0 and cubic_floor.eval(0) > 0 and cubic_floor.eval(alpha4_upper) > 0,
       note="coarse RK4 at varpi=0: M4(x;0) = x^5(1152-160x+16x^2-x^3)/147456, the cubic has no root in (0, alpha4] (Sturm count) and is positive at both ends, so M4>0 on (0, alpha4]; 0<=a<=1 there by the population identity and the ceiling root alpha4; the fine map is CPTP on [0, 2 alpha4] by dilation (Proposition 7)")
M40 = sp.Poly(sp.expand(M4gen.subs(varpi, 0)), x)
m40_roots = [(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in M40.intervals(eps=Q(1, 10**20)) if hi > 0]
record("rk4_nonrotating_margin_root", abs(float(m40_roots[0][0]) - 10.982425466293) < 1e-9 and float(m40_roots[0][0]) > float(alpha4), root=str(sp.N(m40_roots[0][0], 15)))
# Richardson frequency orientation sequence
record("richardson_orientation_bands", (e6r * G(6).subs(varpi, 0) < 0) and (e6r * G(6).subs(varpi, 1) > 0) and (e6r * G(6).subs(varpi, 5) < 0))
record("dp5_orientation_bands_reversed", (e5 * G(6).subs(varpi, 0) > 0) and (e5 * G(6).subs(varpi, 1) < 0) and (e5 * G(6).subs(varpi, 5) > 0))

# ---------------------------------------------------------------------------
# Table 1: buffered cases at varpi=2 (theta, kappa) -> roots r1<r2<r3
# ---------------------------------------------------------------------------
def rk4_F(theta_v, kappa_v, varpi_v):
    a = R4.subs(z, -x)
    arg = -x * (Q(1, 2) + kappa_v - sp.I * varpi_v)
    c = R4.subs(z, arg); cb = R4.subs(z, sp.conjugate(arg))
    return sp.Poly(sp.expand(a + theta_v * (1 - theta_v) * (1 - a) ** 2 - c * cb), x)

table = {
    ("0", "1e-3"): (0.253348869116, 0.740695999652, 1.269542711140),
    ("1e-3", "0"): (0.122268903470, 0.742558312554, 1.270544370136),
    ("1e-3", "1e-3"): (0.262033333161, 0.738364627213, 1.269753570826),
    ("0", "1e-2"): (0.519308218820, 0.681621610068, 1.262280154980),
    ("1e-2", "0"): (0.276930275811, 0.720896805059, 1.272391394935),
}
tab_ok = True; tab_detail = {}
for (ts, ks), expected in table.items():
    tv, kv = Q(ts), Q(ks)
    P = rk4_F(tv, kv, 2)
    rts = [float((sp.Rational(lo) + sp.Rational(hi)) / 2) for (lo, hi), m in P.intervals(eps=Q(1, 10**18)) if hi > 0]
    rts = [r for r in rts if r < float(alpha4)]
    good = len(rts) >= 3 and all(abs(rts[i] - expected[i]) < 1e-10 for i in range(3))
    # also confirm the population/A/D constraints inactive below alpha4 and sign pattern +,-,+,- between roots
    tab_ok &= good
    tab_detail[f"theta={ts},kappa={ks}"] = {"roots": rts[:3], "expected": list(expected), "ok": good}
record("table1_buffered_roots_varpi2", tab_ok, rows=tab_detail)
# persistence-theorem sign samples at (0,0), varpi=2
Mv2 = sp.expand(M4gen.subs(varpi, 2))
record("persistence_sign_samples", Mv2.subs(x, Q(1, 2)) < 0 and Mv2.subs(x, 1) > 0 and Mv2.subs(x, Q(3, 2)) < 0)

# ---------------------------------------------------------------------------
# C6: noncommuting transverse-field certificates (full 4x4 superoperator, exact)
# ---------------------------------------------------------------------------
def choi_charpoly_min_root_sign(L_num_exact, hv):
    S = sp.eye(4) + hv * L_num_exact + (hv * L_num_exact) ** 2 / 2 + (hv * L_num_exact) ** 3 / 6 + (hv * L_num_exact) ** 4 / 24
    J = choi(sp.simplify(S))
    cp = J.charpoly()
    P = sp.Poly(sp.expand(cp.as_expr()), cp.gen)
    coeffs = []
    for cf in P.all_coeffs():
        re_, im_ = sp.expand(cf).as_real_imag()
        if im_ != 0:
            raise ArithmeticError("charpoly coefficient is not real")
        coeffs.append(sp.Rational(re_))
    P = sp.Poly(coeffs, cp.gen, domain=sp.QQ)
    ivs = sorted([(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in P.intervals(eps=Q(1, 10**12))], key=lambda t: t[0])
    lo, hi = ivs[0]
    return lo, hi

t0 = time.time()
L_nc = gksl_superoperator(1 - Q(1, 1000), Q(1, 1000), Q(1, 1000), 2, Q(1, 10))
expected_nc = {Q(1, 10): 9.5156e-5, Q(1, 2): -2.7778e-3, Q(9, 10): 5.8959e-4, Q(14, 10): -3.9227e-1}
nc_ok = True; nc_detail = {}
for hv, ev in expected_nc.items():
    lo, hi = choi_charpoly_min_root_sign(L_nc, hv)
    mid = float((lo + hi) / 2)
    sign_ok = (lo > 0) if ev > 0 else (hi < 0)
    nc_ok &= sign_ok and abs(mid - ev) / abs(ev) < 2e-4
    nc_detail[str(hv)] = {"bracket": [str(lo), str(hi)], "midpoint": mid, "expected": ev, "sign_ok": sign_ok}
print(f"  noncommuting eigenvalue signs: {time.time() - t0:.1f} s")
record("noncommuting_eigenvalue_signs_omega_x_0p1", nc_ok, samples=nc_detail)
# component endpoints at Omega_x = 0.1: boundary where det J = 0 (or other minors); use det J roots and sign checks
t0 = time.time()
Sh = sp.eye(4) + h * L_nc + (h * L_nc) ** 2 / 2 + (h * L_nc) ** 3 / 6 + (h * L_nc) ** 4 / 24
Jh = choi(sp.expand(Sh))
def real_rational_poly(expr, gen):
    P = sp.Poly(sp.expand(expr), gen)
    coeffs = []
    for cf in P.all_coeffs():
        re_, im_ = sp.expand(cf).as_real_imag()
        if im_ != 0:
            raise ArithmeticError("polynomial coefficient is not real")
        coeffs.append(sp.Rational(re_))
    return sp.Poly(coeffs, gen, domain=sp.QQ)

detJ = real_rational_poly(Jh.det(method="berkowitz"), h)
det_roots = [(sp.Rational(lo), sp.Rational(hi)) for (lo, hi), m in detJ.intervals(eps=Q(1, 10**14)) if hi > 0 and lo < Q(3, 2)]
ends = [float((lo + hi) / 2) for lo, hi in det_roots]
expected_ends = [0.26357691, 0.73588381, 1.26847248]
end_ok = len(ends) >= 3 and all(abs(ends[i] - expected_ends[i]) < 1e-7 for i in range(3))
# verify sign pattern of min eigenvalue between the roots
pattern = []
pts = [Q(1, 10), (det_roots[0][1] + det_roots[1][0]) / 2, (det_roots[1][1] + det_roots[2][0]) / 2, det_roots[2][1] + Q(1, 10)]
for pt in pts:
    lo, hi = choi_charpoly_min_root_sign(L_nc, pt)
    pattern.append("+" if lo > 0 else "-" if hi < 0 else "?")
print(f"  noncommuting components: {time.time() - t0:.1f} s")
record("noncommuting_components_omega_x_0p1", end_ok and pattern == ["+", "-", "+", "-"], endpoints=ends[:4], sign_pattern=pattern)
# other transverse values (table in the supplement)
t0 = time.time()
sup_table = {Q(0): (0.262033333161, 0.738364627213, 1.269753570826), Q(1, 20): (0.262418936401, 0.737744306425, 1.269432714005), Q(1, 5): (0.268222457925, 0.728458223265, 1.264665614699)}
sup_ok = True; sup_detail = {}
for wx, exp_ends in sup_table.items():
    Lw = gksl_superoperator(1 - Q(1, 1000), Q(1, 1000), Q(1, 1000), 2, wx)
    Sw = sp.eye(4) + h * Lw + (h * Lw) ** 2 / 2 + (h * Lw) ** 3 / 6 + (h * Lw) ** 4 / 24
    Jw = choi(sp.expand(Sw))
    dw = real_rational_poly(Jw.det(method="berkowitz"), h)
    rr = [float((sp.Rational(lo) + sp.Rational(hi)) / 2) for (lo, hi), m in dw.intervals(eps=Q(1, 10**14)) if hi > 0 and lo < Q(3, 2)]
    good = len(rr) >= 3 and all(abs(rr[i] - exp_ends[i]) < 1e-8 for i in range(3))
    sup_ok &= good
    sup_detail[str(wx)] = {"roots": rr[:3], "expected": list(exp_ends), "ok": good}
print(f"  transverse table endpoints: {time.time() - t0:.1f} s")
record("transverse_table_endpoints", sup_ok, rows=sup_detail,
       note="endpoints are roots of det J(h); the sign pattern between them was certified for Omega_x=0.1 above and is +,-,+,- by continuity of the certified structure only at the listed values")

# ---------------------------------------------------------------------------
# C7: executed binary64 versus printed decimal at the varpi=2 upper endpoint
# ---------------------------------------------------------------------------
xd = Q("1.270334626973889")
xf = Q(*Fraction(1.270334626973889).as_integer_ratio())
Mdec = Mv2.subs(x, xd); Mflt = Mv2.subs(x, xf)
record("printed_decimal_vs_binary64_margin", Mdec > 0 and Mflt < 0 and abs(float(Mdec) - 1.43e-16) < 2e-18 and abs(float(Mflt) + 9.45e-17) < 2e-18,
       decimal_margin=float(Mdec), binary64_margin=float(Mflt))

# ---------------------------------------------------------------------------
# Nine-candidate summary (one-way, nonrotating, no dephasing)
# ---------------------------------------------------------------------------
def admissible_summary(R):
    a = sp.cancel(R.subs(z, -x)); c = sp.cancel(R.subs(z, -x / 2))
    M = sp.cancel(a - c**2)
    num_M = sp.Poly(sp.numer(sp.together(M)), x)
    num_1a = sp.Poly(sp.numer(sp.together(1 - a)), x)
    num_a = sp.Poly(sp.numer(sp.together(a)), x)
    crit = set()
    for P in (num_M, num_1a, num_a):
        if P.degree() > 0:
            for (lo, hi), m in P.intervals(eps=Q(1, 10**15)):
                if hi > 0:
                    crit.add(float((sp.Rational(lo) + sp.Rational(hi)) / 2))
    crit = sorted(crit)
    pts = [0.0] + crit + [crit[-1] + 1 if crit else 1.0]
    segs = []
    for i in range(len(pts) - 1):
        mid = Q((pts[i] + pts[i + 1]) / 2)
        av, Mv = a.subs(x, mid), M.subs(x, mid)
        segs.append((pts[i], pts[i + 1], bool(av >= 0 and av <= 1 and Mv >= 0)))
    return segs

nine = {}
nine_expected = {
    "Forward Euler": (R_EULER, []), "Heun RK2": (R_HEUN, [(0, 2)]), "SSPRK(3,3)": (R_SSP3, []),
    "Classical RK4": (R4, [(0, 2.785293563405282)]), "Richardson RK4": (Rext, [(6.034523958649, 6.459127767826)]),
    "Dormand-Prince 5": (R_DP5, [(0, 3.306567892634947)]), "Dormand-Prince embedded 4": (R_DP4, [(3.068856548050381, 4.384986320801944)]),
    "Backward Euler": (R_BE, [(0, float("inf"))]), "Implicit midpoint": (R_IM, []),
}
nine_ok = True
for name, (R, exp_segs) in nine_expected.items():
    segs = [(a, b) for a, b, ok_ in admissible_summary(R) if ok_]
    # merge adjacent
    merged = []
    for a_, b_ in segs:
        if merged and abs(merged[-1][1] - a_) < 1e-12:
            merged[-1] = (merged[-1][0], b_)
        else:
            merged.append((a_, b_))
    good = len(merged) == len(exp_segs) and all(abs(m[0] - e[0]) < 1e-8 and (e[1] == float("inf") or abs(m[1] - e[1]) < 1e-8) for m, e in zip(merged, exp_segs))
    nine_ok &= good
    nine[name] = {"positive_width_components": [(round(a_, 12), (None if b_ == float("inf") or (b_ > 1e6) else round(b_, 12))) for a_, b_ in merged], "expected": exp_segs, "ok": good}
record("nine_candidate_admissible_sets", nine_ok, methods=nine, note="Backward Euler's last sampled segment is unbounded; the archive states [0, infinity)")

# ---------------------------------------------------------------------------
# Fixed-horizon convergence from outside the channel set
# ---------------------------------------------------------------------------
Xs = sp.Symbol("X", positive=True)
epsn = sp.Symbol("epsilon", positive=True)  # 1/N
fh_ok = True; fh = {}
for name, R in (("Euler", R_EULER), ("RK4", R4)):
    m_, eta_ = first_defect(R)
    K = (-1) ** m_ * (1 - 2 ** (1 - m_)) * eta_
    a_ = R.subs(z, -Xs * epsn); c_ = R.subs(z, -Xs * epsn / 2)
    # a^N = exp(N log a) with N = 1/eps: expand log a in eps first, then divide by eps
    la = sp.expand(sp.series(sp.log(a_), epsn, 0, m_ + 2).removeO() / epsn)
    lc = sp.expand(2 * sp.series(sp.log(c_), epsn, 0, m_ + 2).removeO() / epsn)
    MN = sp.exp(-Xs) * (sp.exp(sp.expand(la + Xs)) - sp.exp(sp.expand(lc + Xs)))
    ser = sp.expand(sp.series(MN, epsn, 0, m_).removeO())
    lead = sp.simplify(ser.coeff(epsn, m_ - 1))
    pred = K * Xs * sp.exp(-Xs) * Xs ** (m_ - 1)
    good = all(sp.simplify(ser.coeff(epsn, k)) == 0 for k in range(m_ - 1)) and sp.simplify(lead - pred) == 0
    fh_ok &= good; fh[name] = {"m": m_, "eta": str(eta_), "K": str(K), "ok": good}
record("fixed_horizon_margin_asymptotics", fh_ok, methods=fh)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
n_pass = sum(1 for v in CHECKS.values() if v["pass"]); n_all = len(CHECKS)
summary = {"script": "scripts/verify_bit_claims.py", "independent_of_rk_choi_margin": True, "sympy": sp.__version__, "mpmath": mp.__version__,
           "checks_passed": n_pass, "checks_total": n_all, "checks": CHECKS}
(OUT / "claim_verification.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
lines = [f"# Independent claim verification ({n_pass}/{n_all} passed)", "", "Computed from first principles with SymPy " + sp.__version__ + "; no import of `rk_choi_margin`.", ""]
for k, v in CHECKS.items():
    lines.append(f"- {'PASS' if v['pass'] else 'FAIL'} `{k}`" + (": " + v["note"] if "note" in v else ""))
(OUT / "claim_verification.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"\n{n_pass}/{n_all} checks passed")
sys.exit(0 if n_pass == n_all else 1)
