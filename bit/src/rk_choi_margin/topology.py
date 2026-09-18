from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import mpmath as mp
import numpy as np
import sympy as sp

from .channels import choi_matrix
from .certification import as_rational
from .methods import all_methods

X = sp.Symbol("x", nonnegative=True, real=True)
Y = sp.Symbol("y", real=True)
QSYM = sp.Symbol("q", real=True)
THETA = sp.Symbol("theta", nonnegative=True, real=True)
KAPPA = sp.Symbol("kappa", nonnegative=True, real=True)
VARPI = sp.Symbol("varpi", real=True)
S = sp.Symbol("s")


@dataclass(frozen=True)
class RootBracket:
    lower: sp.Rational
    upper: sp.Rational

    def midpoint(self) -> sp.Rational:
        return (self.lower + self.upper) / 2


def frequency_sigma(m: int, q: sp.Expr = QSYM) -> sp.Expr:
    """Return sigma_m=2^m(alpha^m+conj(alpha)^m), q=1+4 varpi^2."""
    if m < 0:
        raise ValueError("m must be nonnegative")
    if m in (0, 1):
        return sp.Integer(2)
    s0 = sp.Integer(2)
    s1 = sp.Integer(2)
    for _ in range(2, m + 1):
        s0, s1 = s1, sp.expand(2 * s1 - q * s0)
    return sp.factor(s1)


def frequency_factor_q(m: int, q: sp.Expr = QSYM) -> sp.Expr:
    """Return G_m as a polynomial in q=1+4 varpi^2."""
    return sp.factor(1 - frequency_sigma(m, q) / 2**m)


def frequency_factor(m: int, varpi: sp.Expr = VARPI) -> sp.Expr:
    return sp.factor(frequency_factor_q(m, 1 + 4 * varpi**2))


def degenerate_next_coefficient(
    m: int,
    eta_m: sp.Expr,
    eta_next: sp.Expr,
    varpi: sp.Expr,
) -> sp.Expr:
    """x^(m+1) coefficient when G_m(varpi)=0."""
    alpha = sp.Rational(1, 2) - sp.I * varpi
    term = (
        (-1) ** (m + 1) * eta_next * frequency_factor(m + 1, varpi)
        + 2 * (-1) ** m * eta_m * sp.re(sp.conjugate(alpha) * alpha**m)
    )
    return sp.factor(sp.expand_complex(term))


def rk4_rotating_cubic(q: sp.Expr = QSYM, y: sp.Expr = Y) -> sp.Expr:
    return sp.expand(y**3 - 16 * y**2 - 32 * (q - 6) * y + 384 * (q - 4))


def rk4_rotating_margin(x: sp.Expr = X, varpi: sp.Expr = VARPI) -> sp.Expr:
    q = 1 + 4 * varpi**2
    return sp.factor(-x**5 * q * rk4_rotating_cubic(q, q * x) / sp.Integer(147456))


def rk4_general_margin(
    x: sp.Expr = X,
    theta: sp.Expr = THETA,
    kappa: sp.Expr = KAPPA,
    varpi: sp.Expr = VARPI,
) -> sp.Expr:
    """Exact phase-covariant Choi-block margin for classical RK4 on a ray."""
    R = all_methods()["rk4"].stability_function(x=S)
    a = sp.expand(R.subs(S, -x))
    z = -x * (sp.Rational(1, 2) + kappa - sp.I * varpi)
    c = sp.expand(R.subs(S, z))
    cbar = sp.expand(R.subs(S, sp.conjugate(z)))
    return sp.factor(sp.expand(a + theta * (1 - theta) * (1 - a) ** 2 - c * cbar))


def rk4_population_multiplier(x: sp.Expr = X) -> sp.Expr:
    R = all_methods()["rk4"].stability_function(x=S)
    return sp.factor(R.subs(S, -x))


def exact_phase_covariant_multipliers(x: float, kappa: float, varpi: float) -> tuple[complex, complex]:
    return float(mp.e ** (-x)), complex(mp.e ** (-x / 2 - kappa * x + 1j * varpi * x))


def rk4_multipliers(x: float, kappa: float, varpi: float) -> tuple[complex, complex]:
    R = all_methods()["rk4"].stability_function(x=S)
    f = sp.lambdify(S, R, "mpmath")
    return complex(f(-x)), complex(f(-x / 2 - kappa * x + 1j * varpi * x))


def normalized_choi_defect_and_error(x: float, varpi: float, kappa: float = 0.0) -> tuple[mp.mpf, mp.mpf]:
    """Negative mass and trace distance from the exact normalized Choi state."""
    mp.mp.dps = 100
    def _mp_exact(value):
        rational = as_rational(value)
        return mp.mpf(int(rational.p)) / mp.mpf(int(rational.q))
    xx = _mp_exact(x)
    vv = _mp_exact(varpi)
    kk = _mp_exact(kappa)
    def r4(z):
        return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24
    a = mp.mpf(r4(-xx))
    c = mp.mpc(r4(-xx / 2 - kk * xx + 1j * vv * xx))
    ae = mp.e ** (-xx)
    ce = mp.e ** (-xx / 2 - kk * xx + 1j * vv * xx)

    def eigvals_hermitian_4(a0: mp.mpf, c0: mp.mpc) -> list[mp.mpf]:
        # normalized Choi eigenvalues: 0, (1-a)/2, and block eigenvalues/2
        disc = mp.sqrt((1 - a0) ** 2 + 4 * abs(c0) ** 2)
        return [mp.mpf(0), (1 - a0) / 2, (1 + a0 - disc) / 4, (1 + a0 + disc) / 4]

    eig = eigvals_hermitian_4(a, c)
    defect = sum(-v for v in eig if v < 0)

    # Difference is supported on the same sparse Hermitian pattern. Its eigenvalues
    # consist of 0, -(a-ae)/2, and the two eigenvalues of [[0, dc],[dc*, da]]/2.
    da = a - ae
    dc = c - ce
    block_disc = mp.sqrt(da**2 + 4 * abs(dc) ** 2)
    diff_eigs = [mp.mpf(0), -da / 2, (da - block_disc) / 4, (da + block_disc) / 4]
    trace_distance = mp.mpf("0.5") * sum(abs(v) for v in diff_eigs)
    return defect, trace_distance


def count_positive_roots_exact(poly: sp.Poly) -> int:
    sqf = poly.sqf_part()
    zero = 1 if sqf.eval(0) == 0 else 0
    return int(sqf.count_roots(0, sp.oo)) - zero


def rational_bracket(value: float | str, digits: int = 12, width_units: int = 2) -> RootBracket:
    """Create a decimal rational bracket centered on a numerical value."""
    center = sp.Rational(str(value))
    unit = sp.Rational(1, 10**digits)
    lower = sp.floor(center / unit) * unit - width_units * unit
    upper = sp.ceiling(center / unit) * unit + width_units * unit
    return RootBracket(sp.Rational(lower), sp.Rational(upper))


def isolate_positive_roots(poly: sp.Poly, digits: int = 60) -> list[tuple[sp.Rational, sp.Rational, sp.Expr]]:
    """Return rational isolating intervals and high-precision approximations."""
    roots = []
    for interval in poly.intervals(eps=sp.Rational(1, 10**18)):
        (lo, hi), multiplicity = interval
        if hi <= 0:
            continue
        lo = max(lo, sp.Rational(0))
        if lo == hi == 0:
            continue
        if lo == 0 and hi > 0:
            # skip a bracket that contains the exact zero root
            continue
        approx = sp.N((lo + hi) / 2, digits)
        for _ in range(multiplicity):
            roots.append((sp.Rational(lo), sp.Rational(hi), approx))
    roots.sort(key=lambda row: float(row[2]))
    return roots
