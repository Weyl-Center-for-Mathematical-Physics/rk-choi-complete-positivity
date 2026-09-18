from __future__ import annotations

"""Richardson extrapolation of Runge--Kutta stability functions (BIT revision, extension E1).

For a stability function ``R`` of linear order ``p`` with
``R(z) - exp(z) = eta_{p+1} z^{p+1} + eta_{p+2} z^{p+2} + O(z^{p+3})``, the n-fold Richardson
candidate ``E_n(z) = (n^p R(z/n)^n - R(z)) / (n^p - 1)`` has exponential defect

    E_n(z) - exp(z) = (n-1)(eta_{p+1} - eta_{p+2}) / (n (n^p - 1)) * z^{p+2} + O(z^{p+3}).

Combined with the signed first-defect law, this fixes the local complete-positivity orientation of the
extrapolated candidate on the one-way, zero-dephasing ray for every real-coefficient RK method.
"""

import sympy as sp

from .methods import X


def richardson_extrapolate(R: sp.Expr, p: int, n: int = 2, *, x: sp.Symbol = X) -> sp.Expr:
    """Return the n-fold Richardson extrapolate of the stability function ``R`` of linear order ``p``."""
    p = int(p)
    n = int(n)
    if n < 2:
        raise ValueError("n must be >= 2")
    if p < 1:
        raise ValueError("p must be >= 1")
    weight = sp.Integer(n) ** p
    return sp.cancel((weight * R.subs(x, x / n) ** n - R) / (weight - 1))


def predicted_extrapolation_defect(eta_p1: sp.Expr, eta_p2: sp.Expr, p: int, n: int = 2) -> sp.Expr:
    """Coefficient of z^{p+2} in E_n(z) - exp(z) predicted by the corollary."""
    p = int(p)
    n = int(n)
    if n < 2:
        raise ValueError("n must be >= 2")
    return sp.together(sp.Integer(n - 1) * (eta_p1 - eta_p2) / (sp.Integer(n) * (sp.Integer(n) ** p - 1)))


def frequency_factor(m: int, varpi: sp.Expr) -> sp.Expr:
    """G_m(varpi) = 1 - 2 Re (1/2 - i varpi)^m for real varpi."""
    alpha = sp.Rational(1, 2) - sp.I * varpi
    alphabar = sp.Rational(1, 2) + sp.I * varpi
    return sp.simplify(sp.expand(1 - (alpha ** int(m) + alphabar ** int(m))))


def local_orientation_coefficient(m: int, eta: sp.Expr, varpi: sp.Expr) -> sp.Expr:
    """Leading coefficient (-1)^m eta G_m(varpi) of the one-way Choi margin for a method with first
    exponential defect eta z^m.  Positive means locally completely positive for small steps."""
    return sp.simplify((-1) ** int(m) * eta * frequency_factor(m, varpi))
