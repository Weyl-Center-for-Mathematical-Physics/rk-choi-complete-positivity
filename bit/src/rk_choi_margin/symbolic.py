from __future__ import annotations

from dataclasses import dataclass
import math

import sympy as sp

X = sp.Symbol("x")
Z = sp.Symbol("z", nonnegative=True, real=True)
U = sp.Symbol("u", nonnegative=True, real=True)
NU = sp.Symbol("nu", real=True)
THETA = sp.Symbol("theta", real=True)


@dataclass(frozen=True)
class Defect:
    index: int
    coefficient: sp.Expr


def cp_margin(R: sp.Expr, *, r: sp.Expr = sp.Rational(1, 2), x: sp.Symbol = X, z: sp.Symbol = Z) -> sp.Expr:
    """Return M_{R,r}(z)=R(-z)-R(-r z)^2 on a fixed-rate ray.

    This one-parameter form is recovered from ``two_rate_cp_margin`` by setting
    ``u=(r-1/2)z`` whenever the relaxation rate is nonzero.
    """
    return sp.factor(sp.cancel(R.subs(x, -z) - R.subs(x, -r * z) ** 2))


def two_rate_cp_margin(
    R: sp.Expr,
    *,
    relaxation_step: sp.Expr = Z,
    dephasing_step: sp.Expr = U,
    x: sp.Symbol = X,
) -> sp.Expr:
    """Return M_R(x,u)=R(-x)-R(-x/2-u)^2.

    ``relaxation_step`` is the dimensionless population-relaxation step
    gamma*h and ``dephasing_step`` is the independent pure-dephasing step
    gamma_phi*h.  This coordinate system remains regular at pure dephasing.
    """
    a = R.subs(x, -relaxation_step)
    c = R.subs(x, -sp.Rational(1, 2) * relaxation_step - dephasing_step)
    return sp.factor(sp.cancel(a - c**2))


def rotating_cp_margin(
    R: sp.Expr,
    *,
    r: sp.Expr = sp.Rational(1, 2),
    nu: sp.Symbol | sp.Expr,
    x: sp.Symbol = X,
    z: sp.Symbol = Z,
) -> sp.Expr:
    """Return R(-z)-|R((-r+i nu)z)|^2 for real-coefficient R."""
    q = R.subs(x, (-r + sp.I * nu) * z)
    q_conj = R.subs(x, (-r - sp.I * nu) * z)
    return sp.factor(sp.cancel(R.subs(x, -z) - q * q_conj))



def thermal_cp_margin(
    R: sp.Expr,
    *,
    theta: sp.Expr = THETA,
    relaxation_step: sp.Expr = Z,
    dephasing_step: sp.Expr = U,
    frequency_step: sp.Expr = NU,
    x: sp.Symbol = X,
) -> sp.Expr:
    """Return a+theta(1-theta)(1-a)^2-|c|^2."""
    a = R.subs(x, -relaxation_step)
    arg = -sp.Rational(1, 2) * relaxation_step - dephasing_step + sp.I * frequency_step
    c = R.subs(x, arg)
    c_conj = R.subs(x, sp.conjugate(arg))
    return sp.factor(sp.cancel(a + theta * (1 - theta) * (1 - a) ** 2 - c * c_conj))


def rotating_boundary_factor(m: int, varpi: sp.Expr) -> sp.Expr:
    """Return G_m(varpi)=1-2 Re(1/2-i varpi)^m for real varpi."""
    return sp.factor(sp.expand_complex(1 - 2 * sp.re((sp.Rational(1, 2) - sp.I * varpi) ** m)))

def first_exponential_defect(
    R: sp.Expr,
    *,
    x: sp.Symbol = X,
    max_index: int = 30,
) -> Defect:
    """Find the first nonzero coefficient in R(x)-exp(x)."""
    series = sp.series(R - sp.exp(x), x, 0, max_index + 1).removeO().expand()
    for index in range(max_index + 1):
        coefficient = sp.simplify(series.coeff(x, index))
        if coefficient != 0:
            return Defect(index=index, coefficient=coefficient)
    raise ValueError(f"No nonzero defect found through order {max_index}")


def predicted_boundary_coefficient(defect: Defect) -> sp.Expr:
    """Leading coefficient of M_R at r=1/2."""
    m = defect.index
    return sp.simplify((-1) ** m * (1 - sp.Rational(1, 2) ** (m - 1)) * defect.coefficient)


def zero_multiplicity(expr: sp.Expr, *, z: sp.Symbol = Z) -> int:
    """Multiplicity of z=0 in a rational expression's numerator."""
    numerator = sp.Poly(sp.together(expr).as_numer_denom()[0], z)
    powers = [monomial[0] for monomial, coefficient in numerator.terms() if coefficient != 0]
    return min(powers)


def positive_real_roots(expr: sp.Expr, *, z: sp.Symbol = Z, digits: int = 60) -> list[float]:
    """Positive real roots of a rational expression, excluding poles and z=0."""
    numerator = sp.Poly(sp.together(expr).as_numer_denom()[0], z)
    if numerator.degree() <= 0:
        return []
    # Remove exact zero roots before numerical root finding.
    min_power = min(monomial[0] for monomial, coefficient in numerator.terms() if coefficient != 0)
    reduced = sp.Poly(sp.cancel(numerator.as_expr() / z**min_power), z)
    if reduced.degree() <= 0:
        return []
    roots = sp.nroots(reduced, n=digits, maxsteps=500)
    values: list[float] = []
    for root in roots:
        real, imag = sp.re(root), sp.im(root)
        real_f, imag_f = float(real), float(imag)
        if real_f > 1e-12 and abs(imag_f) < 1e-10:
            values.append(real_f)
    return sorted(values)


def _truth_at(R: sp.Expr, M: sp.Expr, point: float, *, x: sp.Symbol = X, z: sp.Symbol = Z) -> bool:
    a = complex(sp.N(R.subs(x, -point), 50))
    margin = complex(sp.N(M.subs(z, point), 50))
    return abs(a.imag) < 1e-12 and (-1e-12 <= a.real <= 1 + 1e-12) and margin.real >= -1e-12


def cp_admissible_intervals(
    R: sp.Expr,
    *,
    r: sp.Expr = sp.Rational(1, 2),
    x: sp.Symbol = X,
    z: sp.Symbol = Z,
) -> list[tuple[float, float]]:
    """Numerically identify positive-axis CPTP intervals for the scalar test map.

    This routine is intended for low-degree exact rational functions. Endpoints are
    generated from zeros of a, a-1, M, and poles. The singleton z=0 is represented
    as (0.0, 0.0) when no punctured neighborhood is admissible.
    """
    M = cp_margin(R, r=r, x=x, z=z)
    critical: list[float] = []
    for expr in (R.subs(x, -z), R.subs(x, -z) - 1, M):
        critical.extend(positive_real_roots(expr, z=z))
    denominator = sp.together(R.subs(x, -z)).as_numer_denom()[1]
    critical.extend(positive_real_roots(denominator, z=z))
    critical = sorted({round(value, 13) for value in critical})

    points = [0.0] + critical
    intervals: list[tuple[float, float]] = []

    # Determine the punctured-neighborhood behavior symbolically. Numerical
    # sampling near z=0 is unreliable because M can begin at order z^5 or
    # higher and be swallowed by floating-point tolerances.
    local_series = sp.series(M, z, 0, 31).removeO().expand()
    local_coefficient = None
    for power in range(31):
        coefficient = sp.simplify(local_series.coeff(z, power))
        if coefficient != 0:
            local_coefficient = coefficient
            break
    if local_coefficient is None:
        raise ValueError("Unable to determine local CP-margin orientation")
    local_ok = bool(local_coefficient > 0)
    if not local_ok:
        intervals.append((0.0, 0.0))

    for left, right in zip(points, critical):
        if right <= left:
            continue
        sample = (left + right) / 2.0
        if _truth_at(R, M, sample, x=x, z=z):
            intervals.append((left, right))

    last = critical[-1] if critical else 0.0
    samples = [max(1.0, last + 1.0), max(2.0, 2.0 * last + 1.0)]
    if all(_truth_at(R, M, sample, x=x, z=z) for sample in samples):
        intervals.append((last, math.inf))

    # Merge adjacent intervals created by a critical endpoint that remains admissible.
    merged: list[tuple[float, float]] = []
    for interval in intervals:
        if not merged:
            merged.append(interval)
            continue
        prev = merged[-1]
        if prev[1] == interval[0] and prev != (0.0, 0.0):
            merged[-1] = (prev[0], interval[1])
        else:
            merged.append(interval)
    return merged
