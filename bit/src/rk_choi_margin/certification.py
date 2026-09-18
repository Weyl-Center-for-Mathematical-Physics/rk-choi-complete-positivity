from __future__ import annotations

"""Certified complete-positivity predicates and one-dimensional component isolation.

The public acceptance interface is deliberately tri-state.  ``PASS`` and ``FAIL``
are certificates for the decimal/rational inputs supplied to the routine;
``UNCERTAIN`` is never silently converted into acceptance.  Exact rational
arithmetic is used for the supported Runge--Kutta methods and phase-covariant
qubit rays.  Polynomial root decisions use exact real-root isolation rather
than companion-matrix eigenvalues.
"""

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, MutableMapping, Sequence

import sympy as sp

from .methods import RKMethod, all_methods

S = sp.Symbol("s")
X = sp.Symbol("x", nonnegative=True, real=True)


class CPStatus(str, Enum):
    PASS = "PASS"
    FAIL = "FAIL"
    UNCERTAIN = "UNCERTAIN"


@dataclass(frozen=True)
class CertifiedScalar:
    name: str
    value: sp.Expr
    relation: str
    status: CPStatus


@dataclass(frozen=True)
class StepCertificate:
    method: str
    status: CPStatus
    inputs: Mapping[str, sp.Expr]
    constraints: tuple[CertifiedScalar, ...]
    stage_domain: tuple[CertifiedScalar, ...]
    diagnostics: Mapping[str, str]


@dataclass(frozen=True)
class RootIsolation:
    lower: sp.Rational
    upper: sp.Rational
    multiplicity: int
    source: tuple[str, ...]
    exact_root: sp.Expr | None = None
    point_status: CPStatus = CPStatus.UNCERTAIN

    @property
    def midpoint(self) -> sp.Rational:
        return (self.lower + self.upper) / 2

    @property
    def width(self) -> sp.Rational:
        return self.upper - self.lower


@dataclass(frozen=True)
class CertifiedComponent:
    """A certified positive-width interval component.

    ``left`` and ``right`` are root brackets, or ``None`` for 0 / the requested
    upper endpoint.  The open sample interval is certified PASS.  Equality
    endpoints are reported separately and are not selected as operational
    controller targets.
    """

    left: RootIsolation | None
    right: RootIsolation | None
    sample: sp.Rational
    lower_value: sp.Rational
    upper_value: sp.Rational
    status: CPStatus

    @property
    def width_lower_bound(self) -> sp.Rational:
        return self.upper_value - self.lower_value


@dataclass(frozen=True)
class ComponentLocation:
    status: CPStatus
    components: tuple[CertifiedComponent, ...]
    isolated_points: tuple[RootIsolation, ...]
    roots: tuple[RootIsolation, ...]
    upper: sp.Rational
    diagnostics: Mapping[str, object]


def as_rational(value: int | float | str | Fraction | sp.Rational | sp.Integer) -> sp.Rational:
    """Convert a value to the exact rational represented by the supplied object.

    Strings are interpreted as exact decimal rationals.  Python floats are
    interpreted through ``float.as_integer_ratio()``, so a certificate applies
    to the exact IEEE-754 binary value that the numerical workflow received,
    not to a nearby printed decimal surrogate.
    """

    if isinstance(value, sp.Rational):
        return value
    if isinstance(value, sp.Integer):
        return sp.Rational(value)
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    if isinstance(value, int):
        return sp.Rational(value)
    if isinstance(value, float):
        if not sp.Float(value).is_finite:
            raise ValueError("value must be finite")
        numerator, denominator = value.as_integer_ratio()
        return sp.Rational(numerator, denominator)
    return sp.Rational(str(value))


def _sign_status(value: sp.Expr, *, nonnegative: bool = True) -> CPStatus:
    v = sp.factor(sp.simplify(value))
    if v.is_real is False:
        return CPStatus.UNCERTAIN
    if v == 0:
        return CPStatus.PASS if nonnegative else CPStatus.FAIL
    if v.is_positive:
        return CPStatus.PASS if nonnegative else CPStatus.FAIL
    if v.is_negative:
        return CPStatus.FAIL if nonnegative else CPStatus.PASS
    # Algebraic expressions can generally be resolved exactly by ``sign``.
    try:
        sign = sp.sign(v)
        if sign == 1:
            return CPStatus.PASS if nonnegative else CPStatus.FAIL
        if sign == -1:
            return CPStatus.FAIL if nonnegative else CPStatus.PASS
        if sign == 0:
            return CPStatus.PASS if nonnegative else CPStatus.FAIL
    except Exception:
        pass
    return CPStatus.UNCERTAIN




def _to_fraction(value: sp.Rational) -> Fraction:
    return Fraction(int(value.p), int(value.q))


def _pair_mul(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def _rk4_pair(z: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    one = (Fraction(1), Fraction(0))
    z2 = _pair_mul(z, z)
    z3 = _pair_mul(z2, z)
    z4 = _pair_mul(z3, z)
    return (
        one[0] + z[0] + z2[0] / 2 + z3[0] / 6 + z4[0] / 24,
        z[1] + z2[1] / 2 + z3[1] / 6 + z4[1] / 24,
    )


def _rk4_fraction_certificate(x: sp.Rational, theta: sp.Rational, kappa: sp.Rational, varpi: sp.Rational) -> tuple[dict[str, sp.Rational], tuple[CertifiedScalar, ...]]:
    xf, tf, kf, vf = map(_to_fraction, (x, theta, kappa, varpi))
    a_pair = _rk4_pair((-xf, Fraction(0)))
    if a_pair[1] != 0:
        raise ArithmeticError("RK4 population multiplier unexpectedly complex")
    a = a_pair[0]
    c = _rk4_pair((-xf * (Fraction(1, 2) + kf), xf * vf))
    one_minus_a = Fraction(1) - a
    A = Fraction(1) - tf * one_minus_a
    D = tf + (Fraction(1) - tf) * a
    F = a + tf * (Fraction(1) - tf) * one_minus_a**2 - c[0] ** 2 - c[1] ** 2
    vals = {
        "one_minus_a": sp.Rational(one_minus_a.numerator, one_minus_a.denominator),
        "A_theta": sp.Rational(A.numerator, A.denominator),
        "D_theta": sp.Rational(D.numerator, D.denominator),
        "F": sp.Rational(F.numerator, F.denominator),
    }
    constraints = tuple(
        CertifiedScalar(name, value, ">=0", CPStatus.PASS if value >= 0 else CPStatus.FAIL)
        for name, value in vals.items()
    )
    return vals, constraints


def _method(method: str | RKMethod) -> RKMethod:
    if isinstance(method, RKMethod):
        return method
    try:
        return all_methods()[method]
    except KeyError as exc:
        raise ValueError(f"unknown RK method {method!r}") from exc


def exact_phase_covariant_expressions(
    method: str | RKMethod,
    x: sp.Expr,
    theta: sp.Expr,
    kappa: sp.Expr,
    varpi: sp.Expr,
) -> dict[str, sp.Expr]:
    """Return exact phase-covariant one-step expressions on a generator ray."""

    rk = _method(method)
    R = rk.stability_function(x=S)
    a = sp.cancel(R.subs(S, -x))
    z = -x * (sp.Rational(1, 2) + kappa - sp.I * varpi)
    zbar = sp.conjugate(z)
    c = sp.cancel(R.subs(S, z))
    cbar = sp.cancel(R.subs(S, zbar))
    one_minus_a = sp.factor(1 - a)
    A = sp.factor(1 - theta * one_minus_a)
    D = sp.factor(theta + (1 - theta) * a)
    F = sp.factor(sp.cancel(a + theta * (1 - theta) * one_minus_a**2 - c * cbar))
    det = rk.stage_resolvent_determinant(x=S)
    det_pop = sp.factor(det.subs(S, -x))
    det_coh = sp.factor(det.subs(S, z))
    return {
        "a": sp.factor(a),
        "c": sp.factor(c),
        "one_minus_a": one_minus_a,
        "A_theta": A,
        "D_theta": D,
        "F": F,
        "stage_population": det_pop,
        "stage_coherence": det_coh,
    }


def certify_cptp_step(
    method: str | RKMethod,
    *,
    x: int | float | str | Fraction | sp.Rational,
    theta: int | float | str | Fraction | sp.Rational,
    kappa: int | float | str | Fraction | sp.Rational,
    varpi: int | float | str | Fraction | sp.Rational,
) -> StepCertificate:
    """Certify a fixed phase-covariant RK step using exact rational arithmetic.

    The returned PASS certificate proves the CPTP inequalities for the exact
    supplied rational value; a float is certified as its exact IEEE-754 binary
    rational.  No absolute tolerance is used.
    """

    rk = _method(method)
    xr = as_rational(x)
    tr = as_rational(theta)
    kr = as_rational(kappa)
    vr = as_rational(varpi)
    if xr < 0 or tr < 0 or tr > 1 or kr < 0:
        return StepCertificate(
            rk.key,
            CPStatus.FAIL,
            {"x": xr, "theta": tr, "kappa": kr, "varpi": vr},
            tuple(),
            tuple(),
            {"reason": "parameters outside the physical ray domain"},
        )

    if rk.key == "rk4" and all(isinstance(v, sp.Rational) for v in (xr, tr, kr, vr)):
        _vals, constraints = _rk4_fraction_certificate(xr, tr, kr, vr)
        expr = {
            "stage_population": sp.Integer(1),
            "stage_coherence": sp.Integer(1),
        }
    else:
        expr = exact_phase_covariant_expressions(rk, xr, tr, kr, vr)
        constraints = (
            CertifiedScalar("1-a", expr["one_minus_a"], ">=0", _sign_status(expr["one_minus_a"])),
            CertifiedScalar("A_theta", expr["A_theta"], ">=0", _sign_status(expr["A_theta"])),
            CertifiedScalar("D_theta", expr["D_theta"], ">=0", _sign_status(expr["D_theta"])),
            CertifiedScalar("F", expr["F"], ">=0", _sign_status(expr["F"])),
        )
    stage = (
        CertifiedScalar(
            "det(I+xA)",
            expr["stage_population"],
            "!=0",
            CPStatus.FAIL if sp.simplify(expr["stage_population"]) == 0 else CPStatus.PASS,
        ),
        CertifiedScalar(
            "det(I-zA)",
            expr["stage_coherence"],
            "!=0",
            CPStatus.FAIL if sp.simplify(expr["stage_coherence"]) == 0 else CPStatus.PASS,
        ),
    )
    statuses = [item.status for item in (*constraints, *stage)]
    status = CPStatus.FAIL if CPStatus.FAIL in statuses else CPStatus.UNCERTAIN if CPStatus.UNCERTAIN in statuses else CPStatus.PASS
    return StepCertificate(
        rk.key,
        status,
        {"x": xr, "theta": tr, "kappa": kr, "varpi": vr},
        constraints,
        stage,
        {
            "arithmetic": "exact rational/algebraic SymPy arithmetic",
            "input_semantics": "certificate applies to the exact supplied value; floats use exact IEEE-754 rationals",
        },
    )


def _primitive_polynomial(expr: sp.Expr, x: sp.Symbol = X) -> tuple[sp.Poly, sp.Expr]:
    num, den = sp.fraction(sp.cancel(expr))
    poly = sp.Poly(num, x, domain=sp.QQ)
    _, primitive = poly.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
        den = -den
    return primitive, sp.factor(den)


def ray_constraint_polynomials(
    method: str | RKMethod,
    *,
    theta: int | float | str | Fraction | sp.Rational,
    kappa: int | float | str | Fraction | sp.Rational,
    varpi: int | float | str | Fraction | sp.Rational,
) -> dict[str, tuple[sp.Poly, sp.Expr]]:
    """Return exact active-constraint numerators on a fixed ray."""

    tr, kr, vr = map(as_rational, (theta, kappa, varpi))
    expr = exact_phase_covariant_expressions(method, X, tr, kr, vr)
    out: dict[str, tuple[sp.Poly, sp.Expr]] = {}
    for name in ("one_minus_a", "A_theta", "D_theta", "F"):
        out[name] = _primitive_polynomial(expr[name], X)
    # Stage determinants are operational-domain boundaries, not CP inequalities.
    for name in ("stage_population", "stage_coherence"):
        value = sp.factor(expr[name] * sp.conjugate(expr[name]))
        out[name] = _primitive_polynomial(value, X)
    return out


def _algebraic_sign(value: sp.Expr) -> int | None:
    """Return the exact sign of a real algebraic expression when decidable."""

    value = sp.simplify(value)
    if value == 0 or value.is_zero:
        return 0
    if value.is_positive:
        return 1
    if value.is_negative:
        return -1
    try:
        sign = sp.sign(value)
    except Exception:
        return None
    if sign in (-1, 0, 1):
        return int(sign)
    return None


def _same_algebraic_root(left: sp.Expr, right: sp.Expr) -> bool:
    if left == right:
        return True
    try:
        result = left.equals(right)
        if result is not None:
            return bool(result)
    except Exception:
        pass
    return sp.simplify(left - right) == 0


def _bracket_algebraic_root(root: sp.Expr, digits: int) -> tuple[sp.Rational, sp.Rational]:
    """Return an exact rational bracket enclosing a real algebraic root."""

    if isinstance(root, (sp.Integer, sp.Rational)):
        value = sp.Rational(root)
        return value, value
    if hasattr(root, "_get_interval"):
        interval = root._get_interval().refine_size(sp.Rational(1, 10**digits))
        return sp.Rational(interval.a), sp.Rational(interval.b)
    # ``real_roots(..., radicals=False)`` should return rational values or
    # CRootOf objects.  Retain a conservative fallback for future SymPy types.
    approx = sp.N(root, digits + 15)
    center = sp.Rational(str(approx))
    radius = sp.Rational(1, 10**digits)
    return center - radius, center + radius


def _root_in_bracket(root: sp.Expr, lower: sp.Rational, upper: sp.Rational) -> bool:
    left = _algebraic_sign(root - lower)
    right = _algebraic_sign(root - upper)
    return left is not None and right is not None and left >= 0 and right <= 0


def _exact_root_from_bracket(
    poly: sp.Poly, lower: sp.Rational, upper: sp.Rational
) -> tuple[sp.Expr, int] | None:
    """Recover the unique exact real root inside an isolating bracket."""

    matches: list[tuple[sp.Expr, int]] = []
    for root, multiplicity in sp.real_roots(
        poly.as_expr(), multiple=False, radicals=False
    ):
        if _root_in_bracket(root, lower, upper):
            matches.append((root, int(multiplicity)))
    if len(matches) != 1:
        return None
    return matches[0]


def isolate_constraint_roots_exact(
    polynomials: Mapping[str, sp.Poly],
    upper: sp.Rational,
    *,
    digits: int = 30,
    max_digits: int = 240,
    telemetry: MutableMapping[str, int] | None = None,
) -> tuple[list[RootIsolation], bool]:
    """Isolate all named real roots with adaptive separation.

    The fast path uses exact Sturm intervals.  Exact algebraic root objects are
    constructed only for overlapping brackets (to distinguish a shared root
    from a near collision) or later for a potential isolated equality point.
    This keeps the component locator practical for the higher-degree candidate
    polynomials while retaining fail-closed completeness semantics.
    """

    if digits < 8:
        raise ValueError("root isolation digits must be >= 8")
    if telemetry is not None:
        telemetry.setdefault("root_isolation_ops", 0)
        telemetry.setdefault("precision_escalations", 0)

    def _count_root_isolation() -> None:
        if telemetry is not None:
            telemetry["root_isolation_ops"] += 1

    def _count_precision_escalation() -> None:
        if telemetry is not None:
            telemetry["precision_escalations"] += 1

    current_digits = digits
    while True:
        eps = sp.Rational(1, 10**current_digits)
        raw: list[dict[str, object]] = []
        unresolved = False
        for name, poly in polynomials.items():
            if poly.is_zero or poly.degree() <= 0:
                continue
            _count_root_isolation()
            for (lo, hi), multiplicity in poly.intervals(eps=eps):
                lo_r, hi_r = sp.Rational(lo), sp.Rational(hi)
                if hi_r <= 0 or lo_r > upper:
                    continue

                # A Sturm bracket may touch or cross a requested-domain boundary
                # even when its unique root lies outside the domain (for example,
                # lower == upper_requested < root < bracket_upper).  Never clamp
                # such a bracket before proving exact root membership.  Exact
                # Sturm root counts at the rational endpoints decide membership
                # without constructing every algebraic root of a high-degree
                # polynomial; only genuine equality points or later overlap/
                # isolated-point logic need a CRootOf object.
                exact_root = None

                if lo_r <= 0:
                    if lo_r == hi_r == 0:
                        # The public component domain is x>0.  A root at the
                        # origin is represented by the fixed identity point,
                        # not as a positive boundary root.
                        continue
                    _count_root_isolation()
                    restricted = poly.intervals(
                        eps=eps, inf=sp.Rational(0), sup=hi_r
                    )
                    positive = []
                    for (rlo, rhi), rmult in restricted:
                        rlo_r, rhi_r = sp.Rational(rlo), sp.Rational(rhi)
                        if rlo_r == rhi_r == 0:
                            continue
                        if rhi_r <= 0 or rhi_r < lo_r or rlo_r > hi_r:
                            continue
                        positive.append((rlo_r, rhi_r, int(rmult)))
                    if not positive:
                        continue
                    if len(positive) != 1:
                        unresolved = True
                        break
                    lo_r, hi_r, multiplicity = positive[0]
                    if lo_r <= 0:
                        unresolved = True
                        break

                if hi_r >= upper:
                    _count_root_isolation()
                    restricted = poly.intervals(
                        eps=eps, inf=lo_r, sup=upper
                    )
                    in_domain = []
                    for (rlo, rhi), rmult in restricted:
                        rlo_r, rhi_r = sp.Rational(rlo), sp.Rational(rhi)
                        if rhi_r <= 0 or rhi_r < lo_r or rlo_r > hi_r:
                            continue
                        in_domain.append((rlo_r, rhi_r, int(rmult)))
                    if not in_domain:
                        continue
                    if len(in_domain) != 1:
                        unresolved = True
                        break
                    lo_r, hi_r, multiplicity = in_domain[0]
                    if lo_r == hi_r == upper:
                        exact_root = upper
                    elif lo_r <= 0 or hi_r >= upper:
                        # A strict interior root must be separated from both
                        # requested-domain boundaries before interval signs can
                        # certify completeness.
                        unresolved = True
                        break

                if lo_r == hi_r and exact_root is None:
                    exact_root = lo_r

                lo_r = max(lo_r, sp.Rational(0))
                hi_r = min(hi_r, upper)
                if lo_r == hi_r == 0 or hi_r < lo_r:
                    continue
                raw.append(
                    {
                        "name": name,
                        "poly": poly,
                        "lower": lo_r,
                        "upper": hi_r,
                        "multiplicity": int(multiplicity),
                        "root": exact_root,
                    }
                )
            if unresolved:
                break
        raw.sort(key=lambda item: (item["lower"], item["upper"]))
        groups: list[list[dict[str, object]]] = []
        for row in raw:
            if not groups or row["lower"] > max(item["upper"] for item in groups[-1]):
                groups.append([row])
            else:
                groups[-1].append(row)

        roots: list[RootIsolation] = []
        for group in groups:
            if len(group) == 1:
                row = group[0]
                roots.append(
                    RootIsolation(
                        sp.Rational(row["lower"]),
                        sp.Rational(row["upper"]),
                        int(row["multiplicity"]),
                        (str(row["name"]),),
                        exact_root=row.get("root"),
                    )
                )
                continue

            # Overlap can mean either a shared root of several constraints or
            # distinct roots closer than the current bracket width.  Resolve
            # that distinction using exact algebraic roots for this small group.
            exact_rows: list[dict[str, object]] = []
            for row in group:
                match = (
                    (row["root"], int(row["multiplicity"]))
                    if row.get("root") is not None
                    else _exact_root_from_bracket(
                        row["poly"], sp.Rational(row["lower"]), sp.Rational(row["upper"])
                    )
                )
                if match is None:
                    unresolved = True
                    break
                exact_rows.append({**row, "root": match[0], "multiplicity": match[1]})
            if unresolved:
                break

            exact_groups: list[list[dict[str, object]]] = []
            for row in exact_rows:
                placed = False
                for exact_group in exact_groups:
                    if _same_algebraic_root(row["root"], exact_group[0]["root"]):
                        exact_group.append(row)
                        placed = True
                        break
                if not placed:
                    exact_groups.append([row])

            for exact_group in exact_groups:
                root = exact_group[0]["root"]
                lo_r, hi_r = _bracket_algebraic_root(root, current_digits)
                roots.append(
                    RootIsolation(
                        max(sp.Rational(0), lo_r),
                        min(upper, hi_r),
                        max(int(row["multiplicity"]) for row in exact_group),
                        tuple(sorted({str(row["name"]) for row in exact_group})),
                        exact_root=root,
                    )
                )

        if unresolved:
            if current_digits >= max_digits:
                return [], False
            _count_precision_escalation()
            current_digits = min(max_digits, current_digits + 20)
            continue

        roots.sort(key=lambda root: (root.lower, root.upper))
        overlap = any(
            roots[index].upper >= roots[index + 1].lower
            and not (
                roots[index].exact_root is not None
                and roots[index + 1].exact_root is not None
                and _same_algebraic_root(
                    roots[index].exact_root, roots[index + 1].exact_root
                )
            )
            for index in range(len(roots) - 1)
        )
        if overlap:
            if current_digits >= max_digits:
                return [], False
            _count_precision_escalation()
            current_digits = min(max_digits, current_digits + 20)
            continue
        return roots, True


def exact_root_for_isolation(
    root: RootIsolation, polynomials: Mapping[str, sp.Poly]
) -> sp.Expr | None:
    """Return the exact root represented by an isolation record."""

    if root.exact_root is not None:
        return root.exact_root
    for source in root.source:
        poly = polynomials.get(source)
        if poly is None:
            continue
        match = _exact_root_from_bracket(poly, root.lower, root.upper)
        if match is not None:
            return match[0]
    return None


def _certify_phase_covariant_at_expression(
    method: str | RKMethod,
    *,
    x: sp.Expr,
    theta: sp.Rational,
    kappa: sp.Rational,
    varpi: sp.Rational,
) -> StepCertificate:
    """Certify a possibly algebraic step value without rational coercion."""

    rk = _method(method)
    expr = exact_phase_covariant_expressions(rk, x, theta, kappa, varpi)
    constraints = tuple(
        CertifiedScalar(name, sp.factor(expr[name]), ">=0", _sign_status(expr[name]))
        for name in ("one_minus_a", "A_theta", "D_theta", "F")
    )
    stage = tuple(
        CertifiedScalar(
            name,
            sp.factor(expr[name]),
            "!=0",
            CPStatus.FAIL
            if sp.simplify(expr[name]) == 0
            else CPStatus.PASS
            if _sign_status(expr[name] * sp.conjugate(expr[name])) is CPStatus.PASS
            else CPStatus.UNCERTAIN,
        )
        for name in ("stage_population", "stage_coherence")
    )
    statuses = [item.status for item in (*constraints, *stage)]
    status = (
        CPStatus.FAIL
        if CPStatus.FAIL in statuses
        else CPStatus.UNCERTAIN
        if CPStatus.UNCERTAIN in statuses
        else CPStatus.PASS
    )
    return StepCertificate(
        rk.key,
        status,
        {"x": x, "theta": theta, "kappa": kappa, "varpi": varpi},
        constraints,
        stage,
        {
            "arithmetic": "exact real-algebraic SymPy arithmetic",
            "point_semantics": "all active CPTP inequalities and stage-domain predicates evaluated at the exact algebraic root",
        },
    )


def _certify_isolated_phase_covariant_root(
    method: str | RKMethod,
    root: RootIsolation,
    *,
    theta: sp.Rational,
    kappa: sp.Rational,
    varpi: sp.Rational,
) -> CPStatus:
    """Certify every active predicate at an isolated algebraic root.

    Rational equality points are substituted directly.  For a non-rational
    algebraic root, the source constraints vanish by exact root identity and
    every other constraint has no zero in the certified isolating interval;
    its exact rational midpoint sign is therefore its sign at the root.
    """

    if root.exact_root is None:
        return CPStatus.UNCERTAIN
    if isinstance(root.exact_root, (sp.Integer, sp.Rational)):
        return _certify_phase_covariant_at_expression(
            method, x=root.exact_root, theta=theta, kappa=kappa, varpi=varpi
        ).status

    source = set(root.source)
    if "stage_population" in source or "stage_coherence" in source:
        return CPStatus.FAIL
    expr = exact_phase_covariant_expressions(
        method, root.midpoint, theta, kappa, varpi
    )
    statuses: list[CPStatus] = []
    for name in ("one_minus_a", "A_theta", "D_theta", "F"):
        statuses.append(CPStatus.PASS if name in source else _sign_status(expr[name]))
    for name in ("stage_population", "stage_coherence"):
        if name in source:
            statuses.append(CPStatus.FAIL)
        else:
            value = sp.simplify(expr[name])
            statuses.append(CPStatus.PASS if value != 0 else CPStatus.UNCERTAIN)
    if CPStatus.FAIL in statuses:
        return CPStatus.FAIL
    if CPStatus.UNCERTAIN in statuses:
        return CPStatus.UNCERTAIN
    return CPStatus.PASS


def locate_cptp_components(
    method: str | RKMethod,
    *,
    theta: int | float | str | Fraction | sp.Rational,
    kappa: int | float | str | Fraction | sp.Rational,
    varpi: int | float | str | Fraction | sp.Rational,
    upper: int | float | str | Fraction | sp.Rational,
    root_digits: int = 30,
) -> ComponentLocation:
    """Certify the CPTP decomposition on ``0 < x <= upper``.

    Positive-width components use exact Sturm isolation and exact rational
    interval samples.  A zero-width point is listed as CPTP only after every
    active inequality and stage-domain predicate passes at the exact algebraic
    root.  Unresolved root identity or sign information yields ``UNCERTAIN``.
    """

    ur = as_rational(upper)
    if ur <= 0:
        return ComponentLocation(
            CPStatus.FAIL, tuple(), tuple(), tuple(), ur,
            {"reason": "upper must be positive"},
        )
    tr, kr, vr = map(as_rational, (theta, kappa, varpi))
    named = {
        name: poly
        for name, (poly, _denominator) in ray_constraint_polynomials(
            method, theta=tr, kappa=kr, varpi=vr
        ).items()
    }
    telemetry = {
        "root_isolation_ops": 0,
        "certification_ops": 0,
        "precision_escalations": 0,
    }
    roots, separated = isolate_constraint_roots_exact(
        named, ur, digits=root_digits, telemetry=telemetry
    )
    if not separated:
        return ComponentLocation(
            CPStatus.UNCERTAIN, tuple(), tuple(), tuple(), ur,
            {
                "root_method": "exact Sturm isolation with adaptive refinement",
                "reason": "distinct boundary roots could not be certified as separated",
                **telemetry,
            },
        )

    # Represent intervals by root indices first.  This lets exact equality-point
    # statuses be determined after we know which roots are unflanked by a PASS
    # interval, without rebuilding component objects around stale root records.
    interval_rows: list[dict[str, object]] = []
    cursor = sp.Rational(0)
    left_index: int | None = None
    for index, root in enumerate(roots):
        if cursor < root.lower:
            interval_rows.append(
                {"lo": cursor, "hi": root.lower, "left": left_index, "right": index}
            )
        cursor = max(cursor, root.upper)
        left_index = index
    if cursor < ur:
        interval_rows.append(
            {"lo": cursor, "hi": ur, "left": left_index, "right": None}
        )

    interval_uncertain = False
    for row in interval_rows:
        lo, hi = sp.Rational(row["lo"]), sp.Rational(row["hi"])
        sample = (lo + hi) / 2
        row["sample"] = sample
        telemetry["certification_ops"] += 1
        row["status"] = certify_cptp_step(
            method, x=sample, theta=tr, kappa=kr, varpi=vr
        ).status
        if row["status"] is CPStatus.UNCERTAIN:
            interval_uncertain = True

    adjacent_pass = [False] * len(roots)
    for row in interval_rows:
        if row["status"] is not CPStatus.PASS:
            continue
        if row["left"] is not None:
            adjacent_pass[int(row["left"])] = True
        if row["right"] is not None:
            adjacent_pass[int(row["right"])] = True

    final_roots: list[RootIsolation] = []
    unresolved_points = 0
    for index, root in enumerate(roots):
        point_status = CPStatus.UNCERTAIN
        exact_root = root.exact_root
        if not adjacent_pass[index]:
            if exact_root is None:
                exact_root = exact_root_for_isolation(root, named)
            root_with_exact = RootIsolation(
                root.lower, root.upper, root.multiplicity, root.source, exact_root
            )
            telemetry["certification_ops"] += 1
            point_status = _certify_isolated_phase_covariant_root(
                method, root_with_exact, theta=tr, kappa=kr, varpi=vr
            )
            if point_status is CPStatus.UNCERTAIN:
                unresolved_points += 1
        final_roots.append(
            RootIsolation(
                root.lower, root.upper, root.multiplicity, root.source,
                exact_root, point_status,
            )
        )

    components: list[CertifiedComponent] = []
    for row in interval_rows:
        if row["status"] is not CPStatus.PASS:
            continue
        left = None if row["left"] is None else final_roots[int(row["left"])]
        right = None if row["right"] is None else final_roots[int(row["right"])]
        components.append(
            CertifiedComponent(
                left, right, sp.Rational(row["sample"]),
                sp.Rational(row["lo"]), sp.Rational(row["hi"]), CPStatus.PASS,
            )
        )

    isolated = tuple(
        root
        for index, root in enumerate(final_roots)
        if not adjacent_pass[index] and root.point_status is CPStatus.PASS
    )
    status = CPStatus.UNCERTAIN if unresolved_points or interval_uncertain else CPStatus.PASS
    return ComponentLocation(
        status, tuple(components), isolated, tuple(final_roots), ur,
        {
            "root_method": "exact Sturm isolation with adaptive algebraic separation",
            "decision_arithmetic": "exact rational interval samples and exact algebraic equality-point evaluation",
            "isolated_points": "listed only after all active constraints pass at the exact algebraic root",
            "completeness_scope": "rational-coefficient one-dimensional polynomial constraints with exact requested-domain endpoint membership",
            "unresolved_point_count": str(unresolved_points),
            **telemetry,
        },
    )


def choose_certified_interior(
    component: CertifiedComponent,
    *,
    ceiling: int | float | str | Fraction | sp.Rational,
    safety_fraction: sp.Rational = sp.Rational(1, 20),
) -> sp.Rational | None:
    """Choose a rational point strictly inside a certified PASS component."""

    cap = as_rational(ceiling)
    lo = component.lower_value
    hi = min(component.upper_value, cap)
    if hi <= lo:
        return None
    width = hi - lo
    # Bias toward the largest admissible point but retain at least five percent
    # of the certified interval as step-space safety.
    candidate = hi - safety_fraction * width
    if not lo < candidate < hi:
        candidate = (lo + hi) / 2
    return sp.Rational(candidate)


def all_principal_minors(matrix: sp.Matrix) -> dict[tuple[int, ...], sp.Expr]:
    if matrix.rows != matrix.cols:
        raise ValueError("matrix must be square")
    n = matrix.rows
    return {
        idx: sp.factor(matrix.extract(idx, idx).det())
        for r in range(1, n + 1)
        for idx in combinations(range(n), r)
    }


def certify_hermitian_psd_exact(matrix: sp.Matrix) -> tuple[CPStatus, dict[tuple[int, ...], sp.Expr]]:
    if matrix != matrix.conjugate().T:
        return CPStatus.FAIL, {}
    minors = all_principal_minors(matrix)
    statuses = [_sign_status(value) for value in minors.values()]
    if CPStatus.FAIL in statuses:
        return CPStatus.FAIL, minors
    if CPStatus.UNCERTAIN in statuses:
        return CPStatus.UNCERTAIN, minors
    return CPStatus.PASS, minors
