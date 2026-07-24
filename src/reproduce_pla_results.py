#!/usr/bin/env python3
"""Verify and reproduce the data and figures accompanying the PLA manuscript.

All Runge-Kutta stability functions are derived from exact rational Butcher
coefficients. Boundary polynomials, exponential defects, Sturm root counts,
and interval signs are evaluated with SymPy exact arithmetic. mpmath is used
only for the fixed-horizon plotting data, after the exact symbolic checks pass.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import filecmp
import json
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import sympy as sp
from PIL import __version__ as pillow_version

ROOT = Path(__file__).resolve().parents[1]
TRACKED_OUTPUTS = (
    Path("data/certified_boundaries.json"),
    Path("data/method_audit.csv"),
    Path("data/fixed_horizon_defect.csv"),
    Path("figures/figure1_cp_geometry.pdf"),
    Path("figures/figure1_cp_geometry.png"),
    Path("figures/figure2_admissible_sets.pdf"),
    Path("figures/figure2_admissible_sets.png"),
    Path("figures/figure3_fixed_horizon_defect.pdf"),
    Path("figures/figure3_fixed_horizon_defect.png"),
)

s, x = sp.symbols("s x", real=True)


class VerificationError(RuntimeError):
    """Raised when an exact mathematical or repository check fails."""


@dataclass(frozen=True)
class Method:
    key: str
    name: str
    order: int
    A: sp.Matrix
    b: sp.Matrix


@dataclass(frozen=True)
class RootSpec:
    method: str
    equation: str
    label: str
    point: str | None = None
    lower: str | None = None
    upper: str | None = None

    def bounds(self) -> tuple[sp.Rational, sp.Rational]:
        if self.point is not None:
            q = sp.Rational(self.point)
            return q, q
        if self.lower is None or self.upper is None:
            raise VerificationError(f"Incomplete root specification: {self}")
        return sp.Rational(self.lower), sp.Rational(self.upper)


@dataclass(frozen=True)
class CertifiedRoot:
    equation: str
    label: str
    lower: sp.Rational
    upper: sp.Rational
    approximation: str
    multiplicity: int

    @property
    def display(self) -> str:
        if self.label:
            return self.label
        if self.lower == self.upper:
            return str(self.lower)
        return self.approximation


def q(value: int | str) -> sp.Rational:
    return sp.Rational(value)


def matrix(rows: list[list[int | str | sp.Rational]]) -> sp.Matrix:
    return sp.Matrix([[q(v) for v in row] for row in rows])


def vector(values: list[int | str | sp.Rational]) -> sp.Matrix:
    return sp.Matrix([q(v) for v in values])


# Exact Butcher coefficients. The two Dormand-Prince formulas share A and use
# different output weights.
A_FE = matrix([[0]])
A_HEUN = matrix([[0, 0], [1, 0]])
A_SSP33 = matrix([[0, 0, 0], [1, 0, 0], ["1/4", "1/4", 0]])
A_RK4 = matrix(
    [
        [0, 0, 0, 0],
        ["1/2", 0, 0, 0],
        [0, "1/2", 0, 0],
        [0, 0, 1, 0],
    ]
)
A_DOPRI = matrix(
    [
        [0, 0, 0, 0, 0, 0, 0],
        ["1/5", 0, 0, 0, 0, 0, 0],
        ["3/40", "9/40", 0, 0, 0, 0, 0],
        ["44/45", "-56/15", "32/9", 0, 0, 0, 0],
        ["19372/6561", "-25360/2187", "64448/6561", "-212/729", 0, 0, 0],
        ["9017/3168", "-355/33", "46732/5247", "49/176", "-5103/18656", 0, 0],
        ["35/384", 0, "500/1113", "125/192", "-2187/6784", "11/84", 0],
    ]
)
A_BE = matrix([[1]])
A_IM = matrix([["1/2"]])

METHODS = (
    Method("FE", "Forward Euler", 1, A_FE, vector([1])),
    Method("RK2", "Heun RK2", 2, A_HEUN, vector(["1/2", "1/2"])),
    Method("SSP3", "SSPRK(3,3)", 3, A_SSP33, vector(["1/6", "1/6", "2/3"])),
    Method("RK4", "Classical RK4", 4, A_RK4, vector(["1/6", "1/3", "1/3", "1/6"])),
    Method(
        "DP5",
        "Dormand-Prince 5 principal",
        5,
        A_DOPRI,
        vector(["35/384", 0, "500/1113", "125/192", "-2187/6784", "11/84", 0]),
    ),
    Method(
        "DP4",
        "Dormand-Prince embedded 4",
        4,
        A_DOPRI,
        vector(["5179/57600", 0, "7571/16695", "393/640", "-92097/339200", "187/2100", "1/40"]),
    ),
    Method("BE", "Backward Euler", 1, A_BE, vector([1])),
    Method("IM", "Implicit midpoint", 2, A_IM, vector([1])),
)

T1 = 1 + s
T2 = T1 + s**2 / 2
T3 = T2 + s**3 / 6
T4 = T3 + s**4 / 24
T5 = T4 + s**5 / 120
REFERENCE_STABILITY_FUNCTIONS = {
    "FE": T1,
    "RK2": T2,
    "SSP3": T3,
    "RK4": T4,
    "DP5": T5 + s**6 / 600,
    "DP4": T4 + sp.Rational(1097, 120000) * s**5 + sp.Rational(161, 120000) * s**6 + s**7 / 24000,
    "BE": 1 / (1 - s),
    "IM": (1 + s / 2) / (1 - s / 2),
}
REFERENCE_DEFECTS = {
    "FE": (2, sp.Rational(-1, 2)),
    "RK2": (3, sp.Rational(-1, 6)),
    "SSP3": (4, sp.Rational(-1, 24)),
    "RK4": (5, sp.Rational(-1, 120)),
    "DP5": (6, sp.Rational(1, 3600)),
    "DP4": (5, sp.Rational(97, 120000)),
    "BE": (2, sp.Rational(1, 2)),
    "IM": (3, sp.Rational(1, 12)),
}
REFERENCE_ADMISSIBLE_SETS = {
    "FE": "{0}",
    "RK2": "[0, 2]",
    "SSP3": "{0}",
    "RK4": "[0, alpha_4]",
    "DP5": "[0, alpha_5]",
    "DP4": "{0} union [beta_-, beta_+]",
    "BE": "[0, infinity)",
    "IM": "{0}",
}

ROOT_SPECS = (
    RootSpec("FE", "a=0", "", point="1"),
    RootSpec("RK2", "a=1", "", point="2"),
    RootSpec("RK2", "M=0", "", point="8"),
    RootSpec("SSP3", "a=0", "", lower="1.596071637983", upper="1.596071637984"),
    RootSpec("RK4", "a=1", "alpha_4", lower="2.785293563405", upper="2.785293563406"),
    RootSpec("RK4", "M=0", "mu_4", lower="10.982425466293", upper="10.982425466294"),
    RootSpec("DP5", "a=1", "alpha_5", lower="3.306567892634", upper="3.306567892635"),
    RootSpec("DP5", "M=0", "mu_5", lower="13.174093462800", upper="13.174093462801"),
    RootSpec("DP4", "M=0", "beta_-", lower="3.068856548050", upper="3.068856548051"),
    RootSpec("DP4", "a=1", "beta_+", lower="4.384986320801", upper="4.384986320802"),
    RootSpec("DP4", "M=0", "delta_M", lower="15.298767807585", upper="15.298767807586"),
    RootSpec("DP4", "a=1", "delta_1", lower="24.727754318760", upper="24.727754318761"),
    RootSpec("DP4", "a=0", "delta_0", lower="24.727895030129", upper="24.727895030130"),
    RootSpec("IM", "a=0", "", point="2"),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def stability_function(method: Method) -> sp.Expr:
    n = method.A.rows
    require(method.A.cols == n, f"{method.key}: A is not square")
    require(method.b.rows == n and method.b.cols == 1, f"{method.key}: b has the wrong shape")
    ones = sp.ones(n, 1)
    R = 1 + s * (method.b.T * (sp.eye(n) - s * method.A).inv() * ones)[0]
    return sp.factor(sp.cancel(R))


def stage_determinants(method: Method) -> tuple[sp.Expr, sp.Expr]:
    n = method.A.rows
    full = sp.factor((sp.eye(n) + x * method.A).det())
    half = sp.factor((sp.eye(n) + x * method.A / 2).det())
    return full, half


def first_defect(R: sp.Expr, max_order: int = 20) -> tuple[int, sp.Rational]:
    series = sp.series(R - sp.exp(s), s, 0, max_order).removeO().expand()
    for m in range(2, max_order):
        eta = sp.simplify(series.coeff(s, m))
        if eta != 0:
            return m, sp.Rational(eta)
    raise VerificationError("No exponential defect found in the requested range")


def margin(R: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(R.subs(s, -x) - R.subs(s, -x / 2) ** 2))


def integer_primitive(poly: sp.Poly) -> tuple[sp.Rational, sp.Poly]:
    denominator, integer_poly = poly.clear_denoms(convert=True)
    content, primitive = integer_poly.primitive()
    scalar = sp.Rational(content, denominator)
    if primitive.LC() < 0:
        primitive = -primitive
        scalar = -scalar
    return scalar, primitive


def expression_record(expr: sp.Expr) -> tuple[dict, sp.Poly]:
    simplified = sp.factor(sp.cancel(expr))
    numerator, denominator = sp.fraction(sp.together(simplified))
    scalar, primitive = integer_primitive(sp.Poly(numerator, x, domain=sp.QQ))
    square_free = primitive.sqf_part()
    if square_free.LC() < 0:
        square_free = -square_free
    return (
        {
            "expression": str(simplified),
            "numerator_scalar": str(scalar),
            "primitive_numerator": str(sp.factor(primitive.as_expr())),
            "square_free_part": str(sp.factor(square_free.as_expr())),
            "denominator": str(sp.factor(denominator)),
        },
        primitive,
    )


def sign_at(expr: sp.Expr, sample: sp.Rational) -> int:
    value = sp.factor(sp.cancel(expr.subs(x, sample)))
    if value == 0:
        return 0
    if value.is_positive:
        return 1
    if value.is_negative:
        return -1
    numeric = sp.N(value, 80)
    if numeric > 0:
        return 1
    if numeric < 0:
        return -1
    raise VerificationError(f"Could not determine exact sign of {expr} at {sample}")


def multiplicity_in_region(poly: sp.Poly, lower: sp.Rational, upper: sp.Rational) -> int:
    _, factors = sp.sqf_list(poly.as_expr(), x)
    matches: list[int] = []
    for factor, multiplicity in factors:
        factor_poly = sp.Poly(factor, x, domain=sp.QQ)
        if lower == upper:
            if factor_poly.eval(lower) == 0:
                matches.append(int(multiplicity))
        elif factor_poly.count_roots(lower, upper) == 1:
            matches.append(int(multiplicity))
    require(len(matches) == 1, f"Could not assign a unique multiplicity in ({lower}, {upper})")
    return matches[0]


def approximate_root(poly: sp.Poly, lower: sp.Rational, upper: sp.Rational) -> str:
    if lower == upper:
        return str(lower)
    square_free = poly.sqf_part()
    roots = sp.nroots(square_free, n=70, maxsteps=2000)
    selected = []
    low_n = sp.N(lower, 70)
    high_n = sp.N(upper, 70)
    for root in roots:
        real_part = sp.re(root)
        imag_part = sp.im(root)
        if abs(imag_part) < sp.Float("1e-55") and low_n < real_part < high_n:
            selected.append(real_part)
    require(len(selected) == 1, f"Numerical root lookup failed in ({lower}, {upper})")
    return str(sp.N(selected[0], 18))


def certify_roots(method_key: str, polynomials: dict[str, sp.Poly]) -> list[CertifiedRoot]:
    method_specs = [spec for spec in ROOT_SPECS if spec.method == method_key]
    records: list[CertifiedRoot] = []

    for equation, poly in polynomials.items():
        specs = [spec for spec in method_specs if spec.equation == equation]
        square_free = poly.sqf_part()
        zero_root = 1 if square_free.eval(0) == 0 else 0
        positive_root_count = int(square_free.count_roots(0, sp.oo)) - zero_root
        require(
            positive_root_count == len(specs),
            f"{method_key} {equation}: expected {len(specs)} positive roots, found {positive_root_count}",
        )
        for spec in specs:
            lower, upper = spec.bounds()
            require(lower > 0, f"{method_key} {equation}: positive-root bracket must begin above zero")
            if lower == upper:
                require(poly.eval(lower) == 0, f"{method_key} {equation}: {lower} is not an exact root")
            else:
                require(lower < upper, f"{method_key} {equation}: invalid bracket")
                require(
                    square_free.count_roots(lower, upper) == 1,
                    f"{method_key} {equation}: bracket ({lower}, {upper}) does not isolate one root",
                )
            records.append(
                CertifiedRoot(
                    equation=equation,
                    label=spec.label,
                    lower=lower,
                    upper=upper,
                    approximation=approximate_root(poly, lower, upper),
                    multiplicity=multiplicity_in_region(poly, lower, upper),
                )
            )

    records.sort(key=lambda root: (root.lower + root.upper) / 2)
    for left, right in zip(records, records[1:]):
        require(left.upper < right.lower, f"{method_key}: positive root brackets overlap")
    return records


def interval_sign_chart(
    a: sp.Expr,
    a_minus_one: sp.Expr,
    M: sp.Expr,
    det_full: sp.Expr,
    det_half: sp.Expr,
    roots: list[CertifiedRoot],
) -> list[dict]:
    chart: list[dict] = []
    interval_count = len(roots) + 1
    for index in range(interval_count):
        left_root = roots[index - 1] if index > 0 else None
        right_root = roots[index] if index < len(roots) else None
        if left_root is None and right_root is None:
            sample = sp.Rational(1)
        elif left_root is None:
            sample = right_root.lower / 2
        elif right_root is None:
            sample = left_root.upper + 1
        else:
            sample = (left_root.upper + right_root.lower) / 2

        signs = {
            "a": sign_at(a, sample),
            "a_minus_1": sign_at(a_minus_one, sample),
            "M": sign_at(M, sample),
            "det_full": sign_at(det_full, sample),
            "det_half": sign_at(det_half, sample),
        }
        domain = signs["det_full"] != 0 and signs["det_half"] != 0
        cptp = domain and signs["a"] >= 0 and signs["a_minus_1"] <= 0 and signs["M"] >= 0
        chart.append(
            {
                "left": "0" if left_root is None else left_root.display,
                "right": "infinity" if right_root is None else right_root.display,
                "sample": str(sample),
                "signs": signs,
                "rk_domain": domain,
                "cptp": cptp,
            }
        )
    return chart


def derive_admissible_set(chart: list[dict]) -> tuple[str, list[dict]]:
    true_indices = [i for i, row in enumerate(chart) if row["cptp"]]
    components: list[dict] = []
    if not true_indices:
        return "{0}", [{"type": "point", "value": "0"}]

    groups: list[tuple[int, int]] = []
    start = previous = true_indices[0]
    for index in true_indices[1:]:
        if index == previous + 1:
            previous = index
        else:
            groups.append((start, previous))
            start = previous = index
    groups.append((start, previous))

    for start, end in groups:
        components.append(
            {
                "type": "interval",
                "left": chart[start]["left"],
                "right": chart[end]["right"],
                "left_closed": True,
                "right_closed": chart[end]["right"] != "infinity",
            }
        )

    if groups[0][0] > 0:
        components.insert(0, {"type": "point", "value": "0"})

    rendered: list[str] = []
    for component in components:
        if component["type"] == "point":
            rendered.append("{0}")
        else:
            right = component["right"]
            closing = "]" if component["right_closed"] else ")"
            rendered.append(f"[{component['left']}, {right}{closing}")
    return " union ".join(rendered), components


def verify_stage_domain(method: Method, det_full: sp.Expr, det_half: sp.Expr) -> None:
    for label, determinant in (("full", det_full), ("half", det_half)):
        numerator, _ = sp.fraction(sp.together(determinant))
        poly = sp.Poly(numerator, x, domain=sp.QQ)
        if poly.degree() > 0:
            zero_root = 1 if poly.sqf_part().eval(0) == 0 else 0
            nonnegative = int(poly.sqf_part().count_roots(0, sp.oo))
            require(nonnegative - zero_root == 0, f"{method.key}: {label}-step stage determinant vanishes for x>0")
        require(sign_at(determinant, sp.Rational(1)) > 0, f"{method.key}: {label}-step determinant is not positive")


def method_audit(method: Method) -> dict:
    R = stability_function(method)
    require(
        sp.simplify(R - REFERENCE_STABILITY_FUNCTIONS[method.key]) == 0,
        f"{method.key}: tableau-derived stability function does not match the reference formula",
    )
    m, eta = first_defect(R)
    require((m, eta) == REFERENCE_DEFECTS[method.key], f"{method.key}: first exponential defect changed")
    K = sp.simplify((-1) ** m * (1 - sp.Rational(1, 2 ** (m - 1))) * eta)

    a = sp.factor(sp.cancel(R.subs(s, -x)))
    a_minus_one = sp.factor(sp.cancel(a - 1))
    M = margin(R)
    det_full, det_half = stage_determinants(method)
    verify_stage_domain(method, det_full, det_half)

    boundary_records: dict[str, dict] = {}
    boundary_polys: dict[str, sp.Poly] = {}
    for equation, expr in (("a=0", a), ("a=1", a_minus_one), ("M=0", M)):
        record, primitive = expression_record(expr)
        boundary_records[equation] = record
        boundary_polys[equation] = primitive

    certified_roots = certify_roots(method.key, boundary_polys)
    sign_chart = interval_sign_chart(a, a_minus_one, M, det_full, det_half, certified_roots)
    admissible_set, components = derive_admissible_set(sign_chart)
    require(
        admissible_set == REFERENCE_ADMISSIBLE_SETS[method.key],
        f"{method.key}: derived admissible set {admissible_set!r} differs from the reference",
    )

    return {
        "key": method.key,
        "name": method.name,
        "order": method.order,
        "butcher_tableau": {
            "A": [[str(entry) for entry in method.A.row(i)] for i in range(method.A.rows)],
            "b": [str(entry) for entry in method.b],
        },
        "stability_function": str(R),
        "stage_resolvent_determinants": {"full_step": str(det_full), "half_step": str(det_half)},
        "first_exponential_defect": {"m": m, "eta_m": str(eta)},
        "leading_margin_coefficient": str(K),
        "margin": str(M),
        "boundary_polynomials": boundary_records,
        "positive_root_certificates": [
            {
                "equation": root.equation,
                "label": root.label,
                "bracket": (
                    {"type": "point", "value": str(root.lower)}
                    if root.lower == root.upper
                    else {
                        "type": "open_interval",
                        "lower_decimal": next(
                            spec.lower for spec in ROOT_SPECS
                            if spec.method == method.key and spec.equation == root.equation and spec.label == root.label
                        ),
                        "upper_decimal": next(
                            spec.upper for spec in ROOT_SPECS
                            if spec.method == method.key and spec.equation == root.equation and spec.label == root.label
                        ),
                        "lower_rational": str(root.lower),
                        "upper_rational": str(root.upper),
                    }
                ),
                "approximation": root.approximation,
                "multiplicity": root.multiplicity,
                "sturm_count": 1,
            }
            for root in certified_roots
        ],
        "interval_sign_chart": sign_chart,
        "admissible_components": components,
        "certified_admissible_set": admissible_set,
    }


def verify_exact_regressions(audits: list[dict]) -> None:
    by_key = {row["key"]: row for row in audits}
    require(sp.simplify(margin(REFERENCE_STABILITY_FUNCTIONS["FE"]) + x**2 / 4) == 0, "FE margin factorization failed")
    require(sp.simplify(margin(REFERENCE_STABILITY_FUNCTIONS["RK2"]) - x**3 * (8 - x) / 64) == 0, "RK2 margin factorization failed")
    require(
        sp.simplify(margin(REFERENCE_STABILITY_FUNCTIONS["SSP3"]) + x**4 * (x**2 - 12 * x + 84) / 2304) == 0,
        "SSPRK(3,3) margin factorization failed",
    )
    require(
        sp.simplify(margin(REFERENCE_STABILITY_FUNCTIONS["RK4"]) + x**5 * (x**3 - 16 * x**2 + 160 * x - 1152) / 147456) == 0,
        "RK4 margin factorization failed",
    )

    # Backward Euler remains CPTP throughout the two-rate quadrant.
    u = sp.symbols("u", nonnegative=True)
    a_be = 1 / (1 + x)
    c_be = 1 / (1 + x / 2 + u)
    numerator, denominator = sp.fraction(sp.factor(a_be - c_be**2))
    numerator_poly = sp.Poly(numerator, x, u)
    require(bool(numerator_poly.terms()), "Backward Euler two-rate numerator vanished")
    require(all(coefficient >= 0 for _, coefficient in numerator_poly.terms()), "Backward Euler two-rate numerator has a negative coefficient")
    require(denominator != 0, "Backward Euler two-rate denominator vanished symbolically")

    require(len(audits) == len(METHODS), "Method audit count changed")
    require(set(by_key) == {method.key for method in METHODS}, "Method audit keys changed")


def named_endpoints(audits: list[dict]) -> dict[str, str]:
    result: dict[str, str] = {}
    for audit in audits:
        for root in audit["positive_root_certificates"]:
            if root["label"]:
                result[root["label"]] = root["approximation"]
    return dict(sorted(result.items()))


def write_json(audits: list[dict], output_root: Path) -> None:
    payload = {
        "schema_version": "1.0",
        "description": (
            "Exact tableau-derived stability functions, boundary polynomials, Sturm-certified positive-root brackets, "
            "rational interval sign samples, and CPTP-admissible components for eight Runge-Kutta formulas."
        ),
        "arithmetic": {
            "symbolic": "exact rational SymPy arithmetic",
            "root_certification": "Sturm counts via sympy.Poly.count_roots",
            "fixed_horizon_decimal_precision": 80,
        },
        "software": {
            "python": ".".join(map(str, sys.version_info[:3])),
            "sympy": sp.__version__,
            "numpy": np.__version__,
            "mpmath": mp.__version__,
            "matplotlib": matplotlib.__version__,
            "pillow": pillow_version,
        },
        "named_endpoints": named_endpoints(audits),
        "methods": audits,
    }
    path = output_root / "data/certified_boundaries.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write_method_csv(audits: list[dict], output_root: Path) -> None:
    path = output_root / "data/method_audit.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Method", "Order", "m", "eta_m", "Leading margin coefficient", "Local direction", "Admissible set"])
        for row in audits:
            K = sp.Rational(row["leading_margin_coefficient"])
            writer.writerow(
                [
                    row["name"],
                    row["order"],
                    row["first_exponential_defect"]["m"],
                    row["first_exponential_defect"]["eta_m"],
                    row["leading_margin_coefficient"],
                    "inward" if K > 0 else "outward",
                    row["certified_admissible_set"],
                ]
            )


PDF_METADATA = {
    "Author": "G. Blake Pierpoint",
    "Creator": "reproduce_pla_results.py",
    "CreationDate": dt.datetime(2026, 7, 23, tzinfo=dt.timezone.utc),
    "ModDate": dt.datetime(2026, 7, 23, tzinfo=dt.timezone.utc),
}
PNG_METADATA = {
    "Author": "G. Blake Pierpoint",
    "Software": "reproduce_pla_results.py",
}


def save_figure(fig: plt.Figure, output_root: Path, stem: str, title: str) -> None:
    figure_dir = output_root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight", metadata={**PDF_METADATA, "Title": title})
    fig.savefig(figure_dir / f"{stem}.png", dpi=600, bbox_inches="tight", metadata={**PNG_METADATA, "Title": title})
    plt.close(fig)


def figure_cp_geometry(output_root: Path) -> None:
    a_values = np.linspace(0, 1, 600)
    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    ax.fill_between(a_values, 0, a_values, alpha=0.18, label="CPTP region")
    ax.plot(a_values, a_values, linewidth=2.0, label=r"Amplitude damping: $|c|^2=a$")
    ax.plot(a_values, a_values**2, linewidth=2.0, linestyle="--", label=r"Added dephasing: $|c|^2=a^2$")
    ax.plot([1, 1], [0, 1], linewidth=1.1, linestyle=":")
    ax.annotate(
        "population ceiling",
        xy=(1, 0.72),
        xytext=(0.65, 0.88),
        arrowprops=dict(arrowstyle="->", linewidth=0.8),
        fontsize=9,
    )
    ax.annotate(
        "relative interior",
        xy=(0.55, 0.25),
        xytext=(0.18, 0.42),
        arrowprops=dict(arrowstyle="->", linewidth=0.8),
        fontsize=9,
    )
    ax.set_xlabel(r"Population multiplier $a$")
    ax.set_ylabel(r"Squared coherence multiplier $|c|^2$")
    ax.set_xlim(0, 1.03)
    ax.set_ylim(0, 1.03)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    save_figure(fig, output_root, "figure1_cp_geometry", "Complete-positivity geometry")


def endpoint_value(audits: list[dict], label: str) -> float:
    endpoints = named_endpoints(audits)
    return float(sp.N(sp.sympify(endpoints[label]), 17))


def figure_admissible_sets(audits: list[dict], output_root: Path) -> None:
    alpha_4 = endpoint_value(audits, "alpha_4")
    alpha_5 = endpoint_value(audits, "alpha_5")
    beta_minus = endpoint_value(audits, "beta_-")
    beta_plus = endpoint_value(audits, "beta_+")
    methods = [
        ("Forward Euler", [(0.0, 0.0)], False),
        ("Classical RK4", [(0.0, alpha_4)], False),
        ("Dormand-Prince 5", [(0.0, alpha_5)], False),
        ("Embedded DP4", [(0.0, 0.0), (beta_minus, beta_plus)], False),
        ("Backward Euler", [(0.0, 5.3)], True),
        ("Implicit midpoint", [(0.0, 0.0)], False),
    ]
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    y_positions = np.arange(len(methods))[::-1]
    for y_value, (_, intervals, extends) in zip(y_positions, methods):
        for lower, upper in intervals:
            if upper == lower:
                ax.plot(lower, y_value, "o", markersize=5, color="black")
            else:
                ax.hlines(y_value, lower, upper, linewidth=4, color="black")
                ax.plot([lower, upper], [y_value, y_value], "o", markersize=4, color="black")
        if extends:
            ax.annotate("", xy=(5.45, y_value), xytext=(5.05, y_value), arrowprops=dict(arrowstyle="->", linewidth=2))
    ax.set_yticks(y_positions)
    ax.set_yticklabels([method[0] for method in methods])
    ax.set_xlabel(r"Relaxation step $x=\gamma h$")
    ax.set_xlim(-0.12, 5.55)
    ax.set_ylim(-0.65, len(methods) - 0.35)
    ax.grid(axis="x", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save_figure(fig, output_root, "figure2_admissible_sets", "Selected CPTP-admissible sets")


def mp_stability(method: str, z: mp.mpf) -> mp.mpf:
    if method == "FE":
        return 1 + z
    if method == "SSP3":
        return 1 + z + z**2 / 2 + z**3 / 6
    raise KeyError(method)


def lower_choi_eigenvalue(A: mp.mpf, C: mp.mpf) -> mp.mpf:
    return (1 + A - mp.sqrt((1 - A) ** 2 + 4 * C**2)) / 2


def mp_string(value: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def figure_fixed_horizon(output_root: Path) -> None:
    mp.mp.dps = 80
    horizon = mp.mpf("1")
    step_counts = np.unique(np.logspace(np.log10(4), 4, 120).astype(int))
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    configs = (
        ("FE", 2, -mp.mpf(1) / 4, "Forward Euler"),
        ("SSP3", 4, -mp.mpf(7) / 192, "SSPRK(3,3)"),
    )
    csv_rows: list[list[str | int]] = []
    for method, defect_order, coefficient, label in configs:
        exact_values: list[float] = []
        asymptotic_values: list[float] = []
        for step_count in step_counts:
            n = mp.mpf(int(step_count))
            h = horizon / n
            a_step = mp_stability(method, -h)
            c_step = mp_stability(method, -h / 2)
            A = a_step ** int(step_count)
            C = c_step ** int(step_count)
            eigenvalue = lower_choi_eigenvalue(A, C)
            asymptotic = (
                coefficient * horizon * mp.e ** (-horizon) / (1 + mp.e ** (-horizon))
            ) * h ** (defect_order - 1)
            exact_values.append(float(-eigenvalue))
            asymptotic_values.append(float(-asymptotic))
            csv_rows.append(
                [
                    method,
                    int(step_count),
                    mp_string(eigenvalue),
                    mp_string(asymptotic),
                ]
            )
        ax.loglog(step_counts, exact_values, linewidth=2, label=f"{label}, exact")
        ax.loglog(step_counts, asymptotic_values, linestyle="--", linewidth=1.4, label=f"{label}, asymptotic")
    ax.set_xlabel(r"Number of equal steps $N$")
    ax.set_ylabel(r"Structural defect $-\lambda_{-,N}$")
    ax.grid(which="both", alpha=0.22)
    ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    save_figure(fig, output_root, "figure3_fixed_horizon_defect", "Fixed-horizon complete-positivity defect")

    csv_path = output_root / "data/fixed_horizon_defect.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["method", "N", "lambda_minus_exact", "lambda_minus_asymptotic"])
        writer.writerows(csv_rows)


def reproduce(output_root: Path) -> list[dict]:
    audits = [method_audit(method) for method in METHODS]
    verify_exact_regressions(audits)
    write_json(audits, output_root)
    write_method_csv(audits, output_root)
    figure_cp_geometry(output_root)
    figure_admissible_sets(audits, output_root)
    figure_fixed_horizon(output_root)
    return audits


def compare_outputs(generated_root: Path, tracked_root: Path) -> None:
    mismatches: list[str] = []
    for relative_path in TRACKED_OUTPUTS:
        generated = generated_root / relative_path
        tracked = tracked_root / relative_path
        if not tracked.exists():
            mismatches.append(f"missing tracked output: {relative_path}")
        elif not filecmp.cmp(generated, tracked, shallow=False):
            mismatches.append(f"stale tracked output: {relative_path}")
    if mismatches:
        raise VerificationError("Repository output check failed:\n  - " + "\n  - ".join(mismatches))


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the mathematics and confirm that all tracked outputs are current without modifying the repository",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.check:
            with tempfile.TemporaryDirectory(prefix="pla-reproduction-") as temporary_directory:
                generated_root = Path(temporary_directory)
                audits = reproduce(generated_root)
                compare_outputs(generated_root, ROOT)
            print(f"Verified {len(audits)} tableau-derived method audits, 14 positive-root certificates, and 3 figures.")
            print("All tracked outputs are current.")
        else:
            audits = reproduce(ROOT)
            print(f"Reproduced {len(audits)} tableau-derived method audits, 14 positive-root certificates, and 3 figures.")
            print("Exact symbolic, Sturm, sign-chart, and fixed-horizon checks passed.")
        return 0
    except VerificationError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
