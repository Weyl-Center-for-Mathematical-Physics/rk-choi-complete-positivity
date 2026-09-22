"""Exact root and sign certificates of the admissible sets on nonrotating one-way damping.

Writes results/closeout/root_sign_certificate_v2_7.json (records), supplement/root_sign_certificate_v2_7.md,
and supplement/root_sign_certificate_v2_7.tex (a compilable copy of Online Resource 1, Section 9). With
--esm-copy PATH the Online Resource section itself is written to PATH as well.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN, getcontext
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

from rk_choi_margin.methods import X, all_methods
from rk_choi_margin.symbolic import cp_margin, first_exponential_defect, predicted_boundary_coefficient

ROOT = Path(__file__).resolve().parents[1]
OUT_SUPP = ROOT / "supplement"
OUT_RESULTS = ROOT / "results" / "closeout"
OUT_SUPP.mkdir(parents=True, exist_ok=True)
OUT_RESULTS.mkdir(parents=True, exist_ok=True)

x = sp.Symbol("x", nonnegative=True, real=True)
s = X
getcontext().prec = 80


@dataclass(frozen=True)
class RootBracket:
    lo: sp.Rational
    hi: sp.Rational
    approx: str
    multiplicity: int
    exact: sp.Rational | None = None

    @property
    def is_exact(self) -> bool:
        return self.exact is not None

    def display(self) -> str:
        if self.exact is not None:
            return sp.sstr(self.exact)
        return f"({decimal_text(self.lo, digits=15, rounding=ROUND_FLOOR)}, {decimal_text(self.hi, digits=15, rounding=ROUND_CEILING)})"

    def json_obj(self) -> dict[str, Any]:
        return {
            "lo": str(self.lo),
            "hi": str(self.hi),
            "lo_decimal": decimal_text(self.lo, digits=15, rounding=ROUND_FLOOR),
            "hi_decimal": decimal_text(self.hi, digits=15, rounding=ROUND_CEILING),
            "approx": self.approx,
            "multiplicity": self.multiplicity,
            "exact": str(self.exact) if self.exact is not None else None,
        }


def decimal_text(q: sp.Rational, digits: int = 15, rounding: str = ROUND_FLOOR) -> str:
    value = Decimal(int(q.p)) / Decimal(int(q.q))
    quantum = Decimal(1).scaleb(-digits)
    value = value.quantize(quantum, rounding=rounding)
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def decimal_nearest(q: sp.Rational, digits: int = 15) -> str:
    """Render an exact rational to fixed decimal places using nearest rounding."""
    value = Decimal(int(q.p)) / Decimal(int(q.q))
    quantum = Decimal(1).scaleb(-digits)
    value = value.quantize(quantum, rounding=ROUND_HALF_EVEN)
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text



def decimal_inside(lo: sp.Rational, hi: sp.Rational, midpoint: sp.Rational, preferred_digits: int = 16) -> str:
    """Return a rounded decimal approximation certified to lie strictly inside (lo, hi)."""
    for digits in range(preferred_digits, 31):
        text = decimal_nearest(midpoint, digits=digits)
        value = sp.Rational(text)
        if lo < value < hi:
            return text
    # This should be unreachable after exact refinement, but a midpoint fraction
    # remains a valid certified approximation if decimal rounding is pathological.
    return str(midpoint)

def primitive_integer_poly(expr: sp.Expr) -> tuple[sp.Integer, sp.Poly]:
    """Return scalar and primitive integer numerator polynomial preserving sign."""
    num = sp.cancel(sp.together(expr)).as_numer_denom()[0]
    poly_q = sp.Poly(sp.expand(num), x, domain=sp.QQ)
    if poly_q.is_zero:
        return sp.Integer(0), sp.Poly(0, x, domain=sp.ZZ)
    den_lcm, poly_z = poly_q.clear_denoms(convert=True)
    content, primitive = poly_z.primitive()
    scalar = sp.Rational(content, den_lcm)
    # Keep primitive leading coefficient positive and transfer sign to scalar.
    if primitive.LC() < 0:
        primitive = -primitive
        scalar = -scalar
    return sp.simplify(scalar), primitive


def rational_roots(poly: sp.Poly) -> dict[sp.Rational, int]:
    roots = poly.ground_roots()
    return {sp.Rational(root): int(mult) for root, mult in roots.items() if sp.Rational(root) >= 0}


def isolate_nonnegative_roots(poly: sp.Poly, decimal_places: int = 15) -> list[RootBracket]:
    """Isolate every nonnegative real root with exact rational intervals.

    SymPy's real-root isolation uses exact arithmetic and returns certified
    rational endpoints and multiplicities.  This avoids the platform-sensitive
    convergence behavior of high-precision ``nroots`` while retaining narrow,
    readable brackets for the supplemental certificate.
    """
    if poly.is_zero or poly.degree() <= 0:
        return []

    eps = sp.Rational(1, 10**decimal_places)
    brackets: list[RootBracket] = []
    for (lo_raw, hi_raw), multiplicity in poly.intervals(eps=eps):
        lo = sp.Rational(lo_raw)
        hi = sp.Rational(hi_raw)
        if hi < 0:
            continue
        if lo == hi:
            if lo >= 0:
                if poly.eval(lo) != 0:
                    raise RuntimeError(f"Reported exact root {lo} is not a root of {poly.as_expr()}")
                brackets.append(
                    RootBracket(
                        lo=lo,
                        hi=hi,
                        approx=str(lo),
                        multiplicity=int(multiplicity),
                        exact=lo,
                    )
                )
            continue
        if lo < 0 < hi:
            raise RuntimeError(f"Nonzero isolated interval straddles the origin: ({lo}, {hi})")
        if hi <= 0:
            continue
        square_free = poly.sqf_part()
        if int(square_free.count_roots(lo, hi)) != 1:
            raise RuntimeError(f"Interval ({lo}, {hi}) does not isolate exactly one distinct root")
        # Refine far beyond the displayed certified bracket before reporting a
        # decimal approximation.  Earlier versions printed the midpoint of the
        # display bracket to more digits than that bracket justified.
        refined_lo, refined_hi = poly.refine_root(
            lo, hi, eps=sp.Rational(1, 10**20), check_sqf=False
        )
        refined_lo = sp.Rational(refined_lo)
        refined_hi = sp.Rational(refined_hi)
        refined_midpoint = (refined_lo + refined_hi) / 2
        brackets.append(
            RootBracket(
                lo=refined_lo,
                hi=refined_hi,
                approx=decimal_nearest(refined_midpoint, digits=15),
                multiplicity=int(multiplicity),
            )
        )
    return sorted(brackets, key=lambda r: float(r.lo))


def sign_of(expr: sp.Expr, point: sp.Rational) -> int:
    value = sp.cancel(expr).subs(x, point)
    value = sp.simplify(value)
    if value == 0:
        return 0
    sign = sp.sign(value)
    if sign in (-1, 0, 1):
        return int(sign)
    raise RuntimeError(f"Exact sign could not be certified for {value} at {point}")


def bool_cptp(a: sp.Expr, M: sp.Expr, point: sp.Rational, d1: sp.Expr, d2: sp.Expr) -> bool:
    if sp.simplify(d1.subs(x, point)) == 0 or sp.simplify(d2.subs(x, point)) == 0:
        return False
    return sign_of(a - 1, point) <= 0 and sign_of(M, point) >= 0


def format_sign(value: int) -> str:
    return "+" if value > 0 else ("-" if value < 0 else "0")


def bracket_mid(left: RootBracket | None, right: RootBracket | None) -> sp.Rational:
    if left is None and right is None:
        return sp.Rational(1)
    if left is None:
        assert right is not None
        upper = right.lo if not right.is_exact else right.exact
        assert upper is not None
        return upper / 2 if upper > 0 else sp.Rational(1, 1000)
    if right is None:
        lower = left.hi if not left.is_exact else left.exact
        assert lower is not None
        return lower + max(sp.Rational(1), lower / 10)
    lower = left.hi if not left.is_exact else left.exact
    upper = right.lo if not right.is_exact else right.exact
    assert lower is not None and upper is not None
    return (lower + upper) / 2


def merge_unique_events(root_sources: dict[str, list[RootBracket]]) -> list[dict[str, Any]]:
    raw: list[tuple[float, str, RootBracket]] = []
    for source, roots in root_sources.items():
        for root in roots:
            raw.append((float(root.exact if root.exact is not None else (root.lo + root.hi) / 2), source, root))
    raw.sort(key=lambda item: item[0])
    events: list[dict[str, Any]] = []
    for approx, source, root in raw:
        if events and abs(approx - events[-1]["approx_float"]) < 1e-9:
            events[-1]["sources"].append(source)
            events[-1]["multiplicities"][source] = root.multiplicity
            # Prefer exact bracket when available; otherwise widen union safely.
            if root.is_exact:
                events[-1]["root"] = root
            elif not events[-1]["root"].is_exact:
                current: RootBracket = events[-1]["root"]
                events[-1]["root"] = RootBracket(
                    min(current.lo, root.lo),
                    max(current.hi, root.hi),
                    current.approx,
                    max(current.multiplicity, root.multiplicity),
                )
        else:
            events.append(
                {
                    "approx_float": approx,
                    "root": root,
                    "sources": [source],
                    "multiplicities": {source: root.multiplicity},
                }
            )
    return events


def interval_label(left: RootBracket | None, right: RootBracket | None) -> str:
    if left is None:
        ltxt = "0"
    else:
        ltxt = left.approx if left.exact is None else sp.sstr(left.exact)
    if right is None:
        return f"({ltxt}, infinity)"
    rtxt = right.approx if right.exact is None else sp.sstr(right.exact)
    return f"({ltxt}, {rtxt})"


def method_endpoint_name(method_key: str, source: str, positive_index: int) -> str | None:
    mapping = {
        ("rk4", "a-1", 0): r"\alpha_4",
        ("rk4", "M", 0): r"\mu_4",
        ("dp5", "a-1", 0): r"\alpha_5",
        ("dp5", "M", 0): r"\mu_5",
        ("dp4", "M", 0): r"\beta_-",
        ("dp4", "a-1", 0): r"\beta_+",
        ("dp4", "M", 1): r"\delta_M",
        ("dp4", "a-1", 1): r"\delta_1",
        ("dp4", "a", 0): r"\delta_0",
    }
    return mapping.get((method_key, source, positive_index))


def boundary_mechanism(sources: Iterable[str]) -> str:
    source_set = set(sources)
    if "pole" in source_set:
        return "stage-resolvent singularity"
    parts: list[str] = []
    if "M" in source_set:
        parts.append("coherence-population block (M=0)")
    if "a-1" in source_set:
        parts.append("population ceiling (a=1)")
    if "a" in source_set:
        parts.append("population floor (a=0)")
    return "; ".join(parts) if parts else "identity/root bookkeeping"


def event_endpoint_tex(event: dict[str, Any]) -> str:
    if event["labels"]:
        return event["labels"][0]
    root = event["root"]
    if root["exact"] is not None:
        return sp.latex(sp.Rational(root["exact"]))
    return root["approx"]


def derive_admissible_components(sign_rows: list[dict[str, Any]], event_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], str]:
    true_indices = [index for index, row in enumerate(sign_rows) if row["cptp"]]
    components: list[dict[str, Any]] = []
    if not true_indices:
        return [{"type": "point", "value": "0"}], r"\{0\}"

    groups: list[tuple[int, int]] = []
    start = previous = true_indices[0]
    for index in true_indices[1:]:
        if index == previous + 1:
            previous = index
        else:
            groups.append((start, previous))
            start = previous = index
    groups.append((start, previous))

    if groups[0][0] > 0:
        components.append({"type": "point", "value": "0"})

    rendered: list[str] = []
    if components:
        rendered.append(r"\{0\}")
    for start, end in groups:
        left = "0" if start == 0 else event_endpoint_tex(event_rows[start - 1])
        right = r"\infty" if end == len(sign_rows) - 1 else event_endpoint_tex(event_rows[end])
        components.append({"type": "interval", "left": left, "right": right})
        if right == r"\infty":
            rendered.append(rf"[{left},{right})")
        else:
            rendered.append(rf"[{left},{right}]")
    return components, r"\cup".join(rendered)


def analyze_method(key: str) -> dict[str, Any]:
    method = all_methods()[key]
    R = sp.factor(method.stability_function(s))
    a = sp.factor(sp.cancel(R.subs(s, -x)))
    c = sp.factor(sp.cancel(R.subs(s, -x / 2)))
    M = sp.factor(sp.cancel(a - c**2))
    a_minus_1 = sp.factor(sp.cancel(a - 1))
    det = sp.factor(method.stage_resolvent_determinant(s))
    det_full = sp.factor(det.subs(s, -x))
    det_half = sp.factor(det.subs(s, -x / 2))

    expressions = {"a": a, "a-1": a_minus_1, "M": M}
    root_sources: dict[str, list[RootBracket]] = {}
    poly_data: dict[str, Any] = {}
    for name, expr in expressions.items():
        scalar, poly = primitive_integer_poly(expr)
        roots = isolate_nonnegative_roots(poly)
        root_sources[name] = roots
        poly_data[name] = {
            "expr": sp.sstr(expr),
            "numerator_scalar": sp.sstr(scalar),
            "primitive_numerator": sp.sstr(poly.as_expr()),
            "factorization": sp.sstr(sp.factor(poly.as_expr())),
            "square_free": sp.sstr(sp.sqf_part(poly.as_expr(), x)),
            "roots": [r.json_obj() for r in roots],
        }

    pole_roots: list[RootBracket] = []
    for pole_expr in (det_full, det_half):
        _, pole_poly = primitive_integer_poly(pole_expr)
        pole_roots.extend(isolate_nonnegative_roots(pole_poly))
    root_sources["pole"] = pole_roots

    events = merge_unique_events(root_sources)
    # Remove the origin from open-interval boundary list; intervals begin at 0.
    positive_events = [event for event in events if event["approx_float"] > 1e-14]
    roots_ordered = [event["root"] for event in positive_events]

    sign_rows: list[dict[str, Any]] = []
    boundaries: list[RootBracket | None] = [None] + roots_ordered + [None]
    for i in range(len(boundaries) - 1):
        left = boundaries[i]
        right = boundaries[i + 1]
        sample = bracket_mid(left, right)
        row = {
            "interval": interval_label(left, right),
            "left": left.json_obj() if left is not None else None,
            "right": right.json_obj() if right is not None else None,
            "sample": str(sample),
            "sign_a": sign_of(a, sample),
            "sign_a_minus_1": sign_of(a_minus_1, sample),
            "sign_M": sign_of(M, sample),
            "domain": bool(det_full.subs(x, sample) != 0 and det_half.subs(x, sample) != 0),
            "cptp": bool_cptp(a, M, sample, det_full, det_half),
        }
        sign_rows.append(row)

    # Annotate positive roots with labels and active boundary mechanisms.
    source_positive_indices: dict[str, int] = {"a": 0, "a-1": 0, "M": 0, "pole": 0}
    event_rows: list[dict[str, Any]] = []
    for event in positive_events:
        labels: list[str] = []
        for source in event["sources"]:
            name = method_endpoint_name(key, source, source_positive_indices[source])
            source_positive_indices[source] += 1
            if name:
                labels.append(name)
        event_rows.append(
            {
                "root": event["root"].json_obj(),
                "sources": event["sources"],
                "multiplicities": event["multiplicities"],
                "labels": labels,
                "mechanism": boundary_mechanism(event["sources"]),
            }
        )

    defect = first_exponential_defect(R)
    predicted = predicted_boundary_coefficient(defect)

    expected_sets = {
        "euler": r"\{0\}",
        "rk2": r"[0,2]",
        "rk3": r"\{0\}",
        "rk4": r"[0,\alpha_4]",
        "dp5": r"[0,\alpha_5]",
        "dp4": r"\{0\}\cup[\beta_-,\beta_+]",
        "be": r"[0,\infty)",
        "im": r"\{0\}",
    }
    components, derived_set = derive_admissible_components(sign_rows, event_rows)
    if derived_set != expected_sets[key]:
        raise RuntimeError(f"Derived admissible set {derived_set!r} disagrees with expected {expected_sets[key]!r}")

    return {
        "key": key,
        "name": method.name,
        "family": method.family,
        "order": method.order,
        "stages": method.stages,
        "R": sp.sstr(R),
        "R_display": sp.sstr(R.subs(s, sp.Symbol("s"))),
        "R_latex": sp.latex(R.subs(s, sp.Symbol("s"))),
        "a": sp.sstr(a),
        "c": sp.sstr(c),
        "M": sp.sstr(M),
        "stage_determinant": sp.sstr(det),
        "stage_full": sp.sstr(det_full),
        "stage_half": sp.sstr(det_half),
        "polynomials": poly_data,
        "positive_events": event_rows,
        "sign_rows": sign_rows,
        "admissible_components": components,
        "admissible_set_latex": derived_set,
        "defect": {
            "m": defect.index,
            "eta": sp.sstr(defect.coefficient),
            "boundary_coefficient": sp.sstr(predicted),
        },
    }


def sign_table_md(data: dict[str, Any]) -> str:
    lines = ["| Open interval | sgn a | sgn(a-1) | sgn M | RK domain | CPTP |", "|---|---:|---:|---:|:---:|:---:|"]
    for row in data["sign_rows"]:
        lines.append(
            f"| `{row['interval']}` | {format_sign(row['sign_a'])} | {format_sign(row['sign_a_minus_1'])} | {format_sign(row['sign_M'])} | "
            f"{'yes' if row['domain'] else 'no'} | {'yes' if row['cptp'] else 'no'} |"
        )
    return "\n".join(lines)


def roots_table_md(data: dict[str, Any]) -> str:
    if not data["positive_events"]:
        return "No positive real boundary roots or stage-resolvent singularities occur."
    lines = ["| Label | Certified rational bracket | Approximation | Source / active mechanism |", "|---|---|---:|---|"]
    for event in data["positive_events"]:
        root = event["root"]
        label = ", ".join(event["labels"]) if event["labels"] else "-"
        bracket = root["exact"] if root["exact"] is not None else f"({root['lo_decimal']}, {root['hi_decimal']})"
        lines.append(f"| {label} | `{bracket}` | {root['approx']} | {event['mechanism']} |")
    return "\n".join(lines)


def write_markdown(methods_data: list[dict[str, Any]]) -> str:
    lines = [
        "# Supplemental Material: Exact boundary polynomials, root isolation, and sign charts",
        "",
        "This certificate reconstructs the pure-amplitude-damping CPTP-admissible set for each audited Runge-Kutta formula from exact rational data. For every method it records the primitive numerator polynomials associated with `a(x)=R(-x)`, `a(x)-1`, and `M_R(x,0)=a(x)-R(-x/2)^2`; the operational stage-resolvent factors at the full and half step; all nonnegative real roots; rational isolating intervals; and the exact sign on every intervening open interval. Finite decimal endpoints in the brackets are exact rational numbers. Root counts were certified by Sturm sequences (`Poly.count_roots`) and signs were evaluated at exact rational sample points.",
        "",
        "The one-step map is CPTP exactly when it lies in the operational RK domain and satisfies `a(x) <= 1` and `M_R(x,0) >= 0`. The second inequality already forces `a(x) >= 0`.",
        "",
        "## S1.1 Summary",
        "",
        "| Method | Stage domain on x >= 0 | Certified admissible set |",
        "|---|---|---|",
    ]
    for data in methods_data:
        stage_text = "all x >= 0" if not any("pole" in event["sources"] for event in data["positive_events"]) else "see method chart"
        lines.append(f"| {data['name']} | {stage_text} | `${data['admissible_set_latex']}$` |")

    for idx, data in enumerate(methods_data, start=2):
        lines.extend(
            [
                "",
                f"## S1.{idx} {data['name']}",
                "",
                f"Stability function: `${data['R_display']}`.",
                f"Exact population factor: `${data['a']}`.",
                f"Exact Choi margin: `${data['M']}`.",
                "",
                f"Stage-resolvent factors on the full and half steps are `{data['stage_full']}` and `{data['stage_half']}`, respectively.",
                "",
                "Primitive boundary numerators (overall nonzero scalar factors are recorded in the machine-readable JSON):",
                "",
                f"- `P_0(x) = {data['polynomials']['a']['primitive_numerator']}` for `a(x)=0`; square-free part `{data['polynomials']['a']['square_free']}`.",
                f"- `P_1(x) = {data['polynomials']['a-1']['primitive_numerator']}` for `a(x)=1`; square-free part `{data['polynomials']['a-1']['square_free']}`.",
                f"- `P_M(x) = {data['polynomials']['M']['primitive_numerator']}` for `M_R(x,0)=0`; square-free part `{data['polynomials']['M']['square_free']}`.",
                "",
                roots_table_md(data),
                "",
                sign_table_md(data),
                "",
                f"**Certified result:** `${data['admissible_set_latex']}`. The first exponential defect is `(m, eta_m)=({data['defect']['m']}, {data['defect']['eta']})`, with predicted local margin coefficient `{data['defect']['boundary_coefficient']}`.",
            ]
        )

    lines.extend(
        [
            "",
            "## S1.10 Dormand-Prince embedded fourth-order no-omission check",
            "",
            "The embedded fourth-order formula has, in increasing positive order, the first margin zero `beta_-`, the first population-ceiling zero `beta_+`, a second margin zero `delta_M`, a second population-ceiling zero `delta_1`, and finally the population-floor zero `delta_0`. The sign chart shows: `M<0` on `(0,beta_-)`; both inequalities hold only on `(beta_-,beta_+)`; `a>1` on `(beta_+,delta_1)`; `M<0` after `delta_M`; and `a<0` after `delta_0`. Hence no later root creates another admissible component.",
            "",
            "## S1.11 Reproducibility statement",
            "",
            "The companion JSON contains exact expressions, rational brackets, multiplicities, exact rational sign samples, and interval truth values. The generator script reconstructs this document directly from the exact Butcher tableaux in `src/rk_choi_margin/methods.py`.",
        ]
    )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# LaTeX output: Online Resource 1, Section 9, and a standalone copy of it
# ---------------------------------------------------------------------------
# Method names of the article's Table 1, keyed by the method keys of rk_choi_margin.methods.
DISPLAY_NAMES = {
    "euler": "Forward Euler",
    "rk2": "Two-stage order two (Heun)",
    "rk3": "SSPRK(3,3)",
    "rk4": "Classical RK4",
    "dp5": "Dormand--Prince 5 (principal)",
    "dp4": "Dormand--Prince 4 (embedded)",
    "be": "Backward Euler",
    "im": "Implicit midpoint",
}
RICHARDSON_NAME = r"Richardson extrapolate of RK4 ($n=2$)"
RICHARDSON_HEADING = r"Richardson extrapolate of RK4 (\texorpdfstring{$n=2$}{n=2})"
RICHARDSON_SET = r"\{0\}\cup[\rho_-,\rho_+]"

# The Richardson extrapolate is not a Butcher tableau of rk_choi_margin.methods; its brackets and exact sign
# samples are recomputed and checked by scripts/verify_bit_claims.py (claim `richardson_global_sign_chart`).
RICHARDSON_BLOCK = "\n".join([
    rf"\section{{{RICHARDSON_HEADING}}}\label{{s:richardson_chart}}",
    r"Let $R_{\rm ext}$ and $P_{10}$ be the exact polynomials of Section~\ref{s:candidates}, and set $a=R_{\rm ext}(-x)$ and $M=a-R_{\rm ext}(-x/2)^2$. Exact Sturm counts give two positive zeros of $a$, one positive zero of $a-1$, and two positive zeros of $M$; the trivial root $x=0$ is excluded.",
    r"\begin{center}",
    r"\begin{tabular}{@{}lll@{}}",
    r"\toprule",
    r"Root & Certified rational bracket & Boundary \\",
    r"\midrule",
    r"$r_{0a}$ & $(2.919911381617275,2.919911381617276)$ & population floor, $a=0$ \\",
    r"$r_{0b}$ & $(6.034383873449079,6.034383873449080)$ & population floor, $a=0$ \\",
    r"$\rho_-$ & $(6.034523958649170,6.034523958649171)$ & margin, $M=0$ \\",
    r"$\rho_+$ & $(6.459127767825720,6.459127767825721)$ & population ceiling, $a=1$ \\",
    r"$r_M$ & $(22.102437082855815,22.102437082855816)$ & margin, $M=0$ \\",
    r"\bottomrule",
    r"\end{tabular}",
    r"\end{center}",
    r"Each bracket contains one simple root of its indicated polynomial, and these are all positive boundary roots. Exact signs at the rational samples $1,4,6.03445,6.2,10,23$, respectively, give the complete chart:",
    r"\begin{center}",
    r"\begin{tabular}{@{}lcccc@{}}",
    r"\toprule",
    r"Open interval & $\operatorname{sgn}a$ & $\operatorname{sgn}(a-1)$ & $\operatorname{sgn}M$ & CPTP \\",
    r"\midrule",
    r"$(0,r_{0a})$ & $+$ & $-$ & $-$ & no \\",
    r"$(r_{0a},r_{0b})$ & $-$ & $-$ & $-$ & no \\",
    r"$(r_{0b},\rho_-)$ & $+$ & $-$ & $-$ & no \\",
    r"$(\rho_-,\rho_+)$ & $+$ & $-$ & $+$ & yes \\",
    r"$(\rho_+,r_M)$ & $+$ & $+$ & $+$ & no \\",
    r"$(r_M,\infty)$ & $+$ & $+$ & $-$ & no \\",
    r"\bottomrule",
    r"\end{tabular}",
    r"\end{center}",
    r"\textbf{Completeness of the sign chart.} The closed inequalities and continuity include both $\rho_-$ and $\rho_+$; at $r_{0a}$ and $r_{0b}$ the margin is strictly negative, and at $r_M$ the population ceiling is strictly violated. The chart is therefore complete, and the admissible set is exactly $\mathcal A_{\rm ext}=\{0\}\cup[\rho_-,\rho_+]$. First defect $(m,\eta_m)=(6,-1/4320)$; orientation coefficient $-31/138240$.",
])

# Approximate widths in points of math glyphs at 9 pt, used only to break long polynomials across lines.
POLY_LINE_BUDGET = 380.0


def tex_sign(value: int) -> str:
    return "$+$" if value > 0 else ("$-$" if value < 0 else "$0$")


def boundary_tex(sources: Iterable[str]) -> str:
    source_set = set(sources)
    if "pole" in source_set:
        return "stage singularity"
    parts: list[str] = []
    if "a" in source_set:
        parts.append(r"population floor, $a=0$")
    if "a-1" in source_set:
        parts.append(r"population ceiling, $a=1$")
    if "M" in source_set:
        parts.append(r"margin, $M_R=0$")
    return "; ".join(parts)


# Typeset names of positive boundary roots that carry no label in the certificate records: r_0 for a zero of a,
# r_1 for a zero of a - 1, r_M for a zero of M_R (as in the Richardson block), with a, b, ... when a method has
# several.  They exist only in the LaTeX output; the JSON records are unchanged.
FALLBACK_ROOT_NAMES = {"a": "0", "a-1": "1", "M": "M", "pole": "p"}


def root_labels_tex(data: dict[str, Any]) -> list[str]:
    """One label per positive boundary root: the recorded label if there is one, else r_0, r_1, r_M, ..."""
    events = data["positive_events"]
    bases = [None if event["labels"] else FALLBACK_ROOT_NAMES[event["sources"][0]] for event in events]
    seen: dict[str, int] = {}
    labels = []
    for event, base in zip(events, bases):
        if base is None:
            labels.append(", ".join(event["labels"]))
            continue
        if bases.count(base) == 1:
            labels.append(rf"r_{{{base}}}" if len(base) > 1 else rf"r_{base}")
        else:
            seen[base] = seen.get(base, 0) + 1
            labels.append(rf"r_{{{base}{'abcdefgh'[seen[base] - 1]}}}")
    return labels


def event_short_tex(event: dict[str, Any], label: str) -> str:
    """Name of a positive boundary root in a sign chart: its recorded label, else its exact value, else its typeset label."""
    if event["labels"]:
        return event["labels"][0]
    root = event["root"]
    return sp.latex(sp.Rational(root["exact"])) if root["exact"] is not None else label


def _as_x_expr(text: str) -> sp.Expr:
    expr = sp.sympify(text)
    for symbol in expr.free_symbols:
        if symbol.name == "x" and symbol != x:
            expr = expr.subs(symbol, x)
    return expr


def polynomial_lines_tex(expr: sp.Expr) -> list[str]:
    """Terms of an integer polynomial in decreasing degree, grouped into lines of bounded printed width."""
    poly = sp.Poly(sp.expand(expr), x, domain=sp.ZZ)
    lines: list[str] = []
    text, width = "", 0.0
    for index, ((power,), coefficient) in enumerate(poly.terms()):
        magnitude = abs(int(coefficient))
        term = sp.latex(magnitude * x**power)
        sign = "-" if coefficient < 0 else ("" if index == 0 else "+")
        term_width = (9.5 if index else (4.5 if sign else 0.0))
        term_width += 4.5 * len(str(magnitude)) if (magnitude != 1 or power == 0) else 0.0
        term_width += (5.2 if power > 0 else 0.0) + (3.2 * len(str(power)) if power > 1 else 0.0)
        if text and width + term_width > POLY_LINE_BUDGET:
            lines.append(text)
            text, width = "", 0.0
        text += sign + term
        width += term_width
    lines.append(text)
    return lines


def aligned_polynomial(lhs: str, expr: sp.Expr, end: str) -> list[str]:
    chunks = polynomial_lines_tex(expr)
    rows = [rf"{lhs}&={chunks[0]}"] + [rf"&\quad{chunk}" for chunk in chunks[1:]]
    rows[-1] += end
    return rows


def margin_display_tex(data: dict[str, Any]) -> str:
    """M_R(x,0,0) in factored form, or as a rational multiple of P_M(x) when the factored form is very long."""
    M_expr = _as_x_expr(data["M"])
    factored = sp.latex(sp.factor(M_expr))
    if len(factored) <= 160:
        return factored
    ratio = sp.cancel(M_expr / _as_x_expr(data["polynomials"]["M"]["primitive_numerator"]))
    if ratio.free_symbols or not ratio.is_Rational or ratio.p not in (1, -1):
        return factored
    return ("-" if ratio < 0 else "") + rf"\frac{{P_M(x)}}{{{ratio.q}}}"


def method_equations_tex(data: dict[str, Any]) -> str:
    rows = [
        rf"R(s)&={data['R_latex']},",
        rf"a(x)&={sp.latex(sp.factor(_as_x_expr(data['a'])))},",
        rf"M_R(x,0,0)&={margin_display_tex(data)},",
    ]
    numerators = [("a", "P_0"), ("a-1", "P_1"), ("M", "P_M")]
    entries: list[tuple[str, sp.Expr]] = []
    for source, name in numerators:
        entries.append((rf"{name}(x)", _as_x_expr(data["polynomials"][source]["primitive_numerator"])))
    for source, name in numerators:
        primitive = _as_x_expr(data["polynomials"][source]["primitive_numerator"])
        square_free = _as_x_expr(data["polynomials"][source]["square_free"])
        if sp.expand(primitive - square_free) != 0:
            entries.append((rf"\operatorname{{sqf}}{name}(x)", square_free))
    for index, (lhs, expr) in enumerate(entries):
        rows.extend(aligned_polynomial(lhs, expr, "." if index == len(entries) - 1 else ","))
    return "\\begin{align*}\n" + "\\\\\n".join(rows) + "\n\\end{align*}"


def roots_table_tex(data: dict[str, Any]) -> str:
    events = data["positive_events"]
    if not events:
        return "There is no positive boundary root."
    rows = []
    for event, label in zip(events, root_labels_tex(data)):
        root = event["root"]
        if root["exact"] is not None:
            bracket = f"${sp.latex(sp.Rational(root['exact']))}$ (exact)"
        else:
            bracket = f"$({root['lo_decimal']},{root['hi_decimal']})$"
        cells = [f"${label}$", bracket, f"${root['approx']}$", boundary_tex(event["sources"])]
        rows.append(" & ".join(cells) + r" \\")
    header = r"Root & Certified rational bracket & Approx. & Boundary \\"
    return "\n".join([r"\begin{center}", r"\begin{tabular}{@{}llll@{}}", r"\toprule", header, r"\midrule",
                      *rows, r"\bottomrule", r"\end{tabular}", r"\end{center}"])


def sign_table_tex(data: dict[str, Any]) -> str:
    """Sign chart; the stage-domain column is printed only when some interval leaves the stage domain."""
    show_domain = not all(row["domain"] for row in data["sign_rows"])
    names = ["0"] + [event_short_tex(event, label) for event, label in zip(data["positive_events"], root_labels_tex(data))] + [r"\infty"]
    rows = []
    for i, row in enumerate(data["sign_rows"]):
        cells = [f"$({names[i]},{names[i + 1]})$", tex_sign(row["sign_a"]), tex_sign(row["sign_a_minus_1"]), tex_sign(row["sign_M"])]
        if show_domain:
            cells.append("yes" if row["domain"] else "no")
        cells.append("yes" if row["cptp"] else "no")
        rows.append(" & ".join(cells) + r" \\")
    header = r"Open interval & $\operatorname{sgn}a$ & $\operatorname{sgn}(a-1)$ & $\operatorname{sgn}M_R$ & " + ("Stage domain & " if show_domain else "") + r"CPTP \\"
    spec = "lcccc" + ("c" if show_domain else "")
    return "\n".join([r"\begin{center}", rf"\begin{{tabular}}{{@{{}}{spec}@{{}}}}", r"\toprule", header, r"\midrule",
                      *rows, r"\bottomrule", r"\end{tabular}", r"\end{center}"])


def dp4_completeness_tex(data: dict[str, Any]) -> str:
    """Completeness statement for the embedded Dormand--Prince chart, checked against its exact sign rows."""
    labels = [event["labels"] for event in data["positive_events"]]
    rows = data["sign_rows"]
    ok = labels == [[r"\beta_-"], [r"\beta_+"], [r"\delta_M"], [r"\delta_1"], [r"\delta_0"]] and len(rows) == 6
    ok = ok and rows[0]["sign_M"] < 0 and rows[1]["cptp"] and rows[2]["sign_a_minus_1"] > 0 and rows[3]["sign_a_minus_1"] > 0
    ok = ok and all(row["sign_M"] < 0 for row in rows[3:]) and not any(row["cptp"] for row in rows[2:])
    if not ok:
        raise RuntimeError("the embedded Dormand--Prince sign chart no longer matches its completeness statement")
    return (r"\textbf{Completeness of the sign chart.} "
            r"The margin is negative on $(0,\beta_-)$, the population ceiling fails on $(\beta_+,\delta_M]$, and the margin is "
            r"negative on $(\delta_M,\infty)$, which contains $\delta_1$ and $\delta_0$. The inequalities $a\le1$ and $M_R\ge0$ "
            r"therefore hold together only on $[\beta_-,\beta_+]$, and the chart is complete.")


def method_block_tex(data: dict[str, Any]) -> str:
    lines = [rf"\section{{{DISPLAY_NAMES[data['key']]}}}"]
    if data["stage_full"] != "1" or data["stage_half"] != "1":
        full = sp.latex(sp.sympify(data["stage_full"]), fold_short_frac=True)
        half = sp.latex(sp.sympify(data["stage_half"]), fold_short_frac=True)
        lines.append(rf"Stage determinants: ${full}$ (at $-x$) and ${half}$ (at $-x/2$).")
    lines.extend([method_equations_tex(data), roots_table_tex(data), sign_table_tex(data)])
    closing = []
    if data["key"] == "dp4":
        closing.append(dp4_completeness_tex(data))
    defect = data["defect"]
    closing.append(rf"Admissible set: $\mathcal A_R={data['admissible_set_latex']}$. First defect $(m,\eta_m)=({defect['m']},{defect['eta']})$; "
                   rf"orientation coefficient ${defect['boundary_coefficient']}$.")
    lines.append(" ".join(closing))
    return "\n".join(lines)


def esm_certificate_tex(methods_data: list[dict[str, Any]]) -> str:
    """Section 9 of Online Resource 1 (manuscript/nine_candidate_certificate.tex), from the certificate records.

    The fragment depends only on JSON-serializable fields, so it can be rebuilt from
    results/closeout/root_sign_certificate_v2_7.json alone."""
    by_key = {data["key"]: data for data in methods_data}
    order = ["euler", "rk2", "rk3", "rk4", "richardson", "dp5", "dp4", "be", "im"]
    summary = []
    for key in order:
        if key == "richardson":
            summary.append(rf"{RICHARDSON_NAME} & ${RICHARDSON_SET}$ \\")
        else:
            summary.append(rf"{DISPLAY_NAMES[key]} & ${by_key[key]['admissible_set_latex']}$ \\")
    body = [r"\section{Summary}", r"\begin{center}", r"\begin{tabular}{@{}ll@{}}", r"\toprule",
            r"Method & Admissible set on nonrotating one-way damping \\", r"\midrule", *summary, r"\bottomrule",
            r"\end{tabular}", r"\end{center}"]
    for key in order:
        body.append(RICHARDSON_BLOCK if key == "richardson" else method_block_tex(by_key[key]))
    return "\n".join(body) + "\n"


STANDALONE_INTRO = (
    r"This is Section~9 of Online Resource~1 of the article; references to Section~5 are to that document. "
    r"For each method it lists the stability function, the population multiplier $a(x)=R(-x)$ and the margin $M_R(x,0,0)$, "
    r"the numerators $P_0$, $P_1$, $P_M$ whose zeros are those of $a$, $a-1$, and $M_R$ (each $P$ is the primitive integer "
    r"numerator with positive leading coefficient, so its sign can differ from that of $a$, $a-1$, or $M_R$; the sign charts give "
    r"the signs of $a$, $a-1$, and $M_R$ themselves), their square-free parts where these differ, a labelled rational bracket for "
    r"every positive root, and the sign chart. For the explicit methods the stage equations are solvable for every $x$; for "
    r"backward Euler and the implicit midpoint rule the stage determinants $\det(I+xA)$ and $\det(I+xA/2)$ at the population and "
    r"coherence arguments $-x$ and $-x/2$ are positive for $x\ge0$. Decimal bracket endpoints are exact rational numbers, and "
    r"every sign is evaluated exactly at a rational sample point."
)


def standalone_certificate_tex(fragment: str) -> str:
    """A compilable copy of the Online Resource section, kept with the code and data."""
    preamble = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage{lmodern}",
        r"\usepackage{amsmath,amssymb}",
        r"\usepackage{booktabs}",
        r"\usepackage[a4paper,margin=17mm]{geometry}",
        r"\usepackage{microtype}",
        r"\providecommand{\texorpdfstring}[2]{#1}",
        r"\makeatletter\@namedef{r@s:candidates}{{5}{}}\makeatother% Section 5 of Online Resource 1",
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\parskip}{5pt}",
        r"\setlength{\emergencystretch}{3em}",
        r"\title{Exact boundary polynomials, root brackets, and sign charts on nonrotating one-way damping\\[4pt]"
        r"\large Online Resource~1, Section~9, of ``Complete-positivity regions of Runge--Kutta discretizations of phase-covariant qubit dynamics''}",
        r"\author{G. Blake Pierpoint, Olivier Bernard, and Yichen Liu}",
        r"\date{}",
        r"\begin{document}",
        r"\maketitle",
        STANDALONE_INTRO,
        "",
    ]
    return "\n".join(preamble) + "\n" + fragment + r"\end{document}" + "\n"


def write_lf(path: Path, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def main(argv: list[str] | None = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--esm-copy", type=Path, default=None,
                        help="also write the Online Resource section to this path (e.g. ../manuscript/nine_candidate_certificate.tex)")
    args = parser.parse_args(argv)

    keys = ["euler", "rk2", "rk3", "rk4", "dp5", "dp4", "be", "im"]
    methods_data = [analyze_method(key) for key in keys]

    json_path = OUT_RESULTS / "root_sign_certificate_v2_7.json"
    md_path = OUT_SUPP / "root_sign_certificate_v2_7.md"
    tex_path = OUT_SUPP / "root_sign_certificate_v2_7.tex"

    json_path.write_text(json.dumps({"version": "2.3", "certification": "exact rational arithmetic with Sturm root counts and exact sign samples", "methods": methods_data}, indent=2), encoding="utf-8")
    md_path.write_text(write_markdown(methods_data), encoding="utf-8")
    # The LaTeX is rendered from the records as written to the JSON file, so a copy can be rebuilt from it alone.
    fragment = esm_certificate_tex(json.loads(json_path.read_text(encoding="utf-8"))["methods"])
    write_lf(tex_path, standalone_certificate_tex(fragment))

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    print(f"Wrote {tex_path}")
    if args.esm_copy is not None:
        write_lf(args.esm_copy, fragment)
        print(f"Wrote {args.esm_copy}")
    for data in methods_data:
        print(data["key"], data["admissible_set_latex"])
        for event in data["positive_events"]:
            print(" ", event["labels"], event["root"]["approx"], event["sources"])


if __name__ == "__main__":
    main()
