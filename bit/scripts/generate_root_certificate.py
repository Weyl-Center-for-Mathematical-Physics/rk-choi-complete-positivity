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


def latex_interval_label(left: RootBracket | None, right: RootBracket | None) -> str:
    if left is None:
        ltxt = "0"
    else:
        ltxt = root_short_tex(left)
    if right is None:
        return rf"({ltxt},\infty)"
    return rf"({ltxt},{root_short_tex(right)})"


def root_short_tex(root: RootBracket) -> str:
    if root.exact is not None:
        return sp.latex(root.exact)
    return root.approx


def expr_tex(expr: sp.Expr) -> str:
    return sp.latex(sp.factor(expr))


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


def polynomial_align_tex(expr: sp.Expr, lhs: str, *, terms_per_line: int = 3) -> str:
    expr = sp.sympify(expr)
    for symbol in expr.free_symbols:
        if symbol.name == "x" and symbol != x:
            expr = expr.subs(symbol, x)
    poly = sp.Poly(sp.expand(expr), x, domain=sp.ZZ)
    terms = poly.terms()
    chunks = [terms[i:i + terms_per_line] for i in range(0, len(terms), terms_per_line)]
    lines: list[str] = []
    for chunk_index, chunk in enumerate(chunks):
        text = ""
        for term_index, ((power,), coefficient) in enumerate(chunk):
            term = sp.latex(abs(coefficient) * x**power)
            if chunk_index == 0 and term_index == 0:
                text += ("-" if coefficient < 0 else "") + term
            else:
                text += ("-" if coefficient < 0 else "+") + term
        prefix = rf"{lhs}={{}}&" if chunk_index == 0 else r"&"
        suffix = r"\\" if chunk_index < len(chunks) - 1 else ""
        lines.append(prefix + text + suffix)
    return "\\begin{aligned}\n" + "\n".join(lines) + "\n\\end{aligned}"


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


def tex_escape(text: str) -> str:
    return text.replace("&", r"\&").replace("_", r"\_")


def roots_table_tex(data: dict[str, Any]) -> str:
    if not data["positive_events"]:
        return "No positive real boundary roots or stage-resolvent singularities occur."
    rows = []
    for event in data["positive_events"]:
        root = event["root"]
        label = ", ".join(event["labels"]) if event["labels"] else "--"
        if root["exact"] is not None:
            bracket = f"${sp.latex(sp.Rational(root['exact']))}$"
        else:
            bracket = f"$({root['lo_decimal']},{root['hi_decimal']})$"
        rows.append(
            f"${label}$ & {bracket} & ${root['approx']}$ & {tex_escape(event['mechanism'])} \\\\"
        )
    return "\n".join(
        [
            r"\begin{center}",
            r"\begin{tabularx}{\textwidth}{@{}l l l X@{}}",
            r"\toprule",
            r"Label & Certified rational bracket & Approx. & Source / active mechanism \\",
            r"\midrule",
            *rows,
            r"\bottomrule",
            r"\end{tabularx}",
            r"\end{center}",
        ]
    )


def sign_table_tex(data: dict[str, Any]) -> str:
    rows = []
    positive_roots = [
        RootBracket(
            sp.Rational(row["root"]["lo"]),
            sp.Rational(row["root"]["hi"]),
            row["root"]["approx"],
            row["root"]["multiplicity"],
            sp.Rational(row["root"]["exact"]) if row["root"]["exact"] is not None else None,
        )
        for row in data["positive_events"]
    ]
    boundaries: list[RootBracket | None] = [None] + positive_roots + [None]
    for i, row in enumerate(data["sign_rows"]):
        interval = latex_interval_label(boundaries[i], boundaries[i + 1])
        rows.append(
            f"${interval}$ & {format_sign(row['sign_a'])} & {format_sign(row['sign_a_minus_1'])} & {format_sign(row['sign_M'])} & "
            f"{'yes' if row['domain'] else 'no'} & {'yes' if row['cptp'] else 'no'} \\\\"
        )
    return "\n".join(
        [
            r"\begin{center}",
            r"\begin{tabular}{@{}l c c c c c@{}}",
            r"\toprule",
            r"Open interval & $\operatorname{sgn}a$ & $\operatorname{sgn}(a-1)$ & $\operatorname{sgn}M$ & domain & CPTP \\",
            r"\midrule",
            *rows,
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{center}",
        ]
    )


def write_tex(methods_data: list[dict[str, Any]]) -> str:
    summary_rows = []
    for data in methods_data:
        summary_rows.append(
            f"{tex_escape(data['name'])} & all $x\\geq0$ & ${data['admissible_set_latex']}$ \\\\"
        )

    body: list[str] = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage{lmodern}",
        r"\usepackage{amsmath,amssymb,mathtools}",
        r"\usepackage{booktabs,tabularx,array,longtable}",
        r"\usepackage{adjustbox}",
        r"\usepackage[letterpaper,margin=0.72in]{geometry}",
        r"\usepackage{microtype}",
        r"\usepackage[hidelinks]{hyperref}",
        r"\hypersetup{pdftitle={Supplemental Material: Certified complete-positivity stability regions and adaptive guards for Runge--Kutta propagation of phase-covariant qubits},pdfauthor={G. Blake Pierpoint; Olivier Bernard; Yichen Liu},pdfsubject={Exact boundary polynomials, Sturm-certified root isolation, sign charts, and no-omission certificate},pdfkeywords={complete positivity, Runge--Kutta methods, Choi matrix, root isolation, sign charts}}",
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\parskip}{5pt}",
        r"\setlength{\emergencystretch}{3em}",
        r"\title{Supplemental Material for\\[3pt]\large ``Certified complete-positivity stability regions and adaptive guards for Runge--Kutta propagation of phase-covariant qubits''}",
        r"\author{G. Blake Pierpoint$^{1,*}$, Olivier Bernard$^{3,2}$, and Yichen Liu$^{4,2}$\\[6pt]",
        r"\small $^1$Old Dominion University, 5115 Hampton Boulevard,\\",
        r"\small Norfolk, Virginia 23529, USA\\",
        r"\small $^2$Weyl Center for Mathematical Physics, Washington, DC, USA\\",
        r"\small $^3$UFR PhITEM, Universit\'e Grenoble Alpes, Saint-Martin-d'H\`eres, France\\",
        r"\small $^4$Department of Computer Science and Engineering, Shanghai Jiao Tong University, Shanghai, China\\[3pt]",
        r"\small $^*$pierpogb@odu.edu}",
        r"\date{August 2026}",
        r"\begin{document}",
        r"\maketitle",
        r"\section*{Scope of this certificate}",
        r"This document reconstructs every pure-amplitude-damping CPTP-admissible set reported in the article from exact rational Butcher data. It records the boundary numerators for $a=0$, $a=1$, and $M_R=0$; the full- and half-step stage-resolvent factors; rational isolating intervals for every nonnegative real root; and the sign on every intervening interval. Root counts use Sturm sequences, and interval signs use exact rational sample points. Finite decimal bracket endpoints below denote exact rational numbers. A direct one-step map is CPTP exactly when the RK stage equations are defined and $a(x)\leq1$ and $M_R(x,0)\geq0$; the latter already implies $a(x)\geq0$.",
        r"The companion software archive contains the machine-readable JSON certificate, exact sample points and multiplicities, source code, regression tests, and all plotting data. No reported endpoint is inferred from visual sampling; every endpoint is produced by or checked against exact real-root isolation.",
        r"\section{Summary}",
        r"\begin{center}",
        r"\begin{tabularx}{\textwidth}{@{}X l l@{}}",
        r"\toprule",
        r"Method & Stage domain on $x\geq0$ & Certified admissible set \\",
        r"\midrule",
        *summary_rows,
        r"\bottomrule",
        r"\end{tabularx}",
        r"\end{center}",
    ]

    for idx, data in enumerate(methods_data, start=2):
        body.extend(
            [
                rf"\section{{{tex_escape(data['name'])}}}",
                r"The exact stability function is",
                r"\[\begin{adjustbox}{max width=0.98\textwidth}$\displaystyle R(s)=" + data['R_latex'] + r"$\end{adjustbox}\]",
                rf"The full- and half-step stage-resolvent factors are ${sp.latex(sp.sympify(data['stage_full']))}$ and ${sp.latex(sp.sympify(data['stage_half']))}$, respectively.",
                r"The exact scalar functions are",
                r"\[\begin{adjustbox}{max width=0.98\textwidth}$\displaystyle a(x)=" + sp.latex(sp.factor(sp.sympify(data['a']))) + r"$\end{adjustbox}\]",
                r"\[\begin{adjustbox}{max width=0.98\textwidth}$\displaystyle M_R(x,0,0)=" + sp.latex(sp.factor(sp.sympify(data['M']))) + r"$\end{adjustbox}\]",
                r"The primitive boundary numerators below isolate the zeros of $a$, $a-1$, and $M_R$; their signs are retained in the exact expressions above and in every interval row.",
                r"\[" + polynomial_align_tex(sp.sympify(data['polynomials']['a']['primitive_numerator']), r"P_0(x)") + r"\]",
                r"\[" + polynomial_align_tex(sp.sympify(data['polynomials']['a-1']['primitive_numerator']), r"P_1(x)") + r"\]",
                r"\[" + polynomial_align_tex(sp.sympify(data['polynomials']['M']['primitive_numerator']), r"P_M(x)") + r"\]",
                r"These equations correspond respectively to $a(x)=0$, $a(x)=1$, and $M_R(x,0,0)=0$. Their square-free parts are",
                r"\[" + polynomial_align_tex(sp.sympify(data['polynomials']['a']['square_free']), r"\operatorname{sqf}P_0") + r"\]",
                r"\[" + polynomial_align_tex(sp.sympify(data['polynomials']['a-1']['square_free']), r"\operatorname{sqf}P_1") + r"\]",
                r"\[" + polynomial_align_tex(sp.sympify(data['polynomials']['M']['square_free']), r"\operatorname{sqf}P_M") + r"\]",
                roots_table_tex(data),
                sign_table_tex(data),
                rf"The certified admissible set is $\mathcal A_R={data['admissible_set_latex']}$. The first exponential defect is $(m,\eta_m)=({data['defect']['m']},{sp.latex(sp.sympify(data['defect']['eta']))})$, and the signed boundary-law coefficient is ${sp.latex(sp.sympify(data['defect']['boundary_coefficient']))}$.",
            ]
        )

    body.extend(
        [
            r"\section{Embedded Dormand--Prince no-omission check}",
            r"For the embedded fourth-order formula, the ordered positive events are the first margin zero $\beta_-$, the first population-ceiling zero $\beta_+$, the second margin zero $\delta_M$, the second population-ceiling zero $\delta_1$, and the population-floor zero $\delta_0$. The complete sign chart shows that both CPTP inequalities hold only between $\beta_-$ and $\beta_+$. Population overflow excludes the region after $\beta_+$, the margin is again negative after $\delta_M$, and the population multiplier becomes negative after $\delta_0$. Consequently no later root generates an omitted admissible component.",
            r"\section{Reproducibility}",
            r"The companion JSON records exact expressions, rational brackets, root multiplicities, exact rational sample points, and interval truth values. This PDF is generated directly from the exact tableaux in \texttt{src/rk\_choi\_margin/methods.py}; every displayed endpoint is produced by or checked against exact real-root isolation rather than visual sampling.",
            r"\end{document}",
        ]
    )
    return "\n".join(body) + "\n"


def main() -> None:
    keys = ["euler", "rk2", "rk3", "rk4", "dp5", "dp4", "be", "im"]
    methods_data = [analyze_method(key) for key in keys]

    json_path = OUT_RESULTS / "root_sign_certificate_v2_7.json"
    md_path = OUT_SUPP / "root_sign_certificate_v2_7.md"
    tex_path = OUT_SUPP / "root_sign_certificate_v2_7.tex"

    json_path.write_text(json.dumps({"version": "2.3", "certification": "exact rational arithmetic with Sturm root counts and exact sign samples", "methods": methods_data}, indent=2), encoding="utf-8")
    md_path.write_text(write_markdown(methods_data), encoding="utf-8")
    tex_path.write_text(write_tex(methods_data), encoding="utf-8")

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    print(f"Wrote {tex_path}")
    for data in methods_data:
        print(data["key"], data["admissible_set_latex"])
        for event in data["positive_events"]:
            print(" ", event["labels"], event["root"]["approx"], event["sources"])


if __name__ == "__main__":
    main()
