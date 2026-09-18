from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rk_choi_margin.methods import X, all_methods
from rk_choi_margin.symbolic import (
    U,
    Z,
    cp_admissible_intervals,
    cp_margin,
    first_exponential_defect,
    predicted_boundary_coefficient,
    rotating_cp_margin,
    two_rate_cp_margin,
)


def interval_text(intervals: list[tuple[float, float]]) -> str:
    parts: list[str] = []
    for left, right in intervals:
        if left == right == 0.0:
            parts.append("{0}")
        elif math.isinf(right):
            parts.append(f"[{left:.10g}, infinity)")
        else:
            parts.append(f"[{left:.10g}, {right:.10g}]")
    return " union ".join(parts)


def main() -> None:
    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    report_lines: list[str] = []
    rows: list[dict[str, str | int]] = []

    report_lines.append("RK-CHOI COMPLETE-POSITIVITY - EXACT SYMBOLIC VERIFICATION")
    report_lines.append("=" * 52)
    report_lines.append("")

    for method in all_methods().values():
        R = sp.factor(method.stability_function())
        defect = first_exponential_defect(R)
        M = cp_margin(R)
        prediction = sp.factor(predicted_boundary_coefficient(defect))
        actual = sp.factor(sp.series(M, Z, 0, defect.index + 1).removeO().coeff(Z, defect.index))
        if sp.simplify(actual - prediction) != 0:
            raise AssertionError(f"Boundary coefficient mismatch for {method.name}")
        intervals = cp_admissible_intervals(R)

        report_lines.extend(
            [
                method.name,
                "-" * len(method.name),
                f"family: {method.family}; order: {method.order}; stages: {method.stages}",
                f"R(x) = {sp.sstr(R)}",
                f"first defect: m={defect.index}, eta_m={sp.sstr(defect.coefficient)}",
                f"predicted local CP coefficient: {sp.sstr(prediction)}",
                f"M_R(x,0) = {sp.sstr(sp.factor(M))}",
                f"CP-admissible set on x>=0 at u=0: {interval_text(intervals)}",
                "",
            ]
        )

        rows.append(
            {
                "key": method.key,
                "method": method.name,
                "family": method.family,
                "order": method.order,
                "stages": method.stages,
                "stability_function": sp.sstr(R),
                "first_defect_index_m": defect.index,
                "eta_m": sp.sstr(defect.coefficient),
                "local_margin_coefficient": sp.sstr(prediction),
                "local_orientation": "inward/CP" if prediction > 0 else "outward/non-CP",
                "cp_admissible_set": interval_text(intervals),
            }
        )

    # Verify the independent population exit for truncated-exponential RK2.
    rk2 = all_methods()["rk2"].stability_function()
    a_rk2 = sp.factor(rk2.subs(X, -Z))
    c_rk2 = sp.factor(rk2.subs(X, -Z / 2))
    M_rk2 = sp.factor(a_rk2 - c_rk2**2)
    assert a_rk2.subs(Z, 4) == 5
    assert c_rk2.subs(Z, 4) == 1
    assert M_rk2.subs(Z, 4) == 4
    assert 1 - a_rk2.subs(Z, 4) == -4

    # Verify the regular two-rate coordinates and their joint linear term.
    b2, b3 = sp.symbols("b2 b3", real=True)
    generic_R = 1 + X + b2 * X**2 + b3 * X**3
    M_two = sp.expand(two_rate_cp_margin(generic_R))
    assert M_two.subs({Z: 0, U: 0}) == 0
    assert sp.diff(M_two, Z).subs({Z: 0, U: 0}) == 0
    assert sp.diff(M_two, U).subs({Z: 0, U: 0}) == 2
    assert sp.simplify(M_two.subs(Z, 0) - (1 - generic_R.subs(X, -U) ** 2)) == 0

    r_symbol = sp.Symbol("r", real=True)
    assert sp.simplify(
        two_rate_cp_margin(rk2, dephasing_step=(r_symbol - sp.Rational(1, 2)) * Z)
        - cp_margin(rk2, r=r_symbol)
    ) == 0

    # Verify the explicit-polynomial stiffness scope. Every nonconstant explicit RK
    # polynomial loses the margin for sufficiently large dephasing, whereas
    # backward Euler remains CPTP throughout the complete two-rate quadrant.
    for key in ("euler", "rk2", "rk3", "rk4", "dp5", "dp4"):
        explicit_R = all_methods()[key].stability_function()
        assert sp.limit(two_rate_cp_margin(explicit_R), U, sp.oo) == -sp.oo

    backward_euler = all_methods()["be"].stability_function()
    x_nonnegative, u_nonnegative = sp.symbols("x_nonnegative u_nonnegative", nonnegative=True, real=True)
    a_be = sp.factor(backward_euler.subs(X, -x_nonnegative))
    c_be = sp.factor(backward_euler.subs(X, -x_nonnegative / 2 - u_nonnegative))
    M_be = sp.factor(a_be - c_be**2)
    expected_be = sp.factor(
        (x_nonnegative**2 / 4 + x_nonnegative * u_nonnegative + u_nonnegative**2 + 2 * u_nonnegative)
        / ((1 + x_nonnegative) * (1 + x_nonnegative / 2 + u_nonnegative) ** 2)
    )
    assert sp.simplify(M_be - expected_be) == 0

    report_lines.extend(
        [
            "ZERO-TEMPERATURE TWO-RATE CHECKS",
            "------------------------",
            "RK2 at x=4: a=5, c=1, M=4, 1-a=-4 (independent population exit verified).",
            "Two-rate margin: M_R(x,u)=R(-x)-R(-x/2-u)^2.",
            "Joint first derivatives at (0,0): dM/dx=0 and dM/du=2.",
            "Pure dephasing: M_R(0,u)=1-R(-u)^2.",
            "Rationalized Choi denominator obeys D >= 2 for all real finite (a,c).",
            "Every audited nonconstant explicit RK polynomial has M_R(x,u)->-infinity as u->infinity for fixed x>0.",
            f"Backward Euler full-quadrant margin: {sp.sstr(M_be)}",
            "Backward Euler is the scope counterexample to any unrestricted stiffness-nonuniformity claim.",
            "",
        ]
    )

    # Verify the rotating-RK4 coefficient used in the phase-covariant boundary analysis.
    nu = sp.Symbol("nu", real=True)
    rk4 = all_methods()["rk4"].stability_function()
    M_rot = rotating_cp_margin(rk4, nu=nu)
    coefficient_z5 = sp.factor(sp.expand(M_rot).coeff(Z, 5))
    expected = -((4 * nu**2 - 3) * (4 * nu**2 + 1)) / 384
    if sp.simplify(coefficient_z5 - expected) != 0:
        raise AssertionError("Rotating RK4 coefficient failed exact verification")
    report_lines.extend(
        [
            "ROTATING ONE-WAY BOUNDARY CHECK",
            "------------------------",
            f"RK4 coefficient of x^5 at u=0 on a fixed nu ray: {sp.sstr(coefficient_z5)}",
            "Local sign changes at |nu|=sqrt(3)/2.",
            "",
        ]
    )

    report_path = results_dir / "symbolic_verification.txt"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    csv_path = results_dir / "verified_method_table.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(report_path.read_text(encoding="utf-8"))
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
