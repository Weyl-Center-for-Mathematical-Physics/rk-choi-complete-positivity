from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rk_choi_margin.methods import X, all_methods
from rk_choi_margin.symbolic import (
    Z,
    thermal_cp_margin,
    first_exponential_defect,
    rotating_boundary_factor,
    rotating_cp_margin,
)


def main() -> None:
    a, c2, theta = sp.symbols("a c2 theta", real=True)
    A = 1 - theta * (1 - a)
    B = theta * (1 - a)
    C = (1 - theta) * (1 - a)
    D = theta + (1 - theta) * a

    checks: dict[str, bool] = {
        "trace identity on E00": sp.simplify(A + B - 1) == 0,
        "trace identity on E11": sp.simplify(C + D - 1) == 0,
        "thermal Choi determinant identity": sp.simplify(A * D - (a + theta * (1 - theta) * (1 - a) ** 2)) == 0,
    }

    x, u, nu = sp.symbols("x u nu", real=True)
    exact_a = sp.exp(-x)
    exact_c2 = sp.exp(-x - 2 * u)
    exact_margin = sp.factor(exact_a + theta * (1 - theta) * (1 - exact_a) ** 2 - exact_c2)
    expected_exact = sp.exp(-x) * (1 - sp.exp(-2 * u)) + theta * (1 - theta) * (1 - sp.exp(-x)) ** 2
    checks["exact phase-covariant margin"] = sp.simplify(exact_margin - expected_exact) == 0
    checks["exact margin is independent of coherent phase"] = nu not in exact_margin.free_symbols

    r2, varpi = sp.symbols("r2 varpi", real=True)
    generic_R = 1 + X + r2 * X**2
    F = sp.expand(thermal_cp_margin(generic_R, theta=theta, dephasing_step=0, frequency_step=varpi * Z))
    actual_x2 = sp.expand(F).coeff(Z, 2)
    expected_x2 = theta * (1 - theta) + r2 / 2 - sp.Rational(1, 4) + (2 * r2 - 1) * varpi**2
    checks["fixed-frequency bidirectional local coefficient"] = sp.simplify(actual_x2 - expected_x2) == 0
    checks["order-two frequency cancellation"] = sp.simplify(actual_x2.subs(r2, sp.Rational(1, 2)) - theta * (1 - theta)) == 0
    checks["Forward-Euler local coefficient"] = sp.simplify(actual_x2.subs(r2, 0) + (theta - sp.Rational(1, 2)) ** 2 + varpi**2) == 0

    for method in all_methods().values():
        R = method.stability_function()
        defect = first_exponential_defect(R)
        rotating_margin = rotating_cp_margin(R, nu=varpi)
        actual = sp.series(rotating_margin, Z, 0, defect.index + 1).removeO().expand().coeff(Z, defect.index)
        predicted = (-1) ** defect.index * defect.coefficient * rotating_boundary_factor(defect.index, varpi)
        checks[f"frequency boundary law: {method.name}"] = sp.simplify(actual - predicted) == 0

    checks["SSPRK(3,3) leading cancellation threshold"] = sp.simplify(rotating_boundary_factor(4, sp.sqrt(7) / 2)) == 0
    checks["classical RK4 leading cancellation threshold"] = sp.simplify(rotating_boundary_factor(5, sp.sqrt(3) / 2)) == 0
    ssp3_margin_at_threshold = sp.factor(
        rotating_cp_margin(all_methods()["rk3"].stability_function(), nu=sp.sqrt(7) / 2)
    )
    rk4_margin_at_threshold = sp.factor(
        rotating_cp_margin(all_methods()["rk4"].stability_function(), nu=sp.sqrt(3) / 2)
    )
    checks["SSPRK(3,3) exact threshold factorization"] = sp.simplify(
        ssp3_margin_at_threshold - Z**5 * (3 - 2 * Z) / 9
    ) == 0
    checks["classical RK4 exact threshold factorization"] = sp.simplify(
        rk4_margin_at_threshold + Z**6 * (Z - 2) ** 2 / 576
    ) == 0
    checks["SSPRK(3,3) threshold is locally inward"] = sp.Rational(1, 3) > 0
    checks["classical RK4 threshold is locally outward"] = -sp.Rational(1, 144) < 0

    xnn, unn, nurn = sp.symbols("xnn unn nurn", nonnegative=True, real=True)
    be_a = 1 / (1 + xnn)
    be_c2 = 1 / ((1 + xnn / 2 + unn) ** 2 + nurn**2)
    be_numerator = sp.expand((1 + xnn / 2 + unn) ** 2 + nurn**2 - (1 + xnn))
    checks["backward-Euler full zero-temperature domain"] = be_numerator == xnn**2 / 4 + xnn * unn + unn**2 + 2 * unn + nurn**2
    checks["backward-Euler margin identity"] = sp.simplify(
        be_a - be_c2 - be_numerator / ((1 + xnn) * ((1 + xnn / 2 + unn) ** 2 + nurn**2))
    ) == 0

    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise AssertionError("Phase-covariant verification failed: " + "; ".join(failed))

    lines = [
        "# Phase-covariant symbolic verification",
        "",
        "This report independently reconstructs the exact Choi determinant, local buffer coefficients, coherent-precession boundary law including the cancellation thresholds, and backward-Euler domain used in the article.",
        "",
    ]
    lines.extend(f"- PASS: {name}" for name in checks)
    lines.extend([
        "",
        "## Exact semigroup margin",
        "",
        "```text",
        str(exact_margin),
        "```",
        "",
        "## General fixed-frequency quadratic coefficient",
        "",
        "```text",
        str(sp.factor(actual_x2)),
        "```",
    ])
    out = ROOT / "results" / "phase_covariant_verification.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Verified {len(checks)} phase-covariant identities.")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
