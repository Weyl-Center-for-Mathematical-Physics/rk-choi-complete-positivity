from __future__ import annotations

from pathlib import Path

import sympy as sp

from rk_choi_margin.methods import X, all_methods
from rk_choi_margin.symbolic import U, Z, two_rate_cp_margin


def main() -> None:
    a, c, lam = sp.symbols("a c lam", real=True)
    J = sp.Matrix(
        [
            [1, 0, 0, c],
            [0, 0, 0, 0],
            [0, 0, 1 - a, 0],
            [c, 0, 0, a],
        ]
    )
    delta = (1 - a) ** 2 + 4 * c**2
    D = 1 + a + sp.sqrt(delta)
    margin = a - c**2

    checks: dict[str, bool] = {}
    charpoly = J.charpoly()
    spectral_parameter = charpoly.gen
    expected_charpoly = sp.expand(
        spectral_parameter
        * (spectral_parameter - (1 - a))
        * ((spectral_parameter - 1) * (spectral_parameter - a) - c**2)
    )
    checks["Choi characteristic polynomial factors into the scalar and block modes"] = (
        sp.expand(charpoly.as_expr() - expected_charpoly) == 0
    )
    checks["rationalized lower-eigenvalue identity is exact"] = (
        sp.simplify((1 + a - sp.sqrt(delta)) * D - 4 * margin) == 0
    )
    checks["block discriminant vanishes at the stated cusp point"] = delta.subs({a: 1, c: 0}) == 0

    methods = all_methods()
    rk2 = methods["rk2"].stability_function()
    a2 = sp.simplify(rk2.subs(X, -Z))
    c2 = sp.simplify(rk2.subs(X, -Z / 2))
    m2 = sp.factor(a2 - c2**2)
    checks["RK2 population-exit counterexample"] = all(
        [
            sp.simplify(m2 - Z**3 * (8 - Z) / 64) == 0,
            a2.subs(Z, 4) == 5,
            c2.subs(Z, 4) == 1,
            m2.subs(Z, 4) == 4,
            1 - a2.subs(Z, 4) == -4,
        ]
    )

    b2, b3 = sp.symbols("b2 b3", real=True)
    generic_R = 1 + X + b2 * X**2 + b3 * X**3
    generic_margin = sp.expand(two_rate_cp_margin(generic_R))
    checks["two-rate local gradient is (0,2)"] = all(
        [
            generic_margin.subs({Z: 0, U: 0}) == 0,
            sp.diff(generic_margin, Z).subs({Z: 0, U: 0}) == 0,
            sp.diff(generic_margin, U).subs({Z: 0, U: 0}) == 2,
        ]
    )
    checks["pure-dephasing slice is one minus the squared stability factor"] = (
        sp.simplify(generic_margin.subs(Z, 0) - (1 - generic_R.subs(X, -U) ** 2)) == 0
    )

    be = methods["be"].stability_function()
    be_margin = sp.factor(two_rate_cp_margin(be))
    expected_be_margin = sp.factor(
        (Z**2 / 4 + Z * U + U**2 + 2 * U)
        / ((1 + Z) * (1 + Z / 2 + U) ** 2)
    )
    checks["backward Euler margin is nonnegative on the full quadrant"] = (
        sp.simplify(be_margin - expected_be_margin) == 0
    )

    eta = sp.Symbol("eta", nonzero=True, real=True)
    for m in range(2, 9):
        Rm = sum(X**k / sp.factorial(k) for k in range(m)) + (
            sp.Rational(1, sp.factorial(m)) + eta
        ) * X**m
        Mm = sp.expand(two_rate_cp_margin(Rm, dephasing_step=0))
        actual = sp.simplify(Mm.coeff(Z, m))
        predicted = sp.simplify((-1) ** m * (1 - sp.Rational(1, 2) ** (m - 1)) * eta)
        checks[f"signed first-defect coefficient for m={m}"] = sp.simplify(actual - predicted) == 0

    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise AssertionError("Section 3 verification failed: " + "; ".join(failed))

    report = [
        "# Section 3 symbolic verification",
        "",
        "The checks below independently reconstruct the exact Choi reduction and the local asymptotic laws used in Section 3.",
        "",
    ]
    report.extend(f"- PASS: {name}" for name in checks)
    report.extend(
        [
            "",
            "## Exact Choi matrix",
            "",
            "```text",
            str(J),
            "```",
            "",
            "## Backward-Euler two-rate margin",
            "",
            "```text",
            str(be_margin),
            "```",
        ]
    )

    root = Path(__file__).resolve().parents[1]
    out = root / "results" / "section3" / "symbolic_verification.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
