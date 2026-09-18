from __future__ import annotations

from pathlib import Path
import sys

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rk_choi_margin.channels import cp_conditions
from rk_choi_margin.methods import X, all_methods
from rk_choi_margin.symbolic import (
    cp_admissible_intervals,
    first_exponential_defect,
    predicted_boundary_coefficient,
)


OUT = ROOT / "results" / "section5" / "symbolic_verification.md"


def lower_choi_eigenvalue(a: mp.mpf, c: mp.mpf) -> mp.mpf:
    margin = a - c**2
    denominator = 1 + a + mp.sqrt((1 - a) ** 2 + 4 * c**2)
    return 2 * margin / denominator


def main() -> None:
    mp.mp.dps = 100
    methods = all_methods()
    lines: list[str] = [
        "# Section 5 verification",
        "",
        "## Exact repeated-step factorization",
        "",
    ]

    a, c = sp.symbols("a c", real=True)
    for count in (1, 2, 3, 5, 8):
        lhs = a**count - c ** (2 * count)
        rhs = (a - c**2) * sum(a ** (count - 1 - j) * c ** (2 * j) for j in range(count))
        assert sp.expand(lhs - rhs) == 0
        lines.append(f"- N={count}: verified exactly")

    lines.extend(["", "## Fixed-horizon asymptotics at X=1", ""])
    final_x = mp.mpf("1")
    count = 32768
    step = final_x / count
    for key in ("euler", "rk3", "dp4"):
        method = methods[key]
        R = method.stability_function()
        defect = first_exponential_defect(R)
        K = mp.mpf(str(sp.N(predicted_boundary_coefficient(defect), 90)))
        Rf = sp.lambdify(X, R, "mpmath")
        A = mp.mpf(Rf(-step)) ** count
        C = mp.mpf(Rf(-step / 2)) ** count
        margin = A - C**2
        lam = lower_choi_eigenvalue(A, C)
        margin_ratio = margin / step ** (defect.index - 1)
        margin_expected = K * mp.e ** (-1)
        eig_ratio = lam / step ** (defect.index - 1)
        eig_expected = margin_expected / (1 + mp.e ** (-1))
        assert mp.almosteq(margin_ratio, margin_expected, rel_eps=mp.mpf("2e-4"))
        assert mp.almosteq(eig_ratio, eig_expected, rel_eps=mp.mpf("2e-4"))
        lines.extend(
            [
                f"### {method.name}",
                f"- m = {defect.index}",
                f"- K = {mp.nstr(K, 18)}",
                f"- observed M_N/x^(m-1) = {mp.nstr(margin_ratio, 18)}",
                f"- predicted = {mp.nstr(margin_expected, 18)}",
                f"- observed lambda_-/x^(m-1) = {mp.nstr(eig_ratio, 18)}",
                f"- predicted = {mp.nstr(eig_expected, 18)}",
                "",
            ]
        )

    euler_R = methods["euler"].stability_function()
    euler_a = euler_R.subs(X, -1)
    euler_c = euler_R.subs(X, -sp.Rational(1, 2))
    assert abs(euler_a) <= 1 and abs(euler_c) <= 1
    assert euler_a - euler_c**2 == -sp.Rational(1, 4)

    dp4_R = methods["dp4"].stability_function()
    dp4_a = sp.factor(dp4_R.subs(X, -4))
    dp4_c = sp.factor(dp4_R.subs(X, -2))
    dp4_margin = sp.factor(dp4_a - dp4_c**2)
    assert dp4_a == sp.Rational(847, 1875)
    assert dp4_c == sp.Rational(91, 750)
    assert dp4_margin == sp.Rational(245819, 562500)
    assert cp_conditions(complex(float(dp4_a)), complex(float(dp4_c)))

    lines.extend(
        [
            "## Logical-separation checks",
            "",
            "- Forward Euler at x=1: both scalar multipliers satisfy absolute stability, but M=-1/4.",
            f"- SSPRK(3,3) admissible set: {cp_admissible_intervals(methods['rk3'].stability_function())}",
            f"- Classical RK4 admissible set: {cp_admissible_intervals(methods['rk4'].stability_function())}",
            f"- Embedded DP4 at x=4: a={dp4_a}, c={dp4_c}, M={dp4_margin} > 0.",
            "",
        ]
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
