from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rk_choi_margin.methods import X, all_methods
from rk_choi_margin.symbolic import Z, cp_margin, first_exponential_defect, predicted_boundary_coefficient


@dataclass(frozen=True)
class RootCertificate:
    name: str
    polynomial: sp.Expr
    lower: sp.Rational
    upper: sp.Rational

    def verify(self) -> None:
        poly = sp.Poly(self.polynomial, Z)
        count = poly.count_roots(self.lower, self.upper)
        if count != 1:
            raise AssertionError(f"{self.name}: expected one root, found {count}")
        left = sp.sign(poly.eval(self.lower))
        right = sp.sign(poly.eval(self.upper))
        if left == 0 or right == 0 or left == right:
            raise AssertionError(f"{self.name}: bracket does not isolate a simple sign-changing root")

    @property
    def decimal_interval(self) -> str:
        return f"({sp.N(self.lower, 16)}, {sp.N(self.upper, 16)})"


def q(value: str) -> sp.Rational:
    return sp.Rational(value)


# Compact renderings of the sets certified below from exact factors, Sturm
# root counts, isolated rational brackets, and exact interval sign tests.  They
# are output labels, not independent mathematical assumptions.  Separating
# rendering from generic numerical root discovery keeps this closeout check
# deterministic across SymPy and platform versions.
CERTIFIED_AMPLITUDE_SETS = {
    "euler": "{0}",
    "rk2": "[0, 2]",
    "rk3": "{0}",
    "rk4": "[0, 2.78529356341]",
    "dp5": "[0, 3.30656789263]",
    "dp4": "{0} union [3.06885654805, 4.38498632080]",
    "be": "[0, infinity)",
    "im": "{0}",
}


def amplitude_set_text(key: str) -> str:
    return CERTIFIED_AMPLITUDE_SETS[key]


def main() -> None:
    methods = all_methods()
    results_dir = ROOT / "results" / "section4"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Exact low-order factors.
    low_order_expected = {
        "euler": -Z**2 / 4,
        "rk2": Z**3 * (8 - Z) / 64,
        "rk3": -Z**4 * (Z**2 - 12 * Z + 84) / 2304,
        "rk4": -Z**5 * (Z**3 - 16 * Z**2 + 160 * Z - 1152) / 147456,
    }
    for key, expected in low_order_expected.items():
        actual = sp.factor(cp_margin(methods[key].stability_function()))
        if sp.simplify(actual - expected) != 0:
            raise AssertionError(f"Incorrect exact margin for {key}")

    rk4_population = Z**3 - 4 * Z**2 + 12 * Z - 24
    rk4_margin = Z**3 - 16 * Z**2 + 160 * Z - 1152

    dp5_R = methods["dp5"].stability_function()
    dp5_population = sp.factor((dp5_R.subs(X, -Z) - 1) * 600 / Z)
    dp5_margin_numerator = sp.factor(sp.together(cp_margin(dp5_R)).as_numer_denom()[0] / Z**6)

    dp4_R = methods["dp4"].stability_function()
    dp4_population = sp.factor(-(dp4_R.subs(X, -Z) - 1) * 120000 / Z)
    dp4_a_numerator = sp.factor(sp.together(dp4_R.subs(X, -Z)).as_numer_denom()[0])
    dp4_margin_numerator = sp.factor(-sp.together(cp_margin(dp4_R)).as_numer_denom()[0] / Z**5)

    certificates = [
        RootCertificate("RK4 population endpoint", rk4_population, q("2.7852935634"), q("2.7852935635")),
        RootCertificate("RK4 first margin zero", rk4_margin, q("10.9824254662"), q("10.9824254664")),
        RootCertificate("DP5 population endpoint", dp5_population, q("3.3065678926"), q("3.3065678927")),
        RootCertificate("DP5 first margin zero", dp5_margin_numerator, q("13.1740934628"), q("13.1740934629")),
        RootCertificate("DP4 margin re-entry", dp4_margin_numerator, q("3.0688565480"), q("3.0688565481")),
        RootCertificate("DP4 population exit", dp4_population, q("4.3849863208"), q("4.3849863209")),
        RootCertificate("DP4 second margin zero", dp4_margin_numerator, q("15.2987678075"), q("15.2987678077")),
        RootCertificate("DP4 second population crossing", dp4_population, q("24.7277543187"), q("24.7277543189")),
        RootCertificate("DP4 population-multiplier zero", dp4_a_numerator, q("24.7278950301"), q("24.7278950302")),
    ]
    for certificate in certificates:
        certificate.verify()

    # Root counts protect the global interval classification.
    assert sp.Poly(rk4_population, Z).count_roots(0, sp.oo) == 1
    assert sp.Poly(rk4_margin, Z).count_roots(0, sp.oo) == 1
    assert sp.Poly(dp5_population, Z).count_roots(0, sp.oo) == 1
    assert sp.Poly(dp5_margin_numerator, Z).count_roots(0, sp.oo) == 1
    assert sp.Poly(dp4_population, Z).count_roots(0, sp.oo) == 2
    assert sp.Poly(dp4_margin_numerator, Z).count_roots(0, sp.oo) == 2
    assert sp.Poly(dp4_a_numerator, Z).count_roots(0, sp.oo) == 1

    # Exact sign chart for the disconnected embedded fourth-order set.
    a_dp4 = sp.factor(dp4_R.subs(X, -Z))
    M_dp4 = sp.factor(cp_margin(dp4_R))
    sign_samples = {
        q("1"): ("outward neighborhood", False),
        q("4"): ("isolated CPTP component", True),
        q("5"): ("population overflow", False),
        q("20"): ("negative block margin", False),
        q("25"): ("negative population multiplier", False),
    }
    for point, (_, expected_cptp) in sign_samples.items():
        a_value = sp.sign(a_dp4.subs(Z, point))
        one_minus_a = sp.sign((1 - a_dp4).subs(Z, point))
        margin_value = sp.sign(M_dp4.subs(Z, point))
        actual_cptp = a_value >= 0 and one_minus_a >= 0 and margin_value >= 0
        if actual_cptp != expected_cptp:
            raise AssertionError(f"DP4 sign chart failed at x={point}")

    # Implicit counterexamples to a universal order-parity rule.
    be_R = methods["be"].stability_function()
    im_R = methods["im"].stability_function()
    be_margin = sp.factor(cp_margin(be_R))
    im_margin = sp.factor(cp_margin(im_R))
    assert sp.simplify(be_margin - Z**2 / ((Z + 1) * (Z + 2) ** 2)) == 0
    assert sp.simplify(im_margin + 2 * Z**3 / ((Z + 2) * (Z + 4) ** 2)) == 0

    report: list[str] = [
        "SECTION 4 EXACT METHOD AUDIT",
        "=" * 36,
        "",
        "Certified decimal root isolations (exact Sturm counts):",
    ]
    for certificate in certificates:
        report.append(f"- {certificate.name}: {certificate.decimal_interval}")

    report.extend(
        [
            "",
            "Exact low-order margins:",
            *[f"- {methods[key].name}: {sp.sstr(sp.factor(expr))}" for key, expr in low_order_expected.items()],
            "",
            "Pure-amplitude-damping admissible sets:",
        ]
    )
    for key in ("euler", "rk2", "rk3", "rk4", "dp5", "dp4", "be", "im"):
        report.append(f"- {methods[key].name}: {amplitude_set_text(key)}")

    report.extend(
        [
            "",
            "Dormand-Prince embedded fourth-order sign chart:",
        ]
    )
    for point, (description, _) in sign_samples.items():
        report.append(
            f"- x={point}: a={sp.N(a_dp4.subs(Z, point), 16)}, "
            f"M={sp.N(M_dp4.subs(Z, point), 16)} ({description})"
        )

    report.extend(
        [
            "",
            f"Backward Euler margin: {sp.sstr(be_margin)}",
            f"Implicit midpoint margin: {sp.sstr(im_margin)}",
            "",
            "The DP4 finite interval is a structural CPTP component, not an accuracy recommendation.",
        ]
    )

    report_path = results_dir / "exact_method_audit.txt"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    rows: list[dict[str, str | int]] = []
    for key in ("euler", "rk2", "rk3", "rk4", "dp5", "dp4", "be", "im"):
        method = methods[key]
        R = method.stability_function()
        defect = first_exponential_defect(R)
        coefficient = sp.factor(predicted_boundary_coefficient(defect))
        rows.append(
            {
                "key": key,
                "method": method.name,
                "order": method.order,
                "stages": method.stages,
                "m": defect.index,
                "eta_m": sp.sstr(defect.coefficient),
                "local_margin_coefficient": sp.sstr(coefficient),
                "local_orientation": "inward" if coefficient > 0 else "outward",
                "amplitude_damping_admissible_set": amplitude_set_text(key),
            }
        )
    csv_path = results_dir / "method_audit_table.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(report_path.read_text(encoding="utf-8"))
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
