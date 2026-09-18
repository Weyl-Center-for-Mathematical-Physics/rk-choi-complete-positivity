#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

from rk_choi_margin.methods import all_methods
from rk_choi_margin.topology import (
    degenerate_next_coefficient,
    frequency_factor_q,
    normalized_choi_defect_and_error,
    rk4_general_margin,
    rk4_rotating_cubic,
    rk4_rotating_margin,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v26"
OUT.mkdir(parents=True, exist_ok=True)

x, y, q = sp.symbols("x y q", real=True)


def positive_root_records(poly: sp.Poly, eps_power: int = 16) -> list[dict]:
    records: list[dict] = []
    eps = sp.Rational(1, 10**eps_power)
    for (lo, hi), multiplicity in poly.intervals(eps=eps):
        if hi <= 0:
            continue
        if lo == hi == 0:
            continue
        if lo < 0 < hi:
            continue
        midpoint = (lo + hi) / 2
        records.append(
            {
                "lower": str(lo),
                "upper": str(hi),
                "approximation": str(sp.N(midpoint, 18)),
                "multiplicity": int(multiplicity),
            }
        )
    records.sort(key=lambda row: float(sp.N(sp.Rational(row["lower"]), 20)))
    return records


def sign_at(expr: sp.Expr, point: sp.Rational) -> int:
    value = sp.simplify(expr.subs(x, point))
    if value == 0:
        return 0
    return 1 if value > 0 else -1


def exact_case(theta: sp.Rational, kappa: sp.Rational, varpi: sp.Rational, name: str) -> dict:
    margin = sp.factor(rk4_general_margin(x, theta, kappa, varpi))
    numerator = sp.Poly(sp.together(margin).as_numer_denom()[0], x, domain=sp.QQ)
    roots = positive_root_records(numerator)
    endpoints = [sp.Rational(row["upper"]) for row in roots]
    samples: list[dict] = []
    left = sp.Rational(0)
    for root in roots:
        right = sp.Rational(root["lower"])
        sample = (left + right) / 2 if left != right else right / 2
        if sample > 0:
            samples.append({"interval": [str(left), str(right)], "sample": str(sample), "margin_sign": sign_at(margin, sample)})
        left = sp.Rational(root["upper"])
    sample = left + 1
    samples.append({"interval": [str(left), "infinity"], "sample": str(sample), "margin_sign": sign_at(margin, sample)})
    return {
        "name": name,
        "theta": str(theta),
        "kappa": str(kappa),
        "varpi": str(varpi),
        "margin": str(margin),
        "positive_roots": roots,
        "interval_signs": samples,
    }


R4 = all_methods()["rk4"].stability_function(x=sp.Symbol("s"))
s = sp.Symbol("s")
a4 = sp.factor(R4.subs(s, -x))
alpha_poly = sp.Poly(sp.factor((a4 - 1) * 24 / x), x, domain=sp.QQ)
alpha_record = positive_root_records(alpha_poly)[0]

cubic = rk4_rotating_cubic(q, y)
discriminant = sp.factor(sp.discriminant(cubic, y))
q_minus = (sp.Integer(123) - 11 * sp.sqrt(33)) / 16
q_plus = (sp.Integer(123) + 11 * sp.sqrt(33)) / 16
varpi_minus = sp.sqrt(107 - 11 * sp.sqrt(33)) / 8
varpi_plus = sp.sqrt(107 + 11 * sp.sqrt(33)) / 8

# Population/margin active-boundary crossover.
resultant_poly = sp.Poly(
    9 * q**9
    - 48 * q**8
    - 80 * q**7
    + 1536 * q**6
    - 768 * q**5
    - 55296 * q**4
    + 393216 * q**3
    - 1474560 * q**2
    + 2949120 * q
    - 2359296,
    q,
    domain=sp.QQ,
)
q0_candidates = [row for row in positive_root_records(resultant_poly) if 1 < float(sp.N(sp.Rational(row["lower"]), 20)) < float(sp.N(q_minus, 20))]
assert len(q0_candidates) == 1
q0_record = q0_candidates[0]

# Representative exact robustness cases.
robust_cases = [
    exact_case(sp.Rational(0), sp.Rational(1, 1000), sp.Integer(2), "dephasing_1e-3"),
    exact_case(sp.Rational(1, 1000), sp.Rational(0), sp.Integer(2), "thermal_1e-3"),
    exact_case(sp.Rational(1, 1000), sp.Rational(1, 1000), sp.Integer(2), "combined_1e-3"),
    exact_case(sp.Rational(0), sp.Rational(1, 100), sp.Integer(2), "dephasing_1e-2"),
    exact_case(sp.Rational(1, 100), sp.Rational(0), sp.Integer(2), "thermal_1e-2"),
]
for case in robust_cases:
    assert len(case["positive_roots"]) == 3

# Bell-state operational table with norm-consistent error.
mp.mp.dps = 90
bell_rows: list[dict] = []
for varpi in (1.0, 2.0, 5.0):
    for step in (0.1, 0.01, 0.001):
        defect, error = normalized_choi_defect_and_error(step, varpi)
        bell_rows.append(
            {
                "x": f"{step:.3g}",
                "varpi": f"{varpi:g}",
                "negative_born_probability_magnitude": mp.nstr(defect, 18),
                "normalized_choi_trace_distance": mp.nstr(error, 18),
                "fraction_of_error": mp.nstr(defect / error, 12) if error else "nan",
            }
        )
        assert error + mp.mpf("1e-70") >= defect

# High-frequency window samples.
window_rows: list[dict] = []
for varpi in (2, 5, 10, 25, 100, 200):
    qv = sp.Integer(1) + 4 * sp.Integer(varpi) ** 2
    roots = sp.nroots(sp.Poly(cubic.subs(q, qv), y), n=80, maxsteps=2000)
    positive = sorted(sp.N(sp.re(root) / qv, 30) for root in roots if abs(float(sp.im(root))) < 1e-20 and float(sp.re(root)) > 0)
    if len(positive) == 2:
        lo, hi = positive
        window_rows.append(
            {
                "varpi": varpi,
                "x_lower": str(lo),
                "x_upper": str(hi),
                "x_lower_varpi_squared": str(sp.N(lo * varpi**2, 18)),
                "x_upper_varpi": str(sp.N(hi * varpi, 18)),
            }
        )

# Frequency-factor checks.
frequency_factors = {str(m): str(frequency_factor_q(m, q)) for m in range(2, 9)}
assert sp.factor(frequency_factor_q(6, q)) == q * (q**2 - 18 * q + 48) / 32
eta4 = -sp.Rational(1, 24)
eta5 = -sp.Rational(1, 120)
eta6 = -sp.Rational(1, 720)
assert sp.simplify(degenerate_next_coefficient(4, eta4, eta5, sp.sqrt(7) / 2) - sp.Rational(1, 3)) == 0
assert sp.simplify(degenerate_next_coefficient(5, eta5, eta6, sp.sqrt(3) / 2) + sp.Rational(1, 144)) == 0

# CPTP-controller benchmark at varpi=2.
qv2 = sp.Integer(17)
roots2 = sp.nroots(sp.Poly(cubic.subs(q, qv2), y), n=80, maxsteps=2000)
pos2 = sorted(sp.N(sp.re(root) / qv2, 30) for root in roots2 if abs(float(sp.im(root))) < 1e-20 and float(sp.re(root)) > 0)
assert len(pos2) == 2
controller_sequence = [sp.Rational(14,10) / 2**k for k in range(7)]
assert controller_sequence[0] > pos2[1]
assert all(step < pos2[0] for step in controller_sequence[1:])

payload = {
    "schema_version": "2.6",
    "rk4_rotating_margin": str(rk4_rotating_margin(x, sp.Symbol("varpi", real=True))),
    "rk4_cubic": str(cubic),
    "rk4_cubic_discriminant": str(discriminant),
    "critical_parameters": {
        "q_minus": str(q_minus),
        "q_minus_approx": str(sp.N(q_minus, 18)),
        "varpi_minus": str(varpi_minus),
        "varpi_minus_approx": str(sp.N(varpi_minus, 18)),
        "q_unit_circle": "4",
        "varpi_unit_circle": str(sp.sqrt(3) / 2),
        "q_plus": str(q_plus),
        "q_plus_approx": str(sp.N(q_plus, 18)),
        "varpi_plus": str(varpi_plus),
        "varpi_plus_approx": str(sp.N(varpi_plus, 18)),
        "q_population_crossover": q0_record,
    },
    "critical_factorizations": {
        "q_minus": str(sp.factor(cubic.subs(q, q_minus), extension=sp.sqrt(33))),
        "q_4": str(sp.factor(cubic.subs(q, 4))),
        "q_plus": str(sp.factor(cubic.subs(q, q_plus), extension=sp.sqrt(33))),
    },
    "population_ceiling_root": alpha_record,
    "frequency_factors": frequency_factors,
    "dp5_thresholds": [
        str(sp.sqrt(2 - sp.sqrt(33) / 4)),
        str(sp.sqrt(2 + sp.sqrt(33) / 4)),
        str(sp.N(sp.sqrt(2 - sp.sqrt(33) / 4), 18)),
        str(sp.N(sp.sqrt(2 + sp.sqrt(33) / 4), 18)),
    ],
    "robustness_cases": robust_cases,
    "bell_error_table": bell_rows,
    "high_frequency_window": window_rows,
    "controller_example": {"varpi": 2, "interval": [str(pos2[0]), str(pos2[1])], "halving_sequence": [str(sp.N(z,18)) for z in controller_sequence]},
}

(OUT / "v26_verification.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

with (OUT / "bell_error_table.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(bell_rows[0]))
    writer.writeheader(); writer.writerows(bell_rows)
with (OUT / "high_frequency_window.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(window_rows[0]))
    writer.writeheader(); writer.writerows(window_rows)
with (OUT / "robustness_roots.csv").open("w", newline="", encoding="utf-8") as handle:
    fields = ["case", "theta", "kappa", "varpi", "root_index", "lower", "upper", "approximation", "multiplicity"]
    writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
    for case in robust_cases:
        for idx, root in enumerate(case["positive_roots"], start=1):
            writer.writerow({"case": case["name"], "theta": case["theta"], "kappa": case["kappa"], "varpi": case["varpi"], "root_index": idx, **root})

md = f"""# Version 2.6 exact verification\n\n- RK4 rotating cubic: `{cubic}`\n- Discriminant: `{discriminant}`\n- q_- = {sp.N(q_minus, 16)}, varpi_- = {sp.N(varpi_minus, 16)}\n- q = 4, varpi = sqrt(3)/2\n- q_+ = {sp.N(q_plus, 16)}, varpi_+ = {sp.N(varpi_plus, 16)}\n- Correct DP5 factor: G6 = `{frequency_factor_q(6, q)}`\n- Exact off-boundary cases: {len(robust_cases)}, each with three positive margin roots.\n- Bell-state error rows: {len(bell_rows)}; the physicality defect never exceeds the total normalized-Choi trace error.
- Controller counterexample: exact detached interval at varpi=2 with a halving sequence that skips it.\n"""
(OUT / "v26_verification.md").write_text(md, encoding="utf-8")
print(md)
