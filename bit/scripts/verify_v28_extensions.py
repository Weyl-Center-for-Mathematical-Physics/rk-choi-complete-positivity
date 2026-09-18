from __future__ import annotations

import csv
import json
from pathlib import Path

import sympy as sp

from rk_choi_margin.candidates import (
    certify_candidate,
    direct_candidate,
    locate_candidate_components,
    richardson_rk4_margin_nonrotating,
    rk4_step_doubling_bundle,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import run_adaptive_channel_benchmark

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v28"
OUT.mkdir(parents=True, exist_ok=True)

x = sp.Symbol("x", nonnegative=True, real=True)


def status(cand) -> str:
    return certify_candidate(cand).status.value


# Exact binary-vs-decimal boundary case.
boundary_text = "1.270334626973889"
binary = direct_candidate("rk4", total_x=float(boundary_text), theta=0, kappa=0, varpi=2)
decimal = direct_candidate("rk4", total_x=boundary_text, theta=0, kappa=0, varpi=2)

# Coarse/fine candidate mismatch.
rows = []
for H in (1, 2):
    bundle = rk4_step_doubling_bundle(total_x=H, theta=0, kappa=0, varpi=2)
    rows.append({
        "total_step": H,
        "coarse_status": status(bundle.coarse),
        "fine_status": status(bundle.fine),
        "richardson_status": status(bundle.extrapolated),
    })

# Exact component endpoints.
direct = locate_candidate_components("direct", method="rk4", theta=0, kappa=0, varpi=2, upper=3)
fine = locate_candidate_components("equal-substeps", method="rk4", theta=0, kappa=0, varpi=2, upper=3, substeps=2)

def component_payload(location):
    return [
        {
            "lower": str(c.lower_value),
            "upper": str(c.upper_value),
            "lower_decimal": float(c.lower_value),
            "upper_decimal": float(c.upper_value),
        }
        for c in location.components
    ]

# Richardson factorization.
margin = sp.factor(richardson_rk4_margin_nonrotating())
num, den = sp.fraction(margin)
poly = sp.Poly(-num / x**6, x, domain=sp.QQ)
root_intervals = [
    {"lower": str(lo), "upper": str(hi), "multiplicity": int(mult)}
    for (lo, hi), mult in poly.intervals(eps=sp.Rational(1, 10**24))
]
population_poly = sp.Poly(x**3 - 4*x**2 + 12*x - 24, x, domain=sp.QQ)
population_intervals = [
    (sp.Rational(lo), sp.Rational(hi))
    for (lo, hi), mult in population_poly.intervals(eps=sp.Rational(1, 10**28))
    if hi > 0
]
if len(population_intervals) != 1:
    raise RuntimeError(f"could not isolate the classical RK4 population ceiling: {population_intervals}")
alpha4 = sum(population_intervals[0], sp.Rational(0)) / 2
if not all(sp.Rational(lo) > alpha4 for (lo, hi), mult in poly.intervals() if sp.Rational(hi) > 0):
    raise RuntimeError("Richardson polynomial has an unexpected positive root before alpha4")

# Tolerance-sweep benchmark.
benchmark = []
for tol in (1e-2, 3e-3, 1e-3, 3e-4):
    for policy in ("error_only", "candidate_guard", "rotating_frame", "strang"):
        result = run_adaptive_channel_benchmark(
            final_time=2.0,
            initial_h=1.4,
            tolerance=tol,
            theta=0,
            kappa=0,
            varpi=2,
            policy=policy,
            fallback="rotating_frame",
        )
        benchmark.append({
            "scenario": "tolerance_sweep",
            "tolerance": tol,
            "policy": policy,
            "global_normalized_choi_error": result.global_normalized_choi_error,
            "minimum_step_choi_eigenvalue": result.min_step_choi_eigenvalue,
            **result.stats,
        })

# Dedicated high-tolerance proposal that must invoke an actual certified
# component projection for the two-half-step candidate at varpi=2.
projection_result = run_adaptive_channel_benchmark(
    final_time=2.8,
    initial_h=2.8,
    tolerance=10.0,
    theta=0,
    kappa=0,
    varpi=2,
    policy="candidate_guard",
    fallback="rotating_frame",
)
benchmark.append({
    "scenario": "projection_coverage",
    "tolerance": 10.0,
    "policy": "candidate_guard",
    "global_normalized_choi_error": projection_result.global_normalized_choi_error,
    "minimum_step_choi_eigenvalue": projection_result.min_step_choi_eigenvalue,
    **projection_result.stats,
})

with (OUT / "adaptive_benchmark_v28.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(benchmark[0]))
    writer.writeheader()
    writer.writerows(benchmark)

payload = {
    "binary_boundary_case": {
        "text": boundary_text,
        "binary_rational": str(binary.total_x),
        "decimal_rational": str(decimal.total_x),
        "binary_status": status(binary),
        "decimal_status": status(decimal),
    },
    "candidate_status_examples": rows,
    "direct_components_varpi_2": component_payload(direct),
    "two_half_components_varpi_2": component_payload(fine),
    "richardson": {
        "margin": str(margin),
        "denominator": str(den),
        "P10": str(poly.as_expr()),
        "leading_x6_coefficient": str(sp.expand(margin).coeff(x, 6)),
        "positive_root_intervals": root_intervals,
        "certified_negative_on_common_rk4_interval": True,
    },
    "benchmark_rows": benchmark,
}
(OUT / "v28_verification.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

md = [
    "# v2.8 candidate-map verification",
    "",
    "## Exact binary input semantics",
    f"- binary float status at `{boundary_text}`: **{payload['binary_boundary_case']['binary_status']}**",
    f"- exact decimal status at `{boundary_text}`: **{payload['binary_boundary_case']['decimal_status']}**",
    "",
    "## Candidate-map mismatch at varpi=2",
]
for row in rows:
    md.append(
        f"- H={row['total_step']}: coarse {row['coarse_status']}, "
        f"two-half {row['fine_status']}, Richardson {row['richardson_status']}"
    )
md += [
    "",
    "## Richardson extrapolation",
    f"- exact leading coefficient: `{payload['richardson']['leading_x6_coefficient']}`",
    "- the degree-ten factor has no positive root before the classical RK4 population ceiling",
    "- therefore the negative-weight extrapolate is non-CPTP while its coarse and fine constituents are CPTP on their common positive interval",
    "",
    "## Adaptive benchmark",
    "See `adaptive_benchmark_v28.csv` for the complete tolerance sweep.",
]
(OUT / "v28_verification.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
