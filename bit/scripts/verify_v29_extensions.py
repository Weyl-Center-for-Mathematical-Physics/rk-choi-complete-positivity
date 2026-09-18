#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

import sympy as sp

from rk_choi_margin.candidates import (
    CandidateSpec,
    build_candidate_execution,
    certify_candidate,
    locate_candidate_components,
    richardson_rk4_margin_nonrotating,
    richardson_rk4_stability_function,
    rk4_step_doubling_bundle,
    verify_execution_provenance,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import (
    _rk4_step_doubling_evaluation,
    halving_threshold_certificate,
    run_adaptive_channel_benchmark,
)
from rk_choi_margin.symbolic import first_exponential_defect
from rk_choi_margin.topology import frequency_factor_q

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v29"
OUT.mkdir(parents=True, exist_ok=True)

x = sp.Symbol("x", nonnegative=True, real=True)
s = sp.Symbol("s")
q = sp.Symbol("q", real=True)


def payload_components(location):
    return [
        {
            "lower": str(component.lower_value),
            "upper": str(component.upper_value),
            "lower_decimal": float(component.lower_value),
            "upper_decimal": float(component.upper_value),
            "left_sources": [] if component.left is None else list(component.left.source),
            "right_sources": [] if component.right is None else list(component.right.source),
        }
        for component in location.components
    ]


# Full, fine, and Richardson candidate-map regions.
direct = locate_candidate_components(
    "direct", method="rk4", theta=0, kappa=0, varpi=2, upper=3, root_digits=32
)
fine = locate_candidate_components(
    "equal-substeps", method="rk4", theta=0, kappa=0, varpi=2, upper=3, substeps=2, root_digits=32
)
rich_global = locate_candidate_components(
    "richardson-4", method="rk4", theta=0, kappa=0, varpi=0, upper=8, root_digits=32
)

if len(direct.components) != 1 or len(fine.components) != 1 or len(rich_global.components) != 1:
    raise RuntimeError("unexpected candidate component count")

# Richardson as an ordinary stability function.
Rext = sp.expand(richardson_rk4_stability_function())
Rexpected = (
    1 + s + s**2/sp.Integer(2) + s**3/sp.Integer(6) + s**4/sp.Integer(24)
    + s**5/sp.Integer(120) + s**6/sp.Integer(864) + s**7/sp.Integer(8640)
    + s**8/sp.Integer(138240)
)
if sp.expand(Rext - Rexpected) != 0:
    raise RuntimeError("Richardson stability-function expansion mismatch")
defect = first_exponential_defect(Rext, x=s)
if (defect.index, defect.coefficient) != (6, -sp.Rational(1, 4320)):
    raise RuntimeError(f"unexpected Richardson defect {defect}")
G6 = sp.factor(frequency_factor_q(6, q))
if sp.factor(G6 - q*(q**2 - 18*q + 48)/32) != 0:
    raise RuntimeError(f"unexpected G6 factor {G6}")

margin = sp.factor(richardson_rk4_margin_nonrotating())
num, den = sp.fraction(margin)
p10 = sp.Poly(-num/x**6, x, domain=sp.QQ)
p10_intervals = [
    {"lower": str(lo), "upper": str(hi), "multiplicity": int(mult)}
    for (lo, hi), mult in p10.intervals(eps=sp.Rational(1, 10**28))
]
# The actual common constituent interval ends at the ordinary RK4 population ceiling.
population_poly = sp.Poly(x**3 - 4*x**2 + 12*x - 24, x, domain=sp.QQ)
population_positive = [
    (sp.Rational(lo), sp.Rational(hi))
    for (lo, hi), mult in population_poly.intervals(eps=sp.Rational(1, 10**28))
    if hi > 0
]
if len(population_positive) != 1:
    raise RuntimeError("could not isolate alpha4")
alpha4 = sum(population_positive[0], sp.Rational(0))/2
if any(sp.Rational(lo) <= alpha4 for (lo, hi), mult in p10.intervals() if hi > 0):
    raise RuntimeError("P10 has a positive root inside the common constituent CPTP interval")

frequency_rows = []
for varpi in (sp.Integer(0), sp.Integer(1), sp.Rational(3, 2), sp.Integer(2), sp.Integer(5)):
    qv = 1 + 4*varpi**2
    coefficient = sp.factor(defect.coefficient * G6.subs(q, qv))
    frequency_rows.append(
        {
            "varpi": str(varpi),
            "coefficient": str(coefficient),
            "orientation": "inward" if coefficient > 0 else "outward" if coefficient < 0 else "degenerate",
        }
    )

# Disjointness and the n-substep ladder.
dcomp = direct.components[0]
ratio = sp.N(dcomp.upper_value/dcomp.lower_value, 50)
ladder = []
for n in range(1, 9):
    ladder.append(
        {
            "n": n,
            "threshold": str(sp.Rational(n+1, n)),
            "overlap": bool(ratio >= sp.Rational(n+1, n)),
        }
    )
halving = halving_threshold_certificate(45)

# Candidate identity examples.
examples = []
for H in (1.0, 2.0):
    evaluation = _rk4_step_doubling_evaluation(
        h=H, theta=0, kappa=0, varpi=2, frame="lab"
    )
    examples.append(
        {
            "H": H,
            "coarse_status": evaluation.coarse_execution.cp_status.value,
            "fine_status": evaluation.fine_execution.cp_status.value,
            "richardson_status": evaluation.extrapolated_execution.cp_status.value,
            "coarse_hash": evaluation.coarse_execution.provenance_hash,
            "fine_hash": evaluation.fine_execution.provenance_hash,
            "richardson_hash": evaluation.extrapolated_execution.provenance_hash,
        }
    )

# Tolerance sweep with all structural comparators.
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
        benchmark.append(
            {
                "scenario": "boundary_tolerance_sweep",
                "tolerance": tol,
                "policy": policy,
                "global_normalized_choi_error": result.global_normalized_choi_error,
                "minimum_step_choi_eigenvalue": result.min_step_choi_eigenvalue,
                **result.stats,
            }
        )

# Force an actual component projection and record candidate provenance.
projection = run_adaptive_channel_benchmark(
    final_time=2.8,
    initial_h=2.8,
    tolerance=10.0,
    theta=0,
    kappa=0,
    varpi=2,
    policy="candidate_guard",
    fallback="rotating_frame",
)
benchmark.append(
    {
        "scenario": "projection_coverage",
        "tolerance": 10.0,
        "policy": "candidate_guard",
        "global_normalized_choi_error": projection.global_normalized_choi_error,
        "minimum_step_choi_eigenvalue": projection.min_step_choi_eigenvalue,
        **projection.stats,
    }
)
if projection.stats["component_searches"] < 1:
    raise RuntimeError("projection benchmark failed to exercise component search")
if not all(record.provenance_verified for record in projection.records):
    raise RuntimeError("projection benchmark contains an unverified candidate record")

with (OUT / "adaptive_benchmark_v29.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(benchmark[0]))
    writer.writeheader()
    writer.writerows(benchmark)

payload = {
    "candidate_regions_varpi_2": {
        "full": payload_components(direct),
        "two_half": payload_components(fine),
        "disjoint": float(dcomp.upper_value) < 2*float(dcomp.lower_value),
        "ratio_xplus_xminus": str(ratio),
        "ladder": ladder,
        "halving_threshold": halving,
    },
    "candidate_examples": examples,
    "richardson": {
        "stability_function": str(Rext),
        "defect_index": defect.index,
        "eta": str(defect.coefficient),
        "G6": str(G6),
        "frequency_orientations": frequency_rows,
        "margin": str(margin),
        "denominator": str(den),
        "P10": str(p10.as_expr()),
        "P10_root_intervals": p10_intervals,
        "common_constituent_upper_alpha4": str(alpha4),
        "global_components": payload_components(rich_global),
    },
    "benchmark_rows": benchmark,
    "projection_records": [record.__dict__ for record in projection.records],
}
(OUT / "v29_verification.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

md = [
    "# v2.9 candidate-map verification",
    "",
    "## Candidate-map topology",
    f"- full RK4 interval at varpi=2: `{payload_components(direct)}`",
    f"- two-half-step interval: `{payload_components(fine)}`",
    f"- the intervals are disjoint: **{payload['candidate_regions_varpi_2']['disjoint']}**",
    "",
    "## Richardson extrapolation",
    f"- stability function: `{Rext}`",
    f"- first exponential defect: m={defect.index}, eta={defect.coefficient}",
    f"- exact margin leading coefficient: `{sp.expand(margin).coeff(x, 6)}`",
    f"- global remote component: `{payload_components(rich_global)}`",
    "",
    "## Executed-candidate controller",
    f"- projection component searches: `{projection.stats['component_searches']}`",
    "- every accumulated map in the projection benchmark has verified provenance",
]
(OUT / "v29_verification.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
