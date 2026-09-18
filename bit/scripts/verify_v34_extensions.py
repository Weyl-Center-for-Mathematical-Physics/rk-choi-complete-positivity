#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np
import sympy as sp

from rk_choi_margin.candidates import (
    locate_candidate_components,
    richardson_rk4_margin_nonrotating,
    richardson_rk4_stability_function,
    verify_execution_provenance,
)
from rk_choi_margin.controllers import (
    _rk4_step_doubling_evaluation,
    halving_threshold_certificate,
    run_adaptive_channel_benchmark,
)
from rk_choi_margin.liouvillian import (
    exact_gksl_superoperator,
    exact_to_numpy,
    exp_superoperator_numeric,
    min_choi_eigenvalue_numeric,
    normalized_choi_trace_distance,
    rk4_superoperator_numeric,
    strang_gksl_superoperator_numeric,
)
from rk_choi_margin.symbolic import first_exponential_defect
from rk_choi_margin.topology import frequency_factor_q

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v34"
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


def _tex_sci(value: float, significant: int = 4) -> str:
    if value == 0:
        return r"0"
    mantissa, exponent = f"{value:.{significant-1}e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    return rf"{mantissa}\times 10^{{{int(exponent)}}}"


def _run_benchmark(**kwargs):
    start = time.perf_counter()
    result = run_adaptive_channel_benchmark(**kwargs)
    elapsed = time.perf_counter() - start
    return result, elapsed


# ---------------------------------------------------------------------------
# Candidate-map theorems and exact Richardson certificate.
# ---------------------------------------------------------------------------
direct = locate_candidate_components(
    "direct", method="rk4", theta=0, kappa=0, varpi=2, upper=3, root_digits=36
)
fine = locate_candidate_components(
    "equal-substeps", method="rk4", theta=0, kappa=0, varpi=2,
    upper=3, substeps=2, root_digits=36
)
rich_global = locate_candidate_components(
    "richardson-4", method="rk4", theta=0, kappa=0, varpi=0,
    upper=8, root_digits=36
)
if direct.status.value != "PASS" or fine.status.value != "PASS" or rich_global.status.value != "PASS":
    raise RuntimeError("candidate component locator did not return certified PASS")
if len(direct.components) != 1 or len(fine.components) != 1 or len(rich_global.components) != 1:
    raise RuntimeError("unexpected candidate component count")

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
    for (lo, hi), mult in p10.intervals(eps=sp.Rational(1, 10**30))
]
population_poly = sp.Poly(x**3 - 4*x**2 + 12*x - 24, x, domain=sp.QQ)
population_positive = [
    (sp.Rational(lo), sp.Rational(hi))
    for (lo, hi), mult in population_poly.intervals(eps=sp.Rational(1, 10**30))
    if hi > 0
]
if len(population_positive) != 1:
    raise RuntimeError("could not isolate alpha4")
alpha4 = sum(population_positive[0], sp.Rational(0))/2
if any(sp.Rational(lo) <= alpha4 for (lo, hi), mult in p10.intervals() if hi > 0):
    raise RuntimeError("P10 has a positive root inside the common constituent CPTP interval")

frequency_rows = []
for varpi_value in (sp.Integer(0), sp.Integer(1), sp.Rational(3, 2), sp.Integer(2), sp.Integer(5)):
    qv = 1 + 4*varpi_value**2
    coefficient = sp.factor(defect.coefficient * G6.subs(q, qv))
    frequency_rows.append({
        "varpi": str(varpi_value),
        "coefficient": str(coefficient),
        "orientation": "inward" if coefficient > 0 else "outward" if coefficient < 0 else "degenerate",
    })

dcomp = direct.components[0]
ratio = sp.N(dcomp.upper_value/dcomp.lower_value, 50)
ladder = [
    {"n": n, "threshold": str(sp.Rational(n+1, n)),
     "overlap": bool(ratio >= sp.Rational(n+1, n))}
    for n in range(1, 9)
]
halving = halving_threshold_certificate(48)

examples = []
for H in (1.0, 2.0):
    evaluation = _rk4_step_doubling_evaluation(h=H, theta=0, kappa=0, varpi=2, frame="lab")
    examples.append({
        "H": H,
        "coarse_status": evaluation.coarse_execution.cp_status.value,
        "fine_status": evaluation.fine_execution.cp_status.value,
        "richardson_status": evaluation.extrapolated_execution.cp_status.value,
        "coarse_hash": evaluation.coarse_execution.provenance_hash,
        "fine_hash": evaluation.fine_execution.provenance_hash,
        "richardson_hash": evaluation.extrapolated_execution.provenance_hash,
    })

# ---------------------------------------------------------------------------
# Adaptive benchmark matrix.  The boundary sweep demonstrates safe fallback;
# buffered and near-saddle cases exercise direct retention and projection.
# ---------------------------------------------------------------------------
benchmark: list[dict[str, object]] = []

def add_run(*, scenario: str, policy_label: str, role: str, **kwargs):
    result, elapsed = _run_benchmark(**kwargs)
    row: dict[str, object] = {
        "scenario": scenario,
        "scenario_role": role,
        "tolerance": kwargs["tolerance"],
        "policy": policy_label,
        "controller_policy": kwargs["policy"],
        "theta": kwargs["theta"],
        "kappa": kwargs["kappa"],
        "varpi": kwargs["varpi"],
        "final_time": kwargs["final_time"],
        "initial_h": kwargs["initial_h"],
        "global_normalized_choi_error": result.global_normalized_choi_error,
        "minimum_step_choi_eigenvalue": result.min_step_choi_eigenvalue,
        "wall_time_seconds": elapsed,
        **result.stats,
    }
    benchmark.append(row)
    if not all(record.provenance_verified for record in result.records):
        raise RuntimeError(f"{scenario}/{policy_label} contains an unverified execution record")
    return result, row

# Boundary tolerance sweep.
for tol in (1e-2, 3e-3, 1e-3, 3e-4):
    for policy, label in (
        ("error_only", "error_only"),
        ("candidate_guard", "guard_rotating_fallback"),
        ("rotating_frame", "rotating_frame"),
        ("strang", "strang"),
    ):
        add_run(
            scenario="boundary_tolerance_sweep",
            policy_label=label,
            role="boundary safety/fallback comparison",
            final_time=2.0, initial_h=1.4, tolerance=tol,
            theta=0.0, kappa=0.0, varpi=2.0,
            policy=policy, fallback="rotating_frame",
        )

# Scientifically useful candidate-aware cases.
scenario_specs = [
    dict(
        scenario="buffered_projection", policy_label="guard_buffered_projection",
        role="buffered disconnected case with direct retention and projection",
        final_time=2.0, initial_h=0.3, tolerance=1e-3,
        theta=0.0, kappa=1e-3, varpi=2.0,
        policy="candidate_guard", fallback="rotating_frame",
    ),
    dict(
        scenario="buffered_direct", policy_label="guard_buffered_direct",
        role="buffered case dominated by direct accepted candidates",
        final_time=2.0, initial_h=0.3, tolerance=3e-4,
        theta=0.0, kappa=1e-3, varpi=2.0,
        policy="candidate_guard", fallback="rotating_frame",
    ),
    dict(
        scenario="near_saddle_projection", policy_label="guard_near_saddle",
        role="near-saddle case exercising component location and projection",
        final_time=1.0, initial_h=1.0, tolerance=1e-3,
        theta=0.0, kappa=1e-3, varpi=1.65,
        policy="candidate_guard", fallback="rotating_frame",
    ),
]
scenario_rows: list[dict[str, object]] = []
for spec in scenario_specs:
    result, row = add_run(**spec)
    scenario_rows.append(row)

if int(scenario_rows[0]["projected_accepts"]) < 1 or int(scenario_rows[0]["direct_accepts"]) < 1:
    raise RuntimeError("buffered projection case did not exercise both direct and projected accepts")
if int(scenario_rows[2]["locator_calls"]) < 1:
    raise RuntimeError("near-saddle case did not exercise component location")
if float(scenario_rows[0]["global_normalized_choi_error"]) >= 1e-3:
    raise RuntimeError("buffered projection case failed its global-error usefulness criterion")

# Retain the historical large-tolerance projection only as an adversarial path test.
projection, projection_row = add_run(
    scenario="projection_adversarial_path_test",
    policy_label="guard_projection_adversarial",
    role="functionality/adversarial path test; not performance evidence",
    final_time=2.8, initial_h=2.8, tolerance=10.0,
    theta=0.0, kappa=0.0, varpi=2.0,
    policy="candidate_guard", fallback="rotating_frame",
)
if int(projection_row["locator_calls"]) < 1:
    raise RuntimeError("adversarial projection test failed to exercise component location")

# ---------------------------------------------------------------------------
# Noncommuting common dense-qubit benchmark: direct RK4 versus CPTP Strang.
# The exact semigroup is the reference; no broad scalability claim is made.
# ---------------------------------------------------------------------------
L_exact = exact_gksl_superoperator(gamma=1, theta=sp.Rational(1,1000),
                                   kappa=sp.Rational(1,1000), omega_z=2,
                                   omega_x=sp.Rational(1,10))
L_num = exact_to_numpy(L_exact)
T = 1.0
S_ref = exp_superoperator_numeric(L_num, T)
noncommuting_rows: list[dict[str, object]] = []
for N in (1, 2, 4, 8, 16, 32):
    h = T/N
    start = time.perf_counter()
    S_rk = np.linalg.matrix_power(rk4_superoperator_numeric(L_num, h), N)
    rk_time = time.perf_counter() - start
    start = time.perf_counter()
    S_step = strang_gksl_superoperator_numeric(
        h=h, gamma=1, theta=sp.Rational(1,1000), kappa=sp.Rational(1,1000),
        omega_z=2, omega_x=sp.Rational(1,10), substeps=1, dps=80)
    S_strang = np.linalg.matrix_power(S_step, N)
    strang_time = time.perf_counter() - start
    noncommuting_rows.extend([
        {
            "method": "direct_RK4", "N": N, "h": h,
            "global_normalized_choi_error": normalized_choi_trace_distance(S_rk, S_ref),
            "minimum_step_choi_eigenvalue": min_choi_eigenvalue_numeric(rk4_superoperator_numeric(L_num, h)),
            "rhs_stage_evaluations": 4*N, "exponential_actions": 0,
            "wall_time_seconds": rk_time,
        },
        {
            "method": "CPTP_Strang", "N": N, "h": h,
            "global_normalized_choi_error": normalized_choi_trace_distance(S_strang, S_ref),
            "minimum_step_choi_eigenvalue": min_choi_eigenvalue_numeric(S_step),
            "rhs_stage_evaluations": 0, "exponential_actions": 3*N,
            "wall_time_seconds": strang_time,
        },
    ])

# Write current-version result tables.
benchmark_fields = list(benchmark[0].keys())
with (OUT / "adaptive_benchmark_v34.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=benchmark_fields)
    writer.writeheader(); writer.writerows(benchmark)
with (OUT / "noncommuting_method_comparison_v34.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(noncommuting_rows[0].keys()))
    writer.writeheader(); writer.writerows(noncommuting_rows)

# TeX macros mechanically bound to current results.
sweep_rows = [row for row in benchmark if row["scenario"] == "boundary_tolerance_sweep"]
rotating_rows = [row for row in sweep_rows if row["controller_policy"] == "rotating_frame"]
guard_rows = [row for row in sweep_rows if row["controller_policy"] == "candidate_guard"]
rotating_target = next(row for row in rotating_rows if abs(float(row["tolerance"]) - 3e-4) < 1e-16)
rotating_error = float(rotating_target["global_normalized_choi_error"])
nearest_guard = min(guard_rows, key=lambda row: abs(float(row["global_normalized_choi_error"]) - rotating_error))
below_guard = [row for row in guard_rows if float(row["global_normalized_choi_error"]) < rotating_error]
first_guard_below = max(below_guard, key=lambda row: float(row["global_normalized_choi_error"]))
aggregate_guard = {
    name: sum(int(row[name]) for row in guard_rows)
    for name in (
        "accepted_steps", "direct_accepts", "projected_accepts", "rotating_accepts",
        "strang_accepts", "exact_fallback_accepts", "locator_calls",
        "successful_projections", "empty_component_outcomes", "uncertain_outcomes",
        "root_isolation_ops", "certification_ops", "binary_recertifications",
        "precision_escalations", "fallback_events", "rhs_stage_evaluations",
        "exponential_actions",
    )
}
if aggregate_guard["locator_calls"] != (
    aggregate_guard["successful_projections"]
    + aggregate_guard["empty_component_outcomes"]
    + aggregate_guard["uncertain_outcomes"]
):
    raise RuntimeError("aggregate locator telemetry does not partition locator calls")

buffered = scenario_rows[0]
near = scenario_rows[2]
macro_lines = [
    "% Generated mechanically from results/v34/adaptive_benchmark_v34.csv. Do not edit by hand.",
    rf"\providecommand{{\BenchmarkRotatingError}}{{\ensuremath{{{_tex_sci(rotating_error)}}}}}",
    rf"\providecommand{{\BenchmarkRotatingStages}}{{{int(rotating_target['rhs_stage_evaluations'])}}}",
    rf"\providecommand{{\BenchmarkNearestGuardError}}{{\ensuremath{{{_tex_sci(float(nearest_guard['global_normalized_choi_error']))}}}}}",
    rf"\providecommand{{\BenchmarkNearestGuardStages}}{{{int(nearest_guard['rhs_stage_evaluations'])}}}",
    rf"\providecommand{{\BenchmarkFirstGuardBelowError}}{{\ensuremath{{{_tex_sci(float(first_guard_below['global_normalized_choi_error']))}}}}}",
    rf"\providecommand{{\BenchmarkFirstGuardBelowStages}}{{{int(first_guard_below['rhs_stage_evaluations'])}}}",
    rf"\providecommand{{\BenchmarkGuardAccepted}}{{{aggregate_guard['accepted_steps']}}}",
    rf"\providecommand{{\BenchmarkGuardDirectAccepted}}{{{aggregate_guard['direct_accepts']}}}",
    rf"\providecommand{{\BenchmarkGuardProjectedAccepted}}{{{aggregate_guard['projected_accepts']}}}",
    rf"\providecommand{{\BenchmarkGuardRotatingAccepted}}{{{aggregate_guard['rotating_accepts']}}}",
    rf"\providecommand{{\BenchmarkGuardLocatorCalls}}{{{aggregate_guard['locator_calls']}}}",
    rf"\providecommand{{\BenchmarkGuardSuccessfulProjections}}{{{aggregate_guard['successful_projections']}}}",
    rf"\providecommand{{\BenchmarkGuardEmptyComponents}}{{{aggregate_guard['empty_component_outcomes']}}}",
    rf"\providecommand{{\BenchmarkGuardUncertainOutcomes}}{{{aggregate_guard['uncertain_outcomes']}}}",
    rf"\providecommand{{\BenchmarkGuardRootIsolationOps}}{{{aggregate_guard['root_isolation_ops']}}}",
    rf"\providecommand{{\BenchmarkGuardCertificationOps}}{{{aggregate_guard['certification_ops']}}}",
    rf"\providecommand{{\BenchmarkGuardBinaryRecertifications}}{{{aggregate_guard['binary_recertifications']}}}",
    rf"\providecommand{{\BenchmarkGuardPrecisionEscalations}}{{{aggregate_guard['precision_escalations']}}}",
    rf"\providecommand{{\BenchmarkGuardFallbackEvents}}{{{aggregate_guard['fallback_events']}}}",
    rf"\providecommand{{\BenchmarkBufferedError}}{{\ensuremath{{{_tex_sci(float(buffered['global_normalized_choi_error']))}}}}}",
    rf"\providecommand{{\BenchmarkBufferedDirect}}{{{int(buffered['direct_accepts'])}}}",
    rf"\providecommand{{\BenchmarkBufferedProjected}}{{{int(buffered['projected_accepts'])}}}",
    rf"\providecommand{{\BenchmarkBufferedLocatorCalls}}{{{int(buffered['locator_calls'])}}}",
    rf"\providecommand{{\BenchmarkBufferedRootIsolationOps}}{{{int(buffered['root_isolation_ops'])}}}",
    rf"\providecommand{{\BenchmarkBufferedCertificationOps}}{{{int(buffered['certification_ops'])}}}",
    rf"\providecommand{{\BenchmarkBufferedBinaryRecertifications}}{{{int(buffered['binary_recertifications'])}}}",
    rf"\providecommand{{\BenchmarkBufferedStages}}{{{int(buffered['rhs_stage_evaluations'])}}}",
    rf"\providecommand{{\BenchmarkNearSaddleError}}{{\ensuremath{{{_tex_sci(float(near['global_normalized_choi_error']))}}}}}",
    rf"\providecommand{{\BenchmarkNearSaddleProjected}}{{{int(near['projected_accepts'])}}}",
    rf"\providecommand{{\BenchmarkNearSaddleLocatorCalls}}{{{int(near['locator_calls'])}}}",
]
(OUT / "benchmark_macros_v34.tex").write_text("\n".join(macro_lines) + "\n", encoding="utf-8")


def _tex_policy(policy: str) -> str:
    return {
        "guard_rotating_fallback": "guard + rotating fallback",
        "guard_buffered_projection": "guard, buffered",
        "guard_buffered_direct": "guard, buffered/direct",
        "guard_near_saddle": "guard, near saddle",
        "rotating_frame": "rotating frame",
        "strang": "Strang",
        "error_only": "error only",
    }.get(policy, policy.replace("_", " "))


def _tex_value(value: float) -> str:
    if abs(value) < 5e-16:
        return r"$\ge -10^{-15}$"
    return rf"${_tex_sci(value, significant=3)}$"


boundary_repr = [
    next(r for r in sweep_rows if r["policy"] == policy and abs(float(r["tolerance"]) - 1e-3) < 1e-15)
    for policy in ("error_only", "guard_rotating_fallback", "rotating_frame", "strang")
]
table_rows = boundary_repr + scenario_rows
case_names = {
    "boundary_tolerance_sweep": "boundary",
    "buffered_projection": "buffered",
    "buffered_direct": "buffered",
    "near_saddle_projection": "near saddle",
}
lines = [
    "% Generated mechanically from results/v34/adaptive_benchmark_v34.csv. Do not edit by hand.",
    r"\begin{table*}[t]", r"\centering",
    r"\caption{Adaptive safety benchmark. The work proxy counts RK stage evaluations or exponential actions. Algebraic certification costs are excluded and reported separately in Table~\ref{tab:guard_telemetry}.}",
    r"\label{tab:adaptive_matrix}",
    r"\resizebox{0.99\textwidth}{!}{%", r"\begin{tabular}{llrrrrrr}", r"\toprule",
    r"Case & Policy & Choi error & Min. eig. & Direct & Projected & Fallback & RK stages\\",
    r"\midrule",
]
for row in table_rows:
    lines.append(
        rf"{case_names[row['scenario']]} & {_tex_policy(str(row['policy']))} & "
        rf"${_tex_sci(float(row['global_normalized_choi_error']),3)}$ & "
        rf"{_tex_value(float(row['minimum_step_choi_eigenvalue']))} & "
        rf"{int(row['direct_accepts'])} & {int(row['projected_accepts'])} & "
        rf"{int(row['fallback_events'])} & {int(row['rhs_stage_evaluations'])}\\"
    )
lines.extend([r"\bottomrule", r"\end{tabular}%", r"}", r"\end{table*}", ""])

telemetry_rows = guard_rows + scenario_rows
lines.extend([
    r"\begin{table*}[t]", r"\centering",
    r"\caption{Guard telemetry. A locator call is counted on every invocation, including empty and uncertain outcomes. Root-isolation operations count exact Sturm-isolation calls; certification operations count candidate and component certifications, including binary recertification attempts.}",
    r"\label{tab:guard_telemetry}",
    r"\resizebox{0.99\textwidth}{!}{%", r"\begin{tabular}{lrrrrrrrrrr}", r"\toprule",
    r"Case & Tol. & Locator & Proj. found & Empty & Uncertain & Root iso. & Cert. & Binary & Prec. esc. & Fallback\\",
    r"\midrule",
])
for row in telemetry_rows:
    label = case_names.get(str(row["scenario"]), str(row["scenario"]).replace("_", " "))
    lines.append(
        rf"{label} & ${_tex_sci(float(row['tolerance']),2)}$ & {int(row['locator_calls'])} & "
        rf"{int(row['successful_projections'])} & {int(row['empty_component_outcomes'])} & "
        rf"{int(row['uncertain_outcomes'])} & {int(row['root_isolation_ops'])} & "
        rf"{int(row['certification_ops'])} & {int(row['binary_recertifications'])} & "
        rf"{int(row['precision_escalations'])} & {int(row['fallback_events'])}\\"
    )
lines.extend([r"\bottomrule", r"\end{tabular}%", r"}", r"\end{table*}"])
(OUT / "adaptive_benchmark_table_v34.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

payload = {
    "candidate_regions_varpi_2": {
        "full": payload_components(direct), "two_half": payload_components(fine),
        "disjoint_positive_components": float(dcomp.upper_value) < 2*float(dcomp.lower_value),
        "ratio_xplus_xminus": str(ratio), "ladder": ladder,
        "halving_threshold": halving,
    },
    "candidate_examples": examples,
    "richardson": {
        "stability_function": str(Rext), "defect_index": defect.index,
        "eta": str(defect.coefficient), "G6": str(G6),
        "frequency_orientations": frequency_rows, "margin": str(margin),
        "denominator": str(den), "P10": str(p10.as_expr()),
        "P10_root_intervals": p10_intervals,
        "common_constituent_upper_alpha4": str(alpha4),
        "global_components": payload_components(rich_global),
    },
    "benchmark_rows": benchmark,
    "guard_rotating_fallback_action_totals": aggregate_guard,
    "scientifically_useful_candidate_guard_rows": scenario_rows,
    "noncommuting_method_comparison": noncommuting_rows,
    "benchmark_macro_source": "results/v34/benchmark_macros_v34.tex",
}
(OUT / "v34_verification.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
md = [
    "# v3.4 hostile-audit remediation verification", "",
    "## Candidate-map topology",
    f"- full RK4 interval at varpi=2: `{payload_components(direct)}`",
    f"- two-half-step interval: `{payload_components(fine)}`",
    "- nontrivial positive components are disjoint", "",
    "## Richardson extrapolation",
    f"- first defect: m={defect.index}, eta={defect.coefficient}",
    f"- global remote component: `{payload_components(rich_global)}`", "",
    "## Adaptive evidence",
    f"- boundary guard/fallback action totals: `{aggregate_guard}`",
    f"- buffered useful case: direct={buffered['direct_accepts']}, projected={buffered['projected_accepts']}, error={buffered['global_normalized_choi_error']}",
    f"- near-saddle useful case: projected={near['projected_accepts']}, locator calls={near['locator_calls']}, error={near['global_normalized_choi_error']}",
    "- all accumulated records have verified provenance", "",
    "## Comparator scope",
    "- noncommuting dense-qubit RK4/Strang data are a common-representation reference only; no universal scalability claim is made.",
]
(OUT / "v34_verification.md").write_text("\n".join(md)+"\n", encoding="utf-8")
print(json.dumps(payload, indent=2))
