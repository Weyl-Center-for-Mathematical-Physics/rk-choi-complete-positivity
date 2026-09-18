#!/usr/bin/env python3
"""Third-order hostile regressions for the JCP v3.1 submission package.

This script does not modify the submission. It probes two publication-facing
contracts:

1. Typed structural certificates must not return PASS outside the theorem's
   forward-time domain.
2. Benchmark prose must agree with the current v3.1 machine-readable results,
   and the plotted policy labels must describe the path actually exercised.

It also records lower-priority public-API validation gaps. Findings are emitted
as JSON so that a revision agent can convert each one into a failing regression.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


def _json_value(value: Any) -> Any:
    if hasattr(value, "value"):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (float, int, str, bool)) or value is None:
        return value
    return str(value)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--code-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    code_root = args.code_root.resolve()
    source_root = args.source_root.resolve()
    sys.path.insert(0, str(code_root / "src"))

    import numpy as np
    from rk_choi_margin import (
        ExactGKSLExponentialRecipe,
        ExactSubflowStrangRecipe,
        build_exact_gksl_execution,
        build_strang_execution,
        verify_execution_provenance,
    )
    from rk_choi_margin.controllers import (
        _execution_passes_common_gate,
        run_adaptive_channel_benchmark,
    )

    findings: list[dict[str, Any]] = []
    raw: dict[str, Any] = {}

    # ------------------------------------------------------------------
    # P0: structural theorem-domain violation at tiny negative time.
    # ------------------------------------------------------------------
    negative_cases: list[dict[str, Any]] = []
    for h in (-1e-12, -1e-13, -1e-14, -1e-15, -1e-16):
        for builder_name, builder in (
            ("exact-gksl", build_exact_gksl_execution),
            ("strang", build_strang_execution),
        ):
            kwargs = dict(
                total_x=h,
                gamma=1,
                theta=0,
                kappa=0,
                omega_z=2,
                omega_x=0,
                dps=80,
            )
            if builder_name == "strang":
                kwargs["substeps"] = 1
            analytic_population_term = -math.expm1(-h)
            try:
                execution = builder(**kwargs)
            except (TypeError, ValueError) as exc:
                negative_cases.append(
                    {
                        "builder": builder_name,
                        "h": h,
                        "outcome": "rejected",
                        "exception": type(exc).__name__,
                        "message": str(exc),
                        "analytic_population_term_1_minus_exp_minus_h": analytic_population_term,
                    }
                )
                continue
            checks = execution.cp_certificate.realization_checks
            negative_cases.append(
                {
                    "builder": builder_name,
                    "h": h,
                    "outcome": "accepted",
                    "status": execution.cp_status.value,
                    "minimum_realized_choi_eigenvalue": float(
                        checks["minimum_realized_choi_eigenvalue"]
                    ),
                    "realization_tolerance": float(checks["realization_tolerance"]),
                    "analytic_population_term_1_minus_exp_minus_h": analytic_population_term,
                    "provenance_valid": bool(verify_execution_provenance(execution)),
                    "common_gate_accepts_zero_error": bool(
                        _execution_passes_common_gate(execution, 1.0, 0.0)
                    ),
                }
            )

    false_passes = [
        row
        for row in negative_cases
        if row.get("outcome") == "accepted"
        and row.get("status") == "PASS"
        and row["analytic_population_term_1_minus_exp_minus_h"] < 0
    ]
    raw["negative_structural_time_cases"] = negative_cases
    if false_passes:
        findings.append(
            {
                "id": "P0-STRUCTURAL-NEGATIVE-TIME-FALSE-PASS",
                "severity": "P0",
                "status": "CONFIRMED",
                "summary": (
                    "Typed exact-GKSL and Strang builders can return PASS for "
                    "tiny negative dissipative steps, although the forward GKSL "
                    "semigroup theorem does not apply and the realized Choi matrix "
                    "has a negative eigenvalue."
                ),
                "count": len(false_passes),
                "examples": false_passes,
                "required_fix": (
                    "Reject total_x < 0 at the public recipe/builder/certificate "
                    "boundary (or implement an explicit theorem-domain branch); "
                    "bind the nonnegative-time hypothesis into the structural proof "
                    "payload; add exact and Strang negative-time regressions."
                ),
            }
        )

    # Positive control: ordinary forward times should remain accepted.
    positive_controls = []
    for h in (0.0, 1e-16, 0.1, 1.0):
        execution = build_exact_gksl_execution(
            total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0
        )
        positive_controls.append(
            {
                "h": h,
                "status": execution.cp_status.value,
                "provenance_valid": verify_execution_provenance(execution),
            }
        )
    raw["positive_structural_controls"] = positive_controls

    # ------------------------------------------------------------------
    # P0/P1: current manuscript narrative versus current v3.1 CSV.
    # ------------------------------------------------------------------
    csv_path = code_root / "results" / "v31" / "adaptive_benchmark_v31.csv"
    tex_path = source_root / "sections" / "07_step_control.tex"
    rows = list(csv.DictReader(csv_path.open(newline="", encoding="utf-8")))
    tex = tex_path.read_text(encoding="utf-8")

    sweep = [r for r in rows if r["scenario"] == "boundary_tolerance_sweep"]
    rotating = [r for r in sweep if r["policy"] == "rotating_frame"]
    guarded = [r for r in sweep if r.get("controller_policy", r["policy"]) == "candidate_guard"]
    target_rotating = next(r for r in rotating if float(r["tolerance"]) == 3e-4)
    target_error = float(target_rotating["global_normalized_choi_error"])
    below = [r for r in guarded if float(r["global_normalized_choi_error"]) < target_error]
    first_below = below[0] if below else None
    near_equal = min(
        guarded,
        key=lambda r: abs(float(r["global_normalized_choi_error"]) - target_error),
    )

    stale_phrase_present = bool(
        re.search(
            r"3\.7\\times10\^\{-5\}.*96 stage units",
            tex,
            flags=re.DOTALL,
        )
    )
    macro_path = source_root / "generated" / "benchmark_macros_v31.tex"
    macro_present = macro_path.exists() and "BenchmarkRotatingError" in macro_path.read_text(encoding="utf-8")
    current_comparison = {
        "rotating_error": target_error,
        "rotating_stage_units": int(target_rotating["rhs_stage_evaluations"]),
        "first_guarded_below_error": (
            float(first_below["global_normalized_choi_error"])
            if first_below
            else None
        ),
        "first_guarded_below_stage_units": (
            int(first_below["rhs_stage_evaluations"]) if first_below else None
        ),
        "first_guarded_below_tolerance": (
            float(first_below["tolerance"]) if first_below else None
        ),
        "nearest_guarded_error": float(near_equal["global_normalized_choi_error"]),
        "nearest_guarded_stage_units": int(near_equal["rhs_stage_evaluations"]),
        "nearest_guarded_tolerance": float(near_equal["tolerance"]),
        "stale_phrase_present": stale_phrase_present,
        "generated_macro_present": macro_present,
    }
    raw["current_v31_benchmark_comparison"] = current_comparison

    v28_path = code_root / "results" / "v28" / "adaptive_benchmark_v28.csv"
    stale_origin = None
    if v28_path.exists():
        v28_rows = list(csv.DictReader(v28_path.open(newline="", encoding="utf-8")))
        for r in v28_rows:
            if (
                r.get("scenario") in {"boundary_tolerance_sweep", "tolerance_sweep"}
                and r.get("policy") == "candidate_guard"
                and abs(float(r.get("global_normalized_choi_error", "nan")) - 3.7e-5)
                < 1e-5
                and int(r.get("rhs_stage_evaluations", "-1")) == 96
            ):
                stale_origin = r
                break
    raw["probable_v28_origin_of_stale_sentence"] = stale_origin

    if stale_phrase_present or not macro_present:
        findings.append(
            {
                "id": "P0-MANUSCRIPT-BENCHMARK-NUMBERS-STALE",
                "severity": "P0",
                "status": "CONFIRMED",
                "summary": (
                    "Section 6 reports a guarded point at 3.7e-5 in 96 stage "
                    "units or lacks generated v3.1 benchmark macros tied to adaptive_benchmark_v31.csv."
                ),
                "current_v31_values": current_comparison,
                "probable_historical_origin": stale_origin,
                "required_fix": (
                    "Replace prose from current v3.1 results, preferably via "
                    "generated TeX macros, and add a source-to-CSV consistency test."
                ),
            }
        )

    # ------------------------------------------------------------------
    # P1: action composition of the main candidate_guard sweep.
    # ------------------------------------------------------------------
    count_fields = (
        "accepted_steps",
        "direct_accepts",
        "projected_accepts",
        "rotating_accepts",
        "strang_accepts",
        "exact_fallback_accepts",
        "component_searches",
        "fallback_steps",
    )
    composition = {field: sum(int(r[field]) for r in guarded) for field in count_fields}
    raw["candidate_guard_main_sweep_action_composition"] = composition
    guarded_labels = {r["policy"] for r in guarded}
    if (
        composition["accepted_steps"] > 0
        and composition["direct_accepts"] == 0
        and composition["projected_accepts"] == 0
        and composition["rotating_accepts"] == composition["accepted_steps"]
        and guarded_labels != {"guard_rotating_fallback"}
    ):
        findings.append(
            {
                "id": "P1-BENCHMARK-LABEL-HIDES-FALLBACK-COMPOSITION",
                "severity": "P1",
                "status": "CONFIRMED",
                "summary": (
                    "Every accepted step on the main candidate_guard tolerance "
                    "sweep is a rotating-frame fallback; the curve does not measure "
                    "direct or projected laboratory-frame guard performance."
                ),
                "counts": composition,
                "required_fix": (
                    "Relabel the plotted policy as guard + rotating fallback, report "
                    "action counts in the manuscript/caption, and retain the separate "
                    "projection_coverage row as functionality evidence rather than "
                    "performance evidence."
                ),
            }
        )

    # ------------------------------------------------------------------
    # P2: public input validation / silent integer truncation.
    # ------------------------------------------------------------------
    api_cases: list[dict[str, Any]] = []
    for value in (1.5, 2.9):
        try:
            recipe = ExactSubflowStrangRecipe.from_values(substeps=value)
            api_cases.append(
                {
                    "case": "noninteger_substeps",
                    "input": value,
                    "outcome": "accepted",
                    "stored": recipe.substeps,
                }
            )
        except Exception as exc:  # pragma: no cover - expected after repair
            api_cases.append(
                {
                    "case": "noninteger_substeps",
                    "input": value,
                    "outcome": "rejected",
                    "exception": type(exc).__name__,
                }
            )
    for value in (1.5, 0, -5):
        try:
            recipe = ExactGKSLExponentialRecipe.from_values(dps=value)
            execution = build_exact_gksl_execution(total_x=0.1, dps=value)
            api_cases.append(
                {
                    "case": "invalid_precision",
                    "input": value,
                    "outcome": "accepted",
                    "stored": recipe.dps,
                    "certificate_status": execution.cp_status.value,
                }
            )
        except Exception as exc:  # pragma: no cover - expected after repair
            api_cases.append(
                {
                    "case": "invalid_precision",
                    "input": value,
                    "outcome": "rejected",
                    "exception": type(exc).__name__,
                }
            )

    controller_cases = (
        {"final_time": -1.0, "initial_h": 0.1, "tolerance": 0.01, "min_h": 1e-12},
        {"final_time": 1.0, "initial_h": -0.1, "tolerance": 0.01, "min_h": 1e-12},
        {"final_time": 1.0, "initial_h": 0.1, "tolerance": -0.01, "min_h": 1e-12},
        {"final_time": 1.0, "initial_h": 0.1, "tolerance": 0.01, "min_h": -1e-12},
    )
    for case in controller_cases:
        try:
            result = run_adaptive_channel_benchmark(
                theta=0,
                kappa=0,
                varpi=2,
                policy="candidate_guard",
                fallback="exact",
                max_steps=3,
                **case,
            )
            api_cases.append(
                {
                    "case": "invalid_controller_input",
                    "input": case,
                    "outcome": "returned",
                    "accepted_steps": len(result.records),
                    "global_error": result.global_normalized_choi_error,
                }
            )
        except Exception as exc:
            api_cases.append(
                {
                    "case": "invalid_controller_input",
                    "input": case,
                    "outcome": "exception",
                    "exception": type(exc).__name__,
                    "message": str(exc),
                }
            )
    raw["public_api_validation_cases"] = api_cases

    silent_accepts = [r for r in api_cases if r["outcome"] in {"accepted", "returned"}]
    if silent_accepts:
        findings.append(
            {
                "id": "P2-PUBLIC-API-INPUT-VALIDATION",
                "severity": "P2",
                "status": "CONFIRMED",
                "summary": (
                    "Public recipe/controller inputs permit silent integer truncation "
                    "or nonsensical nonpositive domains instead of failing early."
                ),
                "examples": silent_accepts,
                "required_fix": (
                    "Require finite positive final_time, initial_h, tolerance, min_h, "
                    "and dps; require exact positive integers for substeps/dps; reject "
                    "booleans and lossy coercions."
                ),
            }
        )

    report = {
        "script": Path(__file__).name,
        "code_root": str(code_root),
        "source_root": str(source_root),
        "blockers_detected": any(f["severity"] == "P0" for f in findings),
        "finding_count": len(findings),
        "findings": findings,
        "raw_evidence": raw,
    }
    args.json_out.write_text(
        json.dumps(report, indent=2, sort_keys=True, default=_json_value) + "\n",
        encoding="utf-8",
    )

    print(json.dumps({
        "blockers_detected": report["blockers_detected"],
        "finding_count": report["finding_count"],
        "finding_ids": [f["id"] for f in findings],
        "json_out": str(args.json_out),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
