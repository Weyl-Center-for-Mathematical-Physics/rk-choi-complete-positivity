#!/usr/bin/env python3
"""Independent adversarial regressions for JCP v3.1.

Usage:
    python JCP_v2_8_adversarial_regressions.py /path/to/JCP_Code_and_Data_v2_8

Exit status:
    0  all executed-candidate contract checks pass
    1  one or more blocker-class checks fail
    2  audit could not run

The script is intentionally external to the package. It checks the public/controller
behavior at the executable binary values returned by the implementation rather than
trusting internal validation labels.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import traceback
from typing import Any, Callable

import numpy as np
import sympy as sp


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("code_root", type=Path)
    parser.add_argument("--json", type=Path, default=None, help="optional JSON output path")
    args = parser.parse_args()

    root = args.code_root.resolve()
    src = root / "src"
    if not src.is_dir():
        print(f"ERROR: {src} does not exist", file=sys.stderr)
        return 2
    sys.path.insert(0, str(src))

    try:
        from rk_choi_margin.controllers import (
            _rk4_step_doubling_evaluation,
            certified_candidate_guard_step,
            run_adaptive_channel_benchmark,
        )
        from rk_choi_margin.candidates import (
            certify_candidate,
            direct_candidate,
            equal_substep_candidate,
            phase_covariant_superoperator_numeric,
        )
    except Exception:
        traceback.print_exc()
        return 2

    results: list[dict[str, Any]] = []

    def check(name: str, fn: Callable[[], dict[str, Any]]) -> None:
        try:
            payload = fn()
            ok = bool(payload.pop("ok"))
            results.append({"name": name, "status": "PASS" if ok else "FAIL", **payload})
        except Exception as exc:
            results.append(
                {
                    "name": name,
                    "status": "ERROR",
                    "exception": repr(exc),
                    "traceback": traceback.format_exc(),
                }
            )

    def direct_projection_rounding() -> dict[str, Any]:
        rows = []
        all_ok = True
        for sf in (
            sp.Rational(1, 20),
            sp.Rational(1, 10**10),
            sp.Rational(1, 10**20),
            sp.Rational(1, 10**30),
        ):
            d = certified_candidate_guard_step(
                method="rk4",
                candidate_kind="direct",
                proposed_x=1.4,
                theta=0,
                kappa=0,
                varpi=2,
                safety_fraction=sf,
            )
            if d.accepted_x is None:
                recert = None
                row_ok = d.cp_status.value != "PASS"
                f_value = None
            else:
                c = direct_candidate(
                    "rk4", total_x=d.accepted_x, theta=0, kappa=0, varpi=2
                )
                cert = certify_candidate(c)
                recert = cert.status.value
                f_value = str(next(item.value for item in cert.constraints if item.name == "F"))
                row_ok = not (d.cp_status.value == "PASS" and recert != "PASS")
            all_ok &= row_ok
            rows.append(
                {
                    "safety_fraction": str(sf),
                    "accepted_x": d.accepted_x,
                    "guard_status": d.cp_status.value,
                    "executed_binary_recertificate": recert,
                    "executed_binary_F": f_value,
                    "row_ok": row_ok,
                }
            )
        return {"ok": all_ok, "rows": rows}

    def fine_projection_rounding() -> dict[str, Any]:
        d = certified_candidate_guard_step(
            method="rk4",
            candidate_kind="equal-substeps",
            substeps=2,
            proposed_x=2.8,
            theta=0,
            kappa=0,
            varpi=2,
            safety_fraction=sp.Rational(1, 10**20),
        )
        if d.accepted_x is None:
            return {
                "ok": d.cp_status.value != "PASS",
                "accepted_x": None,
                "guard_status": d.cp_status.value,
                "executed_binary_recertificate": None,
            }
        c = equal_substep_candidate(
            "rk4",
            total_x=d.accepted_x,
            substeps=2,
            theta=0,
            kappa=0,
            varpi=2,
        )
        cert = certify_candidate(c)
        f_value = str(next(item.value for item in cert.constraints if item.name == "F"))
        ok = not (d.cp_status.value == "PASS" and cert.status.value != "PASS")
        return {
            "ok": ok,
            "accepted_x": d.accepted_x,
            "guard_status": d.cp_status.value,
            "executed_binary_recertificate": cert.status.value,
            "executed_binary_F": f_value,
        }

    def rotating_fallback_physicality() -> dict[str, Any]:
        r = run_adaptive_channel_benchmark(
            final_time=5.6,
            initial_h=5.6,
            tolerance=2.0,
            theta=0,
            kappa=0,
            varpi=1,
            policy="candidate_guard",
            fallback="rotating_frame",
        )
        accepted_nonphysical = any(
            row.used_h > 0 and row.min_choi_eigenvalue < -1e-12 for row in r.records
        )
        return {
            "ok": not accepted_nonphysical,
            "minimum_step_choi_eigenvalue": r.min_step_choi_eigenvalue,
            "records": [row.__dict__ for row in r.records],
        }

    def rotating_candidate_matrix_identity() -> dict[str, Any]:
        rows = []
        all_ok = True
        for h in (0.5, 1.0, 2.0):
            ev = _rk4_step_doubling_evaluation(
                h=h, theta=0, kappa=0, varpi=2, frame="rotating"
            )
            stored_matrix = phase_covariant_superoperator_numeric(ev.bundle.fine)
            stored_norm = float(np.linalg.norm(stored_matrix - ev.fine_matrix))
            proper_norm = stored_norm
            row_ok = stored_norm <= 1e-12
            all_ok &= row_ok
            rows.append(
                {
                    "h": h,
                    "stored_candidate_matrix_difference": stored_norm,
                    "proper_rotating_candidate_difference": proper_norm,
                    "row_ok": row_ok,
                }
            )
        return {"ok": all_ok, "rows": rows}

    def benchmark_projection_coverage() -> dict[str, Any]:
        csv_path = root / "results" / "v31" / "adaptive_benchmark_v31.csv"
        if not csv_path.exists():
            return {"ok": False, "reason": f"missing {csv_path}"}
        import csv

        with csv_path.open(newline="", encoding="utf-8") as handle:
            guarded = [
                row for row in csv.DictReader(handle) if row.get("controller_policy", row.get("policy")) == "candidate_guard"
            ]
        projection_events = sum(int(row["component_searches"]) for row in guarded)
        fallbacks = sum(int(row["fallback_steps"]) for row in guarded)
        accepted = sum(int(row["accepted_steps"]) for row in guarded)
        # This is an evidence-coverage check: the reported benchmark should include
        # at least one actual component projection if the paper presents it as such.
        return {
            "ok": projection_events > 0,
            "candidate_guard_rows": guarded,
            "projection_events": projection_events,
            "fallback_steps": fallbacks,
            "accepted_steps": accepted,
        }

    def richardson_alpha_static_check() -> dict[str, Any]:
        path = root / "scripts" / "verify_v31_extensions.py"
        text = path.read_text(encoding="utf-8")
        bad = "alpha4 = direct.components[0].upper_value" in text
        true_alpha = float(
            [r for r in sp.nroots(sp.Symbol("x")**3 - 4*sp.Symbol("x")**2 + 12*sp.Symbol("x") - 24, n=40)
             if abs(complex(r).imag) < 1e-30 and complex(r).real > 0][0]
        )
        return {
            "ok": not bad,
            "bad_assignment_present": bad,
            "true_alpha4": true_alpha,
            "script": str(path),
        }

    check("direct_projection_exact_binary_contract", direct_projection_rounding)
    check("fine_projection_exact_binary_contract", fine_projection_rounding)
    check("rotating_fallback_is_physical", rotating_fallback_physicality)
    check("rotating_candidate_matrix_identity", rotating_candidate_matrix_identity)
    check("reported_benchmark_exercises_projection", benchmark_projection_coverage)
    check("richardson_verifier_uses_correct_alpha4", richardson_alpha_static_check)

    failures = [row for row in results if row["status"] != "PASS"]
    report = {
        "code_root": str(root),
        "checks": results,
        "summary": {
            "passed": len(results) - len(failures),
            "failed_or_error": len(failures),
            "overall": "PASS" if not failures else "FAIL",
        },
    }

    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json is not None:
        args.json.write_text(rendered + "\n", encoding="utf-8")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
