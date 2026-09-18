#!/usr/bin/env python3
from __future__ import annotations

"""Independent package-level v3.0 executed-candidate contract checks."""

from dataclasses import replace
import importlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from rk_choi_margin.candidates import (
    CandidateSpec,
    ProofKind,
    build_candidate_execution,
    build_exact_gksl_execution,
    build_strang_execution,
    certify_candidate,
    direct_candidate,
    equal_substep_candidate,
    phase_covariant_superoperator_numeric,
    verify_execution_provenance,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import (
    _execution_passes_common_gate,
    _rk4_step_doubling_evaluation,
    certified_candidate_guard_step,
    run_adaptive_channel_benchmark,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "independent_validation" / "executed_candidate_contract_v30_report.json"
checks: list[dict[str, object]] = []


def add(name: str, condition: bool, detail: str = "") -> None:
    checks.append({"name": name, "pass": bool(condition), "detail": detail})
    if not condition:
        raise AssertionError(f"{name}: {detail}")


# Public attack surface: no arbitrary matrix+theorem structural builder.
module = importlib.import_module("rk_choi_margin.candidates")
add("arbitrary structural builder absent", not hasattr(module, "build_structural_execution"))

# Projection and exact-binary recertification across severe interior fractions.
for kind, proposed, substeps in (("direct", 1.4, 1), ("equal-substeps", 2.8, 2)):
    for exponent in (2, 10, 20, 30):
        decision = certified_candidate_guard_step(
            method="rk4",
            candidate_kind=kind,
            substeps=substeps,
            proposed_x=proposed,
            theta=0,
            kappa=0,
            varpi=2,
            safety_fraction=sp.Rational(1, 10**exponent),
        )
        execution = decision.accepted_execution
        add(f"{kind} projection pass 1e-{exponent}", execution is not None and decision.cp_status is CPStatus.PASS)
        assert execution is not None
        add(f"{kind} provenance 1e-{exponent}", verify_execution_provenance(execution))
        candidate = (
            direct_candidate("rk4", total_x=execution.executed_float, theta=0, kappa=0, varpi=2)
            if kind == "direct"
            else equal_substep_candidate("rk4", total_x=execution.executed_float, substeps=2, theta=0, kappa=0, varpi=2)
        )
        add(f"{kind} exact executed recertificate 1e-{exponent}", certify_candidate(candidate).status is CPStatus.PASS)

# Direct path forwards the estimator execution literally.
for h in (1.5, 2.0, 2.4):
    ev = _rk4_step_doubling_evaluation(h=h, theta=0, kappa=0, varpi=2, frame="lab")
    decision = certified_candidate_guard_step(
        method="rk4", candidate_kind="equal-substeps", substeps=2,
        proposed_x=h, theta=0, kappa=0, varpi=2,
        proposed_execution=ev.fine_execution,
    )
    if ev.fine_execution.cp_status is CPStatus.PASS:
        add(f"direct literal identity h={h}", decision.accepted_execution is ev.fine_execution)

# Projection record is reused by the estimator, not reconstructed.
decision = certified_candidate_guard_step(
    method="rk4", candidate_kind="equal-substeps", substeps=2,
    proposed_x=2.8, theta=0, kappa=0, varpi=2,
)
assert decision.accepted_execution is not None
projected = _rk4_step_doubling_evaluation(
    h=decision.accepted_execution.executed_float,
    theta=0, kappa=0, varpi=2, frame="lab",
    fine_execution=decision.accepted_execution,
)
add("projection literal identity", projected.fine_execution is decision.accepted_execution)

# Coarse/fine/extrapolated identities and distinct hashes.
for frame in ("lab", "rotating"):
    for h in (0.5, 1.0, 2.0):
        ev = _rk4_step_doubling_evaluation(h=h, theta=0, kappa=0, varpi=2, frame=frame)
        add(f"{frame} coarse provenance h={h}", verify_execution_provenance(ev.coarse_execution))
        add(f"{frame} fine provenance h={h}", verify_execution_provenance(ev.fine_execution))
        add(f"{frame} fine matrix identity h={h}", np.array_equal(phase_covariant_superoperator_numeric(ev.bundle.fine), ev.fine_matrix))
        if ev.extrapolated_execution is not None:
            add(f"{frame} extrapolated provenance h={h}", verify_execution_provenance(ev.extrapolated_execution))
            add(f"{frame} distinct execution IDs h={h}", len({ev.coarse_execution.execution_id, ev.fine_execution.execution_id, ev.extrapolated_execution.execution_id}) == 3)

# Typed structural recipes derive their maps and pass semantic rederivation.
for h in (0.1, 1.1, 2.8):
    exact = build_exact_gksl_execution(total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0)
    add(f"exact structural pass h={h}", exact.cp_status is CPStatus.PASS and exact.proof_kind is ProofKind.STRUCTURAL)
    add(f"exact binary time h={h}", exact.exact_input == sp.Rational(*h.as_integer_ratio()))
    add(f"exact provenance h={h}", verify_execution_provenance(exact))

for substeps in (1, 2, 4):
    strang = build_strang_execution(
        total_x=0.5, gamma=1, theta="0.001", kappa="0.001",
        omega_z=2, omega_x="0.1", substeps=substeps,
    )
    add(f"strang structural pass n={substeps}", strang.cp_status is CPStatus.PASS and strang.proof_kind is ProofKind.STRUCTURAL)
    add(f"strang provenance n={substeps}", verify_execution_provenance(strang))

# Deep immutability.
execution = build_candidate_execution(CandidateSpec("equal-substeps", "rk4", 0, 0, 2, 2, "lab"), total_x=1.0)
try:
    execution.diagnostics["attack"] = "x"  # type: ignore[index]
    immutable_mapping = False
except TypeError:
    immutable_mapping = True
add("diagnostics deeply immutable", immutable_mapping)
try:
    execution.numerical_superoperator.setflags(write=True)
    immutable_matrix = False
except ValueError:
    immutable_matrix = True
add("matrix deeply immutable", immutable_matrix)

# Adaptive paths: every accepted map has the same certified/accumulated ID.
for final_time, initial_h, tol, varpi in ((2.0, 1.0, 1e-2, 2), (2.8, 2.8, 10.0, 2), (5.6, 5.6, 2.0, 1), (5.0, 1.0, 10.0, 2)):
    result = run_adaptive_channel_benchmark(
        final_time=final_time, initial_h=initial_h, tolerance=tol,
        theta=0, kappa=0, varpi=varpi,
        policy="candidate_guard", fallback="rotating_frame" if final_time != 5.0 else "exact",
    )
    accepted = [row for row in result.records if row.used_h > 0]
    add(f"benchmark provenance {final_time}", all(row.provenance_verified for row in accepted))
    add(f"benchmark execution continuity {final_time}", all(row.execution_id == row.certified_execution_id for row in accepted))
    add(f"benchmark physicality {final_time}", all(row.min_choi_eigenvalue > -1e-12 for row in accepted))

report = {"checks": checks, "passed": sum(row["pass"] for row in checks), "total": len(checks)}
OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(f"Executed-candidate contract v3.0: {report['passed']}/{report['total']} checks passed.")
print(f"Wrote {OUT}")
