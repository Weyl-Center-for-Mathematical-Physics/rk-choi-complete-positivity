#!/usr/bin/env python3
from __future__ import annotations

"""Second-order adversarial regression suite for JCP v3.0.

Every check targets a blocker found in the hostile v2.9 audit.  The suite exits
nonzero on any unresolved blocker and reads only current-version outputs.
"""

from dataclasses import replace
import importlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from rk_choi_margin.candidates import (
    CandidateExecution,
    CandidateSpec,
    ExactGKSLExponentialRecipe,
    ProofKind,
    StructuralCertificate,
    StructuralTheorem,
    build_exact_gksl_execution,
    verify_execution_provenance,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import (
    _execution_passes_common_gate,
    _rk4_step_doubling_evaluation,
    certified_candidate_guard_step,
    run_adaptive_channel_benchmark,
)
from rk_choi_margin.liouvillian import exp_exact_gksl_superoperator_numeric, min_choi_eigenvalue_numeric

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT.parent / "JCP_LaTeX_Source_v3_0"
OUT = ROOT / "independent_validation" / "second_order_adversarial_v30_report.json"
checks: list[dict[str, object]] = []


def add(name: str, condition: bool, evidence: object = None) -> None:
    row = {"name": name, "pass": bool(condition), "evidence": evidence}
    checks.append(row)
    if not condition:
        raise AssertionError(f"{name}: {evidence}")


# 1. Free-form structural PASS attack surface removed.
mod = importlib.import_module("rk_choi_margin.candidates")
add("no arbitrary matrix theorem builder", not hasattr(mod, "build_structural_execution"))

# 2. A forged structural execution cannot pass semantic rederivation or gate.
recipe = ExactGKSLExponentialRecipe.from_values(gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0)
bad = -np.eye(4, dtype=np.complex128)
cert = StructuralCertificate(
    CPStatus.PASS,
    ProofKind.STRUCTURAL,
    StructuralTheorem.EXACT_GKSL_EXPONENTIAL,
    {"recipe_kind": recipe.recipe_kind.value},
    {"realization_check": "PASS"},
    {"forged": True},
)
forged = CandidateExecution(
    1.0,
    sp.Integer(1),
    CandidateSpec(kind="fallback-exact", method="structural", frame="structural"),
    None,
    recipe,
    bad,
    cert,
    "forged",
    "forged",
    {"forged": True},
)
add("forged structural provenance fails", not verify_execution_provenance(forged))
add("forged structural common gate fails", not _execution_passes_common_gate(forged, 1.0, 0.0), min_choi_eigenvalue_numeric(bad))

# 3. Literal identity on direct acceptance.
evaluation = _rk4_step_doubling_evaluation(h=2.0, theta=0, kappa=0, varpi=2, frame="lab")
decision = certified_candidate_guard_step(
    method="rk4", candidate_kind="equal-substeps", substeps=2,
    proposed_x=2.0, theta=0, kappa=0, varpi=2,
    proposed_execution=evaluation.fine_execution,
)
add("direct estimator guard object identity", decision.accepted_execution is evaluation.fine_execution)

# 4. Literal identity on projection/error recomputation.
projection = certified_candidate_guard_step(
    method="rk4", candidate_kind="equal-substeps", substeps=2,
    proposed_x=2.8, theta=0, kappa=0, varpi=2,
)
assert projection.accepted_execution is not None
recomputed = _rk4_step_doubling_evaluation(
    h=projection.accepted_execution.executed_float,
    theta=0, kappa=0, varpi=2, frame="lab",
    fine_execution=projection.accepted_execution,
)
add("projection guard estimator object identity", recomputed.fine_execution is projection.accepted_execution)

# 5. Exact fallback uses stored IEEE rational.
for h in (0.1, 1.1, 2.8):
    execution = build_exact_gksl_execution(total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0)
    exact = sp.Rational(*h.as_integer_ratio())
    rebuilt = exp_exact_gksl_superoperator_numeric(h=exact, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0)
    add(f"fallback exact input h={h}", execution.exact_input == exact)
    add(f"fallback exact map h={h}", np.array_equal(execution.numerical_superoperator, rebuilt))

# 6. PI family reset persists; no 0.4 collapse.
memory = run_adaptive_channel_benchmark(
    final_time=5.0, initial_h=1.0, tolerance=10.0,
    theta=0, kappa=0, varpi=2, policy="candidate_guard", fallback="exact",
)
add("PI family transition recorded", memory.stats["pi_memory_resets"] >= 1, memory.stats)
add("PI post-transition proposal not minimum shrink", len(memory.records) >= 3 and memory.records[2].proposed_h > 0.4 + 1e-12, [r.proposed_h for r in memory.records[:3]])

# 7. Current-version validation reads current-version output.
current_script = (ROOT / "scripts" / "verify_v30_extensions.py").read_text(encoding="utf-8")
add("current verifier writes v30", 'results" / "v30"' in current_script and "adaptive_benchmark_v30.csv" in current_script)
add("current verifier does not write historical v29", 'results" / "v29"' not in current_script and "adaptive_benchmark_v29.csv" not in current_script)

# 8. Deep immutability.
add("nested diagnostics immutable", type(evaluation.fine_execution.diagnostics).__name__ == "mappingproxy")
try:
    evaluation.fine_execution.diagnostics["attack"] = "x"  # type: ignore[index]
    mapping_blocked = False
except TypeError:
    mapping_blocked = True
add("diagnostics mutation blocked", mapping_blocked)
try:
    evaluation.fine_execution.numerical_superoperator.setflags(write=True)
    matrix_blocked = False
except ValueError:
    matrix_blocked = True
add("matrix mutation blocked", matrix_blocked)

# 9. Every accepted adaptive record has matching certified/accumulated identity.
identity_run = run_adaptive_channel_benchmark(
    final_time=2.8, initial_h=2.8, tolerance=10.0,
    theta=0, kappa=0, varpi=2, policy="candidate_guard", fallback="rotating_frame",
)
accepted = [r for r in identity_run.records if r.used_h > 0]
add("adaptive execution identity", all(r.execution_id == r.certified_execution_id for r in accepted))
add("adaptive provenance", all(r.provenance_verified for r in accepted))

# 10. Manuscript wording and grammar.
if SOURCE_ROOT.exists():
    candidate_text = (SOURCE_ROOT / "sections" / "06_candidate_maps.tex").read_text(encoding="utf-8")
    article_text = (SOURCE_ROOT / "jcp_article.tex").read_text(encoding="utf-8")
    step_text = (SOURCE_ROOT / "sections" / "07_step_control.tex").read_text(encoding="utf-8")
    compare_text = (SOURCE_ROOT / "sections" / "09_comparison.tex").read_text(encoding="utf-8")
    add("factor two scoped to n=1", "For the full-step/two-half-step pair ($n=1$)" in candidate_text)
    add("abstract grammar corrected", "archive support all reported results" in article_text and "archive supports all reported results" not in article_text)
    add("error observation separate from record", "Local-error observations are stored separately" in step_text)
    add("typed structural proof described", "typed structural recipes" in step_text and "\\textsc{structural}" in step_text)
    add("verification distinguishes proof kinds", "\\textsc{algebraic}" in compare_text and "\\textsc{structural}" in compare_text)

report = {
    "version": "3.0",
    "checks": checks,
    "passed": sum(bool(row["pass"]) for row in checks),
    "total": len(checks),
    "technical_decision": "GO" if all(bool(row["pass"]) for row in checks) else "NO-GO",
}
OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(f"Second-order adversarial v3.0: {report['passed']}/{report['total']} checks passed.")
print(f"Wrote {OUT}")
