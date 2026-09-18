from __future__ import annotations

from dataclasses import replace
import importlib
import math

import numpy as np
import pytest
import sympy as sp

from rk_choi_margin.candidates import (
    CandidateExecution,
    CandidateSpec,
    ExactGKSLExponentialRecipe,
    ProofKind,
    StructuralCertificate,
    StructuralTheorem,
    build_candidate_execution,
    build_exact_gksl_execution,
    build_strang_execution,
    verify_execution_provenance,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import (
    PIHistory,
    _execution_passes_common_gate,
    _rk4_step_doubling_evaluation,
    certified_candidate_guard_step,
    run_adaptive_channel_benchmark,
)
from rk_choi_margin.liouvillian import exp_exact_gksl_superoperator_numeric


def test_arbitrary_structural_builder_is_not_public_or_implemented() -> None:
    candidates = importlib.import_module("rk_choi_margin.candidates")
    public = importlib.import_module("rk_choi_margin")
    assert not hasattr(candidates, "build_structural_execution")
    assert not hasattr(public, "build_structural_execution")


def test_forged_structural_execution_fails_semantic_rederivation() -> None:
    recipe = ExactGKSLExponentialRecipe.from_values(
        gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0
    )
    bad = -np.eye(4, dtype=np.complex128)
    cert = StructuralCertificate(
        status=CPStatus.PASS,
        proof_kind=ProofKind.STRUCTURAL,
        theorem_id=StructuralTheorem.EXACT_GKSL_EXPONENTIAL,
        recipe_payload={"recipe_kind": recipe.recipe_kind.value},
        realization_checks={"realization_check": "PASS"},
        diagnostics={"forged": True},
    )
    execution = CandidateExecution(
        executed_float=1.0,
        exact_input=sp.Integer(1),
        spec=CandidateSpec(kind="fallback-exact", method="structural", frame="structural"),
        symbolic_candidate=None,
        structural_recipe=recipe,
        numerical_superoperator=bad,
        cp_certificate=cert,
        execution_id="forged",
        provenance_hash="forged",
        diagnostics={"forged": True},
    )
    assert not verify_execution_provenance(execution)
    assert not _execution_passes_common_gate(execution, 1.0, 0.0)


def test_exact_structural_builder_uses_exact_ieee_rational_time() -> None:
    h = 0.1
    execution = build_exact_gksl_execution(
        total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0
    )
    assert execution.exact_input == sp.Rational(*h.as_integer_ratio())
    rebuilt = exp_exact_gksl_superoperator_numeric(
        h=execution.exact_input,
        gamma=1,
        theta=0,
        kappa=0,
        omega_z=2,
        omega_x=0,
    )
    assert np.array_equal(execution.numerical_superoperator, rebuilt)
    assert execution.proof_kind is ProofKind.STRUCTURAL
    assert execution.cp_status is CPStatus.PASS
    assert verify_execution_provenance(execution)


def test_direct_guard_forwards_literal_estimator_execution() -> None:
    evaluation = _rk4_step_doubling_evaluation(
        h=2.0, theta=0, kappa=0, varpi=2, frame="lab"
    )
    decision = certified_candidate_guard_step(
        method="rk4",
        candidate_kind="equal-substeps",
        substeps=2,
        proposed_x=2.0,
        theta=0,
        kappa=0,
        varpi=2,
        proposed_execution=evaluation.fine_execution,
    )
    assert decision.proposed_execution is evaluation.fine_execution
    assert decision.accepted_execution is evaluation.fine_execution


def test_projection_forwards_literal_guard_execution_into_error_estimator() -> None:
    decision = certified_candidate_guard_step(
        method="rk4",
        candidate_kind="equal-substeps",
        substeps=2,
        proposed_x=2.8,
        theta=0,
        kappa=0,
        varpi=2,
    )
    assert decision.accepted_execution is not None
    evaluation = _rk4_step_doubling_evaluation(
        h=decision.accepted_execution.executed_float,
        theta=0,
        kappa=0,
        varpi=2,
        frame="lab",
        fine_execution=decision.accepted_execution,
    )
    assert evaluation.fine_execution is decision.accepted_execution


def test_execution_record_is_deeply_immutable() -> None:
    execution = build_candidate_execution(
        CandidateSpec("equal-substeps", "rk4", 0, 0, 2, 2, "lab"),
        total_x=1.0,
    )
    with pytest.raises(TypeError):
        execution.diagnostics["new"] = "mutation"  # type: ignore[index]
    with pytest.raises(ValueError):
        execution.numerical_superoperator.setflags(write=True)
    assert verify_execution_provenance(execution)


def test_common_gate_requires_explicit_error_observation() -> None:
    execution = build_exact_gksl_execution(
        total_x=1.0, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0
    )
    assert _execution_passes_common_gate(execution, 1.0, 0.0)
    assert not _execution_passes_common_gate(execution, 1e-3, 1e-2)
    assert not _execution_passes_common_gate(execution, 1.0, math.inf)


def test_pi_history_is_candidate_family_aware() -> None:
    history = PIHistory(10.0)
    assert not history.transition("exact-gksl-exponential")
    history.reset_family("exact-gksl-exponential")
    assert history.previous("exact-gksl-exponential") == 10.0
    assert history.transition("lab:equal-substeps:rk4:2")
    assert history.previous("lab:equal-substeps:rk4:2") == 10.0
    # The exact fallback's zero error never contaminates direct-RK history.
    assert history.previous("lab:equal-substeps:rk4:2") != 1e-300


def test_adaptive_records_bind_certified_and_accumulated_execution_ids() -> None:
    result = run_adaptive_channel_benchmark(
        final_time=2.8,
        initial_h=2.8,
        tolerance=10.0,
        theta=0,
        kappa=0,
        varpi=2,
        policy="candidate_guard",
        fallback="rotating_frame",
    )
    accepted = [row for row in result.records if row.used_h > 0]
    assert accepted
    assert all(row.execution_id == row.certified_execution_id for row in accepted)
    assert all(row.provenance_verified for row in accepted)


def test_family_reset_is_not_overwritten_in_fallback_to_rk_transition() -> None:
    result = run_adaptive_channel_benchmark(
        final_time=5.0,
        initial_h=1.0,
        tolerance=10.0,
        theta=0,
        kappa=0,
        varpi=2,
        policy="candidate_guard",
        fallback="exact",
    )
    assert result.stats["pi_memory_resets"] >= 1
    assert len(result.records) >= 3
    # The v2.9 bug generated a minimum-shrink proposal of 0.4.  A reset family
    # now holds/uses tolerance-scale history instead.
    assert result.records[2].proposed_h > 0.4 + 1e-12


def test_typed_strang_recipe_is_structurally_verified() -> None:
    execution = build_strang_execution(
        total_x=0.5,
        gamma=1,
        theta="0.001",
        kappa="0.001",
        omega_z=2,
        omega_x="0.1",
        substeps=2,
    )
    assert execution.cp_status is CPStatus.PASS
    assert execution.proof_kind is ProofKind.STRUCTURAL
    assert verify_execution_provenance(execution)


def test_structural_certificate_cannot_be_mutated_to_rebind_matrix() -> None:
    execution = build_exact_gksl_execution(
        total_x=0.5, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0
    )
    bad_matrix = -np.eye(4, dtype=np.complex128)
    forged = replace(execution, numerical_superoperator=bad_matrix)
    assert not verify_execution_provenance(forged)
