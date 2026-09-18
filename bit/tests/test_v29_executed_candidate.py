from __future__ import annotations

import math
import numpy as np
import sympy as sp

from rk_choi_margin.candidates import (
    CandidateSpec,
    build_candidate_execution,
    certify_candidate,
    direct_candidate,
    equal_substep_candidate,
    phase_covariant_superoperator_numeric,
    verify_execution_provenance,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import (
    _rk4_step_doubling_evaluation,
    certified_candidate_guard_step,
    run_adaptive_channel_benchmark,
)
from rk_choi_margin.liouvillian import min_choi_eigenvalue_numeric


def test_projection_recopies_exact_executed_binary_value() -> None:
    for safety in (sp.Rational(1, 20), sp.Rational(1, 10**20), sp.Rational(1, 10**30)):
        decision = certified_candidate_guard_step(
            method="rk4",
            candidate_kind="direct",
            proposed_x=1.4,
            theta=0,
            kappa=0,
            varpi=2,
            safety_fraction=safety,
        )
        assert decision.cp_status is CPStatus.PASS
        assert decision.accepted_execution is not None
        assert decision.accepted_execution.cp_status is CPStatus.PASS
        assert verify_execution_provenance(decision.accepted_execution)
        recandidate = direct_candidate(
            "rk4", total_x=decision.accepted_x, theta=0, kappa=0, varpi=2
        )
        assert certify_candidate(recandidate).status is CPStatus.PASS


def test_fine_projection_recopies_exact_executed_binary_value() -> None:
    decision = certified_candidate_guard_step(
        method="rk4",
        candidate_kind="equal-substeps",
        substeps=2,
        proposed_x=2.8,
        theta=0,
        kappa=0,
        varpi=2,
        safety_fraction=sp.Rational(1, 10**20),
    )
    assert decision.accepted_execution is not None
    assert decision.accepted_execution.cp_status is CPStatus.PASS
    assert verify_execution_provenance(decision.accepted_execution)
    recandidate = equal_substep_candidate(
        "rk4",
        total_x=decision.accepted_x,
        substeps=2,
        theta=0,
        kappa=0,
        varpi=2,
    )
    assert certify_candidate(recandidate).status is CPStatus.PASS


def test_rotating_candidate_metadata_matches_executed_matrix() -> None:
    for h in (0.5, 1.0, 2.0):
        evaluation = _rk4_step_doubling_evaluation(
            h=h, theta=0, kappa=0, varpi=2, frame="rotating"
        )
        symbolic_matrix = phase_covariant_superoperator_numeric(evaluation.bundle.fine)
        assert np.linalg.norm(symbolic_matrix - evaluation.fine_matrix) < 1e-13
        assert verify_execution_provenance(evaluation.fine_execution)


def test_candidate_execution_matrix_is_read_only_and_provenance_detects_mutation() -> None:
    execution = build_candidate_execution(
        CandidateSpec("equal-substeps", "rk4", 0, 0, 2, 2, "lab"),
        total_x=1.0,
    )
    assert not execution.numerical_superoperator.flags.writeable
    assert verify_execution_provenance(execution)
    copy = np.array(execution.numerical_superoperator, copy=True)
    copy[0, 0] += 1e-12
    from rk_choi_margin.candidates import candidate_provenance_hash
    changed = candidate_provenance_hash(
        executed_float=execution.executed_float,
        exact_input=execution.exact_input,
        spec=execution.spec,
        symbolic_candidate=execution.symbolic_candidate,
        matrix=copy,
        error_estimate=execution.local_error_estimate,
        certificate=execution.cp_certificate,
    )
    assert changed != execution.provenance_hash


def test_rotating_fallback_must_pass_same_cp_and_error_gate() -> None:
    result = run_adaptive_channel_benchmark(
        final_time=5.6,
        initial_h=5.6,
        tolerance=2.0,
        theta=0,
        kappa=0,
        varpi=1,
        policy="candidate_guard",
        fallback="rotating_frame",
    )
    assert all(row.provenance_verified for row in result.records)
    assert all(row.min_choi_eigenvalue > -1e-12 for row in result.records)
    # At this large step the rotating RK candidate is not physical, so the
    # common gate must switch to the exact fallback rather than accept it.
    assert result.records[0].action == "fallback-exact"


def test_adaptive_benchmark_accumulates_the_certified_execution() -> None:
    result = run_adaptive_channel_benchmark(
        final_time=2.0,
        initial_h=1.0,
        tolerance=0.01,
        theta=0,
        kappa=0,
        varpi=2,
        policy="candidate_guard",
        fallback="rotating_frame",
    )
    assert result.records
    assert all(row.provenance_verified for row in result.records)
    assert all(
        row.candidate_kind == row.certified_kind or row.action.startswith("fallback")
        for row in result.records
    )
    assert result.min_step_choi_eigenvalue > -1e-11


def test_projection_path_is_exercised_and_physical() -> None:
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
    assert result.stats["component_searches"] >= 1
    projected = [row for row in result.records if row.action == "guard-projection"]
    assert projected
    assert all(row.min_choi_eigenvalue > -1e-12 for row in projected)


def test_full_and_fine_candidates_have_opposite_status_at_worked_steps() -> None:
    e1 = _rk4_step_doubling_evaluation(h=1.0, theta=0, kappa=0, varpi=2, frame="lab")
    assert e1.coarse_execution.cp_status is CPStatus.PASS
    assert e1.fine_execution.cp_status is CPStatus.FAIL
    e2 = _rk4_step_doubling_evaluation(h=2.0, theta=0, kappa=0, varpi=2, frame="lab")
    assert e2.coarse_execution.cp_status is CPStatus.FAIL
    assert e2.fine_execution.cp_status is CPStatus.PASS
    assert min_choi_eigenvalue_numeric(e1.coarse_matrix) > -1e-12
    assert min_choi_eigenvalue_numeric(e1.fine_matrix) < 0


def test_rotating_phase_does_not_change_candidate_certificate() -> None:
    lab = build_candidate_execution(
        CandidateSpec("equal-substeps", "rk4", 0, 0, 0, 2, "lab"),
        total_x=1.0,
    )
    rotating = build_candidate_execution(
        CandidateSpec("equal-substeps", "rk4", 0, 0, 2, 2, "rotating"),
        total_x=1.0,
    )
    assert lab.cp_status is rotating.cp_status
    assert np.linalg.norm(lab.numerical_superoperator - rotating.numerical_superoperator) > 0


def test_richardson_execution_binds_non_cptp_candidate() -> None:
    execution = build_candidate_execution(
        CandidateSpec("richardson-4", "rk4", 0, 0, 0, 2, "lab"),
        total_x=1.0,
    )
    assert execution.cp_status is CPStatus.FAIL
    assert min_choi_eigenvalue_numeric(execution.numerical_superoperator) < 0
    assert verify_execution_provenance(execution)


def test_step_doubling_candidates_have_distinct_provenance_hashes() -> None:
    evaluation = _rk4_step_doubling_evaluation(
        h=1.0, theta=0, kappa=0, varpi=2, frame="lab"
    )
    hashes = {
        evaluation.coarse_execution.provenance_hash,
        evaluation.fine_execution.provenance_hash,
        evaluation.extrapolated_execution.provenance_hash,
    }
    assert len(hashes) == 3


def test_exact_fallback_is_certified_structurally() -> None:
    result = run_adaptive_channel_benchmark(
        final_time=1.0,
        initial_h=1.0,
        tolerance=10.0,
        theta=0,
        kappa=0,
        varpi=2,
        policy="candidate_guard",
        fallback="exact",
    )
    assert result.records[0].action == "fallback-exact"
    assert result.records[0].candidate_kind == "fallback-exact"
    assert result.records[0].provenance_verified
