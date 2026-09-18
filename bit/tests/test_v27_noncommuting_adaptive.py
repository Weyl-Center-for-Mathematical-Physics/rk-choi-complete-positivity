from __future__ import annotations

import pytest

from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import run_adaptive_channel_benchmark
from rk_choi_margin.liouvillian import (
    certify_rk4_gksl_step_exact,
    locate_rk4_gksl_cptp_components_exact,
)


@pytest.mark.parametrize("omega_x", [0, 0.05, 0.1, 0.2])
def test_noncommuting_transverse_field_preserves_disconnected_components(omega_x: float) -> None:
    result = locate_rk4_gksl_cptp_components_exact(
        upper=1.5,
        theta="0.001",
        kappa="0.001",
        omega_z=2,
        omega_x=str(omega_x),
        root_digits=18,
    )
    components = result["components"]
    assert len(components) == 2
    first = tuple(float(v) for v in components[0][:2])
    second = tuple(float(v) for v in components[1][:2])
    assert 0 <= first[0] < first[1] < second[0] < second[1] < 1.5


def test_noncommuting_strict_sign_samples_are_exact() -> None:
    expected = {0.1: CPStatus.PASS, 0.5: CPStatus.FAIL, 0.9: CPStatus.PASS, 1.4: CPStatus.FAIL}
    for h, status in expected.items():
        cert = certify_rk4_gksl_step_exact(
            h=str(h), theta="0.001", kappa="0.001", omega_z=2, omega_x="0.1"
        )
        assert cert.status is status


def test_adaptive_error_only_can_accept_non_cptp_steps() -> None:
    result = run_adaptive_channel_benchmark(
        final_time=2.0,
        initial_h=1.4,
        tolerance=0.002,
        theta=0,
        kappa=0,
        varpi=2,
        policy="error_only",
    )
    assert result.stats["accepted_steps"] >= 1
    assert result.min_step_choi_eigenvalue < -1e-5


def test_certified_guard_or_rotating_frame_avoid_false_physical_steps() -> None:
    guarded = run_adaptive_channel_benchmark(
        final_time=2.0,
        initial_h=1.4,
        tolerance=0.002,
        theta=0,
        kappa=0,
        varpi=2,
        policy="certified_guard",
        fallback="rotating_frame",
    )
    rotating = run_adaptive_channel_benchmark(
        final_time=2.0,
        initial_h=1.4,
        tolerance=0.002,
        theta=0,
        kappa=0,
        varpi=2,
        policy="rotating_frame",
    )
    assert guarded.min_step_choi_eigenvalue > -1e-12
    assert rotating.min_step_choi_eigenvalue > -1e-12
    assert guarded.stats["fallback_steps"] >= 1


def test_buffered_guard_recomputes_error_at_actual_candidate() -> None:
    guarded = run_adaptive_channel_benchmark(
        final_time=1.0,
        initial_h=1.4,
        tolerance=0.002,
        theta=0.001,
        kappa=0.001,
        varpi=2,
        policy="certified_guard",
        fallback="rotating_frame",
    )
    projected = [row for row in guarded.records if row.action == "guard-projection"]
    assert projected
    assert all(row.used_h < row.proposed_h for row in projected)
    assert all(row.error_estimate <= 0.002 for row in projected)
    assert guarded.min_step_choi_eigenvalue > -1e-12
