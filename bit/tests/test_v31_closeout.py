from __future__ import annotations

import math

import numpy as np
import pytest

from rk_choi_margin.candidates import (
    ExactGKSLExponentialRecipe,
    ExactSubflowStrangRecipe,
    build_exact_gksl_execution,
    build_strang_execution,
    verify_execution_provenance,
)
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import run_adaptive_channel_benchmark


NEGATIVE_TIMES = (-1.0, -1e-12, -1e-13, -1e-16, np.nextafter(0.0, -math.inf))
NONNEGATIVE_CONTROLS = (0.0, np.nextafter(0.0, math.inf), 1e-13, 0.1, 1.0)


@pytest.mark.parametrize("h", NEGATIVE_TIMES)
@pytest.mark.parametrize("builder", [build_exact_gksl_execution, build_strang_execution])
def test_structural_builders_reject_negative_time(h, builder) -> None:
    with pytest.raises(ValueError, match="total_x must be >= 0"):
        builder(total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0)


@pytest.mark.parametrize("h", NONNEGATIVE_CONTROLS)
@pytest.mark.parametrize("builder", [build_exact_gksl_execution, build_strang_execution])
def test_structural_builders_accept_nonnegative_time(h, builder) -> None:
    execution = builder(total_x=h, gamma=1, theta=0, kappa=0, omega_z=2, omega_x=0)
    assert execution.cp_status is CPStatus.PASS
    assert verify_execution_provenance(execution)
    domain = execution.cp_certificate.recipe_payload["theorem_domain"]
    assert domain["time_nonnegative"] is True
    assert domain["exact_time"] == str(execution.exact_input)


@pytest.mark.parametrize("value", [1.5, 2.0, 2.9, True, False])
def test_strang_substeps_reject_lossy_or_boolean_inputs(value) -> None:
    with pytest.raises((TypeError, ValueError), match="substeps"):
        ExactSubflowStrangRecipe.from_values(substeps=value)


@pytest.mark.parametrize("value", [1.5, 30.0, True, False])
def test_structural_precision_rejects_noninteger_inputs(value) -> None:
    with pytest.raises((TypeError, ValueError), match="dps"):
        ExactGKSLExponentialRecipe.from_values(dps=value)


@pytest.mark.parametrize("value", [0, -5, 1, 29])
def test_structural_precision_enforces_minimum(value) -> None:
    with pytest.raises(ValueError, match="dps must be >= 30"):
        ExactGKSLExponentialRecipe.from_values(dps=value)


def test_valid_exact_integer_inputs_remain_supported() -> None:
    r = ExactSubflowStrangRecipe.from_values(substeps=np.int64(2), dps=np.int64(50))
    assert r.substeps == 2
    assert r.dps == 50


BASE_CONTROLLER = dict(
    final_time=1.0,
    initial_h=0.1,
    tolerance=1e-2,
    theta=0.0,
    kappa=0.0,
    varpi=2.0,
    policy="candidate_guard",
    fallback="exact",
    min_h=1e-12,
    max_steps=100,
)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("final_time", -1.0, "final_time must be > 0"),
        ("final_time", 0.0, "final_time must be > 0"),
        ("initial_h", -0.1, "initial_h must be > 0"),
        ("initial_h", 0.0, "initial_h must be > 0"),
        ("tolerance", -1e-2, "tolerance must be > 0"),
        ("tolerance", 0.0, "tolerance must be > 0"),
        ("min_h", -1e-12, "min_h must be > 0"),
        ("min_h", 0.0, "min_h must be > 0"),
        ("safety", 0.0, "safety must be > 0"),
        ("safety", 1.1, "safety must lie in"),
        ("theta", -1e-3, "theta must be >= 0"),
        ("theta", 1.001, "theta must lie in"),
        ("kappa", -1e-3, "kappa must be >= 0"),
    ],
)
def test_controller_rejects_invalid_numeric_domain(field, value, message) -> None:
    kwargs = dict(BASE_CONTROLLER)
    kwargs[field] = value
    with pytest.raises(ValueError, match=message):
        run_adaptive_channel_benchmark(**kwargs)


@pytest.mark.parametrize("field", ["final_time", "initial_h", "tolerance", "min_h", "theta", "kappa", "varpi"])
def test_controller_rejects_nonfinite_inputs(field) -> None:
    kwargs = dict(BASE_CONTROLLER)
    kwargs[field] = math.nan
    with pytest.raises(ValueError, match=f"{field} must be finite"):
        run_adaptive_channel_benchmark(**kwargs)


@pytest.mark.parametrize("value", [0, -1, 1.5, 3.0, True])
def test_controller_max_steps_requires_exact_positive_integer(value) -> None:
    kwargs = dict(BASE_CONTROLLER)
    kwargs["max_steps"] = value
    with pytest.raises((TypeError, ValueError), match="max_steps"):
        run_adaptive_channel_benchmark(**kwargs)


def test_controller_rejects_minimum_step_larger_than_domain() -> None:
    kwargs = dict(BASE_CONTROLLER)
    kwargs["min_h"] = 0.2
    with pytest.raises(ValueError, match="min_h must not exceed initial_h"):
        run_adaptive_channel_benchmark(**kwargs)


def test_controller_valid_inputs_still_execute() -> None:
    result = run_adaptive_channel_benchmark(**BASE_CONTROLLER)
    assert result.records
    assert all(r.provenance_verified for r in result.records if r.used_h > 0)
