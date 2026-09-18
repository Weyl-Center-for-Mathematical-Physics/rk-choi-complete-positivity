from __future__ import annotations

import math

from rk_choi_margin.v26 import halving_sequence, rk4_cptp_guard_step, rk4_detached_interval
from rk_choi_margin.topology import normalized_choi_defect_and_error


def test_detached_interval_at_varpi_two() -> None:
    interval = rk4_detached_interval(2.0)
    assert interval is not None
    lo, hi = interval
    assert abs(lo - 0.744782427896) < 2e-12
    assert abs(hi - 1.270334626974) < 2e-12


def test_halving_can_skip_physical_window() -> None:
    lo, hi = rk4_detached_interval(2.0)  # type: ignore[misc]
    seq = halving_sequence(1.4, 8)
    assert seq[0] > hi
    assert seq[1] < lo
    assert all(x < lo for x in seq[1:])


def test_guard_projects_to_largest_feasible_step() -> None:
    decision = rk4_cptp_guard_step(1.4, 2.0)
    assert decision.status == "projected-to-certified-interior"
    assert decision.accepted_x is not None
    assert 0.744782427895 < decision.accepted_x < 1.270334626974


def test_guard_detects_accuracy_physicality_incompatibility() -> None:
    decision = rk4_cptp_guard_step(0.7, 2.0)
    assert decision.accepted_x is None
    assert decision.status == "no-cptp-step-below-error-ceiling"


def test_bell_defect_is_order_one_fraction_for_representative_case() -> None:
    defect, error = normalized_choi_defect_and_error(0.01, 2.0)
    ratio = defect / error
    assert 0.90 < float(ratio) < 0.92


def test_high_frequency_guard_scaling_is_consistent() -> None:
    for varpi in (25.0, 100.0):
        lo, hi = rk4_detached_interval(varpi)  # type: ignore[misc]
        assert abs(lo * varpi**2 - 3.0) < 0.01
        assert abs(hi * varpi - 2 * math.sqrt(2)) < 0.03
