from __future__ import annotations

import math
import random

import mpmath as mp
import pytest
import sympy as sp

from rk_choi_margin.certification import CPStatus, certify_cptp_step, locate_cptp_components
from rk_choi_margin.controllers import (
    certified_guard_step,
    every_geometric_sequence_hits_window,
    halving_threshold_certificate,
    reduction_skip_intervals,
    rk4_detached_interval_exact,
)
from rk_choi_margin.topology import rk4_rotating_cubic


@pytest.mark.parametrize("x", [0.004, 0.001, 1e-4, 1e-5, 0.7])
def test_certified_guard_rejects_small_non_cptp_steps(x: float) -> None:
    cert = certify_cptp_step("rk4", x=x, theta=0, kappa=0, varpi=2)
    assert cert.status is CPStatus.FAIL
    decision = certified_guard_step(method="rk4", proposed_x=x, theta=0, kappa=0, varpi=2)
    assert decision.accepted_x is None
    assert decision.policy_status == "BOUNDARY_NO_RECOVERY"


def test_certified_guard_accepts_interior_and_projects_from_above() -> None:
    assert certify_cptp_step("rk4", x=1, theta=0, kappa=0, varpi=2).status is CPStatus.PASS
    decision = certified_guard_step(method="rk4", proposed_x=1.4, theta=0, kappa=0, varpi=2)
    assert decision.cp_status is CPStatus.PASS
    assert decision.accepted_x is not None
    assert 0.744782427895 < decision.accepted_x < 1.270334626974
    assert certify_cptp_step("rk4", x=str(decision.accepted_x), theta=0, kappa=0, varpi=2).status is CPStatus.PASS


def test_points_on_both_sides_of_certified_boundaries() -> None:
    brackets = rk4_detached_interval_exact(2, digits=35)
    assert brackets is not None
    lm, lp, um, up = brackets
    eps = sp.Rational(1, 10**20)
    assert certify_cptp_step("rk4", x=lm - eps, theta=0, kappa=0, varpi=2).status is CPStatus.FAIL
    assert certify_cptp_step("rk4", x=lp + eps, theta=0, kappa=0, varpi=2).status is CPStatus.PASS
    assert certify_cptp_step("rk4", x=um - eps, theta=0, kappa=0, varpi=2).status is CPStatus.PASS
    assert certify_cptp_step("rk4", x=up + eps, theta=0, kappa=0, varpi=2).status is CPStatus.FAIL


def test_component_isolation_survives_near_saddle_node() -> None:
    vp = math.sqrt((float((sp.Integer(123) + 11 * sp.sqrt(33)) / 16) - 1) / 4)
    for delta in (1e-12, 1e-14, 1e-15):
        v = vp + delta
        loc = locate_cptp_components("rk4", theta=0, kappa=0, varpi=str(v), upper=3, root_digits=24)
        # Above the saddle node the exact locator finds a detached component or,
        # at decimal rounding indistinguishable from the algebraic threshold,
        # reports an isolated equality point.  It never infers emptiness from
        # discarded complex roots.
        assert loc.components or loc.isolated_points


def test_companion_matrix_failure_is_reproduced_but_not_used_for_certification() -> None:
    import numpy as np

    vp = math.sqrt((float((sp.Integer(123) + 11 * sp.sqrt(33)) / 16) - 1) / 4)
    failures = 0
    for ulps in range(1, 20):
        v = vp
        for _ in range(ulps):
            v = math.nextafter(v, math.inf)
        q = 1.0 + 4.0 * v**2
        roots = np.roots([1.0, -16.0, -32.0 * (q - 6.0), 384.0 * (q - 4.0)])
        eligible = [r for r in roots if abs(float(r.imag)) < 1e-9 and float(r.real) > 1e-12]
        if len(eligible) != 2:
            failures += 1
    assert failures >= 1
    loc = locate_cptp_components("rk4", theta=0, kappa=0, varpi="1.630712023895", upper=3, root_digits=24)
    assert loc.components


def test_reduction_factor_theorem_and_skip_set() -> None:
    brackets = rk4_detached_interval_exact(2, digits=32)
    assert brackets is not None
    xm = float((brackets[0] + brackets[1]) / 2)
    xp = float((brackets[2] + brackets[3]) / 2)
    threshold = xm / xp
    assert every_geometric_sequence_hits_window(xm, xp, threshold)
    assert every_geometric_sequence_hits_window(xm, xp, min(0.999, threshold + 0.01))
    assert not every_geometric_sequence_hits_window(xm, xp, threshold - 0.01)
    skips = reduction_skip_intervals(xm, xp, 0.5, 4)
    assert skips
    # Choose one point from each skip interval and verify the geometric sequence
    # jumps from above x+ to below x-.
    for lo, hi in skips:
        x0 = (lo + hi) / 2
        seq = [x0 * 0.5**j for j in range(12)]
        assert not any(xm <= x <= xp for x in seq)


def test_halving_threshold_certificate() -> None:
    cert = halving_threshold_certificate(45)
    value = float(cert["varpi_mid"])
    assert abs(value - 2.2482546045149805) < 1e-14
    q = sp.Symbol("q")
    polynomial = 72 * q**3 - 2071 * q**2 + 12552 * q - 21744
    qval = 1 + 4 * sp.Float(cert["varpi_mid"], 60) ** 2
    assert abs(complex(sp.N(polynomial.subs(q, qval), 40))) < 1e-20


def test_boundary_and_buffered_component_policies_differ() -> None:
    boundary = certified_guard_step(method="rk4", proposed_x=0.5, theta=0, kappa=0, varpi=2)
    assert boundary.policy_status == "BOUNDARY_NO_RECOVERY"
    buffered = certified_guard_step(method="rk4", proposed_x=0.5, theta=0.001, kappa=0.001, varpi=2)
    assert buffered.accepted_x is not None
    assert buffered.policy_status == "projected-to-certified-interior"
    assert 0 < buffered.accepted_x < 0.263


def test_randomized_exact_certificates_have_zero_false_positive_passes() -> None:
    rng = random.Random(20260810)
    mp.mp.dps = 90

    def r4(z: mp.mpc) -> mp.mpc:
        return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24

    false_pass = 0
    for _ in range(1000):
        # Decimal grids make the exact certificate and high-precision numerical
        # check refer to precisely the same rational input parameters.
        x = rng.randint(1, 300000) / 100000
        theta = rng.randint(0, 1000) / 1000
        kappa = rng.randint(0, 200) / 1000
        varpi = rng.randint(0, 5000) / 1000
        cert = certify_cptp_step("rk4", x=str(x), theta=str(theta), kappa=str(kappa), varpi=str(varpi))
        xx, tt, kk, vv = map(lambda v: mp.mpf(str(v)), (x, theta, kappa, varpi))
        a = mp.re(r4(-xx))
        c = r4(-xx * (mp.mpf("0.5") + kk - 1j * vv))
        A = 1 - tt * (1 - a)
        D = tt + (1 - tt) * a
        F = a + tt * (1 - tt) * (1 - a) ** 2 - abs(c) ** 2
        high_precision_pass = (1 - a >= 0) and (A >= 0) and (D >= 0) and (F >= 0)
        if cert.status is CPStatus.PASS and not high_precision_pass:
            false_pass += 1
    assert false_pass == 0


def test_logarithmic_small_step_stress() -> None:
    for exponent in range(1, 13):
        x = sp.Rational(1, 10**exponent)
        assert certify_cptp_step("rk4", x=x, theta=0, kappa=0, varpi=2).status is CPStatus.FAIL
