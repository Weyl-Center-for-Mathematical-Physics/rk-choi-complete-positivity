from __future__ import annotations

"""Compatibility wrappers for results introduced in v2.6.

Scientific decisions now delegate to the v2.7 certified implementation.  No
companion-matrix root calculation or fixed absolute CP tolerance is used here.
"""

from dataclasses import dataclass
from typing import Optional

import sympy as sp

from .controllers import certified_guard_step, rk4_detached_interval_exact
from .topology import rk4_rotating_cubic

Y = sp.Symbol("y", real=True)


@dataclass(frozen=True)
class GuardDecision:
    proposed_x: float
    accepted_x: Optional[float]
    status: str
    lower: Optional[float] = None
    upper: Optional[float] = None


def rk4_rotating_positive_roots(varpi: float) -> list[float]:
    """Certified positive roots in x of the exact one-way RK4 Choi cubic."""

    v = sp.Rational(str(varpi))
    q = 1 + 4 * v**2
    poly = sp.Poly(rk4_rotating_cubic(q, Y), Y, domain=sp.QQ)
    rows: list[float] = []
    for (lo, hi), mult in poly.intervals(eps=sp.Rational(1, 10**24)):
        if hi > 0:
            rows.extend([float(((sp.Rational(lo) + sp.Rational(hi)) / 2) / q)] * int(mult))
    return sorted(rows)


def rk4_detached_interval(varpi: float) -> Optional[tuple[float, float]]:
    brackets = rk4_detached_interval_exact(varpi, digits=30)
    if brackets is None:
        return None
    lm, lp, um, up = brackets
    return float((lm + lp) / 2), float((um + up) / 2)


def rk4_cptp_guard_step(error_ceiling_x: float, varpi: float) -> GuardDecision:
    """Compatibility view of the certified one-way RK4 guard.

    Operational steps are selected strictly inside a certified component; an
    equality root is never returned as a supposedly safe floating step.
    """

    result = certified_guard_step(
        method="rk4",
        proposed_x=float(error_ceiling_x),
        theta=0.0,
        kappa=0.0,
        varpi=float(varpi),
    )
    interval = rk4_detached_interval(varpi)
    lower = interval[0] if interval else None
    upper = interval[1] if interval else None
    status_map = {
        "accepted-proposal": "accepted-proposal",
        "projected-to-certified-interior": "projected-to-certified-interior",
        "BOUNDARY_NO_RECOVERY": "no-cptp-step-below-error-ceiling",
        "certification-uncertain": "certification-uncertain",
    }
    return GuardDecision(
        float(error_ceiling_x),
        result.accepted_x,
        status_map.get(result.policy_status, result.policy_status),
        lower,
        upper,
    )


def halving_sequence(x0: float, n: int) -> list[float]:
    return [float(x0) / (2**k) for k in range(n)]
