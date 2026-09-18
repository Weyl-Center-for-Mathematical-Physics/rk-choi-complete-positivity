from __future__ import annotations

"""Candidate-aware certified controller logic and reproducible benchmarks.

The central invariant is literal rather than rhetorical: the immutable
``CandidateExecution`` whose CPTP status is certified is the same execution
record whose numerical superoperator is accumulated.  Local-error observations
are separate from the certified map identity.  Structural candidates are
constructed only from typed recipes that derive their own maps.
"""

from dataclasses import asdict, dataclass, field
import math
import operator
from typing import Literal, Mapping

import numpy as np
import sympy as sp

from .candidates import (
    CandidateBundle,
    CandidateExecution,
    CandidateSpec,
    PhaseCovariantCandidate,
    build_candidate_execution,
    build_exact_gksl_execution,
    build_strang_execution,
    locate_candidate_components,
    verify_execution_provenance,
)
from .certification import CPStatus, as_rational, choose_certified_interior
from .liouvillian import (
    exact_gksl_superoperator,
    exact_to_numpy,
    exp_superoperator_numeric,
    min_choi_eigenvalue_numeric,
    normalized_choi_trace_distance,
)

Q = sp.Symbol("q", real=True)
Y = sp.Symbol("y", real=True)



def _require_finite_float(
    name: str, value: object, *, positive: bool = False, nonnegative: bool = False
) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number, not bool")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and not result > 0.0:
        raise ValueError(f"{name} must be > 0")
    if nonnegative and result < 0.0:
        raise ValueError(f"{name} must be >= 0")
    return result


def _require_exact_positive_integer(name: str, value: object) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact positive integer, not bool")
    try:
        result = int(operator.index(value))
    except TypeError as exc:
        raise TypeError(f"{name} must be an exact positive integer") from exc
    if result < 1:
        raise ValueError(f"{name} must be >= 1")
    return result


@dataclass(frozen=True)
class GuardDecision:
    proposed_x: float
    accepted_x: float | None
    cp_status: CPStatus
    policy_status: str
    candidate_kind: str
    components: tuple[tuple[float, float], ...]
    diagnostics: Mapping[str, object]
    proposed_execution: CandidateExecution | None = None
    accepted_execution: CandidateExecution | None = None
    telemetry: Mapping[str, int] = field(default_factory=dict)


@dataclass
class ControllerStats:
    accepted_steps: int = 0
    direct_accepts: int = 0
    projected_accepts: int = 0
    rotating_accepts: int = 0
    strang_accepts: int = 0
    exact_fallback_accepts: int = 0
    rejected_error: int = 0
    rejected_cp: int = 0
    guard_evaluations: int = 0
    locator_calls: int = 0
    successful_projections: int = 0
    empty_component_outcomes: int = 0
    uncertain_outcomes: int = 0
    root_isolation_ops: int = 0
    certification_ops: int = 0
    binary_recertifications: int = 0
    precision_escalations: int = 0
    fallback_events: int = 0
    pi_memory_resets: int = 0
    rhs_stage_evaluations: int = 0
    exponential_actions: int = 0
    # Backward-compatible aliases retained for older validation scripts.
    candidate_certifications: int = 0
    component_searches: int = 0
    fallback_steps: int = 0


@dataclass(frozen=True)
class AdaptiveStepRecord:
    t: float
    proposed_h: float
    used_h: float
    error_estimate: float
    cp_status: str
    candidate_kind: str
    certified_kind: str
    min_choi_eigenvalue: float
    action: str
    provenance_verified: bool = False
    execution_id: str | None = None
    certified_execution_id: str | None = None
    candidate_family: str | None = None


@dataclass(frozen=True)
class AdaptiveBenchmarkResult:
    policy: str
    final_time: float
    tolerance: float
    global_normalized_choi_error: float
    min_step_choi_eigenvalue: float
    records: tuple[AdaptiveStepRecord, ...]
    stats: dict[str, int]


@dataclass(frozen=True)
class StepDoublingEvaluation:
    error_estimate: float
    coarse_execution: CandidateExecution
    fine_execution: CandidateExecution
    extrapolated_execution: CandidateExecution | None

    @property
    def bundle(self) -> CandidateBundle:
        if self.coarse_execution.symbolic_candidate is None or self.fine_execution.symbolic_candidate is None:
            raise RuntimeError("step-doubling bundle requires symbolic phase-covariant candidates")
        ext = None if self.extrapolated_execution is None else self.extrapolated_execution.symbolic_candidate
        return CandidateBundle(
            self.coarse_execution.symbolic_candidate,
            self.fine_execution.symbolic_candidate,
            ext,
        )

    @property
    def coarse_matrix(self) -> np.ndarray:
        return self.coarse_execution.numerical_superoperator

    @property
    def fine_matrix(self) -> np.ndarray:
        return self.fine_execution.numerical_superoperator

    @property
    def extrapolated_matrix(self) -> np.ndarray | None:
        return None if self.extrapolated_execution is None else self.extrapolated_execution.numerical_superoperator


@dataclass
class PIHistory:
    """Error history owned by a candidate family/representation.

    A CP-only rejection does not update this state.  Entering a new family
    initializes that family's history to the requested tolerance; the reset is
    not subsequently overwritten by a zero-error structural fallback.
    """

    tolerance: float
    active_family: str | None = None
    previous_by_family: dict[str, float] = field(default_factory=dict)

    def previous(self, family: str) -> float:
        return max(self.previous_by_family.get(family, self.tolerance), 1e-300)

    def observe(self, family: str, error: float) -> None:
        self.active_family = family
        self.previous_by_family[family] = max(float(error), 1e-300)

    def transition(self, family: str) -> bool:
        changed = self.active_family is not None and self.active_family != family
        self.active_family = family
        if changed or family not in self.previous_by_family:
            self.previous_by_family[family] = self.tolerance
        return changed

    def reset_family(self, family: str) -> None:
        self.active_family = family
        self.previous_by_family[family] = self.tolerance


# ---------------------------------------------------------------------------
# Exact RK4 detached interval and multiplicative reduction theorem
# ---------------------------------------------------------------------------


def rk4_detached_interval_exact(
    varpi: int | float | str | sp.Rational, digits: int = 40
):
    """Certified outer brackets for the two positive detached RK4 roots."""

    v = as_rational(varpi)
    q = 1 + 4 * v**2
    qplus = (sp.Integer(123) + 11 * sp.sqrt(33)) / 16
    if sp.N(q - qplus, 80) <= 0:
        return None
    C = sp.Poly(
        Y**3 - 16 * Y**2 - 32 * (q - 6) * Y + 384 * (q - 4),
        Y,
        domain=sp.QQ,
    )
    intervals = []
    for (lo, hi), mult in C.intervals(eps=sp.Rational(1, 10**digits)):
        if hi > 0:
            intervals.extend([(sp.Rational(lo) / q, sp.Rational(hi) / q)] * int(mult))
    intervals = sorted(intervals, key=lambda row: row[0])
    if len(intervals) != 2:
        raise RuntimeError(f"expected two positive roots, found {intervals}")
    return intervals[0][0], intervals[0][1], intervals[1][0], intervals[1][1]


def reduction_skip_intervals(
    x_minus: float, x_plus: float, rho: float, max_j: int
) -> list[tuple[float, float]]:
    if not (0 < x_minus < x_plus and 0 < rho < 1):
        raise ValueError("invalid interval or reduction factor")
    if rho >= x_minus / x_plus:
        return []
    return [
        (x_plus * rho ** (-j), x_minus * rho ** (-(j + 1)))
        for j in range(max_j)
    ]


def every_geometric_sequence_hits_window(
    x_minus: float, x_plus: float, rho: float
) -> bool:
    return 0 < rho < 1 and rho >= x_minus / x_plus


def halving_threshold_certificate(digits: int = 40) -> dict[str, str]:
    poly = sp.Poly(
        72 * Q**3 - 2071 * Q**2 + 12552 * Q - 21744,
        Q,
        domain=sp.QQ,
    )
    qplus = (sp.Integer(123) + 11 * sp.sqrt(33)) / 16
    candidates = []
    for (lo, hi), mult in poly.intervals(eps=sp.Rational(1, 10**digits)):
        if sp.N(lo - qplus, 80) > 0:
            candidates.append((sp.Rational(lo), sp.Rational(hi), int(mult)))
    if len(candidates) != 1:
        raise RuntimeError(f"could not isolate unique detached-regime root: {candidates}")
    lo, hi, mult = candidates[0]
    vlo = sp.sqrt((lo - 1) / 4)
    vhi = sp.sqrt((hi - 1) / 4)
    return {
        "polynomial": str(poly.as_expr()),
        "q_lower": str(lo),
        "q_upper": str(hi),
        "multiplicity": str(mult),
        "varpi_lower": str(sp.N(vlo, 30)),
        "varpi_upper": str(sp.N(vhi, 30)),
        "varpi_mid": str(sp.N((vlo + vhi) / 2, 30)),
    }


# ---------------------------------------------------------------------------
# Candidate certification and projection
# ---------------------------------------------------------------------------


def _spec_from_kind(
    kind: str,
    *,
    method: str,
    theta,
    kappa,
    varpi,
    substeps: int,
    frame: str = "lab",
) -> CandidateSpec:
    return CandidateSpec(
        kind=kind,
        method=method,
        theta=theta,
        kappa=kappa,
        varpi=varpi,
        substeps=substeps,
        frame=frame,
    )


def _execution_for_kind(
    kind: str,
    *,
    method: str,
    total_x: float,
    theta,
    kappa,
    varpi,
    substeps: int,
    frame: str = "lab",
) -> CandidateExecution:
    return build_candidate_execution(
        _spec_from_kind(
            kind,
            method=method,
            theta=theta,
            kappa=kappa,
            varpi=varpi,
            substeps=substeps,
            frame=frame,
        ),
        total_x=float(total_x),
    )


def _execution_matches_request(
    execution: CandidateExecution,
    *,
    kind: str,
    method: str,
    total_x: float,
    theta,
    kappa,
    varpi,
    substeps: int,
    frame: str = "lab",
) -> bool:
    expected = _spec_from_kind(
        kind,
        method=method,
        theta=theta,
        kappa=kappa,
        varpi=varpi,
        substeps=substeps,
        frame=frame,
    )
    return (
        execution.spec == expected
        and execution.executed_float == float(total_x)
        and execution.exact_input == as_rational(float(total_x))
        and verify_execution_provenance(execution)
    )


def _empty_guard_telemetry() -> dict[str, int]:
    return {
        "locator_calls": 0,
        "successful_projections": 0,
        "empty_component_outcomes": 0,
        "uncertain_outcomes": 0,
        "root_isolation_ops": 0,
        "certification_ops": 0,
        "binary_recertifications": 0,
        "precision_escalations": 0,
    }


def _locator_telemetry(location) -> dict[str, int]:
    telemetry = _empty_guard_telemetry()
    telemetry["locator_calls"] = 1
    for key in ("root_isolation_ops", "certification_ops", "precision_escalations"):
        telemetry[key] = int(location.diagnostics.get(key, 0))
    return telemetry


def _recertify_float_inside_component(
    *,
    rational_target: sp.Rational,
    component_sample: sp.Rational,
    kind: str,
    method: str,
    theta,
    kappa,
    varpi,
    substeps: int,
) -> tuple[CandidateExecution | None, int]:
    """Convert a rational interior target to binary64 and recertify it.

    The returned integer counts literal binary candidate constructions and exact
    CP certifications.  Failure after the finite search is explicit; the
    theoretical existence of a rational interior point does not imply that this
    search must discover a certifiable binary64 representation.
    """

    value = float(rational_target)
    target = float(component_sample)
    attempts = 0
    for _ in range(256):
        attempts += 1
        execution = _execution_for_kind(
            kind,
            method=method,
            total_x=value,
            theta=theta,
            kappa=kappa,
            varpi=varpi,
            substeps=substeps,
        )
        if execution.cp_status is CPStatus.PASS and verify_execution_provenance(execution):
            return execution, attempts
        next_value = math.nextafter(value, target)
        if next_value == value:
            break
        value = next_value
    return None, attempts


def certified_candidate_guard_step(
    *,
    method: str,
    candidate_kind: Literal["direct", "equal-substeps", "richardson-4"],
    proposed_x: float,
    theta: float,
    kappa: float,
    varpi: float,
    substeps: int = 2,
    safety_fraction: sp.Rational = sp.Rational(1, 20),
    proposed_execution: CandidateExecution | None = None,
) -> GuardDecision:
    """Return the literal executed-binary candidate certified for advancement."""

    if proposed_execution is None:
        proposed_execution = _execution_for_kind(
            candidate_kind,
            method=method,
            total_x=proposed_x,
            theta=theta,
            kappa=kappa,
            varpi=varpi,
            substeps=substeps,
        )
    elif not _execution_matches_request(
        proposed_execution,
        kind=candidate_kind,
        method=method,
        total_x=proposed_x,
        theta=theta,
        kappa=kappa,
        varpi=varpi,
        substeps=substeps,
    ):
        raise ValueError("proposed execution does not match the guarded candidate request")

    cert = proposed_execution.cp_certificate
    if cert.status is CPStatus.PASS:
        return GuardDecision(
            proposed_x,
            proposed_execution.executed_float,
            CPStatus.PASS,
            "accepted-proposal",
            proposed_execution.candidate_kind,
            tuple(),
            cert.diagnostics,
            proposed_execution,
            proposed_execution,
            _empty_guard_telemetry(),
        )
    if cert.status is CPStatus.UNCERTAIN:
        return GuardDecision(
            proposed_x,
            None,
            CPStatus.UNCERTAIN,
            "CERTIFICATION_UNCERTAIN",
            proposed_execution.candidate_kind,
            tuple(),
            cert.diagnostics,
            proposed_execution,
            None,
            _empty_guard_telemetry(),
        )

    location = locate_candidate_components(
        candidate_kind,
        method=method,
        theta=theta,
        kappa=kappa,
        varpi=varpi,
        upper=proposed_x,
        substeps=substeps,
        root_digits=32,
    )
    telemetry = _locator_telemetry(location)
    if location.status is CPStatus.UNCERTAIN:
        telemetry["uncertain_outcomes"] = 1
        return GuardDecision(
            proposed_x,
            None,
            CPStatus.UNCERTAIN,
            "CERTIFICATION_UNCERTAIN",
            proposed_execution.candidate_kind,
            tuple(),
            location.diagnostics,
            proposed_execution,
            None,
            telemetry,
        )
    components = tuple(
        (float(component.lower_value), float(component.upper_value))
        for component in location.components
    )
    exact_proposal = as_rational(proposed_x)
    available = [
        component
        for component in location.components
        if component.lower_value < exact_proposal
        and component.upper_value > component.lower_value
    ]
    if not available:
        telemetry["empty_component_outcomes"] = 1
        boundary = as_rational(theta) in (0, 1) and as_rational(kappa) == 0
        return GuardDecision(
            proposed_x,
            None,
            CPStatus.FAIL,
            "BOUNDARY_NO_RECOVERY" if boundary else "BUFFERED_NO_COMPONENT_BELOW_PROPOSAL",
            proposed_execution.candidate_kind,
            components,
            location.diagnostics,
            proposed_execution,
            None,
            telemetry,
        )

    selected = max(available, key=lambda component: component.upper_value)
    interior = choose_certified_interior(
        selected, ceiling=proposed_x, safety_fraction=safety_fraction
    )
    if interior is None:
        telemetry["uncertain_outcomes"] = 1
        return GuardDecision(
            proposed_x,
            None,
            CPStatus.UNCERTAIN,
            "CERTIFICATION_UNCERTAIN",
            proposed_execution.candidate_kind,
            components,
            location.diagnostics,
            proposed_execution,
            None,
            telemetry,
        )

    accepted_execution, recertification_attempts = _recertify_float_inside_component(
        rational_target=interior,
        component_sample=selected.sample,
        kind=candidate_kind,
        method=method,
        theta=theta,
        kappa=kappa,
        varpi=varpi,
        substeps=substeps,
    )
    telemetry["binary_recertifications"] += recertification_attempts
    telemetry["certification_ops"] += recertification_attempts
    if accepted_execution is None:
        telemetry["uncertain_outcomes"] = 1
        return GuardDecision(
            proposed_x,
            None,
            CPStatus.UNCERTAIN,
            "EXECUTED_BINARY_RECERTIFICATION_FAILED",
            proposed_execution.candidate_kind,
            components,
            {**dict(location.diagnostics), "interior_rational": str(interior)},
            proposed_execution,
            None,
            telemetry,
        )

    telemetry["successful_projections"] = 1
    return GuardDecision(
        proposed_x,
        accepted_execution.executed_float,
        CPStatus.PASS,
        "projected-to-certified-interior",
        accepted_execution.candidate_kind,
        components,
        {
            **dict(location.diagnostics),
            "interior_rational": str(interior),
            "executed_binary_rational": str(accepted_execution.exact_input),
            "execution_id": accepted_execution.execution_id,
        },
        proposed_execution,
        accepted_execution,
        telemetry,
    )


def certified_guard_step(**kwargs) -> GuardDecision:
    """Backward-compatible direct-step guard."""

    return certified_candidate_guard_step(candidate_kind="direct", **kwargs)


# ---------------------------------------------------------------------------
# Step-doubling candidate evaluation
# ---------------------------------------------------------------------------


def _rk4_step_doubling_evaluation(
    *,
    h: float,
    theta: float,
    kappa: float,
    varpi: float,
    frame: Literal["lab", "rotating"],
    fine_execution: CandidateExecution | None = None,
) -> StepDoublingEvaluation:
    if frame == "lab":
        coarse_spec = CandidateSpec("direct", "rk4", theta, kappa, varpi, 1, "lab")
        fine_spec = CandidateSpec("equal-substeps", "rk4", theta, kappa, varpi, 2, "lab")
        ext_spec = CandidateSpec("richardson-4", "rk4", theta, kappa, varpi, 2, "lab")
        extrapolated = build_candidate_execution(ext_spec, total_x=h)
    else:
        coarse_spec = CandidateSpec("equal-substeps", "rk4", theta, kappa, varpi, 1, "rotating")
        fine_spec = CandidateSpec("equal-substeps", "rk4", theta, kappa, varpi, 2, "rotating")
        extrapolated = None

    coarse = build_candidate_execution(coarse_spec, total_x=h)
    if fine_execution is None:
        fine = build_candidate_execution(fine_spec, total_x=h)
    else:
        if not (
            fine_execution.spec == fine_spec
            and fine_execution.executed_float == float(h)
            and fine_execution.exact_input == as_rational(float(h))
            and verify_execution_provenance(fine_execution)
        ):
            raise ValueError("provided fine execution does not match the estimator candidate")
        fine = fine_execution

    error = normalized_choi_trace_distance(
        coarse.numerical_superoperator, fine.numerical_superoperator
    ) / 15.0
    return StepDoublingEvaluation(error, coarse, fine, extrapolated)


def _strang_step_doubling_evaluation(
    *,
    h: float,
    theta: float,
    kappa: float,
    varpi: float,
    omega_x: float = 0.0,
) -> tuple[float, CandidateExecution, CandidateExecution]:
    coarse = build_strang_execution(
        total_x=h,
        gamma=1,
        theta=theta,
        kappa=kappa,
        omega_z=varpi,
        omega_x=omega_x,
        substeps=1,
    )
    fine = build_strang_execution(
        total_x=h,
        gamma=1,
        theta=theta,
        kappa=kappa,
        omega_z=varpi,
        omega_x=omega_x,
        substeps=2,
    )
    error = normalized_choi_trace_distance(
        coarse.numerical_superoperator, fine.numerical_superoperator
    ) / 3.0
    return error, coarse, fine


def _execution_passes_common_gate(
    execution: CandidateExecution,
    tolerance: float,
    error_estimate: float,
) -> bool:
    return (
        execution.cp_status is CPStatus.PASS
        and math.isfinite(error_estimate)
        and error_estimate <= tolerance
        and verify_execution_provenance(execution)
    )


def _candidate_family(execution: CandidateExecution) -> str:
    if execution.structural_recipe is not None:
        return execution.structural_recipe.recipe_kind.value
    return f"{execution.spec.frame}:{execution.spec.kind}:{execution.spec.method}:{execution.spec.substeps}"


# ---------------------------------------------------------------------------
# Adaptive benchmark
# ---------------------------------------------------------------------------


def run_adaptive_channel_benchmark(
    *,
    final_time: float,
    initial_h: float,
    tolerance: float,
    theta: float,
    kappa: float,
    varpi: float,
    policy: Literal[
        "error_only", "candidate_guard", "certified_guard", "rotating_frame", "strang"
    ],
    fallback: Literal["exact", "rotating_frame"] = "exact",
    safety: float = 0.9,
    min_h: float = 1e-12,
    max_steps: int = 10000,
) -> AdaptiveBenchmarkResult:
    """Generator-level benchmark with literal executed-candidate continuity.

    Public inputs are validated before any expensive construction so malformed
    controller domains fail clearly rather than through resource exhaustion.
    """

    final_time = _require_finite_float("final_time", final_time, positive=True)
    initial_h = _require_finite_float("initial_h", initial_h, positive=True)
    tolerance = _require_finite_float("tolerance", tolerance, positive=True)
    min_h = _require_finite_float("min_h", min_h, positive=True)
    safety = _require_finite_float("safety", safety, positive=True)
    if safety > 1.0:
        raise ValueError("safety must lie in (0, 1]")
    max_steps = _require_exact_positive_integer("max_steps", max_steps)
    theta = _require_finite_float("theta", theta, nonnegative=True)
    if theta > 1.0:
        raise ValueError("theta must lie in [0, 1]")
    kappa = _require_finite_float("kappa", kappa, nonnegative=True)
    varpi = _require_finite_float("varpi", varpi)
    if min_h > final_time:
        raise ValueError("min_h must not exceed final_time")
    if min_h > initial_h:
        raise ValueError("min_h must not exceed initial_h")
    if not isinstance(policy, str):
        raise TypeError("policy must be a string")
    if not isinstance(fallback, str):
        raise TypeError("fallback must be a string")
    supported_policies = {"error_only", "candidate_guard", "certified_guard", "rotating_frame", "strang"}
    supported_fallbacks = {"exact", "rotating_frame"}
    if policy not in supported_policies:
        raise ValueError(f"unsupported policy {policy!r}")
    if fallback not in supported_fallbacks:
        raise ValueError(f"unsupported fallback {fallback!r}")

    if policy == "certified_guard":
        policy = "candidate_guard"
    L = exact_to_numpy(
        exact_gksl_superoperator(
            gamma=1,
            theta=theta,
            kappa=kappa,
            omega_z=varpi,
            omega_x=0,
        )
    )
    S_total = np.eye(4, dtype=np.complex128)
    exact_total = exp_superoperator_numeric(L, as_rational(float(final_time)))
    t = 0.0
    h = min(float(initial_h), final_time)
    stats = ControllerStats()
    history = PIHistory(float(tolerance))
    records: list[AdaptiveStepRecord] = []
    min_eig = math.inf
    p = 2 if policy == "strang" else 4
    alpha = 0.7 / (p + 1)
    beta = 0.4 / (p + 1)

    while t < final_time - 1e-15:
        if len(records) >= max_steps or h < min_h:
            raise RuntimeError("adaptive benchmark exceeded resource limits")
        h_prop = min(h, final_time - t)
        used_h = h_prop
        cp_label = "not-checked"
        certified_kind = "none"

        if policy == "strang":
            err, _coarse_execution, proposed_execution = _strang_step_doubling_evaluation(
                h=h_prop, theta=theta, kappa=kappa, varpi=varpi
            )
            stats.exponential_actions += 9
            stats.certification_ops += 2
        else:
            frame = "rotating" if policy == "rotating_frame" else "lab"
            evaluation = _rk4_step_doubling_evaluation(
                h=h_prop,
                theta=theta,
                kappa=kappa,
                varpi=varpi,
                frame=frame,
            )
            proposed_execution = evaluation.fine_execution
            err = evaluation.error_estimate
            stats.rhs_stage_evaluations += 12
            stats.certification_ops += 2 if frame == "rotating" else 3

        proposed_family = _candidate_family(proposed_execution)

        # An LTE rejection uses genuine error information for this family.
        if err > tolerance:
            stats.rejected_error += 1
            previous = history.previous(proposed_family)
            factor = max(
                0.1,
                min(
                    0.8,
                    safety * (tolerance / err) ** alpha * (previous / err) ** beta,
                ),
            )
            records.append(
                AdaptiveStepRecord(
                    t,
                    h_prop,
                    0.0,
                    err,
                    cp_label,
                    proposed_execution.candidate_kind,
                    certified_kind,
                    float("nan"),
                    "error-reject",
                    verify_execution_provenance(proposed_execution),
                    proposed_execution.execution_id,
                    None,
                    proposed_family,
                )
            )
            history.observe(proposed_family, err)
            h = h_prop * factor
            continue

        accepted_execution = proposed_execution
        action = "error-only-accept" if policy == "error_only" else f"{policy}-accept"

        if policy == "candidate_guard":
            stats.guard_evaluations += 1
            decision = certified_candidate_guard_step(
                method="rk4",
                candidate_kind="equal-substeps",
                substeps=2,
                proposed_x=h_prop,
                theta=theta,
                kappa=kappa,
                varpi=varpi,
                proposed_execution=proposed_execution,
            )
            cp_label = decision.policy_status
            certified_kind = decision.candidate_kind
            for key in (
                "locator_calls",
                "successful_projections",
                "empty_component_outcomes",
                "uncertain_outcomes",
                "root_isolation_ops",
                "certification_ops",
                "binary_recertifications",
                "precision_escalations",
            ):
                setattr(stats, key, getattr(stats, key) + int(decision.telemetry.get(key, 0)))

            if decision.accepted_execution is proposed_execution:
                accepted_execution = proposed_execution
                action = "guard-accept"
            elif decision.accepted_execution is not None:
                stats.rejected_cp += 1
                used_h = decision.accepted_execution.executed_float
                projected = _rk4_step_doubling_evaluation(
                    h=used_h,
                    theta=theta,
                    kappa=kappa,
                    varpi=varpi,
                    frame="lab",
                    fine_execution=decision.accepted_execution,
                )
                if projected.fine_execution is not decision.accepted_execution:
                    raise RuntimeError("projected estimator did not retain literal execution identity")
                accepted_execution = decision.accepted_execution
                err = projected.error_estimate
                stats.rhs_stage_evaluations += 12
                # The retained fine execution is not rebuilt; coarse and
                # Richardson estimator candidates are newly constructed.
                stats.certification_ops += 2
                action = "guard-projection"
                if not _execution_passes_common_gate(accepted_execution, tolerance, err):
                    stats.rejected_error += 1
                    records.append(
                        AdaptiveStepRecord(
                            t,
                            h_prop,
                            used_h,
                            err,
                            cp_label,
                            accepted_execution.candidate_kind,
                            certified_kind,
                            float("nan"),
                            "projected-error-reject",
                            verify_execution_provenance(accepted_execution),
                            accepted_execution.execution_id,
                            accepted_execution.execution_id,
                            _candidate_family(accepted_execution),
                        )
                    )
                    family = _candidate_family(accepted_execution)
                    previous = history.previous(family)
                    factor = max(
                        0.1,
                        min(0.8, safety * (tolerance / err) ** alpha * (previous / err) ** beta),
                    )
                    history.observe(family, err)
                    h = max(min_h, used_h * factor)
                    continue
            else:
                stats.rejected_cp += 1
                fallback_execution: CandidateExecution | None = None
                fallback_err = 0.0
                if fallback == "rotating_frame":
                    rotating_eval = _rk4_step_doubling_evaluation(
                        h=h_prop,
                        theta=theta,
                        kappa=kappa,
                        varpi=varpi,
                        frame="rotating",
                    )
                    stats.rhs_stage_evaluations += 12
                    stats.certification_ops += 2
                    if _execution_passes_common_gate(
                        rotating_eval.fine_execution, tolerance, rotating_eval.error_estimate
                    ):
                        fallback_execution = rotating_eval.fine_execution
                        fallback_err = rotating_eval.error_estimate
                        action = "fallback-rotating"
                if fallback_execution is None:
                    fallback_execution = build_exact_gksl_execution(
                        total_x=h_prop,
                        gamma=1,
                        theta=theta,
                        kappa=kappa,
                        omega_z=varpi,
                        omega_x=0,
                    )
                    fallback_err = 0.0
                    stats.exponential_actions += 1
                    stats.certification_ops += 1
                    action = "fallback-exact"
                accepted_execution = fallback_execution
                used_h = accepted_execution.executed_float
                err = fallback_err
                stats.fallback_events += 1
                cp_label = decision.policy_status
                certified_kind = accepted_execution.candidate_kind

        elif policy in {"rotating_frame", "strang"}:
            cp_label = accepted_execution.cp_status.value
            certified_kind = accepted_execution.candidate_kind
            if not _execution_passes_common_gate(accepted_execution, tolerance, err):
                accepted_execution = build_exact_gksl_execution(
                    total_x=h_prop,
                    gamma=1,
                    theta=theta,
                    kappa=kappa,
                    omega_z=varpi,
                    omega_x=0,
                )
                stats.exponential_actions += 1
                stats.certification_ops += 1
                stats.fallback_events += 1
                action = "fallback-exact"
                err = 0.0
                certified_kind = accepted_execution.candidate_kind

        # Error-only intentionally demonstrates physicality failure.  Every
        # other policy must pass one common CP+LTE gate on the accumulated map.
        if policy != "error_only" and not _execution_passes_common_gate(
            accepted_execution, tolerance, err
        ):
            raise RuntimeError("attempt to accumulate a candidate that failed the common CP+LTE gate")
        if not verify_execution_provenance(accepted_execution):
            raise RuntimeError("candidate provenance mismatch before accumulation")

        accepted_family = _candidate_family(accepted_execution)
        family_changed = history.transition(accepted_family)
        if family_changed:
            stats.pi_memory_resets += 1

        candidate_matrix = accepted_execution.numerical_superoperator
        step_min = min_choi_eigenvalue_numeric(candidate_matrix)
        min_eig = min(min_eig, step_min)
        S_total = candidate_matrix @ S_total
        t += used_h
        stats.accepted_steps += 1
        if action == "guard-accept":
            stats.direct_accepts += 1
        elif action == "guard-projection":
            stats.projected_accepts += 1
        elif action == "fallback-rotating" or policy == "rotating_frame":
            stats.rotating_accepts += 1
        elif policy == "strang" and action != "fallback-exact":
            stats.strang_accepts += 1
        elif action == "fallback-exact":
            stats.exact_fallback_accepts += 1

        records.append(
            AdaptiveStepRecord(
                t - used_h,
                h_prop,
                used_h,
                err,
                cp_label,
                accepted_execution.candidate_kind,
                certified_kind,
                step_min,
                action,
                True,
                accepted_execution.execution_id,
                accepted_execution.execution_id,
                accepted_family,
            )
        )

        # On a family transition, hold the accepted step and start the new
        # family's PI memory at tolerance.  A zero-error exact fallback cannot
        # overwrite the reset and force the subsequent 0.4 shrink bug.
        if family_changed:
            history.reset_family(accepted_family)
            factor = 1.0
        elif err == 0:
            # Structural fallbacks supply no LTE observation.  Preserve the
            # reset memory but allow a conservative growth proposal rather
            # than feeding zero into the PI formula.
            history.reset_family(accepted_family)
            factor = 2.0
        else:
            previous = history.previous(accepted_family)
            factor = max(
                0.2,
                min(
                    3.0,
                    safety * (tolerance / err) ** alpha * (previous / err) ** beta,
                ),
            )
            history.observe(accepted_family, err)
        h = used_h * factor

    if stats.locator_calls != (
        stats.successful_projections
        + stats.empty_component_outcomes
        + stats.uncertain_outcomes
    ):
        raise RuntimeError("locator telemetry outcomes do not partition locator calls")
    stats.candidate_certifications = stats.certification_ops
    stats.component_searches = stats.locator_calls
    stats.fallback_steps = stats.fallback_events

    global_error = normalized_choi_trace_distance(S_total, exact_total)
    return AdaptiveBenchmarkResult(
        policy,
        final_time,
        tolerance,
        global_error,
        min_eig,
        tuple(records),
        asdict(stats),
    )
