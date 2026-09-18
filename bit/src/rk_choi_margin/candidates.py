from __future__ import annotations

"""Candidate-map algebra for adaptive Runge--Kutta workflows.

Complete positivity belongs to the map actually advanced.  A step-doubling
workflow may construct a coarse full step, a refined composition, a rotating-
frame candidate, or a negative-weight extrapolate over the same total interval;
these candidates can have different CPTP status.  The execution record in this
module binds the exact IEEE-754 input, symbolic candidate or typed structural
recipe, numerical matrix, and certificate into one immutable object.  Local-
error observations remain separate from certified map identity.
"""

from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from fractions import Fraction
import hashlib
import json
import math
import operator
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np
import sympy as sp

from .certification import (
    CPStatus,
    CertifiedComponent,
    CertifiedScalar,
    ComponentLocation,
    RootIsolation,
    as_rational,
    isolate_constraint_roots_exact,
    exact_root_for_isolation,
)
from .methods import RKMethod, all_methods

S = sp.Symbol("s")
X = sp.Symbol("x", nonnegative=True, real=True)


MIN_STRUCTURAL_DPS = 30


def _require_exact_integer(name: str, value: Any, *, minimum: int) -> int:
    """Return an exact integer input without lossy coercion.

    Booleans and floating-point values such as ``2.0`` or ``2.9`` are rejected
    rather than silently truncated.  NumPy/SymPy integer scalar types are
    accepted through the integer-index protocol.
    """
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer, not bool")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an exact integer") from exc
    result = int(result)
    if result < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return result


def _require_finite_float(
    name: str, value: Any, *, positive: bool = False, nonnegative: bool = False
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


def _validated_rational(name: str, value: Any) -> sp.Rational:
    """Convert a finite real input to the exact rational convention in use."""
    _require_finite_float(name, value)
    try:
        return as_rational(value)
    except Exception as exc:
        raise TypeError(f"{name} must be convertible to an exact finite rational") from exc


def _freeze_mapping(values: Mapping[str, Any] | None = None) -> Mapping[str, Any]:
    """Return an immutable mapping with recursively frozen values."""

    def freeze(value: Any) -> Any:
        if isinstance(value, Mapping):
            return MappingProxyType({str(k): freeze(v) for k, v in value.items()})
        if isinstance(value, list):
            return tuple(freeze(v) for v in value)
        if isinstance(value, tuple):
            return tuple(freeze(v) for v in value)
        if isinstance(value, set):
            return frozenset(freeze(v) for v in value)
        return value

    return MappingProxyType({str(k): freeze(v) for k, v in (values or {}).items()})


class ProofKind(str, Enum):
    ALGEBRAIC = "ALGEBRAIC"
    STRUCTURAL = "STRUCTURAL"


class StructuralTheorem(str, Enum):
    EXACT_GKSL_EXPONENTIAL = "exact-gksl-exponential"
    EXACT_SUBFLOW_STRANG = "exact-subflow-strang"


@dataclass(frozen=True)
class ExactGKSLExponentialRecipe:
    gamma: sp.Rational
    theta: sp.Rational
    kappa: sp.Rational
    omega_z: sp.Rational
    omega_x: sp.Rational
    dps: int = 80
    builder_version: str = "v3.4"
    recipe_kind: StructuralTheorem = StructuralTheorem.EXACT_GKSL_EXPONENTIAL

    @classmethod
    def from_values(cls, *, gamma=1, theta=0, kappa=0, omega_z=0, omega_x=0, dps=80):
        names = ("gamma", "theta", "kappa", "omega_z", "omega_x")
        values = tuple(_validated_rational(name, value) for name, value in zip(names, (gamma, theta, kappa, omega_z, omega_x)))
        # Validate the GKSL parameter domain before a structural proof can be constructed.
        g, th, kap, _, _ = values
        if g < 0:
            raise ValueError("gamma must be >= 0")
        if th < 0 or th > 1:
            raise ValueError("theta must lie in [0, 1]")
        if kap < 0:
            raise ValueError("kappa must be >= 0")
        return cls(*values, _require_exact_integer("dps", dps, minimum=MIN_STRUCTURAL_DPS))


@dataclass(frozen=True)
class ExactSubflowStrangRecipe:
    gamma: sp.Rational
    theta: sp.Rational
    kappa: sp.Rational
    omega_z: sp.Rational
    omega_x: sp.Rational
    substeps: int = 1
    dps: int = 80
    builder_version: str = "v3.4"
    recipe_kind: StructuralTheorem = StructuralTheorem.EXACT_SUBFLOW_STRANG

    @classmethod
    def from_values(
        cls, *, gamma=1, theta=0, kappa=0, omega_z=0, omega_x=0,
        substeps=1, dps=80
    ):
        substeps_i = _require_exact_integer("substeps", substeps, minimum=1)
        dps_i = _require_exact_integer("dps", dps, minimum=MIN_STRUCTURAL_DPS)
        names = ("gamma", "theta", "kappa", "omega_z", "omega_x")
        values = tuple(_validated_rational(name, value) for name, value in zip(names, (gamma, theta, kappa, omega_z, omega_x)))
        g, th, kap, _, _ = values
        if g < 0:
            raise ValueError("gamma must be >= 0")
        if th < 0 or th > 1:
            raise ValueError("theta must lie in [0, 1]")
        if kap < 0:
            raise ValueError("kappa must be >= 0")
        return cls(*values, substeps_i, dps_i)


StructuralRecipe = ExactGKSLExponentialRecipe | ExactSubflowStrangRecipe


@dataclass(frozen=True)
class PhaseCovariantCandidate:
    kind: str
    method: str
    total_x: sp.Rational
    theta: sp.Rational
    kappa: sp.Rational
    varpi: sp.Rational
    a: sp.Expr
    c: sp.Expr
    stage_factors: tuple[sp.Expr, ...]
    constituents: tuple[str, ...]
    diagnostics: Mapping[str, Any]
    phase_angle: sp.Expr = sp.Integer(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "diagnostics", _freeze_mapping(self.diagnostics))


@dataclass(frozen=True)
class CandidateCertificate:
    candidate: PhaseCovariantCandidate
    status: CPStatus
    constraints: tuple[CertifiedScalar, ...]
    stage_domain: tuple[CertifiedScalar, ...]
    diagnostics: Mapping[str, Any]
    proof_kind: ProofKind = ProofKind.ALGEBRAIC

    def __post_init__(self) -> None:
        object.__setattr__(self, "diagnostics", _freeze_mapping(self.diagnostics))


@dataclass(frozen=True)
class StructuralCertificate:
    """Typed structural proof plus numerical-realization defense checks."""

    status: CPStatus
    proof_kind: ProofKind
    theorem_id: StructuralTheorem
    recipe_payload: Mapping[str, Any]
    realization_checks: Mapping[str, Any]
    diagnostics: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "recipe_payload", _freeze_mapping(self.recipe_payload))
        object.__setattr__(self, "realization_checks", _freeze_mapping(self.realization_checks))
        object.__setattr__(self, "diagnostics", _freeze_mapping(self.diagnostics))


@dataclass(frozen=True)
class CandidateBundle:
    coarse: PhaseCovariantCandidate
    fine: PhaseCovariantCandidate
    extrapolated: PhaseCovariantCandidate | None = None


@dataclass(frozen=True)
class CandidateSpec:
    kind: str
    method: str = "rk4"
    theta: float | int | str | sp.Rational = 0
    kappa: float | int | str | sp.Rational = 0
    varpi: float | int | str | sp.Rational = 0
    substeps: int = 1
    frame: str = "lab"


@dataclass(frozen=True, eq=False)
class CandidateExecution:
    """Immutable execution record for the exact map proposed for advancement.

    ``execution_id`` identifies the mathematical/numerical candidate.  The
    estimator, certifier, controller, and accumulator pass this same record;
    local-error observations are deliberately stored outside it.
    """

    executed_float: float
    exact_input: sp.Rational
    spec: CandidateSpec
    symbolic_candidate: PhaseCovariantCandidate | None
    structural_recipe: StructuralRecipe | None
    numerical_superoperator: np.ndarray
    cp_certificate: CandidateCertificate | StructuralCertificate
    execution_id: str
    provenance_hash: str
    diagnostics: Mapping[str, Any]

    def __post_init__(self) -> None:
        source = np.ascontiguousarray(self.numerical_superoperator, dtype=np.complex128)
        # Immutable bytes prevent callers from re-enabling ndarray writes.
        matrix = np.frombuffer(source.tobytes(), dtype=np.complex128).reshape(source.shape)
        object.__setattr__(self, "numerical_superoperator", matrix)
        object.__setattr__(self, "diagnostics", _freeze_mapping(self.diagnostics))

    @property
    def local_error_estimate(self) -> None:
        """Compatibility sentinel; error observations are separate."""
        return None

    @property
    def cp_status(self) -> CPStatus:
        return self.cp_certificate.status

    @property
    def proof_kind(self) -> ProofKind:
        return self.cp_certificate.proof_kind

    @property
    def candidate_kind(self) -> str:
        if self.symbolic_candidate is not None:
            return self.symbolic_candidate.kind
        return self.spec.kind


# ---------------------------------------------------------------------------
# Candidate construction
# ---------------------------------------------------------------------------


def _method(method: str | RKMethod) -> RKMethod:
    if isinstance(method, RKMethod):
        return method
    try:
        return all_methods()[method]
    except KeyError as exc:
        raise ValueError(f"unknown RK method {method!r}") from exc


def _sign_status(value: sp.Expr) -> CPStatus:
    v = sp.factor(sp.simplify(value))
    if v.is_real is False:
        return CPStatus.UNCERTAIN
    if v == 0 or v.is_zero or v.is_positive:
        return CPStatus.PASS
    if v.is_negative:
        return CPStatus.FAIL
    try:
        sign = sp.sign(v)
    except Exception:
        return CPStatus.UNCERTAIN
    if sign in (0, 1):
        return CPStatus.PASS
    if sign == -1:
        return CPStatus.FAIL
    return CPStatus.UNCERTAIN


def _direct_symbolic(
    method: str | RKMethod,
    total_x: sp.Expr,
    kappa: sp.Expr,
    varpi: sp.Expr,
):
    rk = _method(method)
    R = rk.stability_function(x=S)
    z = -total_x * (sp.Rational(1, 2) + kappa - sp.I * varpi)
    a = sp.factor(sp.cancel(R.subs(S, -total_x)))
    c = sp.factor(sp.cancel(R.subs(S, z)))
    det = rk.stage_resolvent_determinant(x=S)
    stage = (sp.factor(det.subs(S, -total_x)), sp.factor(det.subs(S, z)))
    return a, c, stage


def direct_candidate(
    method: str | RKMethod,
    *,
    total_x,
    theta,
    kappa,
    varpi,
) -> PhaseCovariantCandidate:
    rk = _method(method)
    xr, tr, kr, vr = map(as_rational, (total_x, theta, kappa, varpi))
    a, c, stage = _direct_symbolic(rk, xr, kr, vr)
    return PhaseCovariantCandidate(
        "direct",
        rk.key,
        xr,
        tr,
        kr,
        vr,
        a,
        c,
        stage,
        (f"{rk.key}:full",),
        {"input_semantics": "exact rational; floats are exact IEEE-754 values"},
    )


def equal_substep_candidate(
    method: str | RKMethod,
    *,
    total_x,
    substeps: int,
    theta,
    kappa,
    varpi,
) -> PhaseCovariantCandidate:
    substeps = _require_exact_integer("substeps", substeps, minimum=1)
    rk = _method(method)
    xr, tr, kr, vr = map(as_rational, (total_x, theta, kappa, varpi))
    a0, c0, stage = _direct_symbolic(rk, xr / substeps, kr, vr)
    return PhaseCovariantCandidate(
        f"equal-{substeps}-substeps",
        rk.key,
        xr,
        tr,
        kr,
        vr,
        sp.factor(a0**substeps),
        sp.factor(c0**substeps),
        stage,
        tuple(f"{rk.key}:1/{substeps}" for _ in range(substeps)),
        {
            "composition": f"{substeps} identical substeps",
            "input_semantics": "exact rational; floats are exact IEEE-754 values",
        },
    )


def rotating_candidate(
    method: str | RKMethod,
    *,
    total_x,
    substeps: int,
    theta,
    kappa,
    varpi,
) -> PhaseCovariantCandidate:
    """RK on the de-rotated dissipative generator, followed by exact phase.

    The certificate is determined by the de-rotated multipliers.  The exact
    phase is stored separately and applied when the numerical superoperator is
    constructed, so the candidate metadata and executed matrix remain identical.
    """

    base = equal_substep_candidate(
        method,
        total_x=total_x,
        substeps=substeps,
        theta=theta,
        kappa=kappa,
        varpi=0,
    )
    xr = as_rational(total_x)
    vr = as_rational(varpi)
    return PhaseCovariantCandidate(
        f"rotating-frame-{substeps}-substeps",
        base.method,
        base.total_x,
        base.theta,
        base.kappa,
        vr,
        base.a,
        base.c,
        base.stage_factors,
        base.constituents,
        {
            "representation": "de-rotated dissipative RK candidate plus exact lab-frame phase",
            "input_semantics": "exact rational; floats are exact IEEE-754 values",
        },
        phase_angle=sp.factor(vr * xr),
    )


def affine_candidate(
    candidates: Sequence[PhaseCovariantCandidate],
    weights: Sequence,
    *,
    kind: str = "affine",
) -> PhaseCovariantCandidate:
    if not candidates or len(candidates) != len(weights):
        raise ValueError("candidates and weights must be nonempty and have equal length")
    wr = tuple(as_rational(w) for w in weights)
    if sum(wr) != 1:
        raise ValueError("affine weights must sum exactly to one")
    first = candidates[0]
    for cand in candidates[1:]:
        if (
            cand.method,
            cand.total_x,
            cand.theta,
            cand.kappa,
            cand.varpi,
            cand.phase_angle,
        ) != (
            first.method,
            first.total_x,
            first.theta,
            first.kappa,
            first.varpi,
            first.phase_angle,
        ):
            raise ValueError(
                "candidate maps must share method, total step, generator ray, and frame phase"
            )
    return PhaseCovariantCandidate(
        kind,
        first.method,
        first.total_x,
        first.theta,
        first.kappa,
        first.varpi,
        sp.factor(sum(w * c.a for w, c in zip(wr, candidates))),
        sp.factor(sum(w * c.c for w, c in zip(wr, candidates))),
        tuple(s for c in candidates for s in c.stage_factors),
        tuple(c.kind for c in candidates),
        {
            "weights": ",".join(str(w) for w in wr),
            "convex": str(all(w >= 0 for w in wr)),
            "warning": (
                "negative affine weights do not preserve CPTP"
                if any(w < 0 for w in wr)
                else "convex combination"
            ),
        },
        phase_angle=first.phase_angle,
    )


def rk4_step_doubling_bundle(
    *,
    total_x,
    theta,
    kappa,
    varpi,
    include_richardson: bool = True,
) -> CandidateBundle:
    coarse = direct_candidate(
        "rk4", total_x=total_x, theta=theta, kappa=kappa, varpi=varpi
    )
    fine = equal_substep_candidate(
        "rk4",
        total_x=total_x,
        substeps=2,
        theta=theta,
        kappa=kappa,
        varpi=varpi,
    )
    extrapolated = (
        affine_candidate(
            (fine, coarse),
            (sp.Rational(16, 15), -sp.Rational(1, 15)),
            kind="richardson-4",
        )
        if include_richardson
        else None
    )
    return CandidateBundle(coarse, fine, extrapolated)


def candidate_constraints(candidate: PhaseCovariantCandidate) -> dict[str, sp.Expr]:
    a, c, th = candidate.a, candidate.c, candidate.theta
    one_minus_a = sp.factor(1 - a)
    A = sp.factor(1 - th * one_minus_a)
    D = sp.factor(th + (1 - th) * a)
    F = sp.factor(
        sp.cancel(a + th * (1 - th) * one_minus_a**2 - c * sp.conjugate(c))
    )
    return {"1-a": one_minus_a, "A_theta": A, "D_theta": D, "F": F}


def certify_candidate(candidate: PhaseCovariantCandidate) -> CandidateCertificate:
    constraints = tuple(
        CertifiedScalar(name, value, ">=0", _sign_status(value))
        for name, value in candidate_constraints(candidate).items()
    )
    stages = tuple(
        CertifiedScalar(
            f"stage_{i}",
            value,
            "!=0",
            CPStatus.FAIL if sp.simplify(value) == 0 else CPStatus.PASS,
        )
        for i, value in enumerate(candidate.stage_factors)
    )
    statuses = [item.status for item in (*constraints, *stages)]
    status = (
        CPStatus.FAIL
        if CPStatus.FAIL in statuses
        else CPStatus.UNCERTAIN
        if CPStatus.UNCERTAIN in statuses
        else CPStatus.PASS
    )
    return CandidateCertificate(
        candidate,
        status,
        constraints,
        stages,
        {
            **candidate.diagnostics,
            "certified_object": candidate.kind,
            "arithmetic": "exact rational/algebraic",
            "phase_invariance": "post-unitary phase does not alter CPTP status",
        },
    )


def phase_covariant_superoperator_numeric(candidate: PhaseCovariantCandidate) -> np.ndarray:
    a = complex(sp.N(candidate.a, 50))
    c = complex(sp.N(candidate.c, 50))
    phase = complex(sp.N(sp.exp(sp.I * candidate.phase_angle), 50))
    c *= phase
    th = float(candidate.theta)
    A = 1 - th * (1 - a)
    B = th * (1 - a)
    C = (1 - th) * (1 - a)
    D = th + (1 - th) * a
    matrix = np.array(
        [[A, 0, 0, C], [0, np.conjugate(c), 0, 0], [0, 0, c, 0], [B, 0, 0, D]],
        dtype=np.complex128,
    )
    matrix.setflags(write=False)
    return matrix


# ---------------------------------------------------------------------------
# Immutable executed-candidate records
# ---------------------------------------------------------------------------


def _candidate_from_spec(spec: CandidateSpec, total_x) -> PhaseCovariantCandidate:
    if spec.frame == "rotating":
        return rotating_candidate(
            spec.method,
            total_x=total_x,
            substeps=spec.substeps,
            theta=spec.theta,
            kappa=spec.kappa,
            varpi=spec.varpi,
        )
    if spec.kind == "direct":
        return direct_candidate(
            spec.method,
            total_x=total_x,
            theta=spec.theta,
            kappa=spec.kappa,
            varpi=spec.varpi,
        )
    if spec.kind in {"equal-substeps", "fine"}:
        return equal_substep_candidate(
            spec.method,
            total_x=total_x,
            substeps=spec.substeps,
            theta=spec.theta,
            kappa=spec.kappa,
            varpi=spec.varpi,
        )
    if spec.kind == "richardson-4":
        bundle = rk4_step_doubling_bundle(
            total_x=total_x,
            theta=spec.theta,
            kappa=spec.kappa,
            varpi=spec.varpi,
            include_richardson=True,
        )
        if bundle.extrapolated is None:
            raise RuntimeError("Richardson candidate was not constructed")
        return bundle.extrapolated
    raise ValueError(f"unsupported candidate specification {spec}")


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, sp.Basic):
        return str(value)
    if is_dataclass(value):
        return {field.name: _jsonable(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, np.ndarray):
        return {"real": value.real.tolist(), "imag": value.imag.tolist()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(v) for v in value]
    return value


def _certificate_payload(cert: CandidateCertificate | StructuralCertificate) -> dict[str, object]:
    if isinstance(cert, CandidateCertificate):
        return {
            "status": cert.status.value,
            "proof_kind": cert.proof_kind.value,
            "constraints": [(c.name, str(c.value), c.status.value) for c in cert.constraints],
            "stage": [(c.name, str(c.value), c.status.value) for c in cert.stage_domain],
            "diagnostics": _jsonable(cert.diagnostics),
        }
    return {
        "status": cert.status.value,
        "proof_kind": cert.proof_kind.value,
        "theorem_id": cert.theorem_id.value,
        "recipe": _jsonable(cert.recipe_payload),
        "checks": _jsonable(cert.realization_checks),
        "diagnostics": _jsonable(cert.diagnostics),
    }


def _core_payload(
    *, executed_float: float, exact_input: sp.Rational, spec: CandidateSpec,
    symbolic_candidate: PhaseCovariantCandidate | None,
    structural_recipe: StructuralRecipe | None, matrix: np.ndarray,
    certificate: CandidateCertificate | StructuralCertificate,
    diagnostics: Mapping[str, Any] | None,
) -> dict[str, object]:
    return {
        "executed_float_hex": float(executed_float).hex(),
        "exact_input": str(exact_input),
        "spec": _jsonable(spec),
        "symbolic": None if symbolic_candidate is None else _jsonable(symbolic_candidate),
        "structural_recipe": _jsonable(structural_recipe),
        "matrix": _jsonable(np.asarray(matrix)),
        "certificate": _certificate_payload(certificate),
        "diagnostics": _jsonable(diagnostics or {}),
    }


def candidate_execution_id(**kwargs) -> str:
    payload = _core_payload(**kwargs)
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def candidate_provenance_hash(
    *, executed_float: float, exact_input: sp.Rational, spec: CandidateSpec,
    symbolic_candidate: PhaseCovariantCandidate | None, matrix: np.ndarray,
    error_estimate: float | None = None,
    certificate: CandidateCertificate | StructuralCertificate,
    structural_recipe: StructuralRecipe | None = None,
    diagnostics: Mapping[str, Any] | None = None,
    execution_id: str | None = None,
) -> str:
    eid = execution_id or candidate_execution_id(
        executed_float=executed_float, exact_input=exact_input, spec=spec,
        symbolic_candidate=symbolic_candidate, structural_recipe=structural_recipe,
        matrix=matrix, certificate=certificate, diagnostics=diagnostics,
    )
    # error_estimate is accepted only for backward API compatibility; it is not
    # part of the certified map identity.
    return hashlib.sha256(json.dumps({"execution_id": eid}, sort_keys=True).encode("utf-8")).hexdigest()


def _finalize_execution(
    *, executed: float, exact: sp.Rational, spec: CandidateSpec,
    symbolic: PhaseCovariantCandidate | None, structural_recipe: StructuralRecipe | None,
    matrix: np.ndarray, certificate: CandidateCertificate | StructuralCertificate,
    diagnostics: Mapping[str, Any],
) -> CandidateExecution:
    source = np.ascontiguousarray(matrix, dtype=np.complex128)
    mat = np.frombuffer(source.tobytes(), dtype=np.complex128).reshape(source.shape)
    frozen_diag = _freeze_mapping(diagnostics)
    eid = candidate_execution_id(
        executed_float=executed, exact_input=exact, spec=spec,
        symbolic_candidate=symbolic, structural_recipe=structural_recipe,
        matrix=mat, certificate=certificate, diagnostics=frozen_diag,
    )
    digest = candidate_provenance_hash(
        executed_float=executed, exact_input=exact, spec=spec,
        symbolic_candidate=symbolic, structural_recipe=structural_recipe,
        matrix=mat, certificate=certificate, diagnostics=frozen_diag,
        execution_id=eid,
    )
    return CandidateExecution(executed, exact, spec, symbolic, structural_recipe,
                              mat, certificate, eid, digest, frozen_diag)


def build_candidate_execution(spec: CandidateSpec, *, total_x: float) -> CandidateExecution:
    executed = float(total_x)
    exact = as_rational(executed)
    symbolic = _candidate_from_spec(spec, exact)
    matrix = phase_covariant_superoperator_numeric(symbolic)
    cert = certify_candidate(symbolic)
    return _finalize_execution(
        executed=executed, exact=exact, spec=spec, symbolic=symbolic,
        structural_recipe=None, matrix=matrix, certificate=cert,
        diagnostics={
            "executed_float_hex": executed.hex(),
            "exact_binary_rational": str(exact),
            "candidate_identity": symbolic.kind,
            "proof_semantics": ProofKind.ALGEBRAIC.value,
        },
    )


def _structural_realization_checks(matrix: np.ndarray, dim: int = 2) -> dict[str, float | str]:
    from .liouvillian import choi_from_superoperator_numeric
    S_num = np.asarray(matrix, dtype=np.complex128)
    trace_row = np.zeros(dim * dim, dtype=np.complex128)
    for i in range(dim):
        trace_row[i + i * dim] = 1.0
    trace_residual = float(np.linalg.norm(trace_row @ S_num - trace_row))
    hermitian_residual = 0.0
    basis = (
        np.eye(dim, dtype=np.complex128),
        np.array([[0, 1], [1, 0]], dtype=np.complex128),
        np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        np.array([[1, 0], [0, -1]], dtype=np.complex128),
    )
    for B in basis:
        out = (S_num @ B.reshape(-1, order="F")).reshape((dim, dim), order="F")
        hermitian_residual = max(hermitian_residual, float(np.linalg.norm(out-out.conjugate().T)))
    choi = choi_from_superoperator_numeric(S_num, dim)
    eigs = np.linalg.eigvalsh(choi)
    scale = max(1.0, float(np.linalg.norm(S_num, ord=2)))
    tolerance = 2048.0*np.finfo(float).eps*scale
    passed = (trace_residual <= tolerance and hermitian_residual <= tolerance
              and float(eigs[0].real) >= -tolerance)
    return {
        "trace_preservation_residual": trace_residual,
        "hermiticity_preservation_residual": hermitian_residual,
        "minimum_realized_choi_eigenvalue": float(eigs[0].real),
        "realization_tolerance": tolerance,
        "realization_check": "PASS" if passed else "FAIL",
    }


def _structural_certificate(
    recipe: StructuralRecipe, matrix: np.ndarray, *, exact_time: sp.Rational
) -> StructuralCertificate:
    if exact_time < 0:
        raise ValueError("structural CPTP theorem requires exact_time >= 0")
    checks = _structural_realization_checks(matrix)
    status = CPStatus.PASS if checks["realization_check"] == "PASS" else CPStatus.UNCERTAIN
    payload = {
        "recipe": _jsonable(recipe),
        "theorem_domain": {
            "exact_time": str(exact_time),
            "time_nonnegative": True,
            "hypothesis": "forward-time GKSL semigroup/subflow construction, h >= 0",
        },
    }
    return StructuralCertificate(
        status=status, proof_kind=ProofKind.STRUCTURAL,
        theorem_id=recipe.recipe_kind, recipe_payload=payload,
        realization_checks=checks,
        diagnostics={
            "primary_proof": "typed recipe derives an exact forward-time CPTP construction",
            "numerical_realization": "independent trace/Hermiticity/Choi defense checks",
        },
    )


def build_exact_gksl_execution(
    *, total_x: float, gamma=1, theta=0, kappa=0,
    omega_z=0, omega_x=0, dps: int = 80,
) -> CandidateExecution:
    from .liouvillian import exp_exact_gksl_superoperator_numeric
    executed = _require_finite_float("total_x", total_x, nonnegative=True)
    exact = as_rational(executed)
    recipe = ExactGKSLExponentialRecipe.from_values(
        gamma=gamma, theta=theta, kappa=kappa,
        omega_z=omega_z, omega_x=omega_x, dps=dps)
    matrix = exp_exact_gksl_superoperator_numeric(
        h=exact, gamma=recipe.gamma, theta=recipe.theta,
        kappa=recipe.kappa, omega_z=recipe.omega_z,
        omega_x=recipe.omega_x, dps=recipe.dps)
    cert = _structural_certificate(recipe, matrix, exact_time=exact)
    spec = CandidateSpec(kind="fallback-exact", method="structural", frame="structural")
    return _finalize_execution(
        executed=executed, exact=exact, spec=spec, symbolic=None,
        structural_recipe=recipe, matrix=matrix, certificate=cert,
        diagnostics={"executed_float_hex": executed.hex(),
                     "exact_binary_rational": str(exact),
                     "builder": recipe.recipe_kind.value})


def build_strang_execution(
    *, total_x: float, gamma=1, theta=0, kappa=0,
    omega_z=0, omega_x=0, substeps: int = 1, dps: int = 80,
) -> CandidateExecution:
    from .liouvillian import strang_gksl_superoperator_numeric
    executed = _require_finite_float("total_x", total_x, nonnegative=True)
    exact = as_rational(executed)
    recipe = ExactSubflowStrangRecipe.from_values(
        gamma=gamma, theta=theta, kappa=kappa, omega_z=omega_z,
        omega_x=omega_x, substeps=substeps, dps=dps)
    matrix = strang_gksl_superoperator_numeric(
        h=exact, gamma=recipe.gamma, theta=recipe.theta,
        kappa=recipe.kappa, omega_z=recipe.omega_z,
        omega_x=recipe.omega_x, substeps=recipe.substeps, dps=recipe.dps)
    cert = _structural_certificate(recipe, matrix, exact_time=exact)
    spec = CandidateSpec(kind=f"strang-{recipe.substeps}-substeps",
                         method="structural", substeps=recipe.substeps,
                         frame="structural")
    return _finalize_execution(
        executed=executed, exact=exact, spec=spec, symbolic=None,
        structural_recipe=recipe, matrix=matrix, certificate=cert,
        diagnostics={"executed_float_hex": executed.hex(),
                     "exact_binary_rational": str(exact),
                     "builder": recipe.recipe_kind.value,
                     "substeps": recipe.substeps})


def verify_execution_provenance(execution: CandidateExecution) -> bool:
    """Rederive the map/certificate before accepting the provenance binding."""
    if execution.numerical_superoperator.flags.writeable:
        return False
    try:
        if execution.symbolic_candidate is not None:
            if execution.structural_recipe is not None:
                return False
            rebuilt = _candidate_from_spec(execution.spec, execution.exact_input)
            if _jsonable(rebuilt) != _jsonable(execution.symbolic_candidate):
                return False
            rebuilt_matrix = phase_covariant_superoperator_numeric(rebuilt)
            if not np.array_equal(rebuilt_matrix, execution.numerical_superoperator):
                return False
            rebuilt_cert = certify_candidate(rebuilt)
            if _certificate_payload(rebuilt_cert) != _certificate_payload(execution.cp_certificate):
                return False
        else:
            recipe = execution.structural_recipe
            if recipe is None or not isinstance(execution.cp_certificate, StructuralCertificate):
                return False
            if execution.spec.method != "structural" or execution.spec.frame != "structural":
                return False
            if isinstance(recipe, ExactGKSLExponentialRecipe):
                from .liouvillian import exp_exact_gksl_superoperator_numeric
                if execution.spec.kind != "fallback-exact": return False
                rebuilt_matrix = exp_exact_gksl_superoperator_numeric(
                    h=execution.exact_input, gamma=recipe.gamma, theta=recipe.theta,
                    kappa=recipe.kappa, omega_z=recipe.omega_z,
                    omega_x=recipe.omega_x, dps=recipe.dps)
            elif isinstance(recipe, ExactSubflowStrangRecipe):
                from .liouvillian import strang_gksl_superoperator_numeric
                if execution.spec.kind != f"strang-{recipe.substeps}-substeps": return False
                rebuilt_matrix = strang_gksl_superoperator_numeric(
                    h=execution.exact_input, gamma=recipe.gamma, theta=recipe.theta,
                    kappa=recipe.kappa, omega_z=recipe.omega_z,
                    omega_x=recipe.omega_x, substeps=recipe.substeps, dps=recipe.dps)
            else:
                return False
            rebuilt_matrix = np.asarray(rebuilt_matrix, dtype=np.complex128)
            if not np.array_equal(rebuilt_matrix, execution.numerical_superoperator):
                return False
            rebuilt_cert = _structural_certificate(recipe, rebuilt_matrix, exact_time=execution.exact_input)
            if _certificate_payload(rebuilt_cert) != _certificate_payload(execution.cp_certificate):
                return False
        expected_id = candidate_execution_id(
            executed_float=execution.executed_float, exact_input=execution.exact_input,
            spec=execution.spec, symbolic_candidate=execution.symbolic_candidate,
            structural_recipe=execution.structural_recipe,
            matrix=execution.numerical_superoperator,
            certificate=execution.cp_certificate, diagnostics=execution.diagnostics)
        expected_hash = candidate_provenance_hash(
            executed_float=execution.executed_float, exact_input=execution.exact_input,
            spec=execution.spec, symbolic_candidate=execution.symbolic_candidate,
            structural_recipe=execution.structural_recipe,
            matrix=execution.numerical_superoperator,
            certificate=execution.cp_certificate, diagnostics=execution.diagnostics,
            execution_id=expected_id)
        return expected_id == execution.execution_id and expected_hash == execution.provenance_hash
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Symbolic candidate families and exact component location
# ---------------------------------------------------------------------------


def candidate_symbolic_expressions(
    kind: str,
    *,
    method: str = "rk4",
    x: sp.Expr = X,
    kappa: sp.Expr = sp.Integer(0),
    varpi: sp.Expr = sp.Integer(0),
    substeps: int = 2,
):
    substeps = _require_exact_integer("substeps", substeps, minimum=1)
    if kind == "direct":
        return _direct_symbolic(method, x, kappa, varpi)
    if kind == "equal-substeps":
        a0, c0, stage = _direct_symbolic(method, x / substeps, kappa, varpi)
        return sp.factor(a0**substeps), sp.factor(c0**substeps), stage
    if kind == "richardson-4":
        if method != "rk4" or substeps != 2:
            raise ValueError("implemented Richardson candidate is RK4 step doubling")
        ac, cc, sc = _direct_symbolic(method, x, kappa, varpi)
        a0, c0, sf = _direct_symbolic(method, x / 2, kappa, varpi)
        af, cf = sp.factor(a0**2), sp.factor(c0**2)
        return sp.factor((16 * af - ac) / 15), sp.factor((16 * cf - cc) / 15), (*sc, *sf)
    raise ValueError(f"unknown candidate kind {kind!r}")


def candidate_constraint_expressions(
    kind: str,
    *,
    method: str,
    x: sp.Expr,
    theta: sp.Expr,
    kappa: sp.Expr,
    varpi: sp.Expr,
    substeps: int = 2,
):
    a, c, stages = candidate_symbolic_expressions(
        kind,
        method=method,
        x=x,
        kappa=kappa,
        varpi=varpi,
        substeps=substeps,
    )
    oma = sp.factor(1 - a)
    A = sp.factor(1 - theta * oma)
    D = sp.factor(theta + (1 - theta) * a)
    F = sp.factor(sp.cancel(a + theta * (1 - theta) * oma**2 - c * sp.conjugate(c)))
    out = {"1-a": oma, "A_theta": A, "D_theta": D, "F": F}
    for i, stage in enumerate(stages):
        out[f"stage_{i}"] = sp.factor(stage * sp.conjugate(stage))
    return out


def _primitive(expr: sp.Expr) -> sp.Poly:
    num, _ = sp.fraction(sp.cancel(expr))
    poly = sp.Poly(num, X, domain=sp.QQ)
    _, primitive = poly.primitive()
    return primitive if primitive.LC() > 0 else -primitive


def _certify_isolated_candidate_root(
    kind: str,
    root: RootIsolation,
    *,
    method: str,
    theta: sp.Rational,
    kappa: sp.Rational,
    varpi: sp.Rational,
    substeps: int,
) -> CPStatus:
    """Certify every candidate predicate at an isolated algebraic root."""

    if root.exact_root is None:
        return CPStatus.UNCERTAIN
    direct_substitution = isinstance(root.exact_root, (sp.Integer, sp.Rational))
    evaluation_point = root.exact_root if direct_substitution else root.midpoint
    source = set(root.source)
    expr = candidate_constraint_expressions(
        kind,
        method=method,
        x=evaluation_point,
        theta=theta,
        kappa=kappa,
        varpi=varpi,
        substeps=substeps,
    )
    statuses: list[CPStatus] = []
    for name, value in expr.items():
        if name.startswith("stage_"):
            if name in source:
                statuses.append(CPStatus.FAIL)
                continue
            simplified = sp.simplify(value)
            if simplified == 0:
                statuses.append(CPStatus.FAIL if direct_substitution else CPStatus.UNCERTAIN)
            else:
                sign = _sign_status(simplified)
                statuses.append(CPStatus.PASS if sign is CPStatus.PASS else CPStatus.UNCERTAIN)
        else:
            statuses.append(
                _sign_status(value) if direct_substitution or name not in source else CPStatus.PASS
            )
    if CPStatus.FAIL in statuses:
        return CPStatus.FAIL
    if CPStatus.UNCERTAIN in statuses:
        return CPStatus.UNCERTAIN
    return CPStatus.PASS


def locate_candidate_components(
    kind: str,
    *,
    method: str,
    theta,
    kappa,
    varpi,
    upper,
    substeps: int = 2,
    root_digits: int = 28,
) -> ComponentLocation:
    """Certify candidate-map CPTP components and exact equality points."""

    substeps = _require_exact_integer("substeps", substeps, minimum=1)
    root_digits = _require_exact_integer("root_digits", root_digits, minimum=8)
    tr, kr, vr, ur = map(as_rational, (theta, kappa, varpi, upper))
    if ur <= 0:
        return ComponentLocation(
            CPStatus.FAIL, tuple(), tuple(), tuple(), ur,
            {"reason": "upper must be positive"},
        )
    exprs = candidate_constraint_expressions(
        kind, method=method, x=X, theta=tr, kappa=kr, varpi=vr,
        substeps=substeps,
    )
    polys: dict[str, sp.Poly] = {}
    for name, expr in exprs.items():
        poly = _primitive(expr)
        if not poly.is_zero and poly.degree() > 0:
            polys[name] = poly
    telemetry = {
        "root_isolation_ops": 0,
        "certification_ops": 0,
        "precision_escalations": 0,
    }
    roots, separated = isolate_constraint_roots_exact(
        polys, ur, digits=root_digits, telemetry=telemetry
    )
    if not separated:
        return ComponentLocation(
            CPStatus.UNCERTAIN, tuple(), tuple(), tuple(), ur,
            {
                "candidate_kind": kind,
                "root_method": "exact Sturm isolation with adaptive refinement",
                "reason": "distinct boundary roots could not be certified as separated",
                **telemetry,
            },
        )

    interval_rows: list[dict[str, object]] = []
    cursor = sp.Rational(0)
    left_index: int | None = None
    for index, root in enumerate(roots):
        if cursor < root.lower:
            interval_rows.append(
                {"lo": cursor, "hi": root.lower, "left": left_index, "right": index}
            )
        cursor = max(cursor, root.upper)
        left_index = index
    if cursor < ur:
        interval_rows.append(
            {"lo": cursor, "hi": ur, "left": left_index, "right": None}
        )

    interval_uncertain = False
    for row in interval_rows:
        lo, hi = sp.Rational(row["lo"]), sp.Rational(row["hi"])
        sample = (lo + hi) / 2
        a, c, stages = candidate_symbolic_expressions(
            kind, method=method, x=sample, kappa=kr, varpi=vr,
            substeps=substeps,
        )
        candidate = PhaseCovariantCandidate(
            kind, method, sample, tr, kr, vr, a, c, stages,
            (kind,), {"source": "exact rational component sample"},
        )
        row["sample"] = sample
        telemetry["certification_ops"] += 1
        row["status"] = certify_candidate(candidate).status
        if row["status"] is CPStatus.UNCERTAIN:
            interval_uncertain = True

    adjacent_pass = [False] * len(roots)
    for row in interval_rows:
        if row["status"] is not CPStatus.PASS:
            continue
        if row["left"] is not None:
            adjacent_pass[int(row["left"])] = True
        if row["right"] is not None:
            adjacent_pass[int(row["right"])] = True

    final_roots: list[RootIsolation] = []
    unresolved_points = 0
    for index, root in enumerate(roots):
        exact_root = root.exact_root
        point_status = CPStatus.UNCERTAIN
        if not adjacent_pass[index]:
            if exact_root is None:
                exact_root = exact_root_for_isolation(root, polys)
            root_with_exact = RootIsolation(
                root.lower, root.upper, root.multiplicity, root.source, exact_root
            )
            telemetry["certification_ops"] += 1
            point_status = _certify_isolated_candidate_root(
                kind, root_with_exact, method=method, theta=tr, kappa=kr,
                varpi=vr, substeps=substeps
            )
            if point_status is CPStatus.UNCERTAIN:
                unresolved_points += 1
        final_roots.append(
            RootIsolation(
                root.lower, root.upper, root.multiplicity, root.source,
                exact_root, point_status,
            )
        )

    components: list[CertifiedComponent] = []
    for row in interval_rows:
        if row["status"] is not CPStatus.PASS:
            continue
        left = None if row["left"] is None else final_roots[int(row["left"])]
        right = None if row["right"] is None else final_roots[int(row["right"])]
        components.append(
            CertifiedComponent(
                left, right, sp.Rational(row["sample"]),
                sp.Rational(row["lo"]), sp.Rational(row["hi"]), CPStatus.PASS,
            )
        )

    isolated = tuple(
        root
        for index, root in enumerate(final_roots)
        if not adjacent_pass[index] and root.point_status is CPStatus.PASS
    )
    status = CPStatus.UNCERTAIN if unresolved_points or interval_uncertain else CPStatus.PASS
    return ComponentLocation(
        status, tuple(components), isolated, tuple(final_roots), ur,
        {
            "candidate_kind": kind,
            "root_method": "exact Sturm isolation with adaptive algebraic separation",
            "sample_arithmetic": "exact rational",
            "isolated_points": "listed only after exact joint constraint evaluation",
            "completeness_scope": "rational-coefficient one-dimensional candidate constraints with exact requested-domain endpoint membership",
            "unresolved_point_count": str(unresolved_points),
            **telemetry,
        },
    )


def richardson_rk4_margin_nonrotating() -> sp.Expr:
    a, c, _ = candidate_symbolic_expressions(
        "richardson-4", method="rk4", x=X, kappa=0, varpi=0
    )
    return sp.factor(a - c**2)


def richardson_rk4_stability_function() -> sp.Expr:
    R = _method("rk4").stability_function(x=S)
    return sp.factor((16 * R.subs(S, S / 2) ** 2 - R) / 15)
