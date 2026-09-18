from __future__ import annotations

"""Exact and numerical superoperator utilities for qubit GKSL benchmarks."""

from dataclasses import dataclass
from itertools import product

import mpmath as mp
import numpy as np
import sympy as sp

from .certification import (
    CPStatus,
    CertifiedComponent,
    RootIsolation,
    as_rational,
    certify_hermitian_psd_exact,
    isolate_constraint_roots_exact,
)


@dataclass(frozen=True)
class ChoiCertificate:
    status: CPStatus
    choi: sp.Matrix
    principal_minors: dict[tuple[int, ...], sp.Expr]


def _exact_qubit_operators() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    I = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])
    sm = sp.Matrix([[0, 1], [0, 0]])
    return I, sx, sz, sm


def exact_gksl_superoperator(
    *,
    gamma: int | float | str | sp.Rational = 1,
    theta: int | float | str | sp.Rational = 0,
    kappa: int | float | str | sp.Rational = 0,
    omega_z: int | float | str | sp.Rational = 0,
    omega_x: int | float | str | sp.Rational = 0,
) -> sp.Matrix:
    """Column-vectorized qubit GKSL generator in the basis vec_F.

    ``vec_F`` stacks columns, so ``vec(A X B)=(B^T \\otimes A)vec(X)``.
    The dissipative rates are ``gamma_down=gamma*(1-theta)``,
    ``gamma_up=gamma*theta``, and ``gamma_phi=gamma*kappa``.
    """

    g, th, kap, wz, wx = map(as_rational, (gamma, theta, kappa, omega_z, omega_x))
    if g < 0 or th < 0 or th > 1 or kap < 0:
        raise ValueError("parameters outside GKSL domain")
    I, sx, sz, sm = _exact_qubit_operators()
    spm = sm.T
    H = -wz * sz / 2 + wx * sx / 2
    L = -sp.I * (sp.kronecker_product(I, H) - sp.kronecker_product(H.T, I))

    def dissipator(op: sp.Matrix) -> sp.Matrix:
        q = op.conjugate().T * op
        return (
            sp.kronecker_product(sp.conjugate(op), op)
            - sp.kronecker_product(I, q) / 2
            - sp.kronecker_product(q.T, I) / 2
        )

    L += g * (1 - th) * dissipator(sm)
    L += g * th * dissipator(spm)
    # (gamma_phi/2)(sigma_z rho sigma_z-rho) = gamma_phi/2 * D[sigma_z]
    L += g * kap / 2 * dissipator(sz)
    return sp.simplify(L)


def rk4_superoperator_exact(L: sp.Matrix, h: int | float | str | sp.Rational) -> sp.Matrix:
    hr = h if isinstance(h, sp.Basic) else as_rational(h)
    Z = hr * L
    return sp.eye(L.rows) + Z + Z**2 / 2 + Z**3 / 6 + Z**4 / 24


def _mp_from_rational(value: int | float | str | sp.Rational) -> mp.mpf:
    """Convert the exact rational represented by ``value`` to mpmath."""
    rational = as_rational(value)
    return mp.mpf(int(rational.p)) / mp.mpf(int(rational.q))


def _mp_from_sympy(value: sp.Expr) -> mp.mpc:
    value = sp.expand_complex(value)
    return mp.mpc(_mp_from_rational(sp.re(value)), _mp_from_rational(sp.im(value)))


def exp_superoperator_numeric(
    L: np.ndarray, h: int | float | str | sp.Rational, dps: int = 80
) -> np.ndarray:
    """Exponentiate a supplied numerical matrix at the exact supplied time."""
    mp.mp.dps = dps
    M = mp.matrix([[mp.mpc(complex(v)) for v in row] for row in np.asarray(L)])
    E = mp.expm(_mp_from_rational(h) * M)
    return np.array([[complex(E[i, j]) for j in range(E.cols)] for i in range(E.rows)], dtype=np.complex128)


def exp_exact_gksl_superoperator_numeric(
    *, h: int | float | str | sp.Rational, gamma=1, theta=0, kappa=0,
    omega_z=0, omega_x=0, dps: int = 80,
) -> np.ndarray:
    """Derive and exponentiate a GKSL generator from exact recipe fields."""
    mp.mp.dps = dps
    L_exact = exact_gksl_superoperator(
        gamma=gamma, theta=theta, kappa=kappa, omega_z=omega_z, omega_x=omega_x)
    M = mp.matrix([[_mp_from_sympy(L_exact[i,j]) for j in range(L_exact.cols)]
                   for i in range(L_exact.rows)])
    E = mp.expm(_mp_from_rational(h) * M)
    return np.array([[complex(E[i,j]) for j in range(E.cols)] for i in range(E.rows)], dtype=np.complex128)


def matrix_unit(j: int, k: int, dim: int = 2) -> sp.Matrix:
    E = sp.zeros(dim)
    E[j, k] = 1
    return E


def vec_col(matrix: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([matrix[i, j] for j in range(matrix.cols) for i in range(matrix.rows)])


def unvec_col(vector: sp.Matrix, dim: int) -> sp.Matrix:
    return sp.Matrix(dim, dim, lambda i, j: vector[j * dim + i])


def choi_from_superoperator_exact(S: sp.Matrix, dim: int = 2) -> sp.Matrix:
    """Return J=sum_jk |j><k| tensor Phi(|j><k|)."""

    J = sp.zeros(dim * dim)
    for j, k in product(range(dim), repeat=2):
        out = unvec_col(S * vec_col(matrix_unit(j, k, dim)), dim)
        J += sp.kronecker_product(matrix_unit(j, k, dim), out)
    return sp.simplify(J)


def choi_from_superoperator_numeric(S: np.ndarray, dim: int = 2) -> np.ndarray:
    J = np.zeros((dim * dim, dim * dim), dtype=np.complex128)
    for j, k in product(range(dim), repeat=2):
        E = np.zeros((dim, dim), dtype=np.complex128)
        E[j, k] = 1.0
        out = (S @ E.reshape(-1, order="F")).reshape((dim, dim), order="F")
        block = np.zeros_like(J)
        block[j * dim : (j + 1) * dim, k * dim : (k + 1) * dim] = out
        J += block
    return (J + J.conjugate().T) / 2


def certify_rk4_gksl_step_exact(
    *,
    h: int | float | str | sp.Rational,
    gamma: int | float | str | sp.Rational = 1,
    theta: int | float | str | sp.Rational = 0,
    kappa: int | float | str | sp.Rational = 0,
    omega_z: int | float | str | sp.Rational = 0,
    omega_x: int | float | str | sp.Rational = 0,
) -> ChoiCertificate:
    L = exact_gksl_superoperator(
        gamma=gamma,
        theta=theta,
        kappa=kappa,
        omega_z=omega_z,
        omega_x=omega_x,
    )
    S = rk4_superoperator_exact(L, h)
    J = choi_from_superoperator_exact(S)
    status, minors = certify_hermitian_psd_exact(J)
    return ChoiCertificate(status, J, minors)


def min_choi_eigenvalue_numeric(S: np.ndarray, dim: int = 2) -> float:
    return float(np.linalg.eigvalsh(choi_from_superoperator_numeric(S, dim))[0].real)


def normalized_choi_trace_distance(S1: np.ndarray, S2: np.ndarray, dim: int = 2) -> float:
    J1 = choi_from_superoperator_numeric(S1, dim) / dim
    J2 = choi_from_superoperator_numeric(S2, dim) / dim
    return float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(J1 - J2))))


def rk4_superoperator_numeric(L: np.ndarray, h: float) -> np.ndarray:
    Z = float(h) * np.asarray(L, dtype=np.complex128)
    I = np.eye(Z.shape[0], dtype=np.complex128)
    return I + Z + Z @ Z / 2 + Z @ Z @ Z / 6 + Z @ Z @ Z @ Z / 24


def exact_to_numpy(matrix: sp.Matrix, digits: int = 50) -> np.ndarray:
    return np.array([[complex(sp.N(matrix[i, j], digits)) for j in range(matrix.cols)] for i in range(matrix.rows)], dtype=np.complex128)


def rotating_frame_rk4_superoperator_numeric(
    *,
    h: float,
    gamma: float = 1.0,
    theta: float = 0.0,
    kappa: float = 0.0,
    omega_z: float = 0.0,
) -> np.ndarray:
    """RK4 on the exactly de-rotated dissipative generator, then exact lab rotation."""

    Ld = exact_to_numpy(exact_gksl_superoperator(gamma=gamma, theta=theta, kappa=kappa, omega_z=0, omega_x=0))
    Sd = rk4_superoperator_numeric(Ld, h)
    # U=exp(-i H h), H=-(omega_z/2)sigma_z; vec(U rho U^dag)=conj(U) kron U.
    U = np.diag([np.exp(1j * omega_z * h / 2), np.exp(-1j * omega_z * h / 2)]).astype(np.complex128)
    Su = np.kron(np.conjugate(U), U)
    return Su @ Sd


def rk4_gksl_principal_minor_polynomials(
    *,
    gamma: int | float | str | sp.Rational = 1,
    theta: int | float | str | sp.Rational = 0,
    kappa: int | float | str | sp.Rational = 0,
    omega_z: int | float | str | sp.Rational = 0,
    omega_x: int | float | str | sp.Rational = 0,
    h_symbol: sp.Symbol | None = None,
) -> tuple[sp.Symbol, sp.Matrix, dict[tuple[int, ...], sp.Poly]]:
    """Exact Choi matrix and primitive principal-minor polynomials in h."""

    h = h_symbol or sp.Symbol("h", nonnegative=True, real=True)
    L = exact_gksl_superoperator(
        gamma=gamma,
        theta=theta,
        kappa=kappa,
        omega_z=omega_z,
        omega_x=omega_x,
    )
    J = choi_from_superoperator_exact(rk4_superoperator_exact(L, h))
    from .certification import all_principal_minors

    polys: dict[tuple[int, ...], sp.Poly] = {}
    for idx, expr in all_principal_minors(J).items():
        num, _den = sp.fraction(sp.cancel(expr))
        poly = sp.Poly(num, h, domain=sp.QQ)
        _, primitive = poly.primitive()
        polys[idx] = primitive
    return h, J, polys


def locate_rk4_gksl_cptp_components_exact(
    *,
    upper: int | float | str | sp.Rational,
    gamma: int | float | str | sp.Rational = 1,
    theta: int | float | str | sp.Rational = 0,
    kappa: int | float | str | sp.Rational = 0,
    omega_z: int | float | str | sp.Rational = 0,
    omega_x: int | float | str | sp.Rational = 0,
    root_digits: int = 24,
) -> dict[str, object]:
    """Certify map-level CPTP components on the closed domain ``[0, upper]``.

    Raw Sturm brackets are first restricted by exact root membership in the
    requested domain.  They are never clipped into the domain before that
    decision.  A root at ``upper`` is retained as an exact equality point;
    roots strictly above ``upper`` are discarded.  Any unresolved separation or
    equality-point sign decision returns ``UNCERTAIN`` and no usable component.

    ``components`` retains the historical ``(lower, upper, sample)`` tuples.
    ``component_records`` and ``roots`` expose the complete certified endpoint
    semantics used by the phase-covariant and candidate-map locators.
    """

    ur = as_rational(upper)
    if ur <= 0:
        return {
            "status": CPStatus.FAIL,
            "roots": tuple(),
            "components": tuple(),
            "component_records": tuple(),
            "isolated_points": tuple(),
            "excluded_points": tuple(),
            "diagnostics": {"reason": "upper must be positive"},
            "root_digits": root_digits,
        }

    h, J, indexed_polys = rk4_gksl_principal_minor_polynomials(
        gamma=gamma,
        theta=theta,
        kappa=kappa,
        omega_z=omega_z,
        omega_x=omega_x,
    )
    named_polys = {str(index): poly for index, poly in indexed_polys.items()}
    telemetry = {
        "root_isolation_ops": 0,
        "certification_ops": 0,
        "precision_escalations": 0,
    }
    roots, separated = isolate_constraint_roots_exact(
        named_polys,
        ur,
        digits=root_digits,
        telemetry=telemetry,
    )
    if not separated:
        return {
            "status": CPStatus.UNCERTAIN,
            "roots": tuple(),
            "components": tuple(),
            "component_records": tuple(),
            "isolated_points": tuple(),
            "excluded_points": tuple(),
            "choi": J,
            "principal_minor_polynomials": indexed_polys,
            "root_digits": root_digits,
            "diagnostics": {
                "root_method": "exact Sturm isolation with adaptive algebraic separation",
                "reason": "distinct boundary roots could not be certified as separated",
                **telemetry,
            },
        }

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
        telemetry["certification_ops"] += 1
        status, _minors = certify_hermitian_psd_exact(sp.simplify(J.subs(h, sample)))
        row["sample"] = sample
        row["status"] = status
        if status is CPStatus.UNCERTAIN:
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
        telemetry["certification_ops"] += 1
        if adjacent_pass[index]:
            # The PSD cone is closed.  A separated boundary root in the closure
            # of a certified PASS interval is therefore an included equality
            # point, without constructing a costly CRootOf representation.
            point_status = CPStatus.PASS
        else:
            # Each source minor vanishes at the unique root in this Sturm
            # bracket.  Every nonsource minor has no root in the bracket after
            # global separation, so its rational-midpoint sign is its sign at
            # the algebraic equality point.
            source = set(root.source)
            statuses: list[CPStatus] = []
            for name, poly in named_polys.items():
                if name in source:
                    statuses.append(CPStatus.PASS)
                    continue
                value = sp.Rational(poly.eval(root.midpoint))
                if value > 0:
                    statuses.append(CPStatus.PASS)
                elif value < 0:
                    statuses.append(CPStatus.FAIL)
                else:
                    statuses.append(CPStatus.UNCERTAIN)
            point_status = (
                CPStatus.FAIL
                if CPStatus.FAIL in statuses
                else CPStatus.UNCERTAIN
                if CPStatus.UNCERTAIN in statuses
                else CPStatus.PASS
            )
        if point_status is CPStatus.UNCERTAIN:
            unresolved_points += 1
        final_roots.append(
            RootIsolation(
                root.lower,
                root.upper,
                root.multiplicity,
                root.source,
                root.exact_root,
                point_status,
            )
        )

    component_records: list[CertifiedComponent] = []
    for row in interval_rows:
        if row["status"] is not CPStatus.PASS:
            continue
        left = None if row["left"] is None else final_roots[int(row["left"])]
        right = None if row["right"] is None else final_roots[int(row["right"])]
        component_records.append(
            CertifiedComponent(
                left,
                right,
                sp.Rational(row["sample"]),
                sp.Rational(row["lo"]),
                sp.Rational(row["hi"]),
                CPStatus.PASS,
            )
        )

    components = tuple(
        (component.lower_value, component.upper_value, component.sample)
        for component in component_records
    )
    isolated_points = tuple(
        root
        for index, root in enumerate(final_roots)
        if not adjacent_pass[index] and root.point_status is CPStatus.PASS
    )
    excluded_points = tuple(
        root
        for index, root in enumerate(final_roots)
        if not adjacent_pass[index] and root.point_status is CPStatus.FAIL
    )
    status = CPStatus.UNCERTAIN if interval_uncertain or unresolved_points else CPStatus.PASS
    return {
        "status": status,
        "roots": tuple(final_roots),
        "components": components,
        "component_records": tuple(component_records),
        "isolated_points": isolated_points,
        "excluded_points": excluded_points,
        "choi": J,
        "principal_minor_polynomials": indexed_polys,
        "root_digits": root_digits,
        "diagnostics": {
            "root_method": "exact Sturm isolation with adaptive algebraic separation",
            "decision_arithmetic": "exact rational interval samples and exact algebraic equality-point signs",
            "endpoint_semantics": "exact membership in [0, upper] precedes any bracket restriction; upper-boundary equality is retained",
            "unresolved_point_count": unresolved_points,
            **telemetry,
        },
    }

def hermitian_eigenvalue_brackets_exact(matrix: sp.Matrix, digits: int = 30) -> list[tuple[sp.Rational, sp.Rational, int]]:
    """Isolate all real eigenvalues of an exact Hermitian matrix."""

    cp = matrix.charpoly()
    poly = sp.Poly(cp.as_expr(), cp.gen, domain=sp.QQ)
    rows: list[tuple[sp.Rational, sp.Rational, int]] = []
    for (lo, hi), mult in poly.intervals(eps=sp.Rational(1, 10**digits)):
        rows.append((sp.Rational(lo), sp.Rational(hi), int(mult)))
    rows.sort(key=lambda row: row[0])
    return rows


def strang_gksl_superoperator_numeric(
    *, h: int | float | str | sp.Rational, gamma=1, theta=0, kappa=0,
    omega_z=0, omega_x=0, substeps: int = 1, dps: int = 80,
) -> np.ndarray:
    """Strang composition derived from exact Hamiltonian/dissipative recipes."""
    if substeps < 1:
        raise ValueError("substeps must be positive")
    total = as_rational(h); dt = total/substeps; half = dt/2
    UH = exp_exact_gksl_superoperator_numeric(
        h=half, gamma=0, theta=0, kappa=0,
        omega_z=omega_z, omega_x=omega_x, dps=dps)
    UD = exp_exact_gksl_superoperator_numeric(
        h=dt, gamma=gamma, theta=theta, kappa=kappa,
        omega_z=0, omega_x=0, dps=dps)
    one = UH @ UD @ UH
    total_map = np.eye(one.shape[0], dtype=np.complex128)
    for _ in range(substeps): total_map = one @ total_map
    return total_map

