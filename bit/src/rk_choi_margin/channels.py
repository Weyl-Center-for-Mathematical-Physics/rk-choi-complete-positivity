from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def choi_matrix(a: complex, c: complex, *, normalized: bool = False) -> NDArray[np.complex128]:
    """Unnormalized Choi matrix of the one-way excitation-damping map."""
    J = np.array(
        [
            [1.0, 0.0, 0.0, c],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0 - a, 0.0],
            [np.conjugate(c), 0.0, 0.0, a],
        ],
        dtype=np.complex128,
    )
    return J / 2.0 if normalized else J


def thermal_choi_matrix(
    a: complex,
    c: complex,
    theta: float,
    *,
    normalized: bool = False,
) -> NDArray[np.complex128]:
    """Unnormalized Choi matrix of the phase-covariant thermal qubit map."""
    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must lie in [0, 1]")
    if abs(complex(a).imag) > 1e-12:
        raise ValueError("a must be real")
    ar = float(complex(a).real)
    A = 1.0 - theta * (1.0 - ar)
    B = theta * (1.0 - ar)
    C = (1.0 - theta) * (1.0 - ar)
    D = theta + (1.0 - theta) * ar
    J = np.array(
        [[A, 0.0, 0.0, c], [0.0, B, 0.0, 0.0], [0.0, 0.0, C, 0.0], [np.conjugate(c), 0.0, 0.0, D]],
        dtype=np.complex128,
    )
    return J / 2.0 if normalized else J


def cp_conditions(a: complex, c: complex, *, atol: float = 1e-12) -> bool:
    """Necessary-and-sufficient CPTP test for the one-way map."""
    if abs(a.imag) > atol:
        return False
    ar = float(a.real)
    return (-atol <= ar <= 1.0 + atol) and (abs(c) ** 2 <= ar + atol)


def thermal_cp_conditions(a: complex, c: complex, theta: float, *, atol: float = 1e-12) -> bool:
    """Necessary-and-sufficient CPTP test for the thermal phase-covariant map."""
    if not -atol <= theta <= 1.0 + atol or abs(a.imag) > atol:
        return False
    ar = float(a.real)
    A = 1.0 - theta * (1.0 - ar)
    D = theta + (1.0 - theta) * ar
    margin = ar + theta * (1.0 - theta) * (1.0 - ar) ** 2 - abs(c) ** 2
    return ar <= 1.0 + atol and A >= -atol and D >= -atol and margin >= -atol


def minimum_choi_eigenvalue(a: complex, c: complex, *, normalized: bool = False) -> float:
    vals = np.linalg.eigvalsh(choi_matrix(a, c, normalized=normalized))
    return float(vals[0].real)


def apply_phase_covariant_map(matrix: NDArray[np.complex128], a: complex, c: complex) -> NDArray[np.complex128]:
    """Apply the one-way trace-preserving phase-covariant map."""
    X = np.asarray(matrix, dtype=np.complex128)
    if X.shape != (2, 2):
        raise ValueError("matrix must have shape (2, 2)")
    return np.array(
        [[X[0, 0] + (1.0 - a) * X[1, 1], c * X[0, 1]], [np.conjugate(c) * X[1, 0], a * X[1, 1]]],
        dtype=np.complex128,
    )


def apply_thermal_phase_covariant_map(
    matrix: NDArray[np.complex128],
    a: complex,
    c: complex,
    theta: float,
) -> NDArray[np.complex128]:
    """Apply the bidirectional phase-covariant map used in the article."""
    X = np.asarray(matrix, dtype=np.complex128)
    if X.shape != (2, 2):
        raise ValueError("matrix must have shape (2, 2)")
    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must lie in [0, 1]")
    A = 1.0 - theta * (1.0 - a)
    B = theta * (1.0 - a)
    C = (1.0 - theta) * (1.0 - a)
    D = theta + (1.0 - theta) * a
    return np.array(
        [[A * X[0, 0] + C * X[1, 1], c * X[0, 1]], [np.conjugate(c) * X[1, 0], B * X[0, 0] + D * X[1, 1]]],
        dtype=np.complex128,
    )


def repeated_margin(a: float, c: complex, steps: int) -> float:
    if steps < 1:
        raise ValueError("steps must be positive")
    return float(a**steps - abs(c) ** (2 * steps))
