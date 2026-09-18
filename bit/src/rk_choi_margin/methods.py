from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

X = sp.Symbol("x")
Q = sp.Rational


@dataclass(frozen=True)
class RKMethod:
    """Exact Butcher data and metadata for a Runge-Kutta formula."""

    key: str
    name: str
    A: tuple[tuple[sp.Expr, ...], ...]
    b: tuple[sp.Expr, ...]
    order: int
    family: str
    notes: str = ""

    @property
    def stages(self) -> int:
        return len(self.b)

    def stability_function(self, x: sp.Symbol = X) -> sp.Expr:
        return stability_function(self.A, self.b, x=x)

    def stage_resolvent_determinant(self, x: sp.Symbol = X) -> sp.Expr:
        """Return det(I-xA), whose nonzero set is the RK stage domain."""
        return stage_resolvent_determinant(self.A, x=x)


def stage_resolvent_determinant(
    A: tuple[tuple[sp.Expr, ...], ...] | list[list[sp.Expr]],
    *,
    x: sp.Symbol = X,
) -> sp.Expr:
    """Return det(I-xA), before any cancellation in the stability function."""
    A_m = sp.Matrix(A)
    if A_m.rows != A_m.cols:
        raise ValueError("A must be square")
    return sp.factor((sp.eye(A_m.rows) - x * A_m).det())


def stability_function(
    A: tuple[tuple[sp.Expr, ...], ...] | list[list[sp.Expr]],
    b: tuple[sp.Expr, ...] | list[sp.Expr],
    *,
    x: sp.Symbol = X,
) -> sp.Expr:
    """Return the reduced expression R(x)=1+x b^T(I-xA)^(-1)1.

    The operational Runge--Kutta evaluation domain remains the nonzero set of
    ``stage_resolvent_determinant(A, x)``.  A removable cancellation in this
    reduced scalar expression does not regularize singular stage equations.
    """
    A_m = sp.Matrix(A)
    b_m = sp.Matrix(b)
    if A_m.rows != A_m.cols:
        raise ValueError("A must be square")
    if b_m.shape != (A_m.rows, 1):
        raise ValueError("b must contain one weight per stage")
    ones = sp.ones(A_m.rows, 1)
    R = 1 + x * (b_m.T * (sp.eye(A_m.rows) - x * A_m).inv() * ones)[0]
    return sp.factor(sp.cancel(R))


def all_methods() -> dict[str, RKMethod]:
    """Selected exact methods used in the theorem audit."""
    euler = RKMethod(
        key="euler",
        name="Forward Euler",
        A=((Q(0),),),
        b=(Q(1),),
        order=1,
        family="explicit",
    )

    heun2 = RKMethod(
        key="rk2",
        name="Heun RK2",
        A=((Q(0), Q(0)), (Q(1), Q(0))),
        b=(Q(1, 2), Q(1, 2)),
        order=2,
        family="explicit",
        notes="Any two-stage, second-order explicit RK formula has the same scalar stability polynomial.",
    )

    ssprk3 = RKMethod(
        key="rk3",
        name="SSPRK(3,3)",
        A=(
            (Q(0), Q(0), Q(0)),
            (Q(1), Q(0), Q(0)),
            (Q(1, 4), Q(1, 4), Q(0)),
        ),
        b=(Q(1, 6), Q(1, 6), Q(2, 3)),
        order=3,
        family="explicit",
        notes="Representative minimal-stage third-order method; scalar stability polynomial is T_3.",
    )

    rk4 = RKMethod(
        key="rk4",
        name="Classical RK4",
        A=(
            (Q(0), Q(0), Q(0), Q(0)),
            (Q(1, 2), Q(0), Q(0), Q(0)),
            (Q(0), Q(1, 2), Q(0), Q(0)),
            (Q(0), Q(0), Q(1), Q(0)),
        ),
        b=(Q(1, 6), Q(1, 3), Q(1, 3), Q(1, 6)),
        order=4,
        family="explicit",
    )

    A_dp = (
        (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0), Q(0)),
        (Q(1, 5), Q(0), Q(0), Q(0), Q(0), Q(0), Q(0)),
        (Q(3, 40), Q(9, 40), Q(0), Q(0), Q(0), Q(0), Q(0)),
        (Q(44, 45), -Q(56, 15), Q(32, 9), Q(0), Q(0), Q(0), Q(0)),
        (
            Q(19372, 6561),
            -Q(25360, 2187),
            Q(64448, 6561),
            -Q(212, 729),
            Q(0),
            Q(0),
            Q(0),
        ),
        (
            Q(9017, 3168),
            -Q(355, 33),
            Q(46732, 5247),
            Q(49, 176),
            -Q(5103, 18656),
            Q(0),
            Q(0),
        ),
        (
            Q(35, 384),
            Q(0),
            Q(500, 1113),
            Q(125, 192),
            -Q(2187, 6784),
            Q(11, 84),
            Q(0),
        ),
    )

    dp5 = RKMethod(
        key="dp5",
        name="Dormand-Prince 5 principal formula",
        A=A_dp,
        b=(Q(35, 384), Q(0), Q(500, 1113), Q(125, 192), -Q(2187, 6784), Q(11, 84), Q(0)),
        order=5,
        family="explicit embedded",
        notes="Principal fifth-order member of the Dormand-Prince 5(4) pair.",
    )

    dp4 = RKMethod(
        key="dp4",
        name="Dormand-Prince embedded 4 formula",
        A=A_dp,
        b=(Q(5179, 57600), Q(0), Q(7571, 16695), Q(393, 640), -Q(92097, 339200), Q(187, 2100), Q(1, 40)),
        order=4,
        family="explicit embedded",
        notes="Embedded fourth-order member; useful as a disconnected-admissible-set example.",
    )

    backward_euler = RKMethod(
        key="be",
        name="Backward Euler",
        A=((Q(1),),),
        b=(Q(1),),
        order=1,
        family="implicit",
        notes="Simple non-parity counterexample: locally and globally CP-admissible for this generator.",
    )

    implicit_midpoint = RKMethod(
        key="im",
        name="Implicit midpoint",
        A=((Q(1, 2),),),
        b=(Q(1),),
        order=2,
        family="implicit",
        notes="Even-order non-parity counterexample: no nonzero local CP-admissible window at r=1/2.",
    )

    methods = [euler, heun2, ssprk3, rk4, dp5, dp4, backward_euler, implicit_midpoint]
    return {method.key: method for method in methods}
