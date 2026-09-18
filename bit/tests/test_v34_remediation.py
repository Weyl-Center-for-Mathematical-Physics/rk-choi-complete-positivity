from __future__ import annotations

import sympy as sp

from rk_choi_margin import liouvillian
from rk_choi_margin.certification import CPStatus
from rk_choi_margin.controllers import certified_candidate_guard_step


def test_general_qubit_locator_discards_root_strictly_above_upper_before_clipping(monkeypatch) -> None:
    h = sp.Symbol("h", nonnegative=True, real=True)
    upper = sp.Rational(1414213562, 10**9)  # strictly below sqrt(2), within 1e-8

    def fake_polynomials(**_kwargs):
        return h, sp.Matrix([[1]]), {(0,): sp.Poly(h**2 - 2, h, domain=sp.QQ)}

    monkeypatch.setattr(liouvillian, "rk4_gksl_principal_minor_polynomials", fake_polynomials)
    result = liouvillian.locate_rk4_gksl_cptp_components_exact(
        upper=upper, root_digits=8
    )

    assert result["status"] is CPStatus.PASS
    assert result["roots"] == tuple()
    assert result["components"] == ((sp.Rational(0), upper, upper / 2),)


def test_general_qubit_locator_retains_and_certifies_exact_upper_boundary_root(monkeypatch) -> None:
    h = sp.Symbol("h", nonnegative=True, real=True)

    def fake_polynomials(**_kwargs):
        return h, sp.Matrix([[1 - h]]), {(0,): sp.Poly(1 - h, h, domain=sp.QQ)}

    monkeypatch.setattr(liouvillian, "rk4_gksl_principal_minor_polynomials", fake_polynomials)
    result = liouvillian.locate_rk4_gksl_cptp_components_exact(
        upper=1, root_digits=8
    )

    assert result["status"] is CPStatus.PASS
    assert len(result["roots"]) == 1
    root = result["roots"][0]
    assert root.lower == root.upper == root.exact_root == 1
    assert root.point_status is CPStatus.PASS
    assert result["components"] == ((sp.Rational(0), sp.Rational(1), sp.Rational(1, 2)),)


def test_general_qubit_locator_fails_closed_on_unresolved_root_separation(monkeypatch) -> None:
    h = sp.Symbol("h", nonnegative=True, real=True)

    def fake_polynomials(**_kwargs):
        return h, sp.Matrix([[1]]), {(0,): sp.Poly(h - 1, h, domain=sp.QQ)}

    def unresolved(*_args, **_kwargs):
        return [], False

    monkeypatch.setattr(liouvillian, "rk4_gksl_principal_minor_polynomials", fake_polynomials)
    monkeypatch.setattr(liouvillian, "isolate_constraint_roots_exact", unresolved)
    result = liouvillian.locate_rk4_gksl_cptp_components_exact(upper=2)

    assert result["status"] is CPStatus.UNCERTAIN
    assert result["components"] == tuple()


def test_controller_locator_telemetry_partitions_every_call() -> None:
    decision = certified_candidate_guard_step(
        method="rk4",
        candidate_kind="equal-substeps",
        substeps=2,
        proposed_x=1.4,
        theta=0.0,
        kappa=0.0,
        varpi=2.0,
    )
    telemetry = decision.telemetry

    assert telemetry["locator_calls"] == 1
    assert telemetry["locator_calls"] == (
        telemetry["successful_projections"]
        + telemetry["empty_component_outcomes"]
        + telemetry["uncertain_outcomes"]
    )
    assert telemetry["root_isolation_ops"] >= telemetry["locator_calls"]
    assert telemetry["certification_ops"] >= 1
