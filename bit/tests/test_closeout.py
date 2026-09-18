from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sp

from rk_choi_margin.channels import choi_matrix, cp_conditions
from rk_choi_margin.methods import X, all_methods
from rk_choi_margin.symbolic import cp_margin, first_exponential_defect, predicted_boundary_coefficient

ROOT = Path(__file__).resolve().parents[1]
CERT_PATH = ROOT / "results" / "closeout" / "root_sign_certificate_v2_7.json"
TABLE_PATH = ROOT / "results" / "verified_method_table.csv"
FIGURE_DATA_PATH = ROOT / "results" / "benchmark_figures" / "figure3_fixed_horizon.csv"

x = sp.Symbol("x", nonnegative=True, real=True)


def _certificate() -> dict:
    return json.loads(CERT_PATH.read_text(encoding="utf-8"))


def _sign(expr: sp.Expr, point: sp.Rational) -> int:
    value = sp.simplify(expr.subs(x, point))
    if value == 0:
        return 0
    return 1 if value > 0 else -1


def _lower_choi_eigenvalue(a: mp.mpf, c: mp.mpf) -> mp.mpf:
    margin = a - c**2
    denominator = 1 + a + mp.sqrt((1 - a) ** 2 + 4 * c**2)
    return 2 * margin / denominator


def test_closeout_certificate_root_brackets_are_sturm_certified() -> None:
    data = _certificate()
    for method in data["methods"]:
        for source in ("a", "a-1", "M"):
            poly = sp.Poly(sp.sympify(method["polynomials"][source]["primitive_numerator"], locals={"x": x}), x, domain=sp.QQ)
            for root in method["polynomials"][source]["roots"]:
                mult = int(root["multiplicity"])
                if root["exact"] is not None:
                    exact = sp.Rational(root["exact"])
                    assert poly.eval(exact) == 0
                    derivative = poly
                    observed_mult = 0
                    while derivative.degree() >= 0 and derivative.eval(exact) == 0:
                        observed_mult += 1
                        derivative = derivative.diff()
                    assert observed_mult == mult
                else:
                    lo = sp.Rational(root["lo"])
                    hi = sp.Rational(root["hi"])
                    assert poly.eval(lo) != 0 and poly.eval(hi) != 0
                    assert int(poly.count_roots(lo, hi)) == mult




def test_closeout_reported_root_approximations_match_refined_certified_roots() -> None:
    data = _certificate()
    for method in data["methods"]:
        for source in ("a", "a-1", "M"):
            poly = sp.Poly(sp.sympify(method["polynomials"][source]["primitive_numerator"], locals={"x": x}), x, domain=sp.QQ)
            for root in method["polynomials"][source]["roots"]:
                if root["exact"] is not None:
                    continue
                lo = sp.Rational(root["lo"])
                hi = sp.Rational(root["hi"])
                refined_lo, refined_hi = poly.refine_root(lo, hi, eps=sp.Rational(1, 10**28), check_sqf=False)
                midpoint = (sp.Rational(refined_lo) + sp.Rational(refined_hi)) / 2
                reported = sp.Rational(root["approx"])
                assert abs(reported - midpoint) < sp.Rational(6, 10**16)
                displayed_lo = sp.Rational(root["lo_decimal"])
                displayed_hi = sp.Rational(root["hi_decimal"])
                assert displayed_lo <= reported <= displayed_hi


def test_closeout_certificate_sign_rows_recompute_exactly() -> None:
    data = _certificate()
    methods = all_methods()
    for method_data in data["methods"]:
        R = methods[method_data["key"]].stability_function()
        a = sp.factor(R.subs(X, -x))
        M = sp.factor(a - R.subs(X, -x / 2) ** 2)
        a_minus_1 = sp.factor(a - 1)
        det = methods[method_data["key"]].stage_resolvent_determinant()
        d1 = sp.factor(det.subs(X, -x))
        d2 = sp.factor(det.subs(X, -x / 2))
        for row in method_data["sign_rows"]:
            sample = sp.Rational(row["sample"])
            assert _sign(a, sample) == int(row["sign_a"])
            assert _sign(a_minus_1, sample) == int(row["sign_a_minus_1"])
            assert _sign(M, sample) == int(row["sign_M"])
            domain = d1.subs(x, sample) != 0 and d2.subs(x, sample) != 0
            cptp = domain and a.subs(x, sample) <= 1 and M.subs(x, sample) >= 0
            assert bool(domain) == bool(row["domain"])
            assert bool(cptp) == bool(row["cptp"])


def test_closeout_dp4_certificate_excludes_all_later_components() -> None:
    dp4 = next(item for item in _certificate()["methods"] if item["key"] == "dp4")
    assert [row["cptp"] for row in dp4["sign_rows"]] == [False, True, False, False, False, False]
    events = dp4["positive_events"]
    assert [event["labels"][0] for event in events] == [
        r"\beta_-",
        r"\beta_+",
        r"\delta_M",
        r"\delta_1",
        r"\delta_0",
    ]
    assert events[0]["sources"] == ["M"]
    assert events[1]["sources"] == ["a-1"]
    assert events[-1]["sources"] == ["a"]


def test_closeout_stage_resolvent_domain_has_no_nonnegative_poles_for_audited_methods() -> None:
    for method in all_methods().values():
        det = method.stage_resolvent_determinant()
        for scaled in (sp.factor(det.subs(X, -x)), sp.factor(det.subs(X, -x / 2))):
            poly = sp.Poly(sp.together(scaled).as_numer_denom()[0], x, domain=sp.QQ)
            for interval, _mult in poly.intervals():
                lo, hi = interval
                if hi < 0:
                    continue
                if lo == hi and lo < 0:
                    continue
                # No audited full- or half-step stage singularity lies on x >= 0.
                assert hi < 0


def test_closeout_table_is_regenerated_from_exact_method_data() -> None:
    rows = {row["key"]: row for row in csv.DictReader(TABLE_PATH.open(encoding="utf-8"))}
    methods = all_methods()
    assert set(rows) == set(methods)
    for key, method in methods.items():
        defect = first_exponential_defect(method.stability_function())
        coefficient = predicted_boundary_coefficient(defect)
        row = rows[key]
        assert int(row["order"]) == method.order
        assert int(row["stages"]) == method.stages
        assert int(row["first_defect_index_m"]) == defect.index
        assert sp.Rational(row["eta_m"]) == defect.coefficient
        assert sp.Rational(row["local_margin_coefficient"]) == coefficient
        assert row["local_orientation"].split("/")[0] == ("inward" if coefficient > 0 else "outward")


def test_closeout_full_two_rate_no_repair_for_actual_rk_maps() -> None:
    samples = {
        "euler": [(sp.Rational(1, 10), sp.Rational(0)), (sp.Rational(1, 10), sp.Rational(1, 10))],
        "rk2": [(sp.Rational(1), sp.Rational(0)), (sp.Rational(1), sp.Rational(1, 5))],
        "dp4": [(sp.Rational(1), sp.Rational(0)), (sp.Rational(4), sp.Rational(0))],
        "be": [(sp.Rational(10), sp.Rational(100))],
        "im": [(sp.Rational(1, 10), sp.Rational(0))],
    }
    for key, points in samples.items():
        R = all_methods()[key].stability_function()
        for xv, uv in points:
            a = sp.N(R.subs(X, -xv), 80)
            c = sp.N(R.subs(X, -xv / 2 - uv), 80)
            assert a >= 0
            one = cp_conditions(complex(a), complex(c))
            for N in (1, 2, 3, 8):
                repeated = cp_conditions(complex(a**N), complex(c**N))
                assert repeated == one

    assert not cp_conditions(complex(-0.5), complex(0.0))
    assert cp_conditions(complex(0.25), complex(0.0))


def test_closeout_logarithmic_fixed_horizon_coefficients_are_exact() -> None:
    t, final_x = sp.symbols("t final_x", positive=True, real=True)
    for method in all_methods().values():
        R = method.stability_function()
        defect = first_exponential_defect(R)
        m = defect.index
        eta = defect.coefficient
        full = sp.series((final_x / t) * sp.log(R.subs(X, -t)), t, 0, m + 1).removeO().expand()
        half = sp.series((2 * final_x / t) * sp.log(R.subs(X, -t / 2)), t, 0, m + 1).removeO().expand()
        assert sp.simplify(full.coeff(t, m - 1) - (-1) ** m * eta * final_x) == 0
        assert sp.simplify(half.coeff(t, m - 1) - (-1) ** m * 2 ** (1 - m) * eta * final_x) == 0


def test_closeout_fixed_horizon_prefactor_at_multiple_horizons() -> None:
    mp.mp.dps = 100
    for final_x in (mp.mpf("0.4"), mp.mpf("1.0"), mp.mpf("2.0")):
        N = 65536
        step = final_x / N
        for key in ("euler", "rk3", "dp4"):
            R = all_methods()[key].stability_function()
            defect = first_exponential_defect(R)
            K = mp.mpf(str(sp.N(predicted_boundary_coefficient(defect), 90)))
            Rf = sp.lambdify(X, R, "mpmath")
            A = mp.mpf(Rf(-step)) ** N
            C = mp.mpf(Rf(-step / 2)) ** N
            lam = _lower_choi_eigenvalue(A, C)
            observed = lam / step ** (defect.index - 1)
            expected = K * final_x * mp.e ** (-final_x) / (1 + mp.e ** (-final_x))
            assert mp.almosteq(observed, expected, rel_eps=mp.mpf("4e-4"))


def test_closeout_high_precision_near_origin_margin_signs_match_boundary_law() -> None:
    mp.mp.dps = 150
    tiny = mp.mpf("1e-20")
    for method in all_methods().values():
        R = method.stability_function()
        K = predicted_boundary_coefficient(first_exponential_defect(R))
        Rf = sp.lambdify(X, R, "mpmath")
        margin = mp.mpf(Rf(-tiny)) - mp.mpf(Rf(-tiny / 2)) ** 2
        assert mp.sign(margin) == (1 if K > 0 else -1)


def test_closeout_choi_matrix_and_scalar_criterion_agree_on_audit_samples() -> None:
    samples = [
        (0.8, 0.6),
        (0.8, 0.95),
        (1.2, 0.2),
        (0.0, 0.0),
        (1.0, 1.0),
    ]
    for a, c in samples:
        eigs = sp.Matrix(choi_matrix(a, c, normalized=False)).eigenvals()
        matrix_cptp = all(complex(ev).real >= -1e-12 for ev in eigs)
        scalar_cptp = cp_conditions(complex(a), complex(c))
        assert matrix_cptp == scalar_cptp


def test_closeout_fixed_horizon_figure_data_regenerate_exactly() -> None:
    mp.mp.dps = 80
    methods_by_name = {method.name: method for method in all_methods().values()}
    rows = list(csv.DictReader(FIGURE_DATA_PATH.open(encoding="utf-8")))
    assert rows
    for row in rows:
        method = methods_by_name[row["method"]]
        Rf = sp.lambdify(X, method.stability_function(), "mpmath")
        N = int(row["N"])
        step = mp.mpf(row["x_N"])
        A = mp.mpf(Rf(-step)) ** N
        C = mp.mpf(Rf(-step / 2)) ** N
        lam = _lower_choi_eigenvalue(A, C)
        stored = mp.mpf(row["lambda_minus_exact"])
        assert mp.almosteq(lam, stored, rel_eps=mp.mpf("1e-28"), abs_eps=mp.mpf("1e-45"))
        assert first_exponential_defect(method.stability_function()).index - 1 in (1, 3)


def test_benchmark_figure1_data_match_defining_curves_and_resolved_orientations() -> None:
    path = ROOT / "results" / "benchmark_figures" / "figure1_channel_geometry.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    panel_a = [row for row in rows if row["panel"] == "a" and row["admissible"] == "1"]
    assert panel_a
    for row in panel_a[::137]:
        th = float(row["series"].split("=")[1])
        aa = float(row["x_left"])
        yy = float(row["y_value"])
        expected = aa + th * (1 - th) * (1 - aa) ** 2
        assert math.isclose(yy, expected, rel_tol=2e-14, abs_tol=2e-14)
    panel_b = [row for row in rows if row["panel"] == "b"]
    assert len(panel_b) == 6
    thresholds = {(row["series"], row["orientation"]): float(row["x_left"])
                  for row in panel_b if row["kind"] == "threshold"}
    assert math.isclose(thresholds[("SSPRK(3,3)", "inward")], math.sqrt(7) / 2, rel_tol=1e-15)
    assert math.isclose(thresholds[("classical RK4", "outward")], math.sqrt(3) / 2, rel_tol=1e-15)
    intervals = {(row["series"], row["orientation"]): (float(row["x_left"]), float(row["x_right"]))
                 for row in panel_b if row["kind"] == "interval"}
    assert intervals[("SSPRK(3,3)", "outward")][0] == 0.0
    assert intervals[("SSPRK(3,3)", "inward")][1] == 2.0
    assert intervals[("classical RK4", "inward")][0] == 0.0
    assert intervals[("classical RK4", "outward")][1] == 2.0


def test_benchmark_figure2_components_match_certificate() -> None:
    path = ROOT / "results" / "benchmark_figures" / "figure2_admissible_sets.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    by_method = {}
    for row in rows:
        by_method.setdefault(row["method"], []).append((row["component_type"], row["left"], row["right"]))
    assert by_method["Forward Euler"] == [("point", "0.0", "0.0")]
    assert by_method["Heun RK2"] == [("interval", "0.0", "2.0")]
    assert by_method["Dormand--Prince 4 (embedded)"] == [("point", "0.0", "0.0"), ("isolated_interval", "3.068856548050381", "4.384986320801944")]
    assert by_method["Backward Euler"] == [("unbounded", "0.0", "infinity")]

