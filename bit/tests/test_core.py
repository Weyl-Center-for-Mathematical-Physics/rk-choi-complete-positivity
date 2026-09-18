from __future__ import annotations

import math

import mpmath as mp
import numpy as np
import sympy as sp

from rk_choi_margin.channels import (
    apply_phase_covariant_map,
    apply_thermal_phase_covariant_map,
    cp_conditions,
    minimum_choi_eigenvalue,
    thermal_choi_matrix,
    thermal_cp_conditions,
)
from rk_choi_margin.methods import X, all_methods, stability_function, stage_resolvent_determinant
from rk_choi_margin.symbolic import (
    U,
    Z,
    cp_admissible_intervals,
    cp_margin,
    first_exponential_defect,
    predicted_boundary_coefficient,
    rotating_boundary_factor,
    rotating_cp_margin,
    thermal_cp_margin,
    two_rate_cp_margin,
)


def test_minimal_stage_stability_polynomials() -> None:
    methods = all_methods()
    for key, degree in [("euler", 1), ("rk2", 2), ("rk3", 3), ("rk4", 4)]:
        expected = sum(X**k / sp.factorial(k) for k in range(degree + 1))
        assert sp.simplify(methods[key].stability_function() - expected) == 0




def test_reduced_stability_expression_does_not_extend_stage_domain() -> None:
    # A redundant implicit stage is ignored by b, so the reduced scalar
    # stability function cancels the pole.  The direct RK stage equations
    # remain singular at x=1 and therefore are outside the evaluation domain.
    A = ((sp.Integer(0), sp.Integer(0)), (sp.Integer(0), sp.Integer(1)))
    b = (sp.Integer(1), sp.Integer(0))
    R = stability_function(A, b)
    determinant = stage_resolvent_determinant(A)

    assert sp.simplify(R - (1 + X)) == 0
    assert R.subs(X, 1) == 2
    assert sp.factor(determinant) == 1 - X
    assert determinant.subs(X, 1) == 0


def test_exact_low_order_cp_margins() -> None:
    methods = all_methods()
    expected = {
        "euler": -Z**2 / 4,
        "rk2": -Z**3 * (Z - 8) / 64,
        "rk3": -Z**4 * (Z**2 - 12 * Z + 84) / 2304,
        "rk4": -Z**5 * (Z**3 - 16 * Z**2 + 160 * Z - 1152) / 147456,
    }
    for key, formula in expected.items():
        assert sp.simplify(cp_margin(methods[key].stability_function()) - formula) == 0


def test_general_boundary_coefficient_identity() -> None:
    for method in all_methods().values():
        R = method.stability_function()
        defect = first_exponential_defect(R)
        predicted = predicted_boundary_coefficient(defect)
        actual = sp.series(cp_margin(R), Z, 0, defect.index + 1).removeO().coeff(Z, defect.index)
        assert sp.simplify(actual - predicted) == 0


def test_choi_scalar_criterion() -> None:
    samples = [(0.8, 0.7), (0.2, 0.5), (1.0, 1.0), (0.5, -0.6)]
    for a, c in samples:
        scalar = cp_conditions(complex(a), complex(c))
        spectral = minimum_choi_eigenvalue(complex(a), complex(c)) >= -1e-12
        assert scalar == spectral


def test_repeated_step_factorization() -> None:
    a, c = sp.symbols("a c", positive=True)
    for count in range(1, 8):
        lhs = a**count - c ** (2 * count)
        rhs = (a - c**2) * sum(a ** (count - 1 - j) * c ** (2 * j) for j in range(count))
        assert sp.expand(lhs - rhs) == 0


def test_selected_cp_sets() -> None:
    methods = all_methods()
    expected = {
        "euler": [(0.0, 0.0)],
        "rk2": [(0.0, 2.0)],
        "rk3": [(0.0, 0.0)],
        "be": [(0.0, math.inf)],
        "im": [(0.0, 0.0)],
    }
    for key, target in expected.items():
        got = cp_admissible_intervals(methods[key].stability_function())
        assert len(got) == len(target)
        for (left_g, right_g), (left_t, right_t) in zip(got, target):
            assert abs(left_g - left_t) < 1e-9
            if math.isinf(right_t):
                assert math.isinf(right_g)
            else:
                assert abs(right_g - right_t) < 1e-8


def test_dp4_has_disconnected_cp_admissible_set() -> None:
    intervals = cp_admissible_intervals(all_methods()["dp4"].stability_function())
    assert intervals[0] == (0.0, 0.0)
    assert len(intervals) >= 2
    assert 3.0 < intervals[1][0] < 3.2
    assert 4.3 < intervals[1][1] < 4.5


def test_rotating_rk4_leading_coefficient() -> None:
    nu = sp.Symbol("nu", real=True)
    R = all_methods()["rk4"].stability_function()
    M = rotating_cp_margin(R, nu=nu)
    coefficient = sp.factor(sp.expand(M).coeff(Z, 5))
    expected = -((4 * nu**2 - 3) * (4 * nu**2 + 1)) / 384
    assert sp.simplify(coefficient - expected) == 0



def test_rk2_population_exit_is_independent_of_cp_block_margin() -> None:
    R = all_methods()["rk2"].stability_function()
    a = sp.simplify(R.subs(X, -Z))
    c = sp.simplify(R.subs(X, -Z / 2))
    margin = sp.factor(a - c**2)

    assert sp.simplify(margin - Z**3 * (8 - Z) / 64) == 0
    assert a.subs(Z, 4) == 5
    assert c.subs(Z, 4) == 1
    assert margin.subs(Z, 4) == 4
    assert 1 - a.subs(Z, 4) == -4

    # The CP block remains positive while the independent population
    # eigenvalue 1-a has already become negative.
    block = sp.Matrix([[1, c.subs(Z, 4)], [c.subs(Z, 4), a.subs(Z, 4)]])
    assert block.det() == 4
    assert all(value > 0 for value in block.eigenvals())


def test_two_rate_coordinates_recover_fixed_ratio_margin() -> None:
    R = all_methods()["rk4"].stability_function()
    r = sp.Symbol("r", real=True)
    two_rate = two_rate_cp_margin(R, relaxation_step=Z, dephasing_step=(r - sp.Rational(1, 2)) * Z)
    one_rate = cp_margin(R, r=r)
    assert sp.simplify(two_rate - one_rate) == 0


def test_two_rate_joint_linearization_and_pure_dephasing_limit() -> None:
    b2, b3 = sp.symbols("b2 b3", real=True)
    R = 1 + X + b2 * X**2 + b3 * X**3
    margin = sp.expand(two_rate_cp_margin(R))

    assert margin.subs({Z: 0, U: 0}) == 0
    assert sp.diff(margin, Z).subs({Z: 0, U: 0}) == 0
    assert sp.diff(margin, U).subs({Z: 0, U: 0}) == 2
    assert sp.simplify(margin.subs(Z, 0) - (1 - R.subs(X, -U) ** 2)) == 0


def test_choi_rationalized_denominator_has_global_lower_bound() -> None:
    a, c = sp.symbols("a c", real=True)
    delta = (1 - a) ** 2 + 4 * c**2
    denominator = 1 + a + sp.sqrt(delta)

    # The analytic proof uses sqrt(delta) >= |1-a|, yielding
    # D >= 1+a+|1-a| = 2*max(1,a) >= 2.  Sample both branches and
    # the degenerate cusp point to protect the implementation-level claim.
    samples = [(-7, 0), (-2, 3), (sp.Rational(1, 2), 0), (1, 0), (4, 0), (4, 2)]
    for a_value, c_value in samples:
        assert sp.N(denominator.subs({a: a_value, c: c_value}), 50) >= 2

    assert delta.subs({a: 1, c: 0}) == 0
    assert denominator.subs({a: 1, c: 0}) == 2


def test_backward_euler_is_cptp_on_full_two_rate_quadrant() -> None:
    R = all_methods()["be"].stability_function()
    x, u = sp.symbols("x u", nonnegative=True, real=True)
    a = sp.factor(R.subs(X, -x))
    c = sp.factor(R.subs(X, -x / 2 - u))
    margin = sp.factor(a - c**2)

    expected_margin = sp.factor(
        (x**2 / 4 + x * u + u**2 + 2 * u)
        / ((1 + x) * (1 + x / 2 + u) ** 2)
    )
    assert sp.simplify(margin - expected_margin) == 0
    assert sp.simplify(1 - a - x / (1 + x)) == 0

    samples = [(0, 0), (1, 0), (0, 5), (2, 3), (100, 1000)]
    for x_value, u_value in samples:
        a_value = float(a.subs({x: x_value, u: u_value}))
        c_value = float(c.subs({x: x_value, u: u_value}))
        assert cp_conditions(complex(a_value), complex(c_value))


def test_nonconstant_explicit_rk_polynomials_have_stiffness_obstruction() -> None:
    methods = all_methods()
    for key in ("euler", "rk2", "rk3", "rk4", "dp5", "dp4"):
        R = methods[key].stability_function()
        margin = two_rate_cp_margin(R)
        assert sp.limit(margin, U, sp.oo) == -sp.oo


def test_unnormalized_choi_convention_and_normalized_state() -> None:
    from rk_choi_margin.channels import choi_matrix

    a, c = 0.7, 0.4
    J = choi_matrix(a, c, normalized=False)
    rho_J = choi_matrix(a, c, normalized=True)
    assert np.allclose(rho_J, J / 2.0)
    assert abs(np.trace(J) - 2.0) < 1e-12
    assert abs(np.trace(rho_J) - 1.0) < 1e-12


def test_numeric_exact_numerical_map_preserves_trace_and_hermiticity() -> None:
    X = np.array(
        [[0.35, 0.12 + 0.08j], [0.12 - 0.08j, 0.65]],
        dtype=np.complex128,
    )
    a = 0.73
    c = 0.61
    Y = apply_phase_covariant_map(X, a, c)

    assert abs(np.trace(Y) - np.trace(X)) < 1e-14
    assert np.allclose(Y, Y.conjugate().T)
    assert np.allclose(
        Y,
        np.array(
            [
                [X[0, 0] + (1 - a) * X[1, 1], c * X[0, 1]],
                [c * X[1, 0], a * X[1, 1]],
            ],
            dtype=np.complex128,
        ),
    )


def test_exact_numerical_map_rejects_wrong_matrix_shape() -> None:
    with np.testing.assert_raises(ValueError):
        apply_phase_covariant_map(np.zeros((4,), dtype=np.complex128), 1.0, 1.0)


def test_lindblad_action_on_operator_basis_matches_section2() -> None:
    gamma, gamma_phi = sp.symbols("gamma gamma_phi", nonnegative=True, real=True)
    zero = sp.Integer(0)
    one = sp.Integer(1)

    e00 = sp.Matrix([[one, zero], [zero, zero]])
    e01 = sp.Matrix([[zero, one], [zero, zero]])
    e10 = sp.Matrix([[zero, zero], [one, zero]])
    e11 = sp.Matrix([[zero, zero], [zero, one]])
    sigma_minus = e01
    sigma_plus = e10
    sigma_z = sp.diag(one, -one)

    def lindblad(matrix: sp.Matrix) -> sp.Matrix:
        number = sigma_plus * sigma_minus
        relaxation = gamma * (
            sigma_minus * matrix * sigma_plus
            - sp.Rational(1, 2) * (number * matrix + matrix * number)
        )
        dephasing = gamma_phi / 2 * (sigma_z * matrix * sigma_z - matrix)
        return sp.simplify(relaxation + dephasing)

    assert lindblad(e00) == sp.zeros(2)
    assert sp.simplify(lindblad(e11) - gamma * (e00 - e11)) == sp.zeros(2)
    decay = gamma / 2 + gamma_phi
    assert sp.simplify(lindblad(e01) + decay * e01) == sp.zeros(2)
    assert sp.simplify(lindblad(e10) + decay * e10) == sp.zeros(2)


def test_exact_numerical_map_preserves_trace_and_hermiticity() -> None:
    a, c = sp.symbols("a c", real=True)
    x00, x01, x10, x11 = sp.symbols("x00 x01 x10 x11")

    def numerical_map(matrix: sp.Matrix) -> sp.Matrix:
        return sp.Matrix(
            [
                [matrix[0, 0] + (1 - a) * matrix[1, 1], c * matrix[0, 1]],
                [c * matrix[1, 0], a * matrix[1, 1]],
            ]
        )

    X = sp.Matrix([[x00, x01], [x10, x11]])
    assert sp.simplify(sp.trace(numerical_map(X)) - sp.trace(X)) == 0
    assert sp.simplify(numerical_map(X.conjugate().T) - numerical_map(X).conjugate().T) == sp.zeros(2)


def test_section4_rk4_population_boundary_precedes_margin_zero() -> None:
    methods = all_methods()
    R = methods["rk4"].stability_function()
    a = sp.factor(R.subs(X, -Z))
    M = sp.factor(cp_margin(R))

    population_poly = Z**3 - 4 * Z**2 + 12 * Z - 24
    margin_poly = Z**3 - 16 * Z**2 + 160 * Z - 1152
    assert sp.Poly(population_poly, Z).count_roots(sp.Rational("2.7852935634"), sp.Rational("2.7852935635")) == 1
    assert sp.Poly(margin_poly, Z).count_roots(sp.Rational("10.9824254662"), sp.Rational("10.9824254664")) == 1
    assert sp.sign((a - 1).subs(Z, sp.Rational(2))) < 0
    assert sp.sign((a - 1).subs(Z, sp.Rational(3))) > 0
    assert sp.sign(M.subs(Z, sp.Rational(3))) > 0


def test_section4_dp5_population_boundary_precedes_margin_zero() -> None:
    R = all_methods()["dp5"].stability_function()
    a = sp.factor(R.subs(X, -Z))
    M = sp.factor(cp_margin(R))
    population_poly = sp.factor((a - 1) * 600 / Z)
    margin_poly = sp.factor(sp.together(M).as_numer_denom()[0] / Z**6)

    assert sp.Poly(population_poly, Z).count_roots(sp.Rational("3.3065678926"), sp.Rational("3.3065678927")) == 1
    assert sp.Poly(margin_poly, Z).count_roots(sp.Rational("13.1740934628"), sp.Rational("13.1740934629")) == 1
    assert sp.sign((a - 1).subs(Z, sp.Rational(3))) < 0
    assert sp.sign((a - 1).subs(Z, sp.Rational(4))) > 0
    assert sp.sign(M.subs(Z, sp.Rational(4))) > 0


def test_section4_dp4_exact_disconnected_component_boundaries() -> None:
    R = all_methods()["dp4"].stability_function()
    a = sp.factor(R.subs(X, -Z))
    M = sp.factor(cp_margin(R))
    population_poly = sp.factor(-(a - 1) * 120000 / Z)
    margin_poly = sp.factor(-sp.together(M).as_numer_denom()[0] / Z**5)

    assert sp.Poly(margin_poly, Z).count_roots(sp.Rational("3.0688565480"), sp.Rational("3.0688565481")) == 1
    assert sp.Poly(population_poly, Z).count_roots(sp.Rational("4.3849863208"), sp.Rational("4.3849863209")) == 1
    assert sp.sign(M.subs(Z, sp.Rational(1))) < 0
    assert 0 < a.subs(Z, sp.Rational(4)) < 1
    assert sp.sign(M.subs(Z, sp.Rational(4))) > 0
    assert a.subs(Z, sp.Rational(5)) > 1
    assert sp.sign(M.subs(Z, sp.Rational(5))) > 0


def test_section4_implicit_methods_reverse_truncated_exponential_parity() -> None:
    methods = all_methods()
    be_R = methods["be"].stability_function()
    im_R = methods["im"].stability_function()
    be_defect = first_exponential_defect(be_R)
    im_defect = first_exponential_defect(im_R)

    assert methods["be"].order == 1
    assert predicted_boundary_coefficient(be_defect) > 0
    assert cp_admissible_intervals(be_R) == [(0.0, math.inf)]

    assert methods["im"].order == 2
    assert predicted_boundary_coefficient(im_defect) < 0
    assert cp_admissible_intervals(im_R) == [(0.0, 0.0)]


def test_section4_all_reported_admissible_sets() -> None:
    methods = all_methods()
    expected = {
        "euler": [(0.0, 0.0)],
        "rk2": [(0.0, 2.0)],
        "rk3": [(0.0, 0.0)],
        "rk4": [(0.0, 2.7852935634053)],
        "dp5": [(0.0, 3.3065678926349)],
        "dp4": [(0.0, 0.0), (3.0688565480504, 4.3849863208019)],
        "be": [(0.0, math.inf)],
        "im": [(0.0, 0.0)],
    }
    for key, target in expected.items():
        got = cp_admissible_intervals(methods[key].stability_function())
        assert len(got) == len(target)
        for (left_g, right_g), (left_t, right_t) in zip(got, target):
            assert abs(left_g - left_t) < 2e-10
            if math.isinf(right_t):
                assert math.isinf(right_g)
            else:
                assert abs(right_g - right_t) < 2e-10



def test_section5_no_repair_equivalence_on_nonnegative_population_branch() -> None:
    samples = [
        (0.0, 0.0),
        (0.2, 0.3),
        (0.2, 0.6),
        (0.8, 0.7),
        (1.0, 1.0),
        (1.2, 0.5),
    ]
    for a, c in samples:
        one_step = cp_conditions(complex(a), complex(c))
        for count in (1, 2, 3, 5, 8):
            repeated = cp_conditions(complex(a**count), complex(c**count))
            assert repeated == one_step

    # The a >= 0 scope is necessary: an even composition can repair a
    # negative population multiplier outside the local convergence branch.
    assert not cp_conditions(complex(-0.5), complex(0.0))
    assert cp_conditions(complex((-0.5) ** 2), complex(0.0))


def test_section5_fixed_horizon_margin_and_eigenvalue_asymptotics() -> None:
    mp.mp.dps = 100
    final_x = mp.mpf("1.0")
    count = 32768
    step = final_x / count

    for key in ("euler", "rk3", "dp4"):
        R = all_methods()[key].stability_function()
        defect = first_exponential_defect(R)
        K = mp.mpf(str(sp.N(predicted_boundary_coefficient(defect), 90)))
        Rf = sp.lambdify(X, R, "mpmath")

        a_final = mp.mpf(Rf(-step)) ** count
        c_final = mp.mpf(Rf(-step / 2)) ** count
        margin = a_final - c_final**2
        denominator = 1 + a_final + mp.sqrt((1 - a_final) ** 2 + 4 * c_final**2)
        lower_eigenvalue = 2 * margin / denominator

        expected_margin_coefficient = K * final_x * mp.e ** (-final_x)
        observed_margin_coefficient = margin / step ** (defect.index - 1)
        assert mp.almosteq(
            observed_margin_coefficient,
            expected_margin_coefficient,
            rel_eps=mp.mpf("2e-4"),
        )

        expected_eigen_coefficient = expected_margin_coefficient / (1 + mp.e ** (-final_x))
        observed_eigen_coefficient = lower_eigenvalue / step ** (defect.index - 1)
        assert mp.almosteq(
            observed_eigen_coefficient,
            expected_eigen_coefficient,
            rel_eps=mp.mpf("2e-4"),
        )


def test_section5_absolute_stability_of_modes_does_not_imply_cptp() -> None:
    R = all_methods()["euler"].stability_function()
    a = sp.simplify(R.subs(X, -1))
    c = sp.simplify(R.subs(X, -sp.Rational(1, 2)))
    margin = sp.simplify(a - c**2)

    assert abs(a) <= 1
    assert abs(c) <= 1
    assert margin == -sp.Rational(1, 4)
    assert not cp_conditions(complex(float(a)), complex(float(c)))


def test_section5_ssp_and_generator_conditioned_cptp_are_independent() -> None:
    methods = all_methods()
    # SSPRK(3,3) has a positive SSP coefficient in the standard theory but
    # no punctured CPTP neighborhood for this generator.
    assert cp_admissible_intervals(methods["rk3"].stability_function()) == [(0.0, 0.0)]

    # Classical RK4 has zero nonlinear SSP coefficient in the standard
    # theory but a nonzero generator-conditioned CPTP interval here.
    rk4_intervals = cp_admissible_intervals(methods["rk4"].stability_function())
    assert rk4_intervals[0][0] == 0.0
    assert 2.78 < rk4_intervals[0][1] < 2.79


def test_section5_cptp_does_not_imply_accuracy_dp4_at_x4() -> None:
    R = all_methods()["dp4"].stability_function()
    a = sp.factor(R.subs(X, -4))
    c = sp.factor(R.subs(X, -2))
    margin = sp.factor(a - c**2)

    assert a == sp.Rational(847, 1875)
    assert c == sp.Rational(91, 750)
    assert margin == sp.Rational(245819, 562500)
    assert cp_conditions(complex(float(a)), complex(float(c)))
    assert abs(float(a) - math.exp(-4)) > 0.43



def test_thermal_map_and_choi_matrix_are_consistent() -> None:
    a = 0.63
    c = 0.44 + 0.17j
    theta = 0.31
    units = [
        np.array([[1, 0], [0, 0]], dtype=np.complex128),
        np.array([[0, 1], [0, 0]], dtype=np.complex128),
        np.array([[0, 0], [1, 0]], dtype=np.complex128),
        np.array([[0, 0], [0, 1]], dtype=np.complex128),
    ]
    mapped = [apply_thermal_phase_covariant_map(E, a, c, theta) for E in units]
    J = np.block([[mapped[0], mapped[1]], [mapped[2], mapped[3]]])
    assert np.allclose(J, thermal_choi_matrix(a, c, theta))
    X = np.array([[0.4, 0.1 + 0.07j], [0.1 - 0.07j, 0.6]], dtype=np.complex128)
    Y = apply_thermal_phase_covariant_map(X, a, c, theta)
    assert np.allclose(Y, Y.conjugate().T)
    assert abs(np.trace(Y) - np.trace(X)) < 1e-14


def test_thermal_scalar_criterion_agrees_with_choi_spectrum() -> None:
    rng = np.random.default_rng(20260809)
    for _ in range(250):
        a = rng.uniform(-1.2, 1.4)
        c = rng.normal() + 1j * rng.normal()
        theta = rng.uniform(0.0, 1.0)
        scalar = thermal_cp_conditions(complex(a), complex(c), float(theta))
        spectral = np.linalg.eigvalsh(thermal_choi_matrix(a, c, float(theta)))[0] >= -1e-12
        assert scalar == spectral


def test_symbolic_thermal_margin_reduces_to_zero_temperature_margin() -> None:
    R = all_methods()["rk4"].stability_function()
    theta, nu = sp.symbols("theta nu", real=True)
    F = thermal_cp_margin(R, theta=theta, frequency_step=nu)
    zero = sp.factor(F.subs({theta: 0, nu: 0}))
    assert sp.simplify(zero - two_rate_cp_margin(R)) == 0


def test_exact_phase_covariant_margin_identity() -> None:
    x, u, nu, theta = sp.symbols("x u nu theta", real=True)
    a = sp.exp(-x)
    c2 = sp.exp(-x - 2 * u)
    F = sp.factor(a + theta * (1 - theta) * (1 - a) ** 2 - c2)
    expected = sp.exp(-x) * (1 - sp.exp(-2 * u)) + theta * (1 - theta) * (1 - sp.exp(-x)) ** 2
    assert sp.simplify(F - expected) == 0
    assert nu not in F.free_symbols


def test_local_dephasing_and_bidirectional_buffers() -> None:
    r2 = sp.symbols("r2", real=True)
    theta = sp.symbols("theta", real=True)
    R = 1 + X + r2 * X**2
    F = sp.expand(thermal_cp_margin(R, theta=theta, frequency_step=0))
    assert F.subs({Z: 0, U: 0}) == 0
    assert sp.diff(F, U).subs({Z: 0, U: 0}) == 2
    coeff = sp.expand(F.subs(U, 0)).coeff(Z, 2)
    assert sp.simplify(coeff - (theta * (1 - theta) + r2 / 2 - sp.Rational(1, 4))) == 0
    assert sp.simplify(coeff.subs(r2, sp.Rational(1, 2)) - theta * (1 - theta)) == 0


def test_bidirectional_buffer_is_frequency_independent_at_order_two() -> None:
    r2, theta, varpi = sp.symbols("r2 theta varpi", real=True)
    R = 1 + X + r2 * X**2
    F = sp.expand(thermal_cp_margin(R, theta=theta, dephasing_step=0, frequency_step=varpi * Z))
    coeff = sp.expand(F).coeff(Z, 2)
    expected = theta * (1 - theta) + r2 / 2 - sp.Rational(1, 4) + (2 * r2 - 1) * varpi**2
    assert sp.simplify(coeff - expected) == 0
    assert sp.simplify(coeff.subs(r2, sp.Rational(1, 2)) - theta * (1 - theta)) == 0


def test_forward_euler_symmetric_bidirectional_interval() -> None:
    R = all_methods()["euler"].stability_function()
    F = sp.factor(thermal_cp_margin(R, theta=sp.Rational(1, 2), dephasing_step=0, frequency_step=0))
    assert F == 0
    a = sp.factor(R.subs(X, -Z))
    A = sp.factor(1 - sp.Rational(1, 2) * (1 - a))
    D = sp.factor(sp.Rational(1, 2) + sp.Rational(1, 2) * a)
    assert sp.simplify(A - (1 - Z / 2)) == 0
    assert sp.simplify(D - (1 - Z / 2)) == 0


def test_frequency_dependent_boundary_law_for_all_audited_methods() -> None:
    varpi = sp.symbols("varpi", real=True)
    for method in all_methods().values():
        R = method.stability_function()
        defect = first_exponential_defect(R)
        M = rotating_cp_margin(R, nu=varpi)
        actual = sp.series(M, Z, 0, defect.index + 1).removeO().expand().coeff(Z, defect.index)
        predicted = sp.factor((-1) ** defect.index * defect.coefficient * rotating_boundary_factor(defect.index, varpi))
        assert sp.simplify(actual - predicted) == 0


def test_frequency_flip_thresholds_for_ssprk3_and_rk4() -> None:
    v = sp.symbols("v", real=True)
    assert sp.factor(rotating_boundary_factor(4, v) - (sp.Rational(7, 8) + 3 * v**2 - 2 * v**4)) == 0
    assert sp.factor(rotating_boundary_factor(5, v) - (sp.Rational(15, 16) + sp.Rational(5, 2) * v**2 - 5 * v**4)) == 0
    assert sp.simplify(rotating_boundary_factor(4, sp.sqrt(7) / 2)) == 0
    assert sp.simplify(rotating_boundary_factor(5, sp.sqrt(3) / 2)) == 0
    assert sp.simplify(rotating_boundary_factor(2, v) - (sp.Rational(1, 2) + 2 * v**2)) == 0
    assert sp.simplify(rotating_boundary_factor(3, v) - (sp.Rational(3, 4) + 3 * v**2)) == 0


def test_backward_euler_is_cptp_with_arbitrary_commuting_precession() -> None:
    x, u, nu = sp.symbols("x u nu", nonnegative=True, real=True)
    a = 1 / (1 + x)
    c2 = 1 / ((1 + x / 2 + u) ** 2 + nu**2)
    assert sp.simplify(a - c2 - (
        ((1 + x / 2 + u) ** 2 + nu**2 - (1 + x)) /
        ((1 + x) * ((1 + x / 2 + u) ** 2 + nu**2))
    )) == 0
    numerator = sp.expand((1 + x / 2 + u) ** 2 + nu**2 - (1 + x))
    assert numerator == x**2 / 4 + x * u + u**2 + 2 * u + nu**2


def test_frequency_reversal_thresholds_are_resolved_beyond_leading_order() -> None:
    z = sp.Symbol("z", nonnegative=True, real=True)
    methods = all_methods()
    ssp = rotating_cp_margin(methods["rk3"].stability_function(), nu=sp.sqrt(7) / 2, z=z)
    rk4 = rotating_cp_margin(methods["rk4"].stability_function(), nu=sp.sqrt(3) / 2, z=z)
    assert sp.simplify(ssp - z**5 * (3 - 2 * z) / 9) == 0
    assert sp.simplify(rk4 + z**6 * (z - 2) ** 2 / 576) == 0
    assert sp.simplify(sp.limit(ssp / z**5, z, 0, dir="+")) == sp.Rational(1, 3)
    assert sp.simplify(sp.limit(rk4 / z**6, z, 0, dir="+")) == -sp.Rational(1, 144)
