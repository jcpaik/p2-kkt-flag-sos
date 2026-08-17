from fractions import Fraction

import numpy as np

import sos_search as search


def assert_vector_close(actual, expected, tolerance=1e-10):
    assert set(actual) == set(expected)
    for label, value in expected.items():
        assert abs(actual[label] - value) <= tolerance


def test_triangle_reducer_matches_graph_reducer():
    for i in range(7):
        for j in range(7):
            for k in range(7):
                triangle = search.expectation_label((i, j, k))
                graph = search.graph_expectation_label(3, (i, k, j))
                assert triangle == graph


def test_second_moment_is_retained():
    label, coefficient = search.expectation_label((0, 2, 0))
    assert label == ("pair", 2)
    assert coefficient == 1
    label, coefficient = search.graph_expectation_label(2, (2,))
    assert label == ("pair", 2)
    assert coefficient == 1


def test_general_three_flag_identities():
    s2 = [
        (Fraction(1), (2, 0, 0)),
        (Fraction(1), (0, 2, 0)),
        (Fraction(1), (0, 0, 2)),
    ]
    s4 = [
        (Fraction(1), (4, 0, 0)),
        (Fraction(1), (0, 4, 0)),
        (Fraction(1), (0, 0, 4)),
    ]
    first = search.expectation_vector(
        search.multiply_polynomials(search.GRAM_DETERMINANT, s2)
    )
    assert_vector_close(
        first,
        {
            ("pair", 2): 3,
            ("triangle", 1, 1, 3): 6,
            ("pair", 4): -3,
            ("triangle", 0, 2, 2): -6,
        },
    )
    second = search.expectation_vector(
        search.multiply_polynomials(search.GRAM_DETERMINANT, s4)
    )
    assert_vector_close(
        second,
        {
            ("pair", 4): 3,
            ("triangle", 1, 1, 5): 6,
            ("pair", 6): -3,
            ("triangle", 0, 2, 4): -6,
        },
    )
    # On isotropic reference measures the general vectors must agree with
    # the classical isotropic contractions 1/3 - p4 and p4 - p6.
    for evaluate, p4, p6 in (
        (search.onb_label_value, Fraction(1, 3), Fraction(1, 3)),
        (
            lambda label: search.pole_equator_label_value(label, 0),
            Fraction(5, 18),
            Fraction(1, 4),
        ),
    ):
        first_value = sum(
            Fraction(value).limit_denominator(10**6) * evaluate(label)
            for label, value in first.items()
        )
        assert first_value == Fraction(1, 3) - p4
        second_value = sum(
            Fraction(value).limit_denominator(10**6) * evaluate(label)
            for label, value in second.items()
        )
        assert second_value == p4 - p6


def test_unweighted_hessian_identity():
    # The isotropic contraction of both averaged Hessian kernels is
    # -64 + 640 p4 - 448 p6.  The general (non-isotropic) label vectors must
    # reproduce those values on the two exact isotropic reference measures.
    _, hessian, perpendicular_hessian = search.kernel_polynomials()
    for polynomial in (hessian, perpendicular_hessian):
        vector = search.expectation_vector(polynomial)
        for evaluate, p4, p6 in (
            (search.onb_label_value, Fraction(1, 3), Fraction(1, 3)),
            (
                lambda label: search.pole_equator_label_value(label, 0),
                Fraction(5, 18),
                Fraction(1, 4),
            ),
        ):
            value = sum(
                Fraction(coefficient).limit_denominator(10**6)
                * evaluate(label)
                for label, coefficient in vector.items()
            )
            assert value == -64 + 640 * p4 - 448 * p6


def test_orientation_pairing_polynomial():
    rng = np.random.default_rng(20260725)
    for _ in range(20):
        x, z, y, w = rng.normal(size=(4, 3))
        x /= np.linalg.norm(x)
        z /= np.linalg.norm(z)
        y /= np.linalg.norm(y)
        w /= np.linalg.norm(w)
        edges = np.array(
            [
                x @ z,
                x @ y,
                x @ w,
                z @ y,
                z @ w,
                y @ w,
            ]
        )
        polynomial_value = sum(
            float(coefficient)
            * np.prod(edges ** np.array(exponent))
            for coefficient, exponent in search.ORIENTATION_PAIRING
        )
        determinant_value = np.linalg.det(
            np.stack([x, z, y])
        ) * np.linalg.det(np.stack([x, z, w]))
        assert abs(polynomial_value - determinant_value) <= 1e-12


def test_rank_relations_vanish_on_onb():
    relations = search.four_point_rank_relations(8)
    assert relations
    for relation in relations:
        value = sum(
            coefficient * float(search.onb_label_value(label))
            for label, coefficient in relation.items()
        )
        assert abs(value) <= 1e-10


def test_general_target_vanishes_on_onb():
    target = (
        Fraction(-4, 3)
        + 20 * search.onb_label_value(("pair", 2))
        - 48 * search.onb_label_value(("pair", 4))
        + 32 * search.onb_label_value(("pair", 6))
    )
    assert target == 0


def test_pole_equator_equality_moments_are_exact():
    for regular_order in (0, 4, 5):
        p2 = search.pole_equator_label_value(
            ("pair", 2),
            regular_order,
        )
        p4 = search.pole_equator_label_value(
            ("pair", 4),
            regular_order,
        )
        p6 = search.pole_equator_label_value(
            ("pair", 6),
            regular_order,
        )
        assert p2 == Fraction(1, 3)
        assert p4 == Fraction(5, 18)
        assert p6 == Fraction(1, 4)
        assert (
            Fraction(-1, 4) + Fraction(15, 4) * p2 - 9 * p4 + 6 * p6 == 0
        )


def test_lifted_hessian_reproduces_four_point_block():
    parallel, perpendicular = search.four_point_hessian_polynomials()
    basis = [(1,), (3,)]
    for polynomial in (parallel, perpendicular):
        lifted = search.lifted_hessian_expectation_matrix(
            1,
            basis,
            polynomial,
        )
        original = search.four_point_hessian_expectation_matrix(
            [1, 3],
            polynomial,
        )
        original = {
            label: matrix
            for label, matrix in original.items()
            if np.max(np.abs(matrix)) > 1e-12
        }
        assert set(lifted) == set(original)
        for label in lifted:
            assert np.max(np.abs(lifted[label] - original[label])) <= 1e-12


def test_five_point_star_flags_are_psd_on_onb():
    sectors = search.rooted_star_flag_sectors(
        root_count=3,
        maximum_leaf_degree=6,
        leaf_parity=0,
    )
    assert sectors
    for basis in sectors.values():
        label_matrices = search.rooted_star_flag_expectation_matrix(3, basis)
        moment_matrix = sum(
            float(search.onb_label_value(label)) * matrix
            for label, matrix in label_matrices.items()
        )
        assert np.linalg.eigvalsh(moment_matrix)[0] >= -1e-12


def test_five_point_root_weighted_flags_are_psd_on_onb():
    sectors = search.rooted_weighted_flag_sectors(
        root_count=3,
        maximum_total_degree=5,
        maximum_root_degree=2,
    )
    assert sectors
    for basis in sectors.values():
        label_matrices = search.rooted_weighted_flag_expectation_matrix(
            3,
            basis,
        )
        moment_matrix = sum(
            float(search.onb_label_value(label)) * matrix
            for label, matrix in label_matrices.items()
        )
        assert np.linalg.eigvalsh(moment_matrix)[0] >= -1e-12
        exact_moment_matrix = search.exact_onb_moment_matrix(label_matrices)
        exact_nullspace = search.exact_onb_nullspace(label_matrices)
        assert (
            exact_moment_matrix * exact_nullspace
            == search.sp.zeros(
                exact_moment_matrix.rows,
                exact_nullspace.cols,
            )
        )


def test_five_point_hessian_blocks_vanish_on_onb():
    parallel, perpendicular = search.four_point_hessian_polynomials()
    sectors = search.rooted_star_flag_sectors(
        root_count=2,
        maximum_leaf_degree=3,
        leaf_parity=1,
    )
    assert sectors
    for basis in sectors.values():
        for polynomial in (parallel, perpendicular):
            label_matrices = search.lifted_hessian_expectation_matrix(
                2,
                basis,
                polynomial,
            )
            moment_matrix = sum(
                float(search.onb_label_value(label)) * matrix
                for label, matrix in label_matrices.items()
            )
            assert np.max(np.abs(moment_matrix)) <= 1e-10


def test_harmonic_flag_blocks_are_psd_on_reference_measures():
    for order in (2, 4, 6):
        label_matrices = search.harmonic_flag_expectation_matrix(
            order,
            [0, 2, 4],
        )
        assert label_matrices
        for evaluate in (
            search.onb_label_value,
            lambda label: search.pole_equator_label_value(label, 0),
        ):
            moment_matrix = sum(
                float(evaluate(label)) * matrix
                for label, matrix in label_matrices.items()
            )
            assert np.linalg.eigvalsh(moment_matrix)[0] >= -1e-12


def test_harmonic_flag_corner_matches_scalar_harmonic():
    for order in (2, 4, 6):
        label_matrices = search.harmonic_flag_expectation_matrix(
            order,
            [0, 2],
        )
        scalar = search.harmonic_pair_vector(order)
        for label, matrix in scalar.items():
            assert abs(label_matrices[label][0, 0] - matrix[0, 0]) <= 1e-12
        for label, matrix in label_matrices.items():
            if label not in scalar:
                assert abs(matrix[0, 0]) <= 1e-12


def test_spin2_flag_block_is_psd_on_reference_measures():
    basis = search.spin2_flag_basis(4)
    assert ((), 0) in basis
    assert ((2,), 0) in basis
    assert ((1, 1), 1) in basis
    label_matrices = search.spin2_flag_expectation_matrix(basis)
    assert label_matrices
    for evaluate in (
        search.onb_label_value,
        lambda label: search.pole_equator_label_value(label, 0),
    ):
        moment_matrix = sum(
            float(evaluate(label)) * matrix
            for label, matrix in label_matrices.items()
        )
        assert np.linalg.eigvalsh(moment_matrix)[0] >= -1e-12


def test_spin2_empty_flag_squared_is_harmonic_energy():
    # ||D||^2-entry: E[P_2(X.Z)] = h_2 = (3 p2 - 1)/2.
    label_matrices = search.spin2_flag_expectation_matrix([((), 0)])
    assert_vector_close(
        {label: matrix[0, 0] for label, matrix in label_matrices.items()},
        {
            ("pair", 2): 1.5,
            ("constant",): -0.5,
        },
    )


def test_h2_localized_flag_block_vanishes_on_isotropic_measures():
    harmonics = search.tangent_harmonic_polynomials(2)
    for tangent_harmonic in harmonics:
        label_matrices = search.h2_localized_flag_expectation_matrix(
            [0, 2],
            tangent_harmonic,
        )
        # Fresh vertices are independent, so the block factors as
        # h_2 times the plain flag block; both reference measures are
        # isotropic (h_2 = 0), hence the block must vanish there.
        for evaluate in (
            search.onb_label_value,
            lambda label: search.pole_equator_label_value(label, 0),
        ):
            moment_matrix = sum(
                float(evaluate(label)) * matrix
                for label, matrix in label_matrices.items()
            )
            if isinstance(moment_matrix, np.ndarray):
                assert np.max(np.abs(moment_matrix)) <= 1e-12


def test_exact_equality_quotient_annihilates_free_generators():
    constant = ("constant",)
    first = ("pair", 4)
    second = ("pair", 6)
    quotient, generator_count, rank = search.exact_equality_quotient_rows(
        [constant, first, second],
        [{constant: np.array([[1.0]])}],
        [{first: 1.0, second: 1.0}],
    )
    assert generator_count == 2
    assert rank == 2
    assert quotient == [
        {
            first: Fraction(-1),
            second: Fraction(1),
        }
    ]
