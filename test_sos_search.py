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


def test_isotropic_three_flag_identities():
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
    assert_vector_close(
        search.expectation_vector(
            search.multiply_polynomials(search.GRAM_DETERMINANT, s2)
        ),
        {
            ("constant",): 1 / 3,
            ("pair", 4): -1,
        },
    )
    assert_vector_close(
        search.expectation_vector(
            search.multiply_polynomials(search.GRAM_DETERMINANT, s4)
        ),
        {
            ("pair", 4): 1,
            ("pair", 6): -1,
        },
    )


def test_unweighted_hessian_identity():
    _, hessian, perpendicular_hessian = search.kernel_polynomials()
    expected = {
        ("constant",): -64,
        ("pair", 4): 640,
        ("pair", 6): -448,
    }
    assert_vector_close(search.expectation_vector(hessian), expected)
    assert_vector_close(
        search.expectation_vector(perpendicular_hessian),
        expected,
    )


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


def test_isotropic_target_vanishes_on_onb():
    target = (
        1
        - 9 * float(search.onb_label_value(("pair", 4)))
        + 6 * float(search.onb_label_value(("pair", 6)))
    )
    assert target == 0


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
