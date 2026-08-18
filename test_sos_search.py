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


def test_weighted_e1_two_root_classification_small_degree():
    """Weighted-(E1): admissible two-root leaves are the deviatoric span.

    At flag degree 4 the closed form (docs/SHARP_STRUCTURE.md) must
    equal the exact nullspace of the isotropic necessary rows.
    """

    import solve_e1 as e1

    pairs = e1.sample_root_pairs(9)
    for parity, det_sector in ((0, False), (1, True)):
        basis = e1.two_root_basis(4, parity)
        rows = e1.weighted_two_root_rows(basis, det_sector, pairs)
        space = e1.nullspace(rows, len(basis))
        predicted = e1.weighted_two_root_predicted(basis, det_sector)
        assert e1.same_span(space, predicted)


def test_weighted_e1_sufficiency_identities():
    import solve_e1 as e1

    assert e1.weighted_sufficiency_identities(3)


def test_weighted_projection_loader_round_trip(tmp_path):
    """The weighted export is consumed by load_e1_projection with the
    weighted extras (one-root, pair, harmonic-flag bases)."""

    import solve_e1 as e1

    path = tmp_path / "e1w_projection_deg10.json"
    e1.export_projection(str(path), 10, weighted=True)
    degree, sectors, spin2, weighted = search.load_e1_projection(str(path))
    assert degree == 10
    assert weighted is not None
    assert weighted["one_root"][0] == [
        {0: Fraction(-1, 3), 2: Fraction(1)}
    ]
    assert weighted["pair_flags"] == [{0: Fraction(-1), 2: Fraction(3)}]
    assert set(sectors) == {
        "two_root_even_00", "two_root_even_00_minor",
        "two_root_even_11", "two_root_even_11_minor",
        "two_root_odd_01", "two_root_odd_01_minor",
        "two_root_odd_10", "two_root_odd_10_minor",
    }
    # degree-10 caps: flag degree 5, minors 4
    basis, vectors = sectors["two_root_even_00"]
    assert len(vectors) == 6  # alpha,beta in {1,s^2}; gamma in {s,s^3}
    assert spin2 is not None


# ---------------------------------------------------------------------------
# e5(I - A2) machinery (docs/ENRICHMENTS.md)


def test_gap_power_trace_small():
    # tr(A2) = 1, tr(A2^2) = 1 - 12 p2 + 16 p4 (chi2 of rho_x rho_y).
    assert search.gap_power_trace_vector(1) == {
        ("constant",): Fraction(1)
    }
    assert search.gap_power_trace_vector(2) == {
        ("constant",): Fraction(1),
        ("pair", 2): Fraction(-12),
        ("pair", 4): Fraction(16),
    }


def test_gap_elementary_matches_hardcoded_cuts():
    """gap_elementary_vector reproduces the audited e2/e3/e4 cut labels."""

    assert search.gap_elementary_vector(2) == {
        ("constant",): Fraction(6),
        ("pair", 2): Fraction(6),
        ("pair", 4): Fraction(-8),
    }
    assert search.gap_elementary_vector(3) == {
        ("constant",): Fraction(8, 3),
        ("pair", 2): Fraction(32),
        ("pair", 4): Fraction(-32),
        ("triangle", 1, 1, 1): Fraction(-40, 3),
        ("triangle", 0, 2, 2): Fraction(-32),
        ("triangle", 1, 1, 3): Fraction(64),
        ("triangle", 2, 2, 2): Fraction(-64, 3),
    }
    e4 = search.gap_elementary_vector(4)
    assert e4[("constant",)] == Fraction(-22, 3)
    assert e4[("graph_4", 0, 1, 1, 2, 2, 1)] == Fraction(256)
    assert e4[
        search.multiply_labels(("pair", 2), ("pair", 4))
    ] == Fraction(-48)
    assert len(e4) == 16


def test_e5_structure_and_face_values():
    e5 = search.gap_elementary_vector(5)
    assert len(e5) == 35
    cycle_labels = [label for label in e5 if label[0] == "graph_5"]
    assert len(cycle_labels) == 11
    # per-vertex degree <= 4: e5 fits any degree budget >= 10; the
    # obstruction was label arity only.
    for label in cycle_labels:
        matrix = search.graph_matrix(5, tuple(label[1:]))
        assert max(sum(row) for row in matrix) <= 4
    # Face-vanishing (Fact 2): every pole-equator measure kills e5.
    assert sum(
        coefficient * search.pole_equator_label_value(label)
        for label, coefficient in e5.items()
    ) == 0
    for order in (3, 4, 5, 6):
        assert sum(
            coefficient * search.pole_equator_label_value(label, order)
            for label, coefficient in e5.items()
        ) == 0
    # ONB (spectrum {1,1,-1/3,-1/3,-1/3}): e5 = 0 there too.
    assert sum(
        coefficient * search.onb_label_value(label)
        for label, coefficient in e5.items()
    ) == 0
    # Dense positivity anchor (Fact 3): uniform measure gives (4/5)^5.
    assert sum(
        coefficient * search.uniform_label_value(label)
        for label, coefficient in e5.items()
    ) == Fraction(1024, 3125)


def test_uniform_label_values():
    assert search.uniform_label_value(("pair", 2)) == Fraction(1, 3)
    assert search.uniform_label_value(("pair", 4)) == Fraction(1, 5)
    assert search.uniform_label_value(("pair", 3)) == 0
    assert search.uniform_label_value(("triangle", 1, 1, 1)) == Fraction(
        1, 9
    )
    # 4-cycle: integrate one vertex -> (1/3) triangle(1,1,1).
    assert search.uniform_label_value(
        ("graph_4", 1, 0, 1, 1, 0, 1)
    ) == Fraction(1, 27)
    assert search.uniform_label_value(
        search.multiply_labels(("pair", 2), ("pair", 4))
    ) == Fraction(1, 15)


def test_e5_expansion_matches_direct_determinant_on_atoms():
    """Exact rational atomic measure: label expansion == det(I - A2)."""

    import itertools

    import sympy as sp

    atoms = [
        (Fraction(3, 5), Fraction(4, 5), Fraction(0)),
        (Fraction(0), Fraction(3, 5), Fraction(4, 5)),
        (Fraction(4, 9), Fraction(4, 9), Fraction(7, 9)),
    ]
    weights = [Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)]
    for atom in atoms:
        assert sum(value * value for value in atom) == 1

    def dot(left, right):
        return sum(a * b for a, b in zip(left, right))

    def label_value(label):
        if label == ("constant",):
            return Fraction(1)
        if label[0] == "product":
            value = Fraction(1)
            for factor in label[1:]:
                value *= label_value(factor)
            return value
        if label[0] == "pair":
            vertex_count = 2
        elif label[0] == "triangle":
            vertex_count = 3
        else:
            vertex_count = int(label[0].split("_")[1])
        exponents = tuple(label[1:])
        edges = search.graph_edges(vertex_count)
        total = Fraction(0)
        for assignment in itertools.product(
            range(len(atoms)), repeat=vertex_count
        ):
            term = Fraction(1)
            for vertex in assignment:
                term *= weights[vertex]
            for index, (left, right) in enumerate(edges):
                power = exponents[index]
                if power:
                    term *= dot(
                        atoms[assignment[left]], atoms[assignment[right]]
                    ) ** power
            total += term
        return total

    e5 = search.gap_elementary_vector(5)
    expansion_value = sum(
        coefficient * label_value(label)
        for label, coefficient in e5.items()
    )

    basis = [
        sp.diag(1, -1, 0),
        sp.diag(1, 1, -2),
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
    ]
    coordinate_matrix = sp.Matrix(
        [
            [element[0, 0] for element in basis],
            [element[1, 1] for element in basis],
            [element[0, 1] for element in basis],
            [element[0, 2] for element in basis],
            [element[1, 2] for element in basis],
        ]
    )
    columns = []
    for element in basis:
        image = sp.zeros(3, 3)
        for weight, atom in zip(weights, atoms):
            vector = sp.Matrix(
                [sp.Rational(c.numerator, c.denominator) for c in atom]
            )
            rho = 2 * vector * vector.T - sp.eye(3)
            image += (
                sp.Rational(weight.numerator, weight.denominator)
                * rho
                * element
                * rho
            )
        columns.append(
            coordinate_matrix.solve(
                sp.Matrix(
                    [
                        image[0, 0],
                        image[1, 1],
                        image[0, 1],
                        image[0, 2],
                        image[1, 2],
                    ]
                )
            )
        )
    determinant = (sp.eye(5) - sp.Matrix.hstack(*columns)).det()
    assert (
        Fraction(int(determinant.p), int(determinant.q)) == expansion_value
    )


def test_e5_weighted_target_labels_are_covered_and_rationalizable():
    """Every label of kappa*e5*E lies in the e5 cut/localized blocks, and
    every produced float coefficient survives exact rationalization."""

    energy = {
        ("constant",): Fraction(-4, 3),
        ("pair", 2): Fraction(20),
        ("pair", 4): Fraction(-48),
        ("pair", 6): Fraction(32),
    }
    e5 = search.gap_elementary_vector(5)
    product = search.multiply_label_vectors(e5, energy)
    block_labels = set(e5)
    for degree in (2, 4, 6):
        block_labels.update(
            search.multiply_label_vectors(
                e5, search.exact_harmonic_pair_vector(degree)
            )
        )
    assert set(product) <= block_labels
    for kappa in (Fraction(1, 4), Fraction(1), Fraction(4)):
        for vector in (e5, product):
            for coefficient in vector.values():
                scaled = kappa * coefficient
                assert search.rationalize_float(float(scaled)) == scaled


def test_exact_harmonic_pair_vector():
    assert search.exact_harmonic_pair_vector(2) == {
        ("constant",): Fraction(-1, 2),
        ("pair", 2): Fraction(3, 2),
    }
    # h2-weighted target coefficients (the target refactor guard).
    energy = {
        ("constant",): Fraction(-4, 3),
        ("pair", 2): Fraction(20),
        ("pair", 4): Fraction(-48),
        ("pair", 6): Fraction(32),
    }
    h2 = {("constant",): Fraction(-1, 2), ("pair", 2): Fraction(3, 2)}
    weighted = search.multiply_label_vectors(h2, energy)
    assert weighted == {
        ("constant",): Fraction(2, 3),
        ("pair", 2): Fraction(-12),
        ("pair", 4): Fraction(24),
        ("pair", 6): Fraction(-16),
        search.multiply_labels(("pair", 2), ("pair", 2)): Fraction(30),
        search.multiply_labels(("pair", 2), ("pair", 4)): Fraction(-72),
        search.multiply_labels(("pair", 2), ("pair", 6)): Fraction(48),
    }


def test_e5_coverage_families_on_reference_measures():
    """Upper AM-GM cut and e5-Hankel localizations at exact reference
    measures: equality at uniform, vanishing on the pole-equator face."""

    e5 = search.gap_elementary_vector(5)
    complement = {label: -value for label, value in e5.items()}
    complement[("constant",)] = complement.get(
        ("constant",), Fraction(0)
    ) + Fraction(1024, 3125)

    def value(vector, evaluator):
        return sum(
            coefficient * evaluator(label)
            for label, coefficient in vector.items()
        )

    # AM-GM cut: equality exactly at the uniform measure.
    assert value(complement, search.uniform_label_value) == 0
    assert value(complement, search.onb_label_value) == Fraction(
        1024, 3125
    )
    assert value(complement, search.pole_equator_label_value) == Fraction(
        1024, 3125
    )
    # e5-Hankel entries at uniform: e5 * [[1, 1/3], [1/3, 1/5]] and
    # e5 * [[1/3, 1/5], [1/5, 1/7]]; both PSD, dets > 0.
    scale = Fraction(1024, 3125)
    for powers in ((0, 2), (1, 3)):
        entries = {}
        for row_power in powers:
            for column_power in powers:
                total = row_power + column_power
                pair_vector = (
                    {("constant",): Fraction(1)}
                    if total == 0
                    else {("pair", total): Fraction(1)}
                )
                entries[(row_power, column_power)] = value(
                    search.multiply_label_vectors(e5, pair_vector),
                    search.uniform_label_value,
                )
        a = entries[(powers[0], powers[0])]
        b = entries[(powers[0], powers[1])]
        d = entries[(powers[1], powers[1])]
        assert a == scale * Fraction(1, 2 * powers[0] + 1)
        assert d == scale * Fraction(1, 2 * powers[1] + 1)
        assert a > 0 and a * d - b * b > 0
        # and the whole matrix vanishes on the pole-equator face
        face = value(
            search.multiply_label_vectors(
                e5, {("pair", 2 * powers[1]): Fraction(1)}
            ),
            search.pole_equator_label_value,
        )
        assert face == 0
