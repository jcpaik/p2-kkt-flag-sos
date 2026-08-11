from __future__ import annotations

import numpy as np

from collective_transport import (
    collective_transport_degrees,
    collective_transport_expectation_matrix,
    collective_transport_polynomials,
)
from sos_search import exact_onb_moment_matrix


def _evaluate(polynomial, x, y, z, w):
    vectors = (x, y, z, w)
    edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    grams = [float(vectors[i] @ vectors[j]) for i, j in edges]
    return sum(
        float(coefficient)
        * np.prod([value**power for value, power in zip(grams, exponent)])
        for coefficient, exponent in polynomial
    )


def _normalize(vector):
    return vector / np.linalg.norm(vector)


def test_collective_transport_polynomials_match_vector_formula():
    x = _normalize(np.array([1.0, 2.0, -1.0]))
    y = _normalize(np.array([-2.0, 1.0, 3.0]))
    z = _normalize(np.array([2.0, -3.0, 1.0]))
    w = _normalize(np.array([1.0, 4.0, 2.0]))

    a = float(x @ y)
    kp = 192 * a**5 - 192 * a**3 + 40 * a
    kpp = 960 * a**4 - 576 * a**2 + 40

    vxz = z - float(x @ z) * x
    vxw = w - float(x @ w) * x
    vyw = w - float(y @ w) * y

    expected_local = kpp * float(vxz @ y) * float(vxw @ y) - a * kp * float(vxz @ vxw)
    expected_cross = kpp * float(vxz @ y) * float(x @ vyw) + kp * float(vxz @ vyw)

    local, cross = collective_transport_polynomials()
    assert np.isclose(_evaluate(local, x, y, z, w), expected_local, atol=1e-10)
    assert np.isclose(_evaluate(cross, x, y, z, w), expected_cross, atol=1e-10)


def test_collective_transport_matrix_is_symmetric():
    matrices = collective_transport_expectation_matrix([1, 3])
    assert matrices
    for matrix in matrices.values():
        assert np.allclose(matrix, matrix.T)


def test_collective_transport_onb_is_psd():
    matrices = collective_transport_expectation_matrix([1, 3])
    moment_matrix = np.array(exact_onb_moment_matrix(matrices), dtype=float)
    assert np.linalg.eigvalsh(moment_matrix)[0] >= -1e-10


def test_collective_transport_degree_schedule():
    assert collective_transport_degrees(8) == []
    assert collective_transport_degrees(10) == [1]
    assert collective_transport_degrees(14) == [1, 3]
