"""Exact canonical tangent-star certificate.

This proves the Motzkin star inequality when the top marginal eigenvector is
E0=diag(1,-1,0)/sqrt(2) and the remaining orbital eigenbasis is the natural
diagonal/off-diagonal completion displayed below.  It is a sharp subcase,
not an unrestricted proof.
"""

import itertools

import sympy as sp


rt = sp.sqrt
ORB = [sp.diag(1, -1, 0) / rt(2), sp.diag(1, 1, -2) / rt(6)]
for i, j in ((0, 1), (0, 2), (1, 2)):
    matrix = sp.zeros(3)
    matrix[i, j] = matrix[j, i] = 1 / rt(2)
    ORB.append(matrix)

PAIRS = list(itertools.combinations(range(5), 2))
R = sp.Matrix(sp.symbols("x y z"))
PLUECKER = [
    sp.expand(2 * R.dot((ORB[i] * R).cross(ORB[j] * R))) for i, j in PAIRS
]


def tangent_coefficient_span():
    monomials = [
        (a, b, 6 - a - b) for a in range(7) for b in range(7 - a)
    ]
    coordinates = [(a, b) for a in range(10) for b in range(a, 10)]
    out = sp.zeros(55, 28)
    for row, (a, b) in enumerate(coordinates):
        polynomial = sp.Poly(
            sp.expand(PLUECKER[a] * PLUECKER[b] * (rt(2) if a < b else 1)),
            *R,
        )
        for column, (i, j, k) in enumerate(monomials):
            out[row, column] = polynomial.coeff_monomial(
                R[0] ** i * R[1] ** j * R[2] ** k
            )
    return coordinates, out


def pair_with_sign(i, j):
    return ((i, j), 1) if i < j else ((j, i), -1)


def exact_least_norm_form():
    """Derive compatibility and the exact forced-offdiagonal Gram form."""
    coordinates, span = tangent_coefficient_span()
    null = sp.Matrix.hstack(*span.T.nullspace())
    diagonal_rows = [coordinates.index((i, i)) for i in range(10)]
    off_pairs = [(i, j) for i in range(10) for j in range(i + 1, 10)]
    off_rows = [coordinates.index(pair) for pair in off_pairs]
    A_span = null[off_rows, :].T
    B = null[diagonal_rows, :].T

    edge_index = {pair: i for i, pair in enumerate(PAIRS)}
    off_index = {pair: i for i, pair in enumerate(off_pairs)}
    A_marginal = sp.zeros(10, 45)
    row = 0
    for i in range(5):
        for k in range(i + 1, 5):
            for j in range(5):
                if j == i or j == k:
                    continue
                pi, si = pair_with_sign(i, j)
                pk, sk = pair_with_sign(k, j)
                a, b = sorted((edge_index[pi], edge_index[pk]))
                A_marginal[row, off_index[(a, b)]] += sp.Rational(si * sk, 1) / rt(2)
            row += 1

    A = A_span.col_join(A_marginal)
    B = B.col_join(sp.zeros(10, 10))
    compatibility = (A.T.nullspace()[0].T * B)
    _, pivots = A.T.rref()
    Ar = A[list(pivots), :]
    Br = B[list(pivots), :]
    K = -Ar.T * (Ar * Ar.T).inv() * Br
    return compatibility, sp.simplify(K.T * K)


def graph_form():
    p = sp.symbols("p0:10")
    a = sp.Matrix(p[:4])
    b = {
        (i, j): p[PAIRS.index((i, j))]
        for i in range(1, 5)
        for j in range(i + 1, 5)
    }
    degree = sp.Matrix(
        [
            sum(b[tuple(sorted((i, j)))] for j in range(1, 5) if j != i)
            for i in range(1, 5)
        ]
    )
    disjoint = b[1, 2] * b[3, 4] + b[1, 3] * b[2, 4] + b[1, 4] * b[2, 3]
    total = sum(p)
    center = sum(p[:4])
    delta = 3 * center - 2 * total
    u0 = a - degree - delta * sp.ones(4, 1) / 4
    value = 2 * disjoint + u0.dot(u0) / 2 - delta**2 / 24
    return sp.hessian(sp.expand(value), p) / 2


def verify():
    compatibility, forced = exact_least_norm_form()
    # Normalize [1,-3,0,0,1,...] to p02=(p01+p12)/3.
    assert compatibility.rank() == 1
    assert [sp.simplify(v) for v in compatibility] == [1, -3, 0, 0, 1, 0, 0, 0, 0, 0]

    # Eliminate p02.  The independent y coordinates are
    # (p01,p03,p04,p12,p13,p14,p23,p24,p34).
    indices = [0, 2, 3, 4, 5, 6, 7, 8, 9]
    substitution = sp.zeros(10, 9)
    for column, index in enumerate(indices):
        substitution[index, column] = 1
    substitution[1, 0] = substitution[1, 3] = sp.Rational(1, 3)
    Q = sp.simplify(substitution.T * (forced + graph_form()) * substitution)

    one = sp.ones(1, 10) * substitution
    center = sp.Matrix([[int(i == 0) for i, _ in PAIRS]]) * substitution
    forms = {f"y{i}": sp.eye(9).row(i) for i in range(9)}
    for k in range(1, 5):
        degree = sp.Matrix(
            [[int(i == k or j == k) for i, j in PAIRS]]
        ) * substitution
        forms[f"cap{k}"] = center - degree

    def product(a, b):
        return (forms[a].T * forms[b] + forms[b].T * forms[a]) / 2

    sos = (
        Q
        - sp.Rational(4, 3) * product("y3", "y8")
        - sp.Rational(1, 3) * product("y3", "cap1")
        - sp.Rational(1, 3) * (product("y6", "cap3") + product("y7", "cap4"))
        - sp.Rational(1, 12) * (product("y6", "cap4") + product("y7", "cap3"))
    )

    null_vector = sp.Matrix([0, 1, 1, 0, 0, 0, 0, 0, 1])
    assert sos * null_vector == sp.zeros(9, 1)
    quotient_columns = [sp.eye(9).col(i) for i in (0, 3, 4, 5, 6, 7)] + [
        sp.eye(9).col(1) - sp.eye(9).col(2),
        sp.eye(9).col(1) + sp.eye(9).col(2) - 2 * sp.eye(9).col(8),
    ]
    quotient = sp.Matrix.hstack(*quotient_columns)
    reduced = quotient.T * sos * quotient
    leading_minors = [sp.factor(reduced[:k, :k].det()) for k in range(1, 9)]
    expected = [
        sp.Rational(25, 27), sp.Rational(86, 81), sp.Rational(6737, 972),
        sp.Rational(7175, 972), sp.Rational(794251, 186624),
        sp.Rational(281119, 124416), sp.Rational(121295, 15552),
        sp.Rational(34255, 10368),
    ]
    assert leading_minors == expected
    print("Pluecker norm:", sp.factor(sum(z * z for z in PLUECKER)))
    print("span rank:", tangent_coefficient_span()[1].rank())
    print("compatibility:", compatibility)
    print("positive quotient leading minors:", leading_minors)


if __name__ == "__main__":
    verify()
