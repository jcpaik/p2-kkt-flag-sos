"""Exact axial SO(2)-zero-mode certificate for the scalar Motzkin bound.

For the axial top orbital S=diag(1,-2,1)/sqrt(6), put u=y^2 and
let a=E[u], b=E[u^2], c=E[u^3] for an SO(2)-invariant measure.  The
strong scalar target is R0 >= 1/96.  After clearing the denominator 96,
its deficit is the polynomial ``TARGET`` below.  This verifier proves it
by three elementary Hausdorff moment determinants and a positive definite
quadratic remainder.

This handles the azimuthal zero mode at the axial spectral endpoint.  The
nonzero SO(2) modes must still be included before it becomes an unrestricted
axial theorem.
"""

import sympy as sp


a, b, c = sp.symbols("a b c", real=True)
y = sp.Matrix([1, a, b, c])

TARGET = (
    5328 * a**2
    - 28512 * a * b
    + 20160 * a * c
    - 960 * a
    + 41616 * b**2
    - 60480 * b * c
    + 1728 * b
    + 22176 * c**2
    - 960 * c
    + 95
)

HAUSDORFF = [
    b - a**2,
    a * c - b**2,
    (1 - a) * (b - c) - (a - b) ** 2,
]


def quadratic_matrix(polynomial: sp.Expr) -> sp.Matrix:
    """Return M with polynomial=[1,a,b,c] M [1,a,b,c]^T."""
    variables = (a, b, c)
    out = sp.zeros(4)
    for monomial, coefficient in sp.Poly(sp.expand(polynomial), *variables).terms():
        indices = []
        for index, power in enumerate(monomial, 1):
            indices.extend([index] * power)
        if not indices:
            out[0, 0] += coefficient
        elif len(indices) == 1:
            out[0, indices[0]] += coefficient / 2
            out[indices[0], 0] += coefficient / 2
        else:
            i, j = indices
            out[i, j] += coefficient / 2
            out[j, i] += coefficient / 2
    return out


REMAINDER = sp.expand(
    TARGET
    - 16 * HAUSDORFF[0]
    - 384 * HAUSDORFF[1]
    - 64 * HAUSDORFF[2]
)
R = quadratic_matrix(REMAINDER)


def verify():
    assert sp.expand((y.T * R * y)[0] - REMAINDER) == 0
    assert sp.expand(
        (y.T * R * y)[0]
        + 16 * HAUSDORFF[0]
        + 384 * HAUSDORFF[1]
        + 64 * HAUSDORFF[2]
        - TARGET
    ) == 0
    leading_minors = [sp.factor(R[:k, :k].det()) for k in range(1, 5)]
    expected = [
        95,
        283360,
        155772672,
        18717138944,
    ]
    assert leading_minors == expected
    print("positive remainder matrix =")
    print(R)
    print("leading principal minors =", leading_minors)


if __name__ == "__main__":
    verify()
