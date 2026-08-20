"""Exact nonzero SO(2)-mode proof for the axial fermionic bridge.

For S=diag(1,-2,1)/sqrt(6), write p=x+i z and q=x-i z.  Homogeneous
sextic moments split into SO(2) weights m=a-b for p^a q^b y^c.  This
script builds the tangent Pluecker moment G exactly in each weight,
imposes FS=cS in weights one and two, and verifies that every nonzero
weight block of

    R0 = 2||U||^2 + 8||C||^2 + 4 tr(B *B*)

is positive semidefinite.  The zero mode is handled separately in
``tensor_star_axial_zero_mode_exact.py``.
"""

from __future__ import annotations

import itertools

import sympy as sp

import fermionic_axial_zero_mode_exact as axial
import tensor_star_canonical_exact as canonical


I = sp.I
p, q, yy = sp.symbols("p q yy")
PAIRS = canonical.PAIRS
ZETA = axial.axial_pluecker()


def complex_polynomials():
    substitution = {
        canonical.R[0]: (p + q) / 2,
        canonical.R[1]: yy,
        canonical.R[2]: (p - q) / (2 * I),
    }
    return [sp.expand(value.subs(substitution)) for value in ZETA]


COMPLEX_ZETA = complex_polynomials()
ENTRY_POLYNOMIALS = [
    [sp.Poly(sp.expand(COMPLEX_ZETA[i] * COMPLEX_ZETA[j]), p, q, yy) for j in range(10)]
    for i in range(10)
]


def monomials(weight: int):
    return [
        (a, b, 6 - a - b)
        for a in range(7)
        for b in range(7 - a)
        if a - b == weight
    ]


def moment_maps(weight: int):
    out = []
    for a, b, c in monomials(weight):
        monomial = p**a * q**b * yy**c
        out.append(
            sp.Matrix(
                10,
                10,
                lambda i, j: ENTRY_POLYNOMIALS[i][j].coeff_monomial(monomial),
            )
        )
    return out


def contraction(G: sp.Matrix) -> sp.Matrix:
    out = sp.zeros(5)
    for a, (i, j) in enumerate(PAIRS):
        for b, (k, ell) in enumerate(PAIRS):
            value = G[a, b]
            if j == ell:
                out[i, k] += value
            if j == k:
                out[i, ell] -= value
            if i == ell:
                out[j, k] -= value
            if i == k:
                out[j, ell] += value
    return sp.simplify((out + out.T) / 2)


def outer_contraction(B: sp.Matrix) -> sp.Matrix:
    pairs = list(itertools.combinations(range(4), 2))
    out = sp.zeros(4)
    for a, (i, j) in enumerate(pairs):
        for b, (k, ell) in enumerate(pairs):
            value = B[a, b]
            if j == ell:
                out[i, k] += value
            if j == k:
                out[i, ell] -= value
            if i == ell:
                out[j, k] -= value
            if i == k:
                out[j, ell] += value
    return sp.simplify((out + out.T) / 2)


HODGE = sp.zeros(6)
for i, j, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
    HODGE[i, j] = HODGE[j, i] = sign


def components(G: sp.Matrix):
    F = contraction(G)
    A, C, B = G[:4, :4], G[:4, 4:], G[4:, 4:]
    mass = sp.trace(G)
    delta = 3 * F[0, 0] - 2 * mass
    U = sp.simplify(A - outer_contraction(B) - delta * sp.eye(4) / 4)
    return U, C, B


def bilinear(G: sp.Matrix, H: sp.Matrix):
    UG, CG, BG = components(G)
    UH, CH, BH = components(H)
    return sp.simplify(
        2 * sp.trace(UG.T * UH)
        + 8 * sp.trace(CG.T * CH)
        + 4 * sp.trace(BG * HODGE * BH * HODGE)
    )


def eigen_slice(weight: int):
    plus = moment_maps(weight)
    constraints = sp.Matrix(
        [[contraction(G)[0, orbital] for G in plus] for orbital in range(1, 5)]
    )
    null = sp.Matrix.hstack(*constraints.nullspace())
    if null.cols == 0:
        return [], constraints
    sliced = [
        sp.simplify(sum((null[i, j] * plus[i] for i in range(len(plus))), sp.zeros(10)))
        for j in range(null.cols)
    ]
    return sliced, constraints


def hermitian_block(weight: int):
    plus, constraints = eigen_slice(weight)
    minus_raw = moment_maps(-weight)
    # Conjugating p^a q^b exchanges the weight signs.  The eigen-slice
    # basis on the negative side is the entrywise conjugate of the positive.
    minus = [G.conjugate() for G in plus]
    # Audit that these matrices lie in the raw negative-weight span.
    raw_columns = sp.Matrix.hstack(*[G.reshape(100, 1) for G in minus_raw])
    for G in minus:
        assert raw_columns.row_join(G.reshape(100, 1)).rank() == raw_columns.rank()
    block = sp.Matrix(
        len(plus),
        len(plus),
        lambda i, j: sp.simplify(2 * bilinear(minus[i], plus[j])),
    )
    assert sp.simplify(block - block.conjugate().T) == sp.zeros(block.rows)
    return block, constraints


def positive_quotient(block: sp.Matrix):
    null = block.nullspace()
    rank = block.rows - len(null)
    null_matrix = sp.Matrix.hstack(*null) if null else sp.zeros(block.rows, 0)
    chosen = None
    for indices in itertools.combinations(range(block.rows), rank):
        columns = sp.Matrix.hstack(*[sp.eye(block.rows).col(i) for i in indices])
        if sp.Matrix.hstack(columns, null_matrix).rank() == block.rows:
            chosen = columns
            break
    assert chosen is not None
    reduced = sp.simplify(chosen.T * block * chosen)
    minors = [sp.factor(reduced[:k, :k].det()) for k in range(1, rank + 1)]
    assert all(value > 0 for value in minors)
    return null, minors


def verify(verbose: bool = True):
    results = {}
    for weight in range(1, 7):
        block, constraints = hermitian_block(weight)
        null, minors = positive_quotient(block)
        results[weight] = (block, constraints, null, minors)
        if verbose:
            print(
                "weight",
                weight,
                "raw/slice/rank/nullity",
                len(monomials(weight)),
                block.rows,
                block.rank(),
                len(null),
            )
            print("eigen constraint rank", constraints.rank())
            print("positive quotient leading minors", minors)
    return results


if __name__ == "__main__":
    verify()
