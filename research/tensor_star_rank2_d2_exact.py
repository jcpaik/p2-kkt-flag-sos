"""Exact D2 character decomposition of the rank-two tangent-star slice.

This is a verifier/discovery aid for the remaining uncapped fermionic
inequality.  Fix the top orbital

    S = diag(1,-1,0)/sqrt(2)

and impose only F S = c S on the 28-dimensional linear span of genuine
tangent Pluecker projectors.  The resulting slice has dimension 24.  The
physical stabilizer D2 of S splits the non-invariant part into three
five-dimensional character sectors.  In every such sector the purity form
Q=||G||^2-||F||^2/2+tr(G)^2/3 has the sharp completion

    Q_chi + (3/22) ||K_chi||^2 >= 0,

where K=c I-F|_{S^perp}.  The completed form has rank four.  This is an
exact sector lemma, not by itself a proof for arbitrary mixtures of the
four D2 sectors: the negative K terms still have to be absorbed using the
PSD Schur constraints of G and K.
"""

from __future__ import annotations

import itertools

import sympy as sp

import tensor_star_canonical_exact as canonical


rt = sp.sqrt
EDGES = canonical.PAIRS


def contraction_matrix() -> sp.Matrix:
    """Map vec(G) (row-major, actual entries) to vec(F)."""
    out = sp.zeros(25, 100)
    for a, (i, j) in enumerate(EDGES):
        for b, (k, ell) in enumerate(EDGES):
            column = 10 * a + b
            if j == ell:
                out[5 * i + k, column] += 1
            if j == k:
                out[5 * i + ell, column] -= 1
            if i == ell:
                out[5 * j + k, column] -= 1
            if i == k:
                out[5 * j + ell, column] += 1
    return out


def exact_rank_two_slice():
    """Return exact actual-entry maps GM, FM, KM on the 24D slice."""
    coordinates, span = canonical.tangent_coefficient_span()
    # ``span`` uses orthonormal symmetric vectorization: sqrt(2) G_ab for
    # a<b.  Expand it to all 100 actual matrix entries.
    gmap = sp.zeros(100, 28)
    for row, (a, b) in enumerate(coordinates):
        value = span[row, :] / (rt(2) if a < b else 1)
        gmap[10 * a + b, :] = value
        gmap[10 * b + a, :] = value

    fmap = contraction_matrix() * gmap
    eigen_rows = [5 * i for i in range(1, 5)]
    eigen_constraint = fmap[eigen_rows, :]
    null = sp.Matrix.hstack(*eigen_constraint.nullspace())
    assert null.shape == (28, 24)
    gm = sp.simplify(gmap * null)
    fm = sp.simplify(fmap * null)

    # K=c I_4-F|S^perp, flattened as an actual 4x4 matrix.
    km = sp.zeros(16, 24)
    c = fm[0, :]
    for i in range(4):
        for j in range(4):
            km[4 * i + j, :] = (c if i == j else sp.zeros(1, 24)) - fm[
                5 * (i + 1) + (j + 1), :
            ]
    return gm, fm, km


def purity_form(gm: sp.Matrix, fm: sp.Matrix) -> sp.Matrix:
    trace_g = sp.zeros(1, gm.cols)
    for i in range(10):
        trace_g += gm[11 * i, :]
    return sp.simplify(
        gm.T * gm - sp.Rational(1, 2) * fm.T * fm
        + sp.Rational(1, 3) * trace_g.T * trace_g
    )


ORBITAL_CHARS = [(0, 0), (0, 0), (1, 1), (1, 0), (0, 1)]
WEDGE_CHARS = [
    ((ORBITAL_CHARS[i][0] + ORBITAL_CHARS[j][0]) % 2,
     (ORBITAL_CHARS[i][1] + ORBITAL_CHARS[j][1]) % 2)
    for i, j in EDGES
]
NONTRIVIAL = [(1, 1), (1, 0), (0, 1)]


def operator_character(a: int, b: int) -> tuple[int, int]:
    return (
        (WEDGE_CHARS[a][0] + WEDGE_CHARS[b][0]) % 2,
        (WEDGE_CHARS[a][1] + WEDGE_CHARS[b][1]) % 2,
    )


def sector_basis(gm: sp.Matrix, character: tuple[int, int]) -> sp.Matrix:
    outside = [
        10 * a + b
        for a in range(10)
        for b in range(10)
        if operator_character(a, b) != character
    ]
    basis = sp.Matrix.hstack(*gm[outside, :].nullspace())
    assert basis.shape == (24, 5)
    return basis


def nonzero_matrix_entries(flat: sp.Matrix, n: int):
    out = []
    for i in range(n):
        for j in range(i, n):
            value = sp.factor(flat[n * i + j])
            if value != 0:
                out.append(((i, j), value))
    return out


def positive_quotient_data(h: sp.Matrix):
    null = h.nullspace()
    assert len(null) == 1
    v = null[0]
    # Select four standard-coordinate columns independent modulo the null.
    chosen = None
    for indices in itertools.combinations(range(5), 4):
        columns = [sp.eye(5).col(i) for i in indices]
        test = sp.Matrix.hstack(*columns, v)
        if test.det() != 0:
            chosen = sp.Matrix.hstack(*columns)
            break
    assert chosen is not None
    reduced = sp.simplify(chosen.T * h * chosen)
    minors = [sp.factor(reduced[:k, :k].det()) for k in range(1, 5)]
    assert all(value > 0 for value in minors)
    return v, chosen, minors


def verify(verbose: bool = True):
    gm, fm, km = exact_rank_two_slice()
    q = purity_form(gm, fm)
    results = {}
    for character in NONTRIVIAL:
        basis = sector_basis(gm, character)
        q_sector = sp.simplify(basis.T * q * basis)
        k_sector = sp.simplify(basis.T * km.T * km * basis)
        completed = sp.simplify(
            q_sector + sp.Rational(3, 22) * k_sector
        )
        null, quotient, minors = positive_quotient_data(completed)
        physical_g = sp.simplify(gm * basis * null)
        physical_k = sp.simplify(km * basis * null)
        results[character] = {
            "basis": basis,
            "q": q_sector,
            "k": k_sector,
            "completed": completed,
            "null": null,
            "minors": minors,
            "G_null": physical_g,
            "K_null": physical_k,
        }
        assert completed.rank() == 4
        if verbose:
            print("character", character)
            print("q =")
            print(q_sector)
            print("k =")
            print(k_sector)
            print("completion null =", null.T)
            print("positive quotient leading minors =", minors)
            print("null physical G entries =", nonzero_matrix_entries(physical_g, 10))
            print("null physical K entries =", nonzero_matrix_entries(physical_k, 4))
    return results


if __name__ == "__main__":
    verify()
