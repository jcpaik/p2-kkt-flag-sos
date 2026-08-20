"""Exact verifier for the rank-two fermionic scalar bridge.

Fix

    S = diag(1,-1,0)/sqrt(2)

and let G be any element of the exact linear span of tangent Pluecker
projectors for which S is an eigenvector of the contraction F.  With the
star/outer split relative to S, define U, C, B as in
``fermionic_motzkin_unrestricted.md`` and

    m = integral x_3^2 dmu.

This file verifies exactly that

    2 ||U||^2 + 8 ||C||^2 - (tr(G)-3m)^2/6 >= 0.       (*)

Together with tr(B *B*) >= 0, (*) proves

    R0 >= (tr(G)-3m)^2/6

in the rank-two chart.  The proof is finite exact linear algebra on the
24-dimensional tangent/eigen slice, split into the four D2 characters.
"""

from __future__ import annotations

import itertools

import sympy as sp

import tensor_star_rank2_d2_exact as d2


def trace_form(gm: sp.Matrix) -> sp.Matrix:
    out = sp.zeros(1, gm.cols)
    for i in range(10):
        out += gm[10 * i + i, :]
    return out


def hodge_form(gm: sp.Matrix) -> sp.Matrix:
    """Polarized tr(B *B*) on the outer 6x6 block."""
    hodge = sp.zeros(6)
    for i, j, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
        hodge[i, j] = hodge[j, i] = sign
    outer_rows = [10 * (a + 4) + (b + 4) for a in range(6) for b in range(6)]
    bmap = gm[outer_rows, :]
    # vec(B)^T (H tensor H) vec(B) = tr(B H B H), with row-major vec.
    return sp.simplify(bmap.T * sp.kronecker_product(hodge, hodge) * bmap)


def scalar_form():
    gm, fm, _ = d2.exact_rank_two_slice()
    q = d2.purity_form(gm, fm)
    mass = trace_form(gm)
    top = fm[0, :]
    delta = 3 * top - 2 * mass

    # In the canonical Pluecker order, pointwise
    # x_3 = 2 z_{02} + z_{34}; hence m is this linear G functional.
    m = 4 * gm[10 * 1 + 1, :] + gm[10 * 9 + 9, :] + 4 * gm[10 * 1 + 9, :]
    anisotropy = mass - 3 * m
    target = sp.simplify(
        4 * q
        + sp.Rational(1, 6) * delta.T * delta
        - sp.Rational(1, 6) * anisotropy.T * anisotropy
    )
    scalar = sp.simplify(target - 4 * hodge_form(gm))
    return gm, fm, mass, m, scalar


def character_basis(gm: sp.Matrix, character: tuple[int, int]) -> sp.Matrix:
    outside = [
        10 * a + b
        for a in range(10)
        for b in range(10)
        if d2.operator_character(a, b) != character
    ]
    return sp.Matrix.hstack(*gm[outside, :].nullspace())


def positive_quotient(form: sp.Matrix):
    """Choose a coordinate complement to ker(form) and verify Sylvester."""
    null = form.nullspace()
    rank = form.rows - len(null)
    null_matrix = sp.Matrix.hstack(*null) if null else sp.zeros(form.rows, 0)
    chosen = None
    for indices in itertools.combinations(range(form.rows), rank):
        columns = sp.Matrix.hstack(*[sp.eye(form.rows).col(i) for i in indices])
        if sp.Matrix.hstack(columns, null_matrix).rank() == form.rows:
            chosen = columns
            break
    assert chosen is not None
    reduced = sp.simplify(chosen.T * form * chosen)
    minors = [sp.factor(reduced[:k, :k].det()) for k in range(1, rank + 1)]
    assert all(value > 0 for value in minors)
    return null, chosen, minors


def verify(verbose: bool = True):
    gm, fm, mass, m, scalar = scalar_form()
    results = {}
    total_rank = 0
    for character in ((0, 0), *d2.NONTRIVIAL):
        basis = character_basis(gm, character)
        block = sp.simplify(basis.T * scalar * basis)
        null, quotient, minors = positive_quotient(block)
        total_rank += block.rank()
        results[character] = {
            "basis": basis,
            "block": block,
            "null": null,
            "quotient": quotient,
            "minors": minors,
        }
        if verbose:
            print(
                "character",
                character,
                "dimension/rank/nullity",
                block.rows,
                block.rank(),
                len(null),
            )
            print("positive quotient leading minors =", minors)
    assert total_rank == scalar.rank()
    assert scalar.rank() == 18
    if verbose:
        print("full scalar rank/nullity", scalar.rank(), scalar.rows - scalar.rank())
        print("mass map", mass)
        print("m map", m)
    return results


if __name__ == "__main__":
    verify()
