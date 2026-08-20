"""Exact SO(2)-invariant reduction for the axial fermionic bridge (h=1).

The top orbital is S=diag(1,-2,1)/sqrt(6), invariant under rotations in
the xz-plane.  Averaging a tangent Pluecker projector over that SO(2)
orbit leaves a matrix polynomial G(u), u=y^2, of degree three.  This file
derives the exact quadratic R0 and the one-body occupations in the four
moments s_k=integral u^k dmu, k=0,...,3.
"""

from __future__ import annotations

import itertools

import sympy as sp

import tensor_star_canonical_exact as canonical


u = sp.symbols("u", real=True)
s = sp.symbols("s0:4", real=True)
x, y, z = canonical.R


def axial_orbitals():
    return [
        (sp.sqrt(3) * canonical.ORB[0] - canonical.ORB[1]) / 2,
        (canonical.ORB[0] + sp.sqrt(3) * canonical.ORB[1]) / 2,
        *canonical.ORB[2:],
    ]


def axial_pluecker():
    orbitals = axial_orbitals()
    return [
        sp.expand(2 * canonical.R.dot((orbitals[i] * canonical.R).cross(orbitals[j] * canonical.R)))
        for i, j in canonical.PAIRS
    ]


def circle_average(polynomial: sp.Expr) -> sp.Expr:
    """Average a homogeneous even degree-six polynomial at y^2=u."""
    poly = sp.Poly(sp.expand(polynomial), x, y, z)
    out = 0
    for (a, b, c), coefficient in poly.terms():
        if a % 2 or c % 2 or b % 2:
            continue
        aa, cc = a // 2, c // 2
        # E[cos^(2aa) sin^(2cc)] = (2aa)!(2cc)! /
        # (4^(aa+cc) aa! cc! (aa+cc)!).
        angular = (
            sp.factorial(2 * aa)
            * sp.factorial(2 * cc)
            / (
                4 ** (aa + cc)
                * sp.factorial(aa)
                * sp.factorial(cc)
                * sp.factorial(aa + cc)
            )
        )
        out += coefficient * u ** (b // 2) * (1 - u) ** (aa + cc) * angular
    return sp.factor(sp.expand(out))


def contraction(G: sp.Matrix) -> sp.Matrix:
    out = sp.zeros(5)
    for a, (i, j) in enumerate(canonical.PAIRS):
        for b, (k, ell) in enumerate(canonical.PAIRS):
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


def moment_substitute(polynomial: sp.Expr) -> sp.Expr:
    poly = sp.Poly(sp.expand(polynomial), u)
    return sp.expand(sum(coefficient * s[degree[0]] for degree, coefficient in poly.terms()))


def data():
    zeta = axial_pluecker()
    Gu = sp.Matrix(10, 10, lambda i, j: circle_average(zeta[i] * zeta[j]))
    G = Gu.applyfunc(moment_substitute)
    F = contraction(G)
    A, C, B = G[:4, :4], G[:4, 4:], G[4:, 4:]
    mass = sp.trace(G)
    delta = 3 * F[0, 0] - 2 * mass
    U = sp.simplify(A - outer_contraction(B) - delta * sp.eye(4) / 4)
    hodge = sp.zeros(6)
    for i, j, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
        hodge[i, j] = hodge[j, i] = sign
    R0 = sp.factor(
        2 * sp.trace(U.T * U)
        + 8 * sp.trace(C.T * C)
        + 4 * sp.trace(B * hodge * B * hodge)
    )
    target = sp.factor(R0 - mass**2 / 96)
    return Gu, G, F, U, C, B, mass, delta, R0, target


def verify(verbose: bool = True):
    Gu, G, F, U, C, B, mass, delta, R0, target = data()
    assert mass == s[0]
    assert F[0, 0] == 3 * (s[1] - s[2])
    if verbose:
        print("G(u) nonzero entries")
        for i in range(10):
            for j in range(i, 10):
                if Gu[i, j] != 0:
                    print((i, j), Gu[i, j])
        print("F diagonal", [sp.factor(F[i, i]) for i in range(5)])
        print("mass", mass, "delta", sp.factor(delta))
        print("R0 =", R0)
        print("R0-mass^2/96 =", target)
    return Gu, G, F, U, C, B, mass, delta, R0, target


if __name__ == "__main__":
    verify()
