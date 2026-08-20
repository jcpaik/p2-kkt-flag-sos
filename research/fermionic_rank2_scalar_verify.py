"""Exact verifier for fermionic_rank2_scalar_proof.md."""

from __future__ import annotations

import sys

sys.path.insert(0, "research")

import sympy as sp

import tensor_star_canonical_exact as canonical


zeta = sp.Matrix(canonical.PLUECKER)
G = zeta * zeta.T


def contraction(matrix: sp.Matrix) -> sp.Matrix:
    out = sp.zeros(5)
    for a, (i, j) in enumerate(canonical.PAIRS):
        for b, (k, ell) in enumerate(canonical.PAIRS):
            value = matrix[a, b]
            if j == ell:
                out[i, k] += value
            if j == k:
                out[i, ell] -= value
            if i == ell:
                out[j, k] -= value
            if i == k:
                out[j, ell] += value
    return out


def outer_contraction(matrix: sp.Matrix) -> sp.Matrix:
    edges = [(i - 1, j - 1) for i, j in canonical.PAIRS[4:]]
    out = sp.zeros(4)
    for a, (i, j) in enumerate(edges):
        for b, (k, ell) in enumerate(edges):
            value = matrix[a, b]
            if j == ell:
                out[i, k] += value
            if j == k:
                out[i, ell] -= value
            if i == ell:
                out[j, k] -= value
            if i == k:
                out[j, ell] += value
    return out


A = G[:4, :4]
B = G[4:, 4:]
C = G[:4, 4:]
R = outer_contraction(B)
delta = sp.trace(A) - 2 * sp.trace(B)
U = A - R - delta * sp.eye(4) / 4

x, y, z = canonical.R
mass = sp.factor(sum(value**2 for value in zeta))
selected_coordinate = sp.factor(2 * zeta[1] + zeta[9])
m = sp.factor(selected_coordinate**2)

assert mass == (x**2 + y**2 + z**2) ** 3
assert selected_coordinate == z * (x**2 + y**2 + z**2)

left = sp.factor(mass - 3 * m)
right = sp.factor(
    -2 * U[1, 1]
    + U[2, 2]
    + U[3, 3]
    + 2 * sp.sqrt(3) * (C[2, 1] - C[3, 2])
)
assert sp.expand(left - right) == 0

L = sp.diag(0, -2, 1, 1)
M = sp.zeros(4, 6)
M[2, 1] = 2 * sp.sqrt(3)
M[3, 2] = -2 * sp.sqrt(3)
assert sp.trace(L.T * L) == 6
assert sp.trace(M.T * M) == 24

alpha = sp.symbols("alpha", nonnegative=True)
r = sp.symbols("r")
assert sp.expand(
    (r**2 - (r - 3 * alpha) ** 2) / 24
    - sp.Rational(3, 8) * alpha * (sp.Rational(2, 3) - 2 * (1 - r) / 3 - alpha)
) == 0

print("Pluecker mass =", mass)
print("2 z_(S,xy) + z_(xz,yz) =", selected_coordinate)
print("1-3m identity =", left)
print("weighted dual norms =", sp.trace(L.T * L), sp.trace(M.T * M))
