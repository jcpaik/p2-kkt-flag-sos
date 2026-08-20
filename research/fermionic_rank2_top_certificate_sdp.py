"""Search a top-eigenvalue two-cone certificate in the rank-two chart.

This is a discovery script.  Fixes

    S = diag(1,-1,0)/sqrt(2)

as an eigenvector of the one-particle marginal F.  On the exact 24D
tangent-moment slice it tests whether the homogeneous purity gap is a sum
of ordinary squares and polarized positive-cone products built from

    G >= 0,  K = <S,FS>I_4 - F|_{S^perp} >= 0.
"""

from __future__ import annotations

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
import sympy as sp
from scipy.linalg import null_space

import tensor_star_canonical_exact as canonical


MONOMIALS = [(a, b, 6 - a - b) for a in range(7) for b in range(7 - a)]


def tangent_basis() -> list[np.ndarray]:
    x, y, z = canonical.R
    polynomial_entries = [
        [sp.Poly(f * g, x, y, z) for g in canonical.PLUECKER]
        for f in canonical.PLUECKER
    ]
    out = []
    for a, b, c in MONOMIALS:
        monomial = x**a * y**b * z**c
        out.append(
            np.array(
                [
                    [float(polynomial_entries[i][j].coeff_monomial(monomial)) for j in range(10)]
                    for i in range(10)
                ]
            )
        )
    return out


def contraction(G: np.ndarray) -> np.ndarray:
    F = np.zeros((5, 5))
    for a, (i, j) in enumerate(canonical.PAIRS):
        for b, (k, ell) in enumerate(canonical.PAIRS):
            value = G[a, b]
            if j == ell:
                F[i, k] += value
            if j == k:
                F[i, ell] -= value
            if i == ell:
                F[j, k] -= value
            if i == k:
                F[j, ell] += value
    return (F + F.T) / 2


RAW_G = tangent_basis()
RAW_F = [contraction(G) for G in RAW_G]
EIGEN_CONSTRAINT = np.array([[F[i, 0] for F in RAW_F] for i in range(1, 5)])
SLICE = null_space(EIGEN_CONSTRAINT)

G_BASIS = [
    sum(SLICE[i, j] * RAW_G[i] for i in range(28))
    for j in range(SLICE.shape[1])
]
F_BASIS = [contraction(G) for G in G_BASIS]
K_BASIS = [F[0, 0] * np.eye(4) - F[1:, 1:] for F in F_BASIS]
DIMENSION = len(G_BASIS)

# Hodge star on Lambda^2(S^perp), whose edge order is the final six entries
# of canonical.PAIRS.  For B >= 0, tr(B *B*) is nonnegative although it is
# not an ordinary square in the moment coordinates.
OUTER_EDGES = [(i - 1, j - 1) for i, j in canonical.PAIRS[4:]]
HODGE = np.zeros((6, 6))
for a, edge in enumerate(OUTER_EDGES):
    for b, other in enumerate(OUTER_EDGES):
        permutation = list(edge + other)
        if len(set(permutation)) == 4:
            inversions = sum(
                permutation[i] > permutation[j]
                for i in range(4)
                for j in range(i + 1, 4)
            )
            HODGE[a, b] = (-1) ** inversions

HODGE_FORM = np.empty((DIMENSION, DIMENSION))
for a in range(DIMENSION):
    for b in range(DIMENSION):
        A = G_BASIS[a][4:, 4:]
        B = G_BASIS[b][4:, 4:]
        HODGE_FORM[a, b] = np.trace(A @ HODGE @ B @ HODGE)
HODGE_FORM = (HODGE_FORM + HODGE_FORM.T) / 2

Q = np.empty((DIMENSION, DIMENSION))
for a in range(DIMENSION):
    for b in range(DIMENSION):
        Q[a, b] = (
            np.sum(G_BASIS[a] * G_BASIS[b])
            - np.sum(F_BASIS[a] * F_BASIS[b]) / 2
            + np.trace(G_BASIS[a]) * np.trace(G_BASIS[b]) / 3
        )
Q = (Q + Q.T) / 2


def linear_form(values):
    return np.array(values, dtype=float)


MASS = linear_form([np.trace(G) for G in G_BASIS])
TOP = linear_form([F[0, 0] for F in F_BASIS])
# Pointwise, z_coordinate = 2 z_02 + z_34 in the canonical Pluecker
# coordinates, so m=E[z_coordinate^2] is this exact linear functional.
Z2_MOMENT = linear_form(
    [4 * G[1, 1] + G[9, 9] + 4 * G[1, 9] for G in G_BASIS]
)
DELTA = 3 * TOP - 2 * MASS
RANK2_ANISOTROPY = MASS - 3 * Z2_MOMENT
T = (
    4 * Q
    + np.outer(DELTA, DELTA) / 6
    - np.outer(RANK2_ANISOTROPY, RANK2_ANISOTROPY) / 6
)
T = (T + T.T) / 2


def pairs(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def polarized_compound(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    edges = pairs(len(A))
    out = np.empty((len(edges), len(edges)))
    for p, (i, j) in enumerate(edges):
        for q, (k, ell) in enumerate(edges):
            out[p, q] = (
                A[i, k] * B[j, ell]
                + B[i, k] * A[j, ell]
                - A[i, ell] * B[j, k]
                - B[i, ell] * A[j, k]
            ) / 2
    return (out + out.T) / 2


def solve(
    use_g_compound=True,
    use_k_compound=True,
    use_cross=True,
    use_hodge=True,
    target="T",
):
    sos = cp.Variable((DIMENSION, DIMENSION), symmetric=True)
    constraints = [sos >> 0]
    variables = [("sos", sos)]

    hodge_weight = None
    if use_hodge:
        hodge_weight = cp.Variable(nonneg=True)
        variables.append(("hodge", hodge_weight))

    wg = wk = wcross = None
    if use_g_compound:
        wg = cp.Variable((45, 45), symmetric=True)
        constraints.append(wg >> 0)
        variables.append(("wg", wg))
    if use_k_compound:
        wk = cp.Variable((6, 6), symmetric=True)
        constraints.append(wk >> 0)
        variables.append(("wk", wk))
    if use_cross:
        wcross = cp.Variable((40, 40), symmetric=True)
        constraints.append(wcross >> 0)
        variables.append(("wcross", wcross))

    target_matrix = T if target == "T" else Q
    for a in range(DIMENSION):
        for b in range(a, DIMENSION):
            value = sos[a, b]
            if hodge_weight is not None:
                value += hodge_weight * HODGE_FORM[a, b]
            if wg is not None:
                value += cp.sum(
                    cp.multiply(wg, polarized_compound(G_BASIS[a], G_BASIS[b]))
                )
            if wk is not None:
                value += cp.sum(
                    cp.multiply(wk, polarized_compound(K_BASIS[a], K_BASIS[b]))
                )
            if wcross is not None:
                cross = (
                    np.kron(K_BASIS[a], G_BASIS[b])
                    + np.kron(K_BASIS[b], G_BASIS[a])
                ) / 2
                value += cp.sum(cp.multiply(wcross, cross))
            constraints.append(value == target_matrix[a, b])

    problem = cp.Problem(
        cp.Minimize(
            sum(
                variable if variable.ndim == 0 else cp.trace(variable)
                for _, variable in variables
            )
        ),
        constraints,
    )
    try:
        value = problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-8,
            tol_feas=2e-8,
            tol_gap_rel=2e-8,
            max_iter=3000,
        )
    except cp.error.SolverError:
        value = problem.solve(solver="SCS", eps=3e-7, max_iters=500_000)

    print(
        "slice dimension",
        DIMENSION,
        "target",
        target,
        "eig",
        np.linalg.eigvalsh(target_matrix),
    )
    print("status/value", problem.status, value)
    for name, variable in variables:
        if variable.value is not None:
            if variable.ndim == 0:
                print(name, variable.value)
            else:
                print(name, "eig", np.linalg.eigvalsh(variable.value))
    return problem.status


if __name__ == "__main__":
    solve()
