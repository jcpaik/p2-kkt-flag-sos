"""Minor-certificate audit for the rank-two top-orbital tangent slice.

Fix E0 = diag(1,-1,0)/sqrt(2) and use the canonical orthonormal STF basis
from fermionic_tangent_bad_plane.md. The ten tangent Pluecker coordinates
are homogeneous cubics. Their 55 symmetric products span the exact
28-dimensional sextic Hankel space.

We additionally impose that E0 is an eigenvector of the one-particle
marginal. On that linear slice this script searches

    12||U||^2 + 48||C||^2 + 24 tr(B * B *) - delta^2
      = square + <W,C2(G)>,       W >= 0,

where G=[[A,C],[C^T,B]] in the star/outer split,
R=gamma(B), delta=3 tr(A)-2 tr(G), and U=A-R-delta I/4.
Feasibility is only a candidate until rationalized; infeasibility rejects
this smallest quadratic-minor ansatz.
"""

from __future__ import annotations

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIRS)}
OUTER = PAIRS[4:]
MINORS = [(i, j) for i in range(10) for j in range(i + 1, 10)]


def pluecker(x):
    """Canonical rank-two-S tangent Pluecker vector."""
    a, b, c = x
    return np.array(
        [
            2 * np.sqrt(3) * a * b * c,
            c * (a * a + b * b),
            -b * (2 * a * a - c * c),
            -a * (2 * b * b - c * c),
            np.sqrt(3) * c * (a * a - b * b),
            -np.sqrt(3) * b * c * c,
            np.sqrt(3) * a * c * c,
            a * (a * a - b * b - c * c),
            b * (a * a - b * b + c * c),
            -c * (a * a + b * b - c * c),
        ]
    )


def moment_basis(samples=300, seed=20260820):
    """Return an orthonormal basis of the exact orbit-linear span."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(samples):
        x = rng.normal(size=3)
        x /= np.linalg.norm(x)
        z = pluecker(x)
        rows.append(np.outer(z, z).reshape(-1))
    _, singular, vh = np.linalg.svd(np.stack(rows), full_matrices=False)
    rank = int(np.sum(singular > 1e-10))
    assert rank == 28, (rank, singular)
    return [row.reshape(10, 10) for row in vh[:rank]], singular


ORBIT_BASIS, ORBIT_SINGULAR = moment_basis()


def contraction(G):
    """One-particle contraction in the canonical E0,...,E4 basis."""
    out = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            for k in range(5):
                if k in (i, j):
                    continue
                left = (i, k) if i < k else (k, i)
                right = (j, k) if j < k else (k, j)
                left_sign = 1 if i < k else -1
                right_sign = 1 if j < k else -1
                out[i, j] += (
                    left_sign
                    * right_sign
                    * G[PAIR_INDEX[left], PAIR_INDEX[right]]
                )
    return out


def eigenvector_slice():
    equations = np.array(
        [[contraction(G)[0, i] for G in ORBIT_BASIS] for i in range(1, 5)]
    )
    kernel = null_space(equations)
    return [
        sum(kernel[a, k] * ORBIT_BASIS[a] for a in range(len(ORBIT_BASIS)))
        for k in range(kernel.shape[1])
    ], equations


BASIS, EIGENVECTOR_EQUATIONS = eigenvector_slice()
DIMENSION = len(BASIS)


HODGE = np.zeros((6, 6))
outer_index = {edge: i for i, edge in enumerate(OUTER)}
for edge, complement, sign in [
    ((1, 2), (3, 4), 1),
    ((1, 3), (2, 4), -1),
    ((1, 4), (2, 3), 1),
]:
    i = outer_index[edge]
    j = outer_index[complement]
    HODGE[i, j] = HODGE[j, i] = sign


def components(G):
    A = G[:4, :4]
    C = G[:4, 4:]
    B = G[4:, 4:]
    physical = np.zeros((10, 10))
    physical[4:, 4:] = B
    R = contraction(physical)[1:, 1:]
    mass = np.trace(G)
    delta = 3 * np.trace(A) - 2 * mass
    U = A - R - np.eye(4) * delta / 4
    return A, C, B, R, delta, U


def target(A, B):
    _, CA, BA, _, delta_a, UA = components(A)
    _, CB, BB, _, delta_b, UB = components(B)
    hodge = np.trace(BA @ HODGE @ BB @ HODGE)
    return (
        12 * np.sum(UA * UB)
        + 48 * np.sum(CA * CB)
        + 24 * hodge
        - delta_a * delta_b
    )


TARGET = np.array([[target(A, B) for B in BASIS] for A in BASIS])
TARGET = (TARGET + TARGET.T) / 2


def polarized_compound(A, B):
    out = np.zeros((45, 45))
    for p, (i, j) in enumerate(MINORS):
        for q, (k, ell) in enumerate(MINORS):
            out[p, q] = 0.5 * (
                A[i, k] * B[j, ell]
                + B[i, k] * A[j, ell]
                - A[i, ell] * B[j, k]
                - B[i, ell] * A[j, k]
            )
    return (out + out.T) / 2


COMPOUND = {
    (i, j): polarized_compound(BASIS[i], BASIS[j])
    for i in range(DIMENSION)
    for j in range(i, DIMENSION)
}


def solve():
    W = cp.Variable((45, 45), symmetric=True)
    S = cp.Variable((DIMENSION, DIMENSION), symmetric=True)
    constraints = [W >> 0, S >> 0]
    for i in range(DIMENSION):
        for j in range(i, DIMENSION):
            constraints.append(
                S[i, j] + cp.sum(cp.multiply(W, COMPOUND[i, j]))
                == TARGET[i, j]
            )
    problem = cp.Problem(cp.Minimize(cp.trace(W) + cp.trace(S)), constraints)
    try:
        value = problem.solve(
            solver="CLARABEL",
            tol_gap_abs=1e-9,
            tol_feas=1e-9,
            tol_gap_rel=1e-9,
            max_iter=3000,
        )
    except cp.error.SolverError:
        value = problem.solve(solver="SCS", eps=3e-7, max_iters=500000)
    print(
        "orbit/eigen dimensions",
        len(ORBIT_BASIS),
        DIMENSION,
        "target eig",
        np.linalg.eigvalsh(TARGET),
        "status/value",
        problem.status,
        value,
    )
    if W.value is not None:
        print("W eig", np.linalg.eigvalsh(W.value))
        print("S eig", np.linalg.eigvalsh(S.value))
    return problem.status, W.value, S.value


if __name__ == "__main__":
    solve()
