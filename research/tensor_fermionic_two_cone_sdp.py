"""Degree-two certificate search using both fermionic PSD cones.

On the homogeneous tangent linear slice, both

    G >= 0,  L(G):=P(5G_11)-F(G) >= 0

hold.  Besides compounds C2(G), a quadratic Positivstellensatz may therefore
use C2(L) and G tensor L.  This script tests the exact coefficient ansatz

    Q = z^T S z + <Wg,C2(G)> + <Wl,C2(L)> + <Wgl,G tensor L>,

with all multiplier matrices PSD.
"""

import sys
import os

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

import tensor_fermionic_general_relaxation_opt as relations
import tensor_fermionic_l4_compound_sdp as old


SYM = old.SYM
EXPANSION = old.EXPANSION

trace_ratio = np.zeros((1, 100))
for i in range(3):
    trace_ratio[0, 10 * i + i] = 5
for i in range(10):
    trace_ratio[0, 10 * i + i] -= 1

# The upper-LMI route only needs the trace ratio and l=4 matching relation.
constraints_full = np.vstack([trace_ratio, relations.L4_REDUCED])
constraint_coordinates = constraints_full @ EXPANSION
basis_coordinates = null_space(constraint_coordinates)
BASIS = [
    sum(basis_coordinates[a, k] * SYM[a] for a in range(len(SYM)))
    for k in range(basis_coordinates.shape[1])
]
DIMENSION = len(BASIS)


def contraction(G):
    return old.contraction(G)


def upper_slack(G):
    F = contraction(G)
    P = (relations.UPPER_NP @ G[:3, :3].reshape(-1)).reshape(5, 5)
    return (P - F + (P - F).T) / 2


Q = np.zeros((DIMENSION, DIMENSION))
WEIGHTED = os.environ.get("FERMIONIC_WEIGHTED", "0") == "1"
for a, A in enumerate(BASIS):
    FA = contraction(A)
    Adev = A[:3, :3] - np.eye(3) * np.trace(A) / 15
    for b, B in enumerate(BASIS):
        FB = contraction(B)
        Bdev = B[:3, :3] - np.eye(3) * np.trace(B) / 15
        Q[a, b] = (
            np.sum(A * B)
            - 0.5 * np.sum(FA * FB)
            + np.trace(A) * np.trace(B) / 3
            - (25 / 12) * np.sum(Adev * Bdev) * WEIGHTED
        )
Q = (Q + Q.T) / 2


def compound_pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def polarized_compound(A, B):
    pairs = compound_pairs(A.shape[0])
    out = np.zeros((len(pairs), len(pairs)))
    for p, (i, j) in enumerate(pairs):
        for q, (k, ell) in enumerate(pairs):
            out[p, q] = 0.5 * (
                A[i, k] * B[j, ell]
                + B[i, k] * A[j, ell]
                - A[i, ell] * B[j, k]
                - B[i, ell] * A[j, k]
            )
    return (out + out.T) / 2


LS = [upper_slack(A) for A in BASIS]
CG = {}
CL = {}
GL = {}
for a in range(DIMENSION):
    for b in range(a, DIMENSION):
        CG[a, b] = polarized_compound(BASIS[a], BASIS[b])
        CL[a, b] = polarized_compound(LS[a], LS[b])
        GL[a, b] = 0.5 * (
            np.kron(BASIS[a], LS[b]) + np.kron(BASIS[b], LS[a])
        )


def solve(use_cg=False, use_cl=True, use_gl=True):
    S = cp.Variable((DIMENSION, DIMENSION), symmetric=True)
    variables = [("S", S)]
    constraints = [S >> 0]
    Wg = Wl = Wgl = None
    if use_cg:
        Wg = cp.Variable((45, 45), symmetric=True)
        variables.append(("Wg", Wg))
        constraints.append(Wg >> 0)
    if use_cl:
        Wl = cp.Variable((10, 10), symmetric=True)
        variables.append(("Wl", Wl))
        constraints.append(Wl >> 0)
    if use_gl:
        Wgl = cp.Variable((50, 50), symmetric=True)
        variables.append(("Wgl", Wgl))
        constraints.append(Wgl >> 0)

    for a in range(DIMENSION):
        for b in range(a, DIMENSION):
            value = S[a, b]
            if Wg is not None:
                value += cp.sum(cp.multiply(Wg, CG[a, b]))
            if Wl is not None:
                value += cp.sum(cp.multiply(Wl, CL[a, b]))
            if Wgl is not None:
                value += cp.sum(cp.multiply(Wgl, GL[a, b]))
            constraints.append(value == Q[a, b])

    problem = cp.Problem(
        cp.Minimize(sum(cp.trace(variable) for _, variable in variables)),
        constraints,
    )
    try:
        value = problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_feas=2e-9,
            tol_gap_rel=2e-9,
            max_iter=3000,
        )
    except cp.error.SolverError:
        value = problem.solve(solver="SCS", eps=3e-7, max_iters=500000)
    print(
        "dimension/rank",
        DIMENSION,
        np.linalg.matrix_rank(constraint_coordinates),
        "variant",
        "weighted",
        WEIGHTED,
        use_cg,
        use_cl,
        use_gl,
    )
    print("Q eig", np.linalg.eigvalsh(Q))
    print("status/value", problem.status, value)
    for name, variable in variables:
        if variable.value is not None:
            print(name, "eig", np.linalg.eigvalsh(variable.value))
    return problem.status


if __name__ == "__main__":
    solve(use_cg=False, use_cl=True, use_gl=True)
