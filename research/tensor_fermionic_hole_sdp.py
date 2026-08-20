"""Degree-two certificate using the valid two-hole PSD cone.

For every PSD two-fermion density G on wedge^2(R^5), its two-hole matrix

    H(G) = tr(G) I - (F wedge I) + G

is PSD.  On the 40-dimensional spin-2/spin-4 slice this script tests

    Q(G) = square(G) + <W, G tensor H(G)>,  W >= 0.

Unlike K=F wedge I-2G, positivity of H is valid for entangled states.
"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

import tensor_fermionic_l4_compound_sdp as old
import tensor_wedge_relation as wedge


BASIS = old.BASIS
DIMENSION = old.DIMENSION
Q = old.Q
T = old.relations.T


def two_hole(G):
    F = old.contraction(G)
    physical = T @ G @ T.T
    Hphysical = np.trace(G) * np.eye(10) - wedge.additive_wedge(F) + physical
    return (T.T @ Hphysical @ T + T.T @ Hphysical.T @ T) / 2


HS = [two_hole(G) for G in BASIS]
GH = {}
for a in range(DIMENSION):
    for b in range(a, DIMENSION):
        GH[a, b] = 0.5 * (
            np.kron(BASIS[a], HS[b]) + np.kron(BASIS[b], HS[a])
        )


def solve():
    S = cp.Variable((DIMENSION, DIMENSION), symmetric=True)
    W = cp.Variable((100, 100), symmetric=True)
    constraints = [S >> 0, W >> 0]
    for a in range(DIMENSION):
        for b in range(a, DIMENSION):
            constraints.append(
                S[a, b] + cp.sum(cp.multiply(W, GH[a, b])) == Q[a, b]
            )
    problem = cp.Problem(cp.Minimize(cp.trace(S) + cp.trace(W)), constraints)
    try:
        value = problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_feas=2e-9,
            tol_gap_rel=2e-9,
            max_iter=3000,
        )
    except cp.error.SolverError:
        value = problem.solve(solver="SCS", eps=5e-7, max_iters=500000)
    print("dimension", DIMENSION, "status/value", problem.status, value)
    if S.value is not None:
        print("S eig", np.linalg.eigvalsh(S.value))
        print("W eig", np.linalg.eigvalsh(W.value))
        np.savez(
            "research/tensor_fermionic_hole_solution.npz",
            S=S.value,
            W=W.value,
            basis=np.stack(BASIS),
            Q=Q,
        )


if __name__ == "__main__":
    solve()
