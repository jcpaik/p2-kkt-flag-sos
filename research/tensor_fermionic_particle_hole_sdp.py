"""Rejected quadratic certificate search using a false extra cone.

For a convex mixture of decomposable Pluecker states G, with one-body
marginal F,

    K(G) = F wedge I - 2 G >= 0.

This implication is false for an arbitrary two-fermion density and is not a
consequence of the 40-dimensional spin-2/spin-4 relaxation.  For example,
G=|omega><omega| with omega=(e_12+e_34)/sqrt(2) gives lambda_min(K)=-1.
Accordingly the SDP below is retained only as an audit of a certificate that
would require the missing Slater separability; it is not a valid proof search
for the abstract relaxation.

On the homogeneous spin-2/spin-4 tangent slice we test

    Q(G) = square(G) + <W, G tensor K(G)>,   W >= 0.

The tensor multiplier is deliberately unrestricted in this diagnostic.
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


def particle_hole(G):
    """Return K=F wedge I-2G in the H1+H3 block basis."""
    F = old.contraction(G)
    physical = T @ G @ T.T
    Kphysical = wedge.additive_wedge(F) - 2 * physical
    return (T.T @ Kphysical @ T + T.T @ Kphysical.T @ T) / 2


KS = [particle_hole(G) for G in BASIS]
GK = {}
for a in range(DIMENSION):
    for b in range(a, DIMENSION):
        GK[a, b] = 0.5 * (
            np.kron(BASIS[a], KS[b]) + np.kron(BASIS[b], KS[a])
        )


def solve():
    S = cp.Variable((DIMENSION, DIMENSION), symmetric=True)
    W = cp.Variable((100, 100), symmetric=True)
    constraints = [S >> 0, W >> 0]
    for a in range(DIMENSION):
        for b in range(a, DIMENSION):
            constraints.append(
                S[a, b] + cp.sum(cp.multiply(W, GK[a, b])) == Q[a, b]
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
            "research/tensor_fermionic_particle_hole_solution.npz",
            S=S.value,
            W=W.value,
            basis=np.stack(BASIS),
            Q=Q,
        )


if __name__ == "__main__":
    solve()
