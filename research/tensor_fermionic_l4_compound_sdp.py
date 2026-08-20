"""Minor/SOS certificate search for the spin-2/spin-4 fermionic relaxation.

The relaxation consists of PSD G on H1+H3, the homogeneous block trace
ratio 5 tr(G11)=tr(G), and the exact tangent-orbit Clebsch relations in spin
2 and spin 4.  Numerically its purity gap is nonnegative.  We search for a
certificate on the resulting linear subspace:

  Q(G)=||G||^2-||2Tr_2G||^2/2+(trG)^2/3
      = z^T S z + tr(W C2(G)),  S,W >= 0.
"""

import sys
import os

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

import tensor_fermionic_general_relaxation_opt as relations
import tensor_fermionic_relaxation_opt as base


def symmetric_orthonormal_basis(n):
    out = []
    for i in range(n):
        E = np.zeros((n, n))
        E[i, i] = 1
        out.append(E)
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n))
            E[i, j] = E[j, i] = 1 / np.sqrt(2)
            out.append(E)
    return out


SYM = symmetric_orthonormal_basis(10)
EXPANSION = np.stack([E.reshape(-1) for E in SYM], axis=1)

# Homogeneous linear constraints.
trace_ratio = np.zeros((1, 100))
for i in range(3):
    trace_ratio[0, 10 * i + i] = 5
for i in range(10):
    trace_ratio[0, 10 * i + i] -= 1
constraints_full = np.vstack(
    [
        trace_ratio,
        relations.LINEAR_CONSTRAINT_REDUCED,
        relations.L4_REDUCED,
    ]
)
constraint_coordinates = constraints_full @ EXPANSION
BASIS_COORDINATES = null_space(constraint_coordinates)
BASIS = [
    sum(BASIS_COORDINATES[a, k] * SYM[a] for a in range(len(SYM)))
    for k in range(BASIS_COORDINATES.shape[1])
]
DIMENSION = len(BASIS)


def contraction(G):
    physical = relations.T @ G @ relations.T.T
    return (base.GAMMA.numpy() @ physical.reshape(-1)).reshape(5, 5)


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

PAIRS = [(i, j) for i in range(10) for j in range(i + 1, 10)]


def polarized_compound(A, B):
    out = np.zeros((45, 45))
    for p, (i, j) in enumerate(PAIRS):
        for q, (k, ell) in enumerate(PAIRS):
            out[p, q] = 0.5 * (
                A[i, k] * B[j, ell]
                + B[i, k] * A[j, ell]
                - A[i, ell] * B[j, k]
                - B[i, ell] * A[j, k]
            )
    return (out + out.T) / 2


COMPOUND = {}
for a in range(DIMENSION):
    for b in range(a, DIMENSION):
        COMPOUND[a, b] = polarized_compound(BASIS[a], BASIS[b])


def solve():
    W = cp.Variable((45, 45), symmetric=True)
    S = cp.Variable((DIMENSION, DIMENSION), symmetric=True)
    constraints = [W >> 0, S >> 0]
    for a in range(DIMENSION):
        for b in range(a, DIMENSION):
            constraints.append(
                S[a, b] + cp.sum(cp.multiply(W, COMPOUND[a, b])) == Q[a, b]
            )
    problem = cp.Problem(cp.Minimize(cp.trace(W) + cp.trace(S)), constraints)
    try:
        value = problem.solve(
            solver="CLARABEL",
            tol_gap_abs=1e-9,
            tol_feas=1e-9,
            tol_gap_rel=1e-9,
            max_iter=2000,
        )
    except cp.error.SolverError:
        value = problem.solve(solver="SCS", eps=2e-7, max_iters=500000)
    print("weighted", WEIGHTED, "dimension", DIMENSION, "constraint rank", np.linalg.matrix_rank(constraint_coordinates))
    print("Q eig", np.linalg.eigvalsh(Q))
    print("status", problem.status, "value", value)
    if W.value is not None:
        print("W eig", np.linalg.eigvalsh(W.value))
        print("S eig", np.linalg.eigvalsh(S.value))
        np.savez(
            "research/tensor_fermionic_l4_compound_solution.npz",
            W=W.value,
            S=S.value,
            basis=np.stack(BASIS),
            Q=Q,
        )
    return problem.status, W.value, S.value


if __name__ == "__main__":
    solve()
