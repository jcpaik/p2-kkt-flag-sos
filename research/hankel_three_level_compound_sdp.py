"""Minor certificate for the filtered-bosonic three-level inequality."""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import null_space

import hankel_three_level_cp as data
import tensor_fermionic_l4_compound_sdp as compound_data


SYM = compound_data.SYM
EXPANSION = compound_data.EXPANSION
trace_ratio = np.array(
    [5 * np.trace(data.P1 @ E) - np.trace(E) for E in SYM]
)[None, :]
coordinates = null_space(trace_ratio)
BASIS = [
    sum(coordinates[a, k] * SYM[a] for a in range(len(SYM)))
    for k in range(coordinates.shape[1])
]
DIMENSION = len(BASIS)


def bilinear(A, B):
    a = data.levels(A)
    b = data.levels(B)
    return (
        np.sum(A * B)
        - 0.5 * np.sum(a[3] * b[3])
        + np.sum(a[4] * b[4]) / 6
        + a[-1] * b[-1] / 18
    )


Q = np.array([[bilinear(A, B) for B in BASIS] for A in BASIS])
Q = (Q + Q.T) / 2

COMPOUND = {}
for a in range(DIMENSION):
    for b in range(a, DIMENSION):
        COMPOUND[a, b] = compound_data.polarized_compound(BASIS[a], BASIS[b])


def solve():
    S = cp.Variable((DIMENSION, DIMENSION), symmetric=True)
    W = cp.Variable((45, 45), symmetric=True)
    constraints = [S >> 0, W >> 0]
    for a in range(DIMENSION):
        for b in range(a, DIMENSION):
            constraints.append(
                S[a, b] + cp.sum(cp.multiply(W, COMPOUND[a, b])) == Q[a, b]
            )
    problem = cp.Problem(cp.Minimize(cp.trace(S) + cp.trace(W)), constraints)
    for solver in ("CLARABEL", "SCS"):
        try:
            if solver == "CLARABEL":
                value = problem.solve(
                    solver=solver,
                    tol_gap_abs=1e-9,
                    tol_feas=1e-9,
                    tol_gap_rel=1e-9,
                    max_iter=3000,
                )
            else:
                value = problem.solve(solver=solver, eps=2e-7, max_iters=500000)
            print("dimension", DIMENSION, "Q eig", np.linalg.eigvalsh(Q))
            print(solver, problem.status, value)
            if W.value is not None:
                print("W eig", np.linalg.eigvalsh(W.value))
                print("S eig", np.linalg.eigvalsh(S.value))
                np.savez(
                    "research/hankel_three_level_compound_solution.npz",
                    W=W.value,
                    S=S.value,
                    basis=np.stack(BASIS),
                    Q=Q,
                )
            if problem.status in ("optimal", "optimal_inaccurate"):
                break
        except cp.error.SolverError as error:
            print(solver, error)


if __name__ == "__main__":
    solve()
