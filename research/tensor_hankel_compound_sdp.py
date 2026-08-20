"""Search for a 2x2-minor certificate on the PSD Hankel relaxation.

Let H(y) be the 10x10 cubic catalecticant of a ternary sextic moment
functional.  The homogeneous target is a quadratic form Q(y).  This script
asks for

    Q(y) = y^T S y + tr(W * C2(H(y))),   S >= 0, W >= 0,

where C2 is the second compound matrix.  Both terms are nonnegative whenever
H(y)>=0, so a feasible solution is a complete proof on the full PSD-Hankel
relaxation (and hence for every measure).
"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

from tensor_hankel_fw import maps


N = maps[3].shape[2]
H_BASIS = [maps[3][:, :, k] for k in range(N)]
PAIR_SETS = {
    degree: [(i, j) for i in range(maps[degree].shape[0]) for j in range(i + 1, maps[degree].shape[0])]
    for degree in (1, 2, 3)
}


def gram_of_map(tensor):
    return np.einsum("abk,abl->kl", tensor, tensor)


NORMALIZER = np.trace(maps[3], axis1=0, axis2=1)
TARGET = (
    32 * gram_of_map(maps[3])
    - 48 * gram_of_map(maps[2])
    + 20 * gram_of_map(maps[1])
    - (4 / 3) * np.outer(NORMALIZER, NORMALIZER)
)
TARGET = (TARGET + TARGET.T) / 2


def polarized_compound(A, B, pairs):
    """Coefficient of st in C2(sA+tB), divided by two for a=b convention.

    The returned symmetric 45x45 matrix C satisfies
      C2(H(y)) = sum_a y_a^2 C_aa + 2 sum_{a<b} y_a y_b C_ab.
    """
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


COMPOUND_COEFFICIENTS = {degree: {} for degree in (1, 2, 3)}
for degree in (1, 2, 3):
    basis = [maps[degree][:, :, k] for k in range(N)]
    for a in range(N):
        for b in range(a, N):
            COMPOUND_COEFFICIENTS[degree][a, b] = polarized_compound(
                basis[a], basis[b], PAIR_SETS[degree]
            )


def solve(solver="CLARABEL"):
    W = {
        degree: cp.Variable((len(PAIR_SETS[degree]), len(PAIR_SETS[degree])), symmetric=True)
        for degree in (1, 2, 3)
    }
    S = cp.Variable((N, N), symmetric=True)
    constraints = [S >> 0] + [W[degree] >> 0 for degree in (1, 2, 3)]
    for a in range(N):
        for b in range(a, N):
            constraints.append(
                S[a, b]
                + sum(
                    cp.sum(
                        cp.multiply(
                            W[degree], COMPOUND_COEFFICIENTS[degree][a, b]
                        )
                    )
                    for degree in (1, 2, 3)
                )
                == TARGET[a, b]
            )
    problem = cp.Problem(
        cp.Minimize(sum(cp.trace(W[d]) for d in (1, 2, 3)) + cp.trace(S)),
        constraints,
    )
    try:
        value = problem.solve(
            solver=solver,
            tol_gap_abs=1e-9,
            tol_feas=1e-9,
            tol_gap_rel=1e-9,
            max_iter=2000,
        )
    except (cp.error.SolverError, TypeError):
        value = problem.solve(solver="SCS", eps=2e-7, max_iters=500000, verbose=True)
    print("status", problem.status, "value", value)
    if W[3].value is not None:
        for degree in (1, 2, 3):
            print("W", degree, "eig", np.linalg.eigvalsh(W[degree].value))
        print("S eig", np.linalg.eigvalsh(S.value))
        residual = np.zeros_like(TARGET)
        for a in range(N):
            for b in range(a, N):
                residual[a, b] = (
                    S.value[a, b]
                    + sum(
                        np.sum(
                            W[degree].value
                            * COMPOUND_COEFFICIENTS[degree][a, b]
                        )
                        for degree in (1, 2, 3)
                    )
                    - TARGET[a, b]
                )
                residual[b, a] = residual[a, b]
        print("max residual", np.max(np.abs(residual)))
        np.savez(
            "research/tensor_hankel_compound_solution.npz",
            W1=W[1].value,
            W2=W[2].value,
            W3=W[3].value,
            S=S.value,
            residual=residual,
        )
    return problem.status, {degree: W[degree].value for degree in W}, S.value


if __name__ == "__main__":
    print("target eig", np.linalg.eigvalsh(TARGET))
    solve()
