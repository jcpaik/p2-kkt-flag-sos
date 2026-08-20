"""Search a structured Schur-dual witness Y=B W B^T, W>=0.

For R=M+A6>=0 with A6 in the spin-6 operator subspace, every Y>=0 gives

  ||A6||^2 >= -<Y,M> - ||Proj_6 Y||^2/4.

This script optimizes that bound over Y=B W B^T, where B:H1->H3 is the
spin-4 cross block and W is an arbitrary 3-by-3 PSD matrix.
"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

import tensor_schur_feasible_scan as data


def symmetric_basis3():
    out = []
    for i in range(3):
        M = np.zeros((3, 3))
        M[i, i] = 1
        out.append(M)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        M = np.zeros((3, 3))
        M[i, j] = M[j, i] = 1
        out.append(M)
    return out


W_BASIS = symmetric_basis3()


def dual_bound(A4, B, b):
    t = np.sqrt(b)
    Bactual = t * B
    Q = Bactual @ Bactual.T
    M = (2 / 35) * np.eye(7) + t * A4 - 5 * Q
    ybasis = [Bactual @ W @ Bactual.T for W in W_BASIS]
    y6basis = [data.proj3(Y, 6) for Y in ybasis]

    w = cp.Variable(6)
    W = sum(w[i] * W_BASIS[i] for i in range(6))
    linear = np.array([-np.sum(Y * M) for Y in ybasis]) @ w
    y6 = sum(w[i] * y6basis[i] for i in range(6))
    problem = cp.Problem(cp.Maximize(linear - cp.sum_squares(y6) / 4), [W >> 0])
    try:
        value = problem.solve(solver="CLARABEL")
    except cp.error.SolverError:
        return None
    return value, W.value


def two_generator_bound(A4, B, b):
    """Optimize a manifestly PSD polynomial witness in L=A4-5BB*."""
    t = np.sqrt(b)
    A4actual = t * A4
    Q = b * (B @ B.T)
    M = (2 / 35) * np.eye(7) + A4actual - 5 * Q
    L = A4actual - 5 * Q
    powers = [L, L @ L, L @ L @ L]
    gram_basis = []
    for W0 in W_BASIS:
        Y0 = np.zeros((7, 7))
        for i in range(3):
            for j in range(3):
                Y0 += W0[i, j] * (powers[i] @ powers[j])
        gram_basis.append((Y0 + Y0.T) / 2)
    ybasis = [Q] + gram_basis
    y6basis = [data.proj3(Y, 6) for Y in ybasis]
    alpha = cp.Variable(nonneg=True)
    gram = cp.Variable((3, 3), symmetric=True)
    coeff = cp.hstack(
        [alpha, gram[0, 0], gram[1, 1], gram[2, 2], gram[0, 1], gram[0, 2], gram[1, 2]]
    )
    # W_BASIS order is diagonal entries followed by (01),(02),(12).
    Y6 = sum(coeff[i] * y6basis[i] for i in range(len(ybasis)))
    linear = np.array([-np.sum(Y * M) for Y in ybasis]) @ coeff
    problem = cp.Problem(
        cp.Maximize(linear - cp.sum_squares(Y6) / 4), [gram >> 0]
    )
    value = problem.solve(solver="CLARABEL")
    return value, np.hstack([alpha.value, gram.value.reshape(-1)])


def main():
    rng = np.random.default_rng(20260820)
    worst = None
    for trial in range(30):
        coeff = rng.normal(size=9)
        A4, B, C4, inv = data.direction(coeff)
        for b in np.linspace(0.0004, 0.0796, 60):
            eig = np.linalg.eigvalsh(np.eye(5) - 5 * np.sqrt(b) * C4)
            if eig[0] < -1e-9:
                continue
            result = two_generator_bound(A4, B, b)
            if result is None:
                continue
            value, W = result
            determinant = np.prod(eig)
            gap = (
                2 / 105
                - 5 * b / 11
                + 2 * value
                - (2 / 105) * determinant**2
            )
            record = (gap, trial, b, value, determinant, eig, W, coeff)
            if worst is None or gap < worst[0]:
                worst = record
                print("worst", worst[:-1], flush=True)


if __name__ == "__main__":
    main()
