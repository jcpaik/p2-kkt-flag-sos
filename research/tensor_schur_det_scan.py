"""Does the isotropic cubic-block Schur relaxation imply the det-gap bound?"""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np

import tensor_schur_feasible_scan as data


def scan_direction(coeff):
    A4, B, C4, inv = data.direction(coeff)
    Q = B @ B.T
    z = cp.Variable(13)
    A6 = sum(z[i] * data.basis6[i] for i in range(13))
    tpar = cp.Parameter(nonneg=True)
    R = (2 / 35) * np.eye(7) + tpar * A4 + A6 - 5 * tpar * tpar * Q
    problem = cp.Problem(cp.Minimize(cp.sum_squares(A6)), [R >> 0])
    best = None
    for b in np.linspace(0, 0.08, 101):
        t = np.sqrt(b)
        tpar.value = t
        try:
            value = problem.solve(solver="CLARABEL", warm_start=True)
        except cp.error.SolverError:
            continue
        if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            continue
        determinant = np.linalg.det(np.eye(5) - 5 * t * C4)
        gap = 2 / 105 - 5 * b / 11 + 2 * value - (2 / 105) * determinant**2
        result = (gap, b, value, determinant, np.linalg.eigvalsh(np.eye(5) - 5*t*C4))
        if best is None or result[0] < best[0]:
            best = result
    return best


def main():
    rng = np.random.default_rng(20260820)
    worst = None
    for i in range(80):
        coeff = rng.normal(size=9)
        result = scan_direction(coeff)
        if result is not None and (worst is None or result[0] < worst[0]):
            worst = (result[0], i, result, coeff)
            print("worst", worst, flush=True)


if __name__ == "__main__":
    main()
