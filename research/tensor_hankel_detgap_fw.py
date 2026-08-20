"""Frank-Wolfe audit of the determinant gap on the PSD Hankel relaxation."""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.optimize import minimize_scalar

from tensor_detgap_opt import Q as Q_PHYSICAL
from tensor_hankel_fw import (
    deg6,
    direction,
    energy,
    gradient,
    maps,
    normalizer,
    uniform_moment,
)


def q_coordinates(Q):
    # maps[2] basis order: z^2,yz,y^2,xz,xy,x^2 with off-diagonal
    # tensors normalized by sqrt(2).
    return np.array(
        [Q[2, 2], np.sqrt(2) * Q[1, 2], Q[1, 1], np.sqrt(2) * Q[0, 2], np.sqrt(2) * Q[0, 1], Q[0, 0]]
    )


QCOORD = np.stack([q_coordinates(Q) for Q in Q_PHYSICAL])
BMAP = np.zeros((5, 5, len(deg6)))
for k in range(len(deg6)):
    M = maps[1][:, :, k]
    C2 = maps[2][:, :, k]
    for i, Qi in enumerate(Q_PHYSICAL):
        for j, Qj in enumerate(Q_PHYSICAL):
            first = np.trace(Qi @ Qj @ M)
            second = QCOORD[i] @ C2 @ QCOORD[j]
            BMAP[i, j, k] = 5 * (first - second)


y_variable = cp.Variable(len(deg6))
direction = cp.Parameter(len(deg6))
hankel = sum(y_variable[k] * maps[3][:, :, k] for k in range(len(deg6)))
bhat_cvx = sum(y_variable[k] * BMAP[:, :, k] for k in range(len(deg6)))
rho1_cvx = sum(y_variable[k] * maps[1][:, :, k] for k in range(len(deg6)))
problem = cp.Problem(
    cp.Minimize(direction @ y_variable),
    [
        hankel >> 0,
        bhat_cvx >> 0,
        rho1_cvx == np.eye(3) / 3,
        normalizer @ y_variable == 1,
    ],
)


def bhat(y):
    return np.einsum("abk,k->ab", BMAP, y)


def determinant_and_gradient(y):
    B = bhat(y)
    determinant = np.linalg.det(B)
    if abs(determinant) > 1e-10 and np.linalg.cond(B) < 1e10:
        adjugate = determinant * np.linalg.inv(B)
    else:
        # Stable cofactor definition, valid at singular matrices too.
        adjugate = np.empty((5, 5))
        for i in range(5):
            for j in range(5):
                minor = np.delete(np.delete(B, j, axis=0), i, axis=1)
                adjugate[i, j] = (-1) ** (i + j) * np.linalg.det(minor)
    grad = np.einsum("ab,abk->k", adjugate.T, BMAP)
    return determinant, grad


def objective(y):
    determinant, _ = determinant_and_gradient(y)
    return energy(y) - (32 / 105) * determinant**2


def objective_gradient(y):
    determinant, det_gradient = determinant_and_gradient(y)
    return gradient(y) - (64 / 105) * determinant * det_gradient


def frank_wolfe(start, iterations=500):
    y = start.copy()
    best = objective(y)
    for iteration in range(iterations):
        direction.value = objective_gradient(y)
        try:
            problem.solve(solver=cp.CLARABEL)
        except cp.error.SolverError:
            problem.solve(solver=cp.SCS, eps=1e-8, max_iters=200000)
        vertex = y_variable.value.copy()
        delta = vertex - y
        line = minimize_scalar(
            lambda step: objective(y + step * delta),
            bounds=(0, 1),
            method="bounded",
            options={"xatol": 1e-12},
        )
        step = line.x if line.fun < best - 1e-13 else 0.0
        if step == 0:
            break
        y += step * delta
        best = objective(y)
    return y, best, iteration


def main():
    rng = np.random.default_rng(20260820)
    uniform = np.array([uniform_moment(alpha) for alpha in deg6])
    starts = [uniform]
    for _ in range(30):
        direction.value = rng.normal(size=len(deg6))
        problem.solve(solver=cp.CLARABEL)
        starts.append(y_variable.value.copy())
    for i, start in enumerate(starts):
        y, value, iteration = frank_wolfe(start)
        B = bhat(y)
        print(
            i,
            iteration,
            value,
            "E",
            energy(y),
            "det",
            np.linalg.det(B),
            "Beig",
            np.linalg.eigvalsh(B),
            "Heig",
            np.linalg.eigvalsh(np.einsum("abk,k->ab", maps[3], y)),
            flush=True,
        )


if __name__ == "__main__":
    main()
