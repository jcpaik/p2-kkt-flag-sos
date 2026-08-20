"""CCCP searches on symmetry-fixed faces of the 40D relaxation."""

import sys

sys.path.insert(0, "research")

import cvxpy as cp
import numpy as np
from scipy.linalg import expm

import tensor_fermionic_general_relaxation_opt as rel
from tensor_weighted_h2_relaxation_audit import fermionic_weighted_value


GENERATORS = [
    np.block(
        [
            [rel.wedge.gw1[a], np.zeros((3, 7))],
            [np.zeros((7, 3)), rel.wedge.gw3[a]],
        ]
    )
    for a in range(3)
]


def rotation(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    return expm(angle * sum(axis[a] * GENERATORS[a] for a in range(3)))


GROUPS = {
    "so2": [GENERATORS[2]],  # commutation with the infinitesimal generator
    "c2": [rotation([0, 0, 1], np.pi)],
    "d2": [rotation([1, 0, 0], np.pi), rotation([0, 1, 0], np.pi)],
    "d4": [rotation([0, 0, 1], np.pi / 2), rotation([1, 0, 0], np.pi)],
    "tetra": [rotation([1, 1, 1], 2 * np.pi / 3), rotation([1, 0, 0], np.pi)],
    "octa": [rotation([0, 0, 1], np.pi / 2), rotation([1, 0, 0], np.pi / 2)],
}


def solve(group_name, weighted=False, restarts=40, iterations=100, seed=23):
    rng = np.random.default_rng(seed + 100 * weighted)
    X = cp.Variable((10, 10), symmetric=True)
    xvec = cp.reshape(X, (100,), order="C")
    constraints = [
        X >> 0,
        cp.trace(X) == 1,
        cp.trace(X[:3, :3]) == 1 / 5,
        rel.LINEAR_CONSTRAINT_REDUCED @ xvec == 0,
        rel.L4_REDUCED @ xvec == 0,
    ]
    for item in GROUPS[group_name]:
        if group_name == "so2":
            constraints.append(item @ X - X @ item == 0)
        else:
            constraints.append(item @ X @ item.T == X)

    direction = cp.Parameter((10, 10), symmetric=True)
    extreme = cp.Problem(cp.Minimize(cp.sum(cp.multiply(direction, X))), constraints)
    best = None
    for restart in range(restarts):
        raw = rng.normal(size=(10, 10))
        direction.value = (raw + raw.T) / 2
        try:
            extreme.solve(solver="CLARABEL")
        except cp.error.SolverError:
            extreme.solve(solver="SCS", eps=1e-8, max_iters=300000)
        current = X.value
        for _ in range(iterations):
            F0 = (rel.F_LINEAR @ current.reshape(-1)).reshape(5, 5)
            linear = (rel.F_LINEAR.T @ F0.reshape(-1)).reshape(10, 10)
            center = (linear + linear.T) / 4
            if weighted:
                center[:3, :3] += (25 / 12) * (
                    current[:3, :3] - np.eye(3) / 15
                )
            problem = cp.Problem(cp.Minimize(cp.sum_squares(X - center)), constraints)
            try:
                problem.solve(solver="CLARABEL")
            except cp.error.SolverError:
                problem.solve(solver="SCS", eps=1e-8, max_iters=300000)
            current = X.value
        F = (rel.F_LINEAR @ current.reshape(-1)).reshape(5, 5)
        Q = np.sum(current * current) - np.sum(F * F) / 2 + 1 / 3
        value = fermionic_weighted_value(current) / 8 if weighted else Q
        if best is None or value < best[0]:
            best = value, current
            print(
                group_name,
                "weighted" if weighted else "plain",
                restart,
                value,
                np.linalg.eigvalsh(current),
                flush=True,
            )
    return best


if __name__ == "__main__":
    for name in GROUPS:
        solve(name)
        solve(name, weighted=True)
