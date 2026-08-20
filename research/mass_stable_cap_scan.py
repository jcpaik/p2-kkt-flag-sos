"""Scan the conjectural weight-stability dichotomy E>=0 or fusion cap."""

import numpy as np

from weighted_certificate_fermion_cap import cap_matrix


def kernel(t):
    return 32 * t**6 - 48 * t**4 + 20 * t**2 - 4 / 3


def stationary_weights(matrix):
    inverse_one = np.linalg.solve(matrix, np.ones(len(matrix)))
    return inverse_one / np.sum(inverse_one)


def scan(restarts=200_000, seed=20260822):
    rng = np.random.default_rng(seed)
    counts = {}
    best = None
    for atoms in range(2, 11):
        accepted = 0
        for _ in range(restarts // 9):
            points = rng.normal(size=(atoms, 3))
            points /= np.linalg.norm(points, axis=1)[:, None]
            matrix = kernel(points @ points.T)
            try:
                weights = stationary_weights(matrix)
            except np.linalg.LinAlgError:
                continue
            if np.min(weights) <= 1e-9:
                continue
            energy = 1 / np.sum(np.linalg.solve(matrix, np.ones(atoms)))
            basis = np.eye(atoms)[:, :-1] - np.eye(atoms)[:, -1, None]
            if np.min(np.linalg.eigvalsh(basis.T @ matrix @ basis)) < -1e-8:
                continue
            accepted += 1
            gap = np.min(np.linalg.eigvalsh(cap_matrix(points, weights)))
            score = max(-energy, 0) * max(-gap, 0)
            if best is None or energy < best[0] or score > best[2]:
                best = energy, gap, score, points, weights
        counts[atoms] = accepted
        print(atoms, accepted, None if best is None else best[:3], flush=True)
    return counts, best


if __name__ == "__main__":
    print(scan())
