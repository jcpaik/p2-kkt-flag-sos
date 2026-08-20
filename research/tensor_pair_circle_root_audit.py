"""Rootwise audit of J_x-108D_x-288A_x."""

import numpy as np

from tensor_pair_circle_projection import circle_square


def rooted_values(points, weights):
    values = []
    for x in points:
        inner = points @ x
        potential_j = np.sum(
            weights * (144 * inner**6 - 216 * inner**4 + 87 * inner**2 - 5)
        )
        seed = np.eye(3)[np.argmin(np.abs(x))]
        u = np.cross(x, seed)
        u /= np.linalg.norm(u)
        v = np.cross(x, u)
        t = points @ x
        z = points @ u + 1j * (points @ v)
        radial = np.sum(weights * t**2 * np.abs(z)**4)
        spin = np.sum(weights * t**2 * z**4)
        determinant = radial**2 - abs(spin)**2
        pair = 0.0
        for y, weight in zip(points, weights):
            a = x @ y
            pair += weight * a**4 * (1 - a*a)**6 * circle_square(
                points, weights, x, y
            )
        values.append((potential_j, determinant, pair,
                       potential_j - 108 * determinant - 288 * pair))
    return np.asarray(values)


if __name__ == "__main__":
    rng = np.random.default_rng(20260822)
    for _ in range(10):
        points = rng.normal(size=(8, 3))
        points /= np.linalg.norm(points, axis=1)[:, None]
        weights = rng.dirichlet(np.ones(8))
        value = rooted_values(points, weights)
        print(np.min(value[:, -1]), weights @ value[:, -1])
